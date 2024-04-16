# ***********************************************
# Title       : QC peptide processing
# Description : Functions for the processing
#               of QC peptides
# Author      : Karolina Krystofova
# Date        : 2024/03/22
# Version     : 1.0.0
# ***********************************************

# Function for getting already calculated files and those that need recalculation
alreadyCalculatedFiles <- function(s, bool = FALSE, dir, suffix = '.csv'){
  if (bool) {
    load <- NULL
    calc <- s
  } else {
    searched <- sapply(s, function(x) {
      paste0(x, suffix) %in% list.files(dir)
    })
    
    load <- names(searched)[which(searched == TRUE)]
    calc <- names(searched)[which(searched == FALSE)]
  }
  
  return(list(load = load, calc = calc))
}

# Load files that have been calculated previously
loadCalculatedTable <- function(x, dir, sample.name = 'sample', suffix = '.tsv') {
  out <- list()
  
  for (xi in x) {
    fn <- paste0(dir, xi, suffix)
    tmp <- read.table(fn, header = TRUE, sep = ',', strip.white = TRUE) |>
      mutate(!!sample.name := xi)
    out[[length(out)+1]] <- tmp
  }
  
  return(do.call(rbind, out))
  
}

# Load matrix for the acquisition of QC peptides
loadQCReferenceMatrix <- function(path) {
  tbl <- read.table(path, header = TRUE, sep = ',', strip.white = TRUE, fill = FALSE)
  m <- unique(tbl$precursor.mz)
  i <- unique(tbl$id)
  
  return(list(matrix = tbl, mz = m, ids = i))
  
}

# Get name of spectra with highest correlation matrix
qcCorr <- function(s, ref.matrix) {
  # all data from a given spectrum
  tmp <- data.frame(mz = mz(s),
                    intensity = intensity(s))
    
  # create a matrix with real values
  comp <- data.frame(mz = ref.matrix$mz,
                     intensity = 0)
  
  for (j in 1:nrow(comp)) {
    mzi <- comp$mz[j]
    
    matches <- mzi > tmp$mz - 0.015 & mzi < tmp$mz + 0.015
    index <- which(matches == TRUE)
    
    if(length(index == 1)) comp$intensity[j] <- tmp$intensity[index]
    if(length(index > 1)) comp$intensity[j] <- max(tmp$intensity[index], na.rm = TRUE)
  }
  
  # compute correlation coefficient
  if (sum(comp$intensity) == 0) {
    return(NA)
  } else {
    corr <- cor.test(as.matrix(comp), as.matrix(ref.matrix))
    return(unname(corr$estimate))
  }
}

# Process single m/z
processMzOfMzMLfile <- function(x, m, ref.matrix) {
  ref <- ref.matrix$matrix
  mzs <- ref.matrix$mz
  ids <- ref.matrix$ids
  
  # Filter out only data for relevant m/z
  px <- x |> filterPrecursorMz(m, ppm = (0.015/m)*1e6)
  refx <- ref |>
    dplyr::filter(precursor.mz == m) |>
    select(mz, intensity)
  
  # Get info about spectra
  spec <- spectra(px)
  ns <- names(spec)
  
  # Create XIC chromatogram
  xc <- chromatogram(x, mz = c(m-0.015, m+0.015))
  xc[,1]@intensity[is.na(xc[,1]@intensity)] <- 0
  xc[,1]@rtime <- xc[,1]@rtime/60
  
  # Get spectrum with highest correlation
  if(length(spec) > 0) {
    res <- sapply(ns, function(nsx){
      s <- spec[[nsx]]
      qcCorr(s, refx)
    })

    if(sum(!is.na(res)) != 0) {
      poi <- res[which.max(res)]
    } else {
      poi <- NA
    }
    
  } else {
    poi <- NA
  }
  
  # Get info about QC peak
  if(!is.na(poi)) {
    mcor <- unname(poi)
    
    # smooth using Golay-Savitzky
    smooth_int <- sgolayfilt(xc[,1]@intensity, p = 5)
    smooth_peaks <- chromatographicPeakDetector(smooth_int)
    
    # at what RT is peak of interest
    poix <- spec[[names(poi)]]@rt/60
    poix_i <- min(which(unname(xc[,1]@rtime) > poix))
    
    poi_row <- smooth_peaks[which.max(smooth_peaks[,1][smooth_peaks[,1] < poix_i]),]
    poi_start <- unname(xc[,1]@rtime)[poi_row[1]]
    poi_end <- unname(xc[,1]@rtime)[poi_row[2]]
    
    if (poix_i %in% smooth_peaks[,1]) {
      scnd_peak <- smooth_peaks[which(smooth_peaks[,1] == poix_i),]
      poi_end <- unname(xc[,1]@rtime)[scnd_peak[2]]
    }
    
  } else {
    mcor <- NA
    poix <- NA
    poi_start <- NA
    poi_end <- NA
  }
  
  # Create output table
  d <- data.frame(QC = ids[mzs == m],
                  rt = xc[,1]@rtime,
                  intensity = xc[,1]@intensity,
                  matrix_correlation = mcor,
                  peak_start = poi_start,
                  peak_end = poi_end,
                  peak_x = poix)
  
  locmin <- which(islocalminimum(d$intensity)==-1)
  d$baseline <- IPA_baselineDeveloper(locmin, d$intensity)
  rownames(d) <- NULL
  
  return(d)
}

# Process single sample file
processQCfileMzML <- function(x, ref) {
  p <- pickPeaks(x)
  mzs <- ref$mz

  l <- lapply(mzs, function(m){
    processMzOfMzMLfile(p, m, ref)
  })
  
  return(do.call(rbind, l))
}

# Save created dataframe
saveDataCsv <- function(x, d, dir, suffix = '.csv') {
  fn <- paste0(dir, x, suffix)
  write.table(d, fn, sep = ',', append = FALSE,
              row.names = FALSE, quote = FALSE)
  say(fn, type = 'output')
}

# Wrapper function for QC data calculation
calculateQCmzML <- function(x, dir = input, output.dir = qc_dir, 
                            ref = loadQCReferenceMatrix(qc_mx_path), sample.name = 'sample') {
  # read .mzML files
  fn <- paste0(dir, x, '.mzML')
  readMS <- readMSData(fn, mode = "onDisk")
  
  # return success message
  say(paste0('.mzML files for QC peptide XIC chromatogram computation read from ', dir, ' for samples:'))
  for (xi in x) {
    say(xi, type = 'input')
  }
  say('QC peptide data calculated and saved as:')
  
  l <- lapply(x, function(xi) {
    tmp <- readMS |>
      filterFile(paste0(xi, '.mzML'))
    d <- processQCfileMzML(tmp, ref)
    
    saveDataCsv(xi, d, output.dir)
    d <- mutate(d, !!sample.name := xi)
    return(d)
  })
  
  return(do.call(rbind, l))
}