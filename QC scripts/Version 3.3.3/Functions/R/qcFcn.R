# ***********************************************
# Title       : QC peptide processing
# Description : Functions for the processing
#               of QC peptides
# Author      : Karolina Krystofova
# Date        : 2024/03/22
# Version     : 1.0.0
# ***********************************************

#------------ Functions for calculation ------------
euclNorm <- function(x) {
  res <- sqrt(sum(x^2))
  return(x/res)
}

cosine <- function(x, y) {
  mgvx <- sqrt(sum(x^2))
  mgvy <- sqrt(sum(y^2))
  dot <- sum(x*y)
  
  sim <- dot/(mgvx*mgvy)
  
  return(sim)
}

# Save created dataframe
saveDataCsv <- function(x, d, dir, suffix = '.csv') {
  fn <- paste0(dir, x, suffix)
  write.table(d, fn, sep = ',', append = FALSE,
              row.names = FALSE, quote = FALSE)
  say(fn, type = 'output')
}

# Load matrix for the acquisition of QC peptides
loadQCReferenceMatrix <- function(path) {
  tbl <- read.table(path, header = TRUE, sep = ',', strip.white = TRUE, fill = FALSE)
  i <- unique(tbl$id)
  m <- tbl |>
    select(id, precursor.mz) |>
    distinct() |>
    pull(precursor.mz)
  
  return(list(matrix = tbl, mz = m, ids = i))
  
}

#------------ Identification of peptides ------------
# Get name of spectra with highest correlation matrix
qcCorr <- function(s, ref.matrix, format) {
  # all data from a given spectrum
  if(format == '.mzML') {
    tmp <- data.frame(mz = mz(s),
                      intensity = intensity(s))
  } else if(format == '.raw'){
    tmp <- data.frame(mz = s$mZ,
                      intensity = s$intensity)
  }
  
  # create a matrix with real values
  comp <- data.frame(mz = ref.matrix$mz,
                     intensity = 0)
  
  for (j in 1:nrow(comp)) {
    mzi <- comp$mz[j]
    
    matches <- mzi > tmp$mz - 0.015 & mzi < tmp$mz + 0.015
    index <- which(matches == TRUE)
    
    if(length(index) == 1) comp$intensity[j] <- tmp$intensity[index]
    if(length(index) > 1) comp$intensity[j] <- sum(tmp$intensity[index], na.rm = TRUE)
  }
  
  # compute correlation coefficient
  if (sum(comp$intensity) == 0) {
    return(NA)
  } else {
    comp$intensity <- euclNorm(comp$intensity)
    ref.matrix$intensity <- euclNorm(ref.matrix$intensity)
    corr <- cosine(comp$intensity, ref.matrix$intensity)
    return(corr)
  }
}

#---------------- .mzML processing -----------------
# Process single m/z
processIdMzML <- function(x, i, ref, identify.bool) {
  mx <- ref$matrix
  mzs <- ref$mz
  ids <- ref$ids
  
  m <- mzs[ids == i]
  
  # Create XIC chromatogram
  xc <- chromatogram(x, mz = c(m-0.015, m+0.015))
  xc[,1]@intensity[is.na(xc[,1]@intensity)] <- 0
  xc[,1]@rtime <- xc[,1]@rtime/60
  
  if(identify.bool) {  
    # Filter out only data for relevant m/z
    px <- x |> filterPrecursorMz(m, ppm = (0.015/m)*1e6)
    mxx <- mx |>
      dplyr::filter(id == i) |>
      select(mz, intensity)
    
    # Get info about spectra
    spec <- spectra(px)
    ns <- names(spec)
    
    # Get spectrum with highest correlation
    if(length(spec) > 0) {
      res <- sapply(ns, function(nsx){
        s <- spec[[nsx]]
        qcCorr(s, mxx, format = '.mzML')
      })
      
      if(sum(!is.na(res)) != 0 & max(res, na.rm = TRUE) > 0) {
        poi <- res[which.max(res)]
      } else {
        poi <- NA
      }
      
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
  d <- data.frame(QC = i,
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
processQCfileMzML <- function(x, ref, identify.bool) {
  p <- pickPeaks(x)
  ids <- ref$ids
  
  l <- lapply(ids, function(id){
    processIdMzML(p, id, ref, identify.bool)
  })
  
  return(do.call(rbind, l))
}

# Wrapper function for calculation from mzML files
calculateQCmzML <- function(x, dir, output.dir, matrix, 
                            identify.bool = config$qc.identify, sample.name = 'sample') {
  
  ref <- loadQCReferenceMatrix(matrix)
  
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
    d <- processQCfileMzML(tmp, ref, identify.bool)
    
    saveDataCsv(xi, d, output.dir)
    d <- mutate(d, !!sample.name := xi)
    return(d)
  })
  
  return(do.call(rbind, l))
}

#----------------- .raw processing -----------------
# Process a single QC peptide
processIdRaw <- function(x, i, ref, identify.bool) {
  mx <- ref$matrix
  mzs <- ref$mz
  ids <- ref$ids
  
  m <- as.numeric(mzs[ids == i])
  xc <- readChromatogram(x, mass = m, tol = (0.015/m)*1e6, type = 'xic')
  
  # Get spectrum with highest correlation
  if(identify.bool){
    mxx <- mx |>
      dplyr::filter(id == i) |>
      select(mz, intensity)
    
    scan.df <- readIndex(x) |>
      dplyr::filter(MSOrder == 'Ms2') |>
      select(scan, precursorMass)
    
    scannum <- scan.df |>
      dplyr::filter(precursorMass >= m-0.015 & precursorMass <= m+0.015) |>
      pull(scan)
    
    if(length(scannum) > 0) {
      spec <- readSpectrum(x, scan = scannum)
      ns <- 1:length(spec)
      
      res <- sapply(ns, function(nsx){
        s <- spec[[nsx]]
        qcCorr(s, mxx, format = '.raw')
      })
      
      if(sum(!is.na(res)) != 0 & max(res, na.rm = TRUE) > 0) {
        poi <- res[which.max(res)]
        names(poi) <- which.max(res)
      } else {
        poi <- NA
      }
      
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
    smooth_int <- sgolayfilt(xc[[1]]$intensities, p = 5)
    smooth_peaks <- chromatographicPeakDetector(smooth_int)
    
    # at what RT is peak of interest
    poix <- spec[[as.numeric(names(poi))]]$StartTime
    poix_i <- min(which(xc[[1]]$times > poix))
    
    poi_row <- smooth_peaks[which.max(smooth_peaks[,1][smooth_peaks[,1] < poix_i]),]
    poi_start <- xc[[1]]$times[poi_row[1]]
    poi_end <- xc[[1]]$times[poi_row[2]]
    
    if (poix_i %in% smooth_peaks[,1]) {
      scnd_peak <- smooth_peaks[which(smooth_peaks[,1] == poix_i),]
      poi_end <- xc[[1]]$times[scnd_peak[2]]
    }
    
  } else {
    mcor <- NA
    poix <- NA
    poi_start <- NA
    poi_end <- NA
  }
  
  # Create output table
  d <- data.frame(QC = i,
                  rt = xc[[1]]$times,
                  intensity = xc[[1]]$intensities,
                  matrix_correlation = mcor,
                  peak_start = poi_start,
                  peak_end = poi_end,
                  peak_x = poix)
  
  locmin <- which(islocalminimum(d$intensity)==-1)
  d$baseline <- IPA_baselineDeveloper(locmin, d$intensity)
  rownames(d) <- NULL
  
  return(d)
  
}

# Wrapper function for calculation from raw files
calculateQCraw <- function(x, dir, matrix, identify.bool,
                           output.dir, sample.name = 'sample') {
  
  ref <- loadQCReferenceMatrix(matrix)
  fn <- paste0(dir, x, '.raw')
  
  l <- lapply(ref$ids, function(id){
    processIdRaw(fn, id, ref, identify.bool)
  })

  say(paste('.raw data read from', dir, 'for sample:'))
  say(x, type = 'input')
  
  d <- do.call(rbind, l)
  
  say(paste('QC peptide data calculated and saved as:'))
  saveDataCsv(x, d, output.dir)
  
  d <- mutate(d, !!sample.name := x)
  return(d)
}

# Final wrapper function
getQCs <- function(x, opts, suffix.csv = '.csv') {
  
  if(opts$qc.reprocess | opts$file.reanalyze) {
    flist <- list(load = NULL,
                  calc = x)
  } else {
    flist <- alreadyCalculatedFiles(x, dir = opts$qc.dir, suffix = suffix.csv)
  }
  
  df <- NULL
  
  if(length(flist$load) > 0) {
    df <- loadCalculatedTable(flist$load, dir = opts$qc.dir, suffix = suffix.csv)
    say(paste('QC data read from', opts$qc.dir))
  }
  
  if(length(flist$calc) > 0) {
    if(opts$file.format == '.mzML') {
      calc <- calculateQCmzML(flist$calc, dir = opts$file.input,
                              identify.bool = opts$qc.identify,
                              output.dir = opts$qc.dir,
                              matrix = opts$qc.reference.matrix)
    } else if(opts$file.format == '.raw'){
      calc <- lapply(flist$calc, function(x) {
        calculateQCraw(x, dir = opts$file.input,
                       matrix = opts$qc.reference.matrix, 
                       identify.bool = opts$qc.identify,
                       output.dir = opts$qc.dir)
      })
      calc <- do.call(rbind, calc)
    }
    
    df <- rbind(df, calc)
    
  }
  
  return(list(qc.df = df,
              new.qc = length(flist$calc)))
}

#-------------------- Plotting ---------------------
# Create variables for plot outputs
initiateQcPlotVariables <- function(d, by) {
  d <- length(unique(d[,by]))
  
  template <- data.frame(QC = numeric(),
                         matrix_correlation = numeric(),
                         idsl_auc = numeric(),
                         sym = numeric(),
                         SN_ratio = numeric(),
                         peak_y = numeric(),
                         width = numeric(),
                         sample = character())
  
  return(list(qc.subplot.number = d,
              qc.stats = template))
}

# Create a dataframe for creating an inset plot
createInputDfForInsetPlotting <- function(x, group) {
  df <- x |>
    mutate(rt = ifelse(rt < peak_start | rt > peak_end, NA, rt)) |>
    dplyr::filter(!is.na(rt)) |>
    group_by(!!as.symbol(group)) |>
    mutate(rt_at_max = ifelse(is.na(matrix_correlation), NA, rt[which.max(intensity)]),
           max_intensity = max(intensity, na.rm = TRUE),
           idsl_auc = ifelse(is.na(matrix_correlation), NA, peakAreaCalculator(rt, intensity)),
           auc_left = ifelse(is.na(matrix_correlation), NA, peakAreaCalculator(rt[rt <=rt_at_max], intensity)),
           auc_right = ifelse(is.na(matrix_correlation), NA, peakAreaCalculator(rt[rt >=rt_at_max], intensity)),
           sym = ifelse(is.na(matrix_correlation), NA, auc_left/auc_right),
           SN_ratio = ifelse(is.na(matrix_correlation), NA, SNRbaseline(intensity, baseline)),
           width = ifelse(is.na(matrix_correlation), NA, peak_end - peak_start)) |>
    as.data.frame()
  
  return(subset(df, !is.na(matrix_correlation)))
}

# Wrapper function for creating a separate QC plot
createSeparateQCplot <- function(x, df, i, mark.peak, config, by, batch, inset.bool = TRUE) {
  group.col <- ifelse(by == 'QC', 'sample', 'QC')
  
  # only get data for a single sample
  tmp <- dplyr::filter(df, !!as.symbol(by) == x) |>
    group_by(!!as.symbol(group.col)) |>
    mutate(rt = as.numeric(rt),
           intensity = as.numeric(intensity)) |>
    mutate(abs_intensity = max(intensity, na.rm = TRUE),
           peak_y = max(intensity, na.rm = TRUE),
           peak_start = ifelse(is.na(peak_start), min(rt, na.rm = TRUE), peak_start),
           peak_end = ifelse(is.na(peak_end), min(rt, na.rm = TRUE), peak_end)) |>
    ungroup() |>
    as.data.frame()
  
  # information for inset plot
  if(inset.bool) {
    inset.df <- createInputDfForInsetPlotting(tmp, group.col)
    inset.coords <- tibble(x = c(rep(max(tmp$rt, na.rm = TRUE), length(unique(inset.df[,group.col])))),
                           y = inset.df |>
                             select(!!as.symbol(group.col), peak_y) |>
                             distinct() |>
                             pull(peak_y),
                           plot = getInset(inset.df),
                           !!group.col := unique(inset.df[,group.col]),
                           !!by := unique(inset.df[,by]))
    
    inset.stats <- inset.df |>
      select(!!as.symbol(group.col), matrix_correlation, idsl_auc, sym, SN_ratio, peak_y, width) |>
      distinct() |>
      mutate(!!by := x)
    
    
    # final info about inset plot ready for plotting
    inset.stats.plot <- inset.stats |>
      mutate(label = ifelse(is.na(matrix_correlation), '',
                            paste0('SYM: ', signif(sym , 3), 
                                   '\nAUC: ', as.integer(idsl_auc),
                                   '\nSNR: ', signif(SN_ratio, 3),
                                   '\nCORR: ', signif(matrix_correlation, 3),
                                   '\nWDTH: ', signif(width, 3))),
             y_label = Inf)
  } else {
    inset.stats <- NULL
    inset.df <- tibble()
  }
  
  if (!nrow(inset.df) > 0) {
    inset.bool <- FALSE
    mark.peak <- FALSE
  }
  
  tmp[,group.col] <- factor(tmp[,group.col], levels = unique(tmp[,group.col]))
  tmp$sample <- plotAdjustNames(tmp$sample, config$plot.prefix, config$plot.suffix)
  
  # plotting
  p <- tmp |>
    ggplot(aes(x = rt, y = intensity, color = !!as.symbol(group.col))) +
    geom_line(na.rm = TRUE) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = 'RT [min]', y = 'Intensity')
  
  if(by == 'QC'){
    col <- rep(qcColors$color[qcColors$QC == x], times = length(unique(tmp[,group.col])))
    
    p <- p + 
      facet_wrap(sample~., scales = 'free_y', ncol = config$plot.max.columns) +
      theme(legend.position = 'none') +
      labs(subtitle = paste('peptide:', x)) +
      scale_color_manual(values = col)
  
    if(i == 1) p <- p + labs(title = 'QC peptides', subtitle = paste(batch, '| level: MS\npeptide:', x))
    #if(i == 2) p <- p + labs(title = ' ', subtitle = paste('\npeptide:', x))
      
  } else {
    colors <- tmp |> dplyr::select(QC, color) |>
      distinct() |>
      pull(color)
    
    p <- p + 
      facet_grid(QC~., scales = 'free', space = 'free_x') +
      labs(subtitle = 'QC peptides') +
      theme(legend.position = 'none') +
      scale_color_manual(values = colors)
  }
  
  if(inset.bool) {
    inset.coords$sample <- plotAdjustNames(inset.coords$sample, config$plot.prefix, config$plot.suffix)
    p <- p + 
      geom_plot_npc(data = inset.coords, 
                    aes(npcx = x, npcy = y, label = plot, vp.width = 0.25, vp.height = 0.9))
  }
  
  if(mark.peak) {
    inset.df$sample <- plotAdjustNames(inset.df$sample, config$plot.prefix, config$plot.suffix)
    inset.stats.plot$sample <- plotAdjustNames(inset.stats.plot$sample, config$plot.prefix, config$plot.suffix)
    
    p <- p +
      geom_point(data = inset.df, aes(x = rt_at_max, y = max_intensity + abs_intensity*0.1), color = 'black', shape = "\u2193", size = 3) +
      geom_text(data = inset.stats.plot, aes(x = max(tmp$rt)*0.75, y = y_label, 
                                             label = label), 
                size = 3, lineheight = 0.75, color = 'black', hjust = 1, vjust = 1.1)
  }
  
  p.info <- list(plot = list(p),
                 sample = x,
                 height = ((150 * length(unique(tmp[,group.col])))/config$plot.max.columns) + ifelse(i == 1 & by == 'QC', 250, 0),
                 stats = inset.stats)
  
  return(info = p.info)
}

createJoinedQCplot <- function(df, plots, current.batch, opts) {
  plots$row <- as.integer(rownames(plots))
  #plots$column <- ((plots$number+1) %% opts$plot.max.columns)+1
  #plots$row <- rep(1:ceiling(nrow(plots)/opts$plot.max.columns), each = opts$plot.max.columns)[1:nrow(plots)]
  plots$column <- rep(1, times = nrow(plots))
  
  height <- sum(plots$height) + 100
  
  p <- plots |>
    ggplot(aes(x = column, y = row)) +
    geom_plot_npc(aes(npcx = 0, npcy = 0, label = plot, vp.width = 1, vp.height = 1)) +
    facet_grid(rows = vars(row), cols = vars(column), scales = 'free') +
    theme_void() +
    theme(strip.text = element_blank(),
          plot.background = element_rect(fill = 'white', color = 'white'))
  
  tab <- tibble(name = 'QC peptides',
                link = 'qc_peptides',
                plot = list(p),
                filename = paste0(removeSpecialCharacters(current.batch), '_QC.png'),
                width = NA, 
                height = height, 
                type = 'QC peptides', 
                level = 'MS',
                full.path = NA,
                include.in.report = TRUE)
  return(tab)
}

# Create inset for iRT plots
getInset <- function(df) {
  out <- map(unique(df$QC), function(x) {
    tmp <- df |>
      dplyr::filter(QC == x)
    
    if(!is.na(tmp$matrix_correlation[1])){
      ggplot(data = tmp |> 
               group_by(QC),
             aes(x = rt, y = intensity)) +
        geom_line() +
        scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
        geom_vline(aes(xintercept = rt_at_max), color = 'red') +
        theme_bw() +  ## makes everything smaller
        expand_limits(y = 0) +
        theme(panel.background = element_rect(fill="white"),  ## white plot background 
              axis.title.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text = element_text(size=rel(0.7)), ## tiny axis text
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank())
    } else {
      NULL
    }
  })
  return(out)
}

# Final wrapper function
plotQc <- function(df, opts, by) {
  current.batch <- opts$file.batch
  
  plot_list <- list()
  plot_stats <- list()
  
  say('Generating QC peptide plots')
  
  for (x in unique(df[,by])) {
    res <- createSeparateQCplot(x, df, which(unique(df[,by]) == x), opts$qc.identify, opts, by, current.batch, opts$qc.identify)
    plot_list[[length(plot_list)+1]] <- tibble(plot = res$plot,
                                               sample = res$sample,
                                               height = res$height)
    plot_stats[[length(plot_stats)+1]] <- res$stats
  }
  
  say('QC peptide plots created')
  
  plot_list <- do.call(rbind, plot_list)
  
  if (length(plot_stats) > 0) {
    plot_stats <- do.call(rbind, plot_stats) |>
      relocate(sample, QC) |>
      mutate(idsl_auc = as.integer(idsl_auc),
             sym = signif(sym, 3),
             matrix_correlation = signif(matrix_correlation, 3),
             idsl_auc = signif(idsl_auc, 3),
             SN_ratio = signif(SN_ratio, 3),
             width = signif(width, 3)) |>
      na.omit()
    
    colnames(plot_stats) <- c('Sample', 'QC peptide', 'Reference similarity', 'AUC',
                              'Peak symmetry', 'SNR', 'Peak height', 'Peak width')
  } else {
    plot_stats <- NULL
  }
  
  plot_final <- createJoinedQCplot(df, plot_list, current.batch, opts)
  say('Joined QC peptide plot created')
  
  return(list(plot = plot_final,
              qc.stats = plot_stats))
  
}

# Function specifically for TIMS
plotQcTims <- function(df, config = opts, by) {
  current.batch <- config$file.batch
  
  if(config$qc.identify == TRUE){
    say('Warning!')
    say('Currently not possible to identify QC peptides, proceeding without identification', type = 'warn')
  }
  
  plot_list <- list()
  
  say('Generating QC peptide plots')
  
  for (x in unique(df[,by])){
    res <- createSeparateQCplot(x, df, which(unique(df[,by]) == x), FALSE, config, by, current.batch, inset.bool = FALSE)
    plot_list[[length(plot_list)+1]] <- tibble(plot = res$plot,
                                               sample = res$sample,
                                               height = res$height)
  }
  
  say('QC peptide plots created')
  
  plot_list <- do.call(rbind, plot_list)
  plot_final <- createJoinedQCplot(df, plot_list, current.batch, config)
  
  say('Joined QC peptide plot created')
  
  return(plot_final)
}
