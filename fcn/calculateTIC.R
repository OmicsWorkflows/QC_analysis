getDateFromMzML <- function(x) {
  r <- runInfo(openMSfile(x))$startTimeStamp
  f <- format(as.POSIXct(r, format = '%Y-%m-%dT%H:%M:%S'),
         '%Y-%m-%d %H:%M:%S')
  return(f)
}

processTICfileMzML <- function(x) {
  d <- data.frame(RT = x@rtime/60,
                  intensity = x@intensity,
                  row.names = NULL)
  return(d)
}

# Calculate TICs from data that haven't been previously calculated
calculateTICmzML <- function(x, ms.level = 1, dir = input, output.dir = tic_dir, 
                             sample.name = 'sample') {
  # read .mzML files
  fn <- paste0(dir, x, '.mzML')
  readMS <- readMSData(fn, mode = "onDisk")
  
  # return success message
  say(paste0('.mzML files for TIC chromatogram computation read from ', dir, ' for samples:'))
  for (xi in x) {
    say(xi, type = 'input')
  }
  
  # calculate chromatogram
  t <- suppressWarnings(chromatogram(readMS, aggregationFun = "sum", msLevel = ms.level,
                                       BPPARAM = BiocParallel::SerialParam()))
  
  say('TIC data calculated and saved as:')
  
  # save as data frame
  l <- lapply(1:ncol(t), function(i) {
    tmp <- t[,i]
    date <- getDateFromMzML(paste0(dir, x[i], '.mzML'))
    d <- processTICfileMzML(tmp) |>
      mutate(time = date)
    
    saveDataCsv(x[i], d, output.dir)
    
    d <- mutate(d, !!sample.name := x[i])
    return(d)
  })
  
  return(do.call(rbind, l))
}

# Wrapper function for mzMLs
getTICmzML <- function(x, ms.level = msn, ms.dict = ms_lvls, output.dir = tic_dir, 
                       suffix.csv = '.csv', counter.var = 'new_tic'){
  ms_val <- ms.dict$value[ms.dict$name == ms.level]
  ms_dir <- addSlash(paste0(output.dir, ms.level))
  
  flist <- alreadyCalculatedFiles(x, dir = ms_dir, suffix = suffix.csv)
  df <- NULL
  
  if(length(flist$load) > 0) {
    df <- loadCalculatedTable(flist$load, ms_dir, suffix = suffix.csv)
    say(paste('TIC data read from', ms_dir, 'for samples:'))
    for (s in flist$load) {
      say(s, type = 'input')
    }
  }
  
  if(length(flist$calc) > 0) {
    calc <- calculateTICmzML(flist$calc, output.dir = ms_dir, ms.level = ms_val)
    df <- rbind(df, calc)
    
    counter <- get(counter.var)
    counter <- counter + length(flist$calc)
    assign(counter.var, counter, envir = parent.env())
  }
  
  return(df)
}

# calculate TICs from .raw files
calculateTICraw <- function(x, ms.level = ms, dir = input, output.dir = tic_dir, 
                            sample.name = 'sample') {
  fn <- paste0(dir, x, '.raw')
  cr <- readChromatogram(fn, filter = ms.level, type = 'tic')
  
  say(paste0('.raw for TIC chromatogram computation read from ', dir, ' for sample:'))
  say(fn, type = 'input')
  
  date_read <- readFileHeader(fn)$`Creation date`
  date <- format(as.POSIXct(date_read, format = '%m/%d/%Y %H:%M:%S'),
                      '%Y-%m-%d %H:%M:%S')
  
  d <- data.frame(RT = as.numeric(cr$times),
                  intensity = cr$intensities,
                  level = ms.level,
                  time = date,
                  row.names = NULL)
  
  say('TIC data calculated and saved as:')
  saveDataCsv(x, d, output.dir)
  d <- mutate(d, !!sample.name := x)

  return(d)
}

# Wrapper function for Thermo .raw
getTICraw <- function(x, ms.level = msn, ms.dict = ms_lvls, 
                      output.dir = tic_dir, suffix.csv = '.csv', counter.var = 'new_tic'){
  ms_val <- ms.dict$value[ms.dict$name == ms.level]
  ms_dir <- addSlash(paste0(output.dir, ms.level))
  flist <- alreadyCalculatedFiles(x, dir = ms_dir, suffix = suffix.csv)
  
  df <- NULL
  
  if(length(flist$load) > 0) {
    df <- loadCalculatedTable(flist$load, ms_dir, suffix = suffix.csv)
    say(paste('TIC data read from', ms_dir, 'for samples:'))
    for (s in flist$load) {
      say(s, type = 'input')
    }
  }
  
  if(length(flist$calc) > 0) {
    calc <- lapply(flist$calc, function(x) calculateTICraw(x, output.dir = ms_dir, ms.level = ms_val))
    calc <- do.call(rbind, calc)
    df <- rbind(df, calc)
    
    counter <- get(counter.var)
    counter <- counter + length(flist$calc)
    assign(counter.var, counter, envir = parent.env())
  }
  
  return(df)
  
}