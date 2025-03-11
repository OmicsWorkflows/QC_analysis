# Wrapper function for .mzML files
getBPCmzML <- function(x, suffix.csv = '.csv', 
                       config = opts) {
  
  output.dir <- config$bpc.dir
  counter <- config$new.bpc
  
  # Check if all files should be reanalyzed
  if(config$file.reanalyze){
    say('Initiating de novo BPC analysis of all samples')
    flist <- list(load = NULL,
                  calc = x)
  } else {
    say('Locating already calculated BPC data')
    flist <- alreadyCalculatedFiles(x, dir = output.dir, suffix = suffix.csv)
  }

  df <- NULL
  
  if(length(flist$load) > 0) {
    df <- loadCalculatedTable(flist$load, output.dir, suffix = suffix.csv)
    say(paste('BPC data loaded from', output.dir))
  }
  
  if(length(flist$calc) > 0) {
    calc <- calculateChromatogramMzML(x = flist$calc, ms.level = 1, dir = config$file.input,
                                      output.dir = output.dir, aggregation = 'max', config)
    df <- rbind(df, calc)
    
    counter <- counter + length(flist$calc)
  }
  
  return(list(bpc.df = df,
              new.bpc = counter))
}

# Wrapper function for .mzML files
getBPCraw <- function(x, suffix.csv = '.csv', 
                      config = opts) {
  
  output.dir <- config$bpc.dir
  counter <- config$new.bpc
  
  # Check if all files should be reanalyzed
  if(config$file.reanalyze){
    say('Initiating de novo BPC analysis of all samples')
    flist <- list(load = NULL,
                  calc = x)
  } else {
    say('Locating already calculated BPC data')
    flist <- alreadyCalculatedFiles(x, dir = output.dir, suffix = suffix.csv)
  }

  df <- NULL
  
  if(length(flist$load) > 0) {
    df <- loadCalculatedTable(flist$load, output.dir, suffix = suffix.csv)
    say(paste('BPC data loaded from', output.dir, 'for samples:'))
    for (s in flist$load) {
      say(s, type = 'input')
    }
  }
  
  if(length(flist$calc) > 0) {
    calc <- lapply(flist$calc, function(x) calculateChromatogramRaw(x, ms.level = 'ms',
                                                                    output.dir = output.dir,
                                                                    cr.type = 'bpc', config = config))
    
    calc <- do.call(rbind, calc)
    df <- rbind(df, calc)
    
    counter <- counter + length(flist$calc)
  }
  
  return(list(bpc.df = df,
              new.bpc = counter))
}