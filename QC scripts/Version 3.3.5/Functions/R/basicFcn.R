# ***********************************************
# Title       : Functions for QC calibration
# Description : plot_tic(), trapezoid_area()
# Author      : Karolina Krystofova
# Date        : 2024/02/22
# ***********************************************

# Function for detecting keyboard hit
kbhit <- function() {
  ch <- .C("do_kbhit", result=as.integer(0))
  return(ch$result)
}

# Function for extracting batch name
getBatchName <- function(config) {
  r <- read.table(config, skip = 1, nrows = 1, sep = '\n')
  batch <- r[1,1]
  batch <- trimws(sub('batch name:', '', batch, fixed = TRUE))
  return(batch)
}

# Add slash if not there
addSlash <- function(x) {
  if(substr(x, nchar(x), nchar(x)) != '/') x <- paste0(x, '/')
  return(x)
}

# Clean path - replace backslashes
cleanPathString <- function(x) {
  x <- gsub('\\', '/', x, fixed = TRUE)
  x <- gsub('//', '/', x, fixed = TRUE)
  return(x)
}

# Remove potential double underscore
removeDoubleUnderscore <- function(x) {
  gsub('_+', '_', x, fixed = FALSE)
}

# Merge all data beyond second column
twoColumns <- function(x) {temp <- unname(unlist(strsplit(x, ':')))
                            temp_df <- data.frame(variable = trimws(temp[1]), 
                                                  value = ifelse(length(temp) > 1,
                                                                 trimws(paste(temp[2:length(temp)], collapse = ':')),
                                                                 ''),
                                                  row.names = NULL)
}

# Read input form
readConfig <- function(x) {
  r <- read.table(x, sep = '\n', blank.lines.skip = FALSE,
                  strip.white = TRUE, fill = TRUE)
  cfg <- lapply(r$V1, twoColumns)
  cfg <- do.call(rbind, cfg)
  cfg <- subset(cfg, !is.na(variable))
  
  # configuration of calibration
  opts <- list()
  
  opts$config <- r
  opts$calibration.path <- cleanPathString(addSlash(cfg[cfg$variable == 'calibration output',2]))
  opts$calibration.file <- paste0(opts$calibration.path, 'calibrations.tsv')
  opts$calibration.bool <- ifelse(cfg[cfg$variable == 'use calibration (Y/N)',2] == 'Y', TRUE, FALSE)
  opts$calibration.tag.read <- cfg[cfg$variable == 'tag', 2]
  opts$calibration.method <- cfg[cfg$variable == 'method', 2]
  opts$calibration.start <- 0
  opts$calibration.end <- 0
  
  # configuration of file processing
  opts$file.input <- cleanPathString(addSlash(cfg[cfg$variable == 'input', 2]))
  opts$file.output <- cleanPathString(addSlash(cfg[cfg$variable == 'output', 2]))
  lvls_str <- trimws(unlist(strsplit(cfg[cfg$variable == 'MS levels',2], ',')))
  lvls_l <- lapply(strsplit(lvls_str, '='), trimws)
  lvls_df <- do.call(rbind, lvls_l) |>
    as.data.frame() |>
    set_names(c('name', 'value'))
  opts$MS.levels <- lvls_df
  opts$file.batch <- cfg[which(cfg[,1] == 'batch name'),2]
  opts$file.format <- cfg[cfg$variable == 'file format', 2]
  opts$file.offset <- as.numeric(cfg[cfg$variable == 'time offset (min)', 2])
  opts$file.date <- cfg[cfg$variable == 'date format', 2]
  opts$file.reanalyze <- ifelse(cfg[cfg$variable == 'reprocess raw files (Y/N)', 2] == 'Y', TRUE, FALSE)
  opts$file.overwrite <- ifelse(cfg[cfg$variable == 'overwrite existing outputs (Y/N)', 2] == 'Y', TRUE, FALSE)
  opts$file.wait.time <- ifelse(cfg[cfg$variable == 'wait time between loops (min)', 2] %in% c('Inf', 'NA', 'N', 'inf', 'na', 'Na'), 
                                NA, as.numeric(cfg[cfg$variable == 'wait time between loops (min)', 2]))
  
  # configuration for QC peptide calculation
  opts$qc.reference.matrix <- cleanPathString(cfg[cfg$variable == 'QC reference matrix', 2])
  opts$qc.bool <- ifelse(cfg[cfg$variable == 'analyze QC peptides (Y/N)', 2] == 'Y', TRUE, FALSE)
  opts$qc.reprocess <- ifelse(cfg[cfg$variable == 'reprocess QC peptide data (Y/N)', 2] == 'Y', TRUE, FALSE)
  opts$qc.identify <- ifelse(cfg[cfg$variable == 'identify QC peptides (Y/N)', 2] == 'Y', TRUE, FALSE)

  # configuration of plots
  opts$plot.max.columns <- as.numeric(cfg[cfg$variable == 'max columns', 2])
  opts$plot.max.samples <- as.numeric(cfg[cfg$variable == 'max samples in overlay plot', 2])
  opts$plot.bar.max.samples <- as.numeric(cfg[cfg$variable == 'max samples in bar plots', 2])
  opts$plot.prefix <- cfg[cfg$variable == 'prefix to remove', 2]
  opts$plot.suffix <- cfg[cfg$variable == 'suffix to remove', 2]
  opts$plot.fluctuation.threshold <- as.numeric(cfg[cfg$variable == 'TIC fluctuation threshold', 2])
  
  # save info about newly processed files
  opts$new.qc <- 0
  opts$new.tic <- 0
  opts$new.bpc <- 0
  
  return(opts)
}

readCalibrationConfig <- function(x) {
  r <- read.table(x, sep = '\n', blank.lines.skip = FALSE,
                  strip.white = TRUE, fill = TRUE)
  cfg <- lapply(r$V1, twoColumns)
  cfg <- do.call(rbind, cfg)
  cfg <- subset(cfg, !is.na(variable))
  
  opts <- list()
  opts$config <- r
  for (v in c('input', 'output')) opts[paste0('file.',v)] <- addSlash(cfg[cfg$variable == v, 2])
  for (v in c('date', 'method', 'tag', 'formula')) opts[v] <- cfg[cfg$variable == v, 2]
  
  opts$calibration.path <- cfg[cfg$variable == 'calibration table', 2]
  opts$file.format <- cfg[cfg$variable == 'file format', 2]
  opts$file.date <- cfg[cfg$variable == 'date format', 2]
  opts$tic.start <- as.numeric(cfg[cfg$variable == 'TIC trace start', 2])
  opts$tic.end <- as.numeric(cfg[cfg$variable == 'TIC trace end', 2])
  opts$amounts <- trimws(unlist(strsplit(cfg[cfg$variable == 'include amounts', 2], ',')))
  opts$file.reanalyze <- ifelse(cfg[cfg$variable == 'reprocess raw files (Y/N)', 2] == 'Y', TRUE, FALSE)
  
  opts$qc.bool <- FALSE
  
  opts$real.amounts <- cfg[(which(cfg[,1] == 'REAL AMOUNTS') + 1):nrow(cfg),]
  rownames(opts$real.amounts) <- NULL
  colnames(opts$real.amounts) <- c('file', 'amount')
  
  opts$plot.width <- 500
  opts$plot.max.columns <- 1
  opts$plot.subtitle <- opts$tag
  opts$plot.prefix <- cfg[cfg$variable == 'prefix to remove', 2]
  opts$plot.suffix <- cfg[cfg$variable == 'suffix to remove', 2]
   
  return(opts)
}

# Get sample names
getSampleNames <- function(batch, config = opts) {
  all.files <- list.files(config$file.input, full.names = TRUE)
  
  # Only files of the appropriate format
  fn <- unlist(sapply(all.files, function(x) {
    s <- substr(x, nchar(x)-nchar(config$file.format)+1, nchar(x))
    if(s == config$file.format) return(x)
  }, USE.NAMES = FALSE))
  
  # Get the names of all the files which match the searched pattern
  if(grepl('?', batch, fixed = TRUE) | grepl('*', batch, fixed = TRUE)){
    lookup <- gsub('?', '.', batch, fixed = TRUE)
    lookup <- gsub('*', '.*', lookup, fixed = TRUE)
    lookup <- paste0('^', lookup)
    fn <- fn[grepl(lookup, sub(config$file.input, '', fn))]
  } else {
    fn <- fn[grepl(batch, fn, fixed = TRUE)]
  }
  
  if(config$file.format == '.d') { 
    fn.times <- sapply(fn, function(x) {
      analysis.file <- paste(x, 'analysis.tdf', sep = '/')
      file.info(analysis.file)$mtime
    })
  } else {
    fn.times <- file.info(fn)$mtime
  }
  
  if(length(fn) > 0){
    fn.times.diff <- as.numeric(difftime(Sys.time(), fn.times, units = 'mins'))
    
    fn <- fn[fn.times.diff > config$file.offset]
  }
  
  # Throw an error when no data with given pattern are found
  if (length(fn) == 0 & is.na(config$file.wait.time)) {
    say('Error!')
    say(paste0('No files containing the pattern "', batch, '" in selected input folder'), type = 'error')
    quit()
  }
  
  # Remove address and suffix, get pure sample names
  samples <- sapply(fn, function(x) sub(config$file.input, '', x, fixed = TRUE), USE.NAMES = FALSE)
  samples <- sapply(samples, function(x) substr(x, 1, nchar(x)-nchar(opts$file.format)), USE.NAMES = FALSE)
  
  return(samples)
}

# Remove "*" and "?" from characters
removeSpecialCharacters <- function(char){
  char <- gsub('?', 'x', char, fixed = TRUE)
  char <- gsub('*', 'x', char, fixed = TRUE)
  return(char)
}

# Get reverse string
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = '')

# Remove prefix and suffix from sample names
plotAdjustNames <- function(x, prefix, suffix) {
  x <- as.character(x)
  x.unique <- unique(x)
  x.unique.original <- x.unique
  
  prefix.bool <- grepl('*', prefix, fixed = TRUE) |  grepl('?', prefix, fixed = TRUE)
  suffix.bool <- grepl('*', suffix, fixed = TRUE) |  grepl('?', suffix, fixed = TRUE)
  
  if(prefix != '') x.unique <- sapply(x.unique, function(xi) {
    # Use regex?
    if(prefix.bool) {
      
      xi.rev <- strReverse(xi)
      
      pre.rev <- strReverse(prefix)
      pre.rev <- gsub('?', '.', pre.rev, fixed = TRUE)
      pre.rev <- gsub('*', '.*', pre.rev, fixed = TRUE)
      pre.rev <- paste0('(.*)', pre.rev)
      
      xi.rev <- sub(pre.rev, '\\1', xi.rev)
      xi <- strReverse(xi.rev)
      return(xi)
      
      # If not, just remove the fixed string  
    } else {
      
      if(grepl(prefix, xi, fixed = TRUE)) {
        return(substr(xi, nchar(prefix)+1, nchar(xi)))
      } else {
        return(xi)
      }
      
    }
  }, USE.NAMES = FALSE)
  
  if(suffix != '') x.unique <- sapply(x.unique, function(xi) {
    # Use regex?
    if(suffix.bool) {
      suf <- gsub('?', '.', suffix, fixed = TRUE)
      suf <- gsub('*', '.*', suf, fixed = TRUE)
      suf <- paste0('(.*)', suf)
      
      xi <- sub(suf, '\\1', xi)
      return(xi)
      
      # If not, just remove the fixed string  
    } else {
      if(grepl(suffix, xi, fixed = TRUE)) {
        return(substr(xi, 1, nchar(xi)-nchar(suffix)))
      } else {
        return(xi)
      }
    }
  }, USE.NAMES = FALSE)
  
  # Check if there are any duplicate sample names
  multiple.i <- which(table(x.unique) > 1)
  
  if(length(multiple.i) > 0){
    multiple.names <- names(multiple.i)
    
    for(mn in multiple.names){
      mn.i <- which(x.unique == mn)
      renum <- 1:length(mn.i)
      x.unique[mn.i] <- paste0(mn, ' (', renum, ')')
    }
  }
  
  replace.df <- data.frame(original = x.unique.original,
                           replace = x.unique)
  
  for(i in 1:nrow(replace.df)) {
    x[which(x == replace.df$original[i])] <- replace.df$replace[i]
  }
  
  x <- factor(x, levels = x.unique)
  
  return(x)
}

# Concatenate logs
concatLogs <- function(output.dir, logfile) {
  r <- readLines(logfile, warn = FALSE)
  r.collapse <- paste0(paste(r, collapse = '\n'), '\n')
  
  log <- paste0(output.dir, 'log.txt')

  if(file.exists(log)) r.collapse <- paste0('\n', r.collapse)
  
  cat(r.collapse, file = log, append = TRUE)
  unlink(logfile, force = TRUE)
}
