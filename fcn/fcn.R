# ***********************************************
# Title       : Functions for QC calibration
# Description : plot_tic(), trapezoid_area()
# Author      : Karolina Krystofova
# Date        : 2024/02/22
# Version     : 1.0.5
# ***********************************************

# Add slash if not there
addSlash <- function(x) {
  if(substr(x, nchar(x), nchar(x)) != '/') x <- paste0(x, '/')
  return(x)
}

# Remove potential double underscore
removeDoubleUnderscore <- function(x) {
  gsub('__', '_', x, fixed = TRUE)
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
  r <- read.table(x, sep = '\n', 
                  strip.white = TRUE, fill = TRUE)
  cfg <- lapply(r$V1, twoColumns)
  cfg <- do.call(rbind, cfg)
  
  # assign individual variables
  for (v in c('tag', 'method')) {
    assign(v, cfg[cfg$variable == v, 2], envir = parent.frame())
  }
  
  for (v in c('input', 'output')) {
    assign(v, addSlash(cfg[cfg$variable == v, 2]), envir = parent.frame())
  }

  assign('qc_mx_path', cfg[cfg$variable == 'QC reference matrix', 2], envir = parent.frame())
  assign('plot_qc_bool', ifelse(cfg[cfg$variable == 'plot QC peptides (Y/N)',2] == 'Y', TRUE, FALSE), envir = parent.frame())
  assign('re_qc_bool', ifelse(cfg[cfg$variable == 'reprocess QC peptide data (Y/N)',2] == 'Y', TRUE, FALSE), envir = parent.frame())
  assign('cal_bool', ifelse(cfg[cfg$variable == 'use calibration (Y/N)',2] == 'Y', TRUE, FALSE), envir = parent.frame())
  
  lvls_str <- trimws(unlist(strsplit(cfg[cfg$variable == 'MS levels',2], ',')))
  lvls_l <- lapply(strsplit(lvls_str, '='), trimws)
  lvls_df <- do.call(rbind, lvls_l) |>
    as.data.frame() |>
    set_names(c('name', 'value'))
  
  assign('ms_lvls', lvls_df, envir = parent.frame())
  
  assign('file_format', cfg[cfg$variable == 'file format', 2], envir = parent.frame())
  assign('plot_colmax', as.numeric(cfg[cfg$variable == 'max columns',2]), envir = parent.frame())
  assign('overlay_maxsamps', as.numeric(cfg[cfg$variable == 'max samples in overlay plot',2]), envir = parent.frame())
  
  calOut <- cfg[cfg$variable == 'calibration output',2]
  calOut <- addSlash(calOut)
  assign('cal_path', calOut, envir = parent.frame())
  assign('cal_file', paste0(calOut, 'calibrations.tsv'), envir = parent.frame())
  
  assign('batches', cfg[(which(cfg[,1] == 'FILENAMES') + 1):nrow(cfg),1], envir = parent.frame())
}

# Create new folders
createFolders <- function(o, ms.levels = ms_lvls$name) {
  o <- addSlash(o)

  t <- paste0(o, "tic_data/")
  q <- paste0(o, "QC_data/")
  
  if(!dir.exists(o)) dir.create(o)
  if(!dir.exists(t)) dir.create(t)
  if(!dir.exists(q)) dir.create(q)
  
  for (n in ms.levels){
    if(!dir.exists(paste0(t, n))) dir.create(paste0(t, n))
  }
  
  assign('tic_dir', t, envir = parent.frame())
  assign('qc_dir', q, envir = parent.frame())
}
