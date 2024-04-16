createExportTables <- function(d, ms.dict){
  # Creates first table with sample information
  samp <- data.frame('No.' = 1:length(unique(d$sample)),
                    Sample = unique(d$sample),
                    'Time of acquisition' = unique(d$time),
                    check.names = FALSE)
  samp_exp <- tableHTML(samp, rownames = FALSE, border = 0)
  samp_txt <- as.character(samp_exp)
  
  # Creates final table with calculated metrics
  sums <- matrix(ncol = 2 + length(ms.dict$name), nrow = 0) |>
    as.data.frame()
  colnames(sums) <- c('Sample', 'Metric', ms.dict$name)

  return(list(samples_tbl_exp = samp_txt,
              sums_tbl = sums))
}

createExportDirectory <- function(name, tag = '', tag.bool = cal_bool) {
  # Create folders for output
  dir_name <- ifelse(tag.bool, paste(name, tag, sep = '_'), name)
  
  if(dir.exists(paste0(output, dir_name))) {
    all_dirs <- grep(dir_name, list.dirs(output, recursive = FALSE, full.names = FALSE))
    dir_count <- sprintf('%02d', length(all_dirs))
    output_dir <- paste0(output, dir_name, '_', dir_count)
    output_dir <- removeDoubleUnderscore(output_dir)
  } else {
    output_dir <- paste0(output, dir_name)
  }
  
  output_dir <- addSlash(output_dir)
  dir.create(output_dir)
  assign('output_dir', output_dir, envir = parent.frame())
  
  # Other folders for storing outputs
  fn <- c('plots', 'tables')

  for (f in fn) {
    dir <- addSlash(paste0(output_dir, f))
    dir.create(dir)
    assign(paste0(f, '_dir'), dir, envir = parent.frame())
  }
}