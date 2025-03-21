# HTML to character string
tableToHtml <- function(x) {
  output <- as.character(tableHTML(x, rownames = FALSE, border = 0))
  return(output)
}

# Create new folders
createFolders <- function(opts, ms.levels) {
  output <- list()
  
  o <- addSlash(opts$file.output)
  t <- paste0(o, "tic_data/")
  b <- paste0(o, "bpc_data/")
  
  if(!dir.exists(o)) dir.create(o)
  if(!dir.exists(t)) dir.create(t)
  if(!dir.exists(b)) dir.create(b)
  
  output$tic.dir <- t
  output$bpc.dir <- b
  
  if(opts$qc.bool){
    q <- paste0(o, "QC_data/")
    if(!dir.exists(q)) dir.create(q)
    output$qc.dir <- q
  }
  
  for (n in ms.levels){
      if(!dir.exists(paste0(t, n))) dir.create(paste0(t, n))
  }
  
  # Change the file.input folder to a tic folder if file.format is .csv
  if(opts$file.format == '.csv') {
    output$file.input <- addSlash(paste0(t, opts$MS.levels$name[1]))
  }
  
  return(output)
}

createExportTables <- function(d, ms.dict, removed = NULL){
 tmp <- d |>
	dplyr::select(sample,time) |>
	distinct()

 # Creates first table with sample information
  samp <- data.frame('No.' = 1:nrow(tmp),
                    Sample = tmp$sample,
                    Code = plotAdjustNames(tmp$sample, opts$plot.prefix, opts$plot.suffix),
                    'Time of acquisition' = tmp$time,
                    check.names = FALSE)
  
  if(length(removed) > 0) {
    samp.removed <- data.frame('No.' = (1:length(removed)) + length(unique(d$sample)),
                               Sample = removed,
                               Code = 'NA',
                               'Time of acquisition' = 'NA',
                               check.names = FALSE)
    samp <- rbind(samp, samp.removed)
  }
  
  if(opts$plot.prefix == '' & opts$plot.suffix == '') {
    samp <- dplyr::select(samp, -Code)
  }
  
  samp_txt <- tableToHtml(samp)
  
  # Creates final table with calculated metrics
  sums <- matrix(ncol = 4, nrow = 0) |>
    as.data.frame()
  colnames(sums) <- c('Sample', 'Metric', 'Level', 'Value')

  return(list(samples_tbl_exp = samp_txt,
              sums_tbl = sums))
}

createExportTablesCalibration <- function(d, fn, summ, curve, config){
  # Creates first table with sample information
  samp_table <- data.frame('No.' = 1:length(unique(d$sample)),
                           'Calibration sample' = levels(d$sample),
                           'Real amount' = paste(unique(d$real_amount), 'ng'),
                           Filename = fn,
                           'Time of acquisition' = unique(d$time),
                           check.names = FALSE)
  
  samp_txt <- tableToHtml(samp_table)
  
  # Creates final table with calculated metrics
  summ_table <- cbind(summ[,c('sample', 'real_amount')], curve$d$resid, summ$AUC)
  colnames(summ_table) <- c('Calibration sample', 'Real amount', 'Residuals', 'AUC')
  summ_txt <- tableToHtml(summ_table)
  
  curve_table <- data.frame(Equation = curve$equation,
                            'TIC start' = paste(config$tic.start, 'min'),
                            'TIC end' = paste(config$tic.end, 'min'),
                            Time = format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
                            check.names = FALSE)
  curve_txt <- tableToHtml(curve_table)
  
  return(list(samples_tbl_exp = samp_txt,
              sums_tbl_exp = paste0(summ_txt, '<br>', curve_txt)))
}

# Creates temporary export directory
createTmpDirectory <- function(batch.name = opts$file.batch, output = opts$file.output,
                               tag = opts$calibration.tag, tag.bool = opts$calibration.bool) {
  dir_name <- ifelse(tag.bool, paste(batch.name, tag, sep = '_'), batch.name)
  dir_name <- removeDoubleUnderscore(dir_name)
  dir_name <- gsub('?', 'x', dir_name, fixed = TRUE)
  dir_name <- gsub('x', 'x', dir_name, fixed = TRUE)
  
  if(substr(dir_name, nchar(dir_name), nchar(dir_name))=='_') {
    dir_name <- substr(dir_name, 1, (nchar(dir_name)-1))
  }
  
  path <- paste0(output, '.', dir_name, '/')
  
  if(dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)

  dir.create(path)
  
  dir_list <- list(sample.dir = path)
  
  # Other folders for storing outputs
  dn <- c('plots', 'tables')
  
  for (d in dn) {
    dir <- addSlash(paste0(path, d))
    dir.create(dir)
    dir_list[paste0(d, '.dir')] <- dir
  }
  
  return(dir_list)
}

# Replace tmp folder with a new folder
replaceTmpFolder <- function(tmp.folder = opts$sample.dir, batch.name = opts$file.batch, 
                             output = opts$file.output, tag = opts$calibration.tag, 
                             tag.bool = opts$calibration.bool, overwrite = opts$file.overwrite) {
  
  dir_name <- ifelse(tag.bool, paste(batch.name, tag, sep = '_'), batch.name)
  dir_name <- removeDoubleUnderscore(dir_name)
  dir_name <- gsub('?', 'x', dir_name, fixed = TRUE)
  dir_name <- gsub('x', 'x', dir_name, fixed = TRUE)
  
  if(substr(dir_name, nchar(dir_name), nchar(dir_name))=='_') {
    dir_name <- substr(dir_name, 1, (nchar(dir_name)-1))
  }
  
  if (dir.exists(paste0(output, dir_name)) & overwrite == FALSE) {
    found_dirs <- grep(dir_name, list.dirs(output, recursive = FALSE, full.names = FALSE), value = TRUE)
    pattern <- paste0("^", dir_name, "(?: \\((\\d+)\\))?$")
    
    found_dirs <- grep(pattern, found_dirs, value = TRUE)
    num <- suppressWarnings(as.numeric(sub(paste0(dir_name, " \\((\\d+)\\)"), "\\1", found_dirs)))
    
    if(all(is.na(num))) {
      ord <- 0
    } else {
      ord <- num[length(num)]
    }
    
    dir_name <- paste0(dir_name, ' (', ord+1, ')')
  } else if (dir.exists(paste0(output, dir_name)) & overwrite == TRUE) {
    unlink(paste0(output, dir_name, '/'), recursive = TRUE, force = TRUE)
  }
  
  full_name <- paste0(output, dir_name, '/')
  dir.create(full_name)
  
  for (f in list.files(tmp.folder, full.names = TRUE)) {
    file.copy(f, full_name, recursive = TRUE)
  }

  say('All output files moved to folder:')
  say(full_name, type = 'output')
  
  unlink(tmp.folder, recursive = TRUE, force = TRUE)
}

# Creates real export directory
createExportDirectory <- function(name, output = opts$file.output, tag = opts$calibration.tag, 
                                  tag.bool = opts$calibration.bool, overwrite = opts$file.overwrite) {
  # Create folders for output
  dir_name <- ifelse(tag.bool, paste(name, tag, sep = '_'), name)
  dir_name <- removeDoubleUnderscore(dir_name)
  
  if(!dir.exists(paste0(output, dir_name))) {
    output_dir <- paste0(output, dir_name)
  } else if (dir.exists(paste0(output, dir_name)) & overwrite == FALSE) {
    all_dirs <- grep(dir_name, list.dirs(output, recursive = FALSE, full.names = FALSE))
    dir_count <- sprintf('%02d', length(all_dirs))
    output_dir <- paste0(output, dir_name, '_', dir_count)
    output_dir <- removeDoubleUnderscore(output_dir)
  } else {
    output_dir <- paste0(output, dir_name)
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
  
  output_dir <- addSlash(output_dir)
  dir.create(output_dir)
  
  dir_list <- list(sample.dir = output_dir)
  
  # Other folders for storing outputs
  dn <- c('plots', 'tables')

  for (d in dn) {
    dir <- addSlash(paste0(output_dir, d))
    dir.create(dir)
    dir_list[paste0(d, '.dir')] <- dir
  }
  
  return(dir_list)
}

# Save all created plots
saveCreatedPlots <- function(plot.df, output = opts$plots.dir) {
  tmp <- subset(plot.df, is.na(full.path))
  
  for(i in 1:nrow(tmp)){
    x <- tmp[i,]
    x$full.path <- paste0(output, x$filename)
    #ggsave(x$full.path, x$plot, scale = 2, width = x$width, dpi = 150, height = x$height,
    #       units = 'px', limitsize = FALSE)
    png(x$full.path, width = x$width*2, height = x$height*2, units = 'px', res = 150)
    if(is.ggplot(x$plot[[1]])) {
      print(x$plot[[1]])
    } else {
      grid.draw(x$plot[[1]])
    }
    dev.off()
    
    #if(file.exists(paste0(current_dir, 'Rplots.pdf'))) {
      #unlink(paste0(current_dir, 'Rplots.pdf'))
    #}
  }
  
  say('All created plots saved to folder:')
  say(output, type = 'output')
  
  return(plot.df |>
    mutate(full.path = ifelse(is.na(full.path), paste0(output, filename), full.path)))
}

# Create .tsv table
saveExportTable <- function(d, cal.tag, output = opts$tables.dir) {
  export <- data.frame('Sample name' = d$Sample,
                       'Calibration tag' = cal.tag,
                       'Integrated AUC' = d$Integral,
                       'Amount (ng)' = d$Amount,
                       'Integral time' = format(Sys.time(), '%Y-%m-%d %H:%M:%S'),
                       'Comment' = d$Note,
                       check.names = FALSE) |>
    mutate(Comment = case_when(Comment == '*' ~ 'extrapolated amount value, integral lies outside calibration curve boundary',
                               Comment == '&dagger;' ~ 'no amount value, integral lies outside permitted calibration curve value range',
                               TRUE ~ ''))
  
  say('Table with calculated amounts saved as:')
  
  fn <- paste0(output, 'calculated_amounts.tsv')
  write.table(export, file = fn, sep = '\t', row.names = FALSE, quote = FALSE, col.names = colnames(export))
  
  say(fn, type = 'output')
}
