# ***********************************************
# Title       : Calibration calculation
# Description : Creates a calibration curve
# Author      : Karolina Krystofova
# Date        : 2025/03/21
# ***********************************************
version <- '3.3.6'

# File paths
args <- commandArgs(trailingOnly=TRUE)

file_config <- args[1]
fcn_src <- args[2]
logfile <- args[3]

file_fcn <- paste0(fcn_src, '/R/')

# Time at the start of analysis
time_start <- Sys.time()

#------------------- ANALYSIS -------------------
tryCatch({
  source(paste0(file_fcn, 'console.R'))
  
  # Opens relevant libraries
  say(paste0('Started script QC_calibration_', version, '.R'))
  
  #---------------- LOAD LIBRARIES ---------------
  # Opens relevant libraries
  for(f in list.files(file_fcn, full.names = TRUE)){
    source(f)
  }
  
  for(l in libs){
    suppressWarnings(suppressMessages(library(l, character.only = TRUE, quietly = TRUE, 
                             warn.conflicts = FALSE, verbose = FALSE), 
                     classes = c('message', 'warning')))
  }
  
  say('Libraries loaded')
  
#----------------- CONFIGURATION ----------------
  # Read config file
  opts <- readCalibrationConfig(file_config)
  dirs <- createFolders(opts, ms.levels = 'calibration')
  opts <- modifyList(opts, dirs)
  
  opts$new.tic <- 0
  
  # Get the names of all the files which match the searched pattern
  fn <- sapply(opts$amounts, function(x) {
    list.files(opts$file.input)[grep(paste0('QC_cal', x, 'ng-', opts$method, '-', opts$date),
                                            list.files(opts$file.input), fixed = TRUE)]})
  
  samples <- gsub(opts$file.format, '', fn)
  
# -------------- CREATE EXPORT DIRS -------------
  # Output folder for TIC data
  if(dir.exists(paste0(opts$file.output, opts$tag,'_calibration'))) unlink(paste0(opts$file.output, opts$tag,'_calibration'), recursive = TRUE)
  opts <- c(opts, createExportDirectory(paste0(opts$tag,'_calibration'), 
                                               tag.bool = FALSE))

#------------------- READ TICs ------------------
  msn <- 1
  
  if(opts$file.format == '.mzML') {
    output <- getTICmzML(unname(samples), 'calibration')
  } else if(opts$file.format == '.raw'){
    output <- getTICraw(unname(samples), 'calibration')
  }
  
  df <- output$tic.df
  opts$new.tic <- output$new.tic
  msn <- ''
  
  # Adjust the output df
  df <- adjustDfCalibration(df)
  
  # Save as a new long df
  say('Full calibration data saved as:')
  saveDataCsv(paste0(opts$tag, '_full_calibration'), df, paste0(opts$tic.dir, 'calibration/'))
  
  # Are there at least 4 valid amounts?
  if (length(unique(df$sample)) < 4) {
    say('Error!')
    say(paste0('At least 4 calibration points needed (', length(unique(df$sample)),' provided)'), type = 'error')
    quit()
  }
  
#---------------- CALCULATE CURVE ---------------
  # Summarization
  df_summary <- summaryCalibration(df, opts)
  
  # Get curve
  curve.config <- getCalibrationCurve(df_summary, opts$formula)
  
#------------------- PLOT TIC -------------------
  # Create table where info about images is stored
  img_tbl <- img_tbl_template
 
  # Create variable for plot creation
  opts$plot.height <- (120*length(unique(df$sample))) + 100
  opts$plot.colors <- calibrationColors$colors[calibrationColors$amounts %in% unique(df$sample)]

  # Create plots
  say('Generating individual TIC plots')
  
  img_tbl[nrow(img_tbl)+1,] <- df |>
    mutate(intensity = ifelse(intensity > max(df_summary$max)*1.05, Inf, intensity)) |>
    plotChromatogramSingle('RT', 'intensity', cr.type = 'TIC', type = 'fixed', opts, 
                  start = opts$tic.start, end = opts$tic.end, colors = opts$plot.colors,
		  strip.position = 'right',
                  filename = paste0(opts$tag, '_calibration_fixedY.png'))
  
  img_tbl[nrow(img_tbl)+1,] <- merge(df, df_summary) |>
    mutate(max = max * 1.05,
           intensity = ifelse(intensity > max, Inf, intensity)) |>
    arrange(sample) |>
    plotChromatogramSingle('RT', 'intensity', cr.type = 'TIC', type = 'free', opts,
                  start = opts$tic.start, end = opts$tic.end, colors = opts$plot.colors,
		  strip.position = 'right',
                  filename = paste0(opts$tag, '_calibration_freeY.png'))
  
  say('Generating overlay TIC plot')
  img_tbl[nrow(img_tbl)+1,] <- plotChromatogramOverlay(df, 'RT', 'intensity', cr.type = 'TIC', 
                                                       start = opts$tic.start, end = opts$tic.end, 
                                                       chunk.no = length(unique(df$sample)), 
                                                       colors = opts$plot.colors, opts, 
                                                       calibration.bool = FALSE, batch = NULL, 
                                                       tag = opts$tag, 
                                                       filename = paste0(opts$tag, '_calibration_overlay.png'))
  say('TIC plots created')

#------------ PLOT CALIBRATION CURVE ------------
  say('Generating curve plot')
  img_tbl[nrow(img_tbl)+1,] <- plotCalibrationCurve(curve.config, opts$tag,
                                                    paste0(opts$tag, '_calibration_curve.png'))
  say('Calibration curve plot created')
  
#------------- CREATE FINAL TABLES --------------
  df.cal <- createCalibrationTable(opts, curve.config)
  new.df.cal <- newCalibrationTable(df.cal, opts$calibration.path)
  
  # Create plot for the comparison of calibration samples
  say('Generating calibration comparison plot')
  img_tbl[nrow(img_tbl)+1,] <- plotCalibrationSamplesComparison(new.df.cal, opts$method,
                                                                paste0(opts$tag, '_comparison.png'))
  say('Calibration comparison plot created')

  # Create output table
  names(fn) <- paste(as.integer(names(fn)), 'ng')
  fn <- fn[levels(df$sample)]
  
  # Create export tables
  export_tbls <- createExportTablesCalibration(df, fn, df_summary,
                                               curve.config, opts)
  
  # Save all created plots into the output folder
  say('Saving plots')
  img_tbl <- saveCreatedPlots(img_tbl)
  
  header <- paste('<b>Tag:</b>', opts$tag,
                  '<br><b>Method:</b>', opts$method,
		  '<br><b>Script version:</b>', version, '<br>')
  
  say('Generating HTML output')
  createHTMLreport(header, title = 'Calibration calculation', img.table = img_tbl,
                   export.tables = export_tbls, output.dir = opts$sample.dir, qc = FALSE, cal = TRUE,
                   type = 'calibration')

},

error = function(err) {
  e(err)
  quit()
})

diff <- Sys.time() - time_start

say(paste('Done!', opts$new.tic, 'new sample data processed in', round(as.numeric(diff), 2), 
          units(diff), '\n'))