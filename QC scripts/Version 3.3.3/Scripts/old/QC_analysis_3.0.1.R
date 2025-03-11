# ***********************************************
# Title       : Sample injection calculation
# Description : Calculates the injection volume
#               based on previously acquired
#               calibration curve
# Author      : Karolina Krystofova
# Date        : 2024/08/12
# Version     : 3.0.1
# ***********************************************

# File paths
args <- commandArgs(trailingOnly=TRUE)

file_config <- args[1]
fcn_src <- args[2]
logfile <- args[3]
python_path <- args[4]

file_fcn <- paste0(fcn_src, 'R/')

# Time at the start of analysis
time_start <- Sys.time()

#------------------- ANALYSIS -------------------
tryCatch({
  source(paste0(file_fcn, 'console.R'))
  
  # Opens relevant libraries
  say('Started script sample_calculation.R')

#---------------- LOAD LIBRARIES ---------------
  # Opens relevant libraries
  for(f in list.files(file_fcn, full.names = TRUE)){
    source(f)
  }
  
  for(l in libs){
    suppressMessages(library(l, character.only = TRUE, quietly = TRUE, 
                             warn.conflicts = FALSE, verbose = FALSE), 
                     classes = c('message', 'warning'))
  }
  
  say('Libraries loaded')
  
#----------------- CONFIGURATION ----------------
  # Read config file and initiate startup variables
  opts <- readConfig(file_config)
  opts$src <- fcn_src
  opts$python <- python_path

  # Output folder for TIC data
  dirs <- createFolders(opts$file.output,
                        qc.bool = opts$qc.bool,
                        ms.levels = opts$MS.levels$name)
  opts <- modifyList(opts, dirs)

#------------------ LOOP START ------------------
  # Loop over batches of files
  for(b in opts$file.batches){
    say(paste0('Processing sample batch ', b, ' [',
                which(opts$file.batches == b), '/', length(opts$file.batches), ']'))
    
#------------------ RAW FILE LOOKUP ----------------
    # Get the names of all the files which match the searched pattern
    fn <- list.files(opts$file.input)[grep(b, list.files(opts$file.input), fixed = TRUE)]
    fn <- fn[grep(opts$file.format, fn, fixed = TRUE)]
    
    # Throw an error when no data with given pattern are found
    if (length(fn) == 0) {
      say('Warning!')
      say(paste0('No files containing the pattern "', b, '" in selected input folder'), type = 'warn')
      quit()
    }
    
    # Remove suffix, get pure sample names
    samples <- sapply(fn, function(x) substr(x, 1, nchar(x)-nchar(opts$file.format)), USE.NAMES = FALSE)

#--------------- CREATE EXPORT DIRS -----------------
    opts <- modifyList(opts, createExportDirectory(b))
    
    # Create table where info about images is stored
    img_tbl <- img_tbl_template
    
    # QC peptide identification and plot creation
    if (opts$qc.bool) {
      out <- getQCs(x = samples, opts = opts)
      qc_df <- out$df
      opts$new.qc <- out$new.qc
      
      opts <- modifyList(opts, initiateQcPlotVariables(qc_df))
      qc_output <- plotQc(qc_df, opts)
      
      img_tbl[nrow(img_tbl)+1,] <- qc_output$plot
    }
    
#-------------------- MS LEVEL LOOP -----------------
    samp_tbl_bool <- TRUE
    
    for (msn in opts$MS.levels$name) {
      say(paste0('Processing MS level ', msn, ' [',
                 which(opts$MS.levels == msn), '/', length(opts$MS.levels), ']'))
    
#------------------- CALCULATE TICS -----------------
      if(opts$file.format == '.mzML') {
        out <- getTICmzML(samples, type = 'MS', ms.level = msn)
      } else if(opts$file.format == '.raw') {
        out <- getTICraw(samples, type = 'MS', ms.level = msn)
      } else if(opts$file.format == '.d') {
        out <- getTICtimsD(samples, type = 'MS', ms.level = msn)
      }

      df <- out$df
      opts$new.tic <- out$new.tic
      
      # Create a factor column out of the sample column
      # Order levels based on datetime
      df <- arrange(df, time, RT)
      df$sample <- factor(df$sample, levels = unique(df$sample))

#------------- GET CALIBRATION TAG --------------
      if(tolower(msn) == 'ms' & opts$calibration.bool) {
        opts <- modifyList(opts, getCalibrationInfo(samples, df$time, opts))
        img_tbl[nrow(img_tbl)+1,] <- getCalibrationCurveInfo(opts)
      }
      
#----------------- SAMPLE TABLE -----------------
      if(samp_tbl_bool){
        export_tbls <- createExportTables(df, opts$MS.levels)
        samp_tbl_bool <- FALSE
        
        if(opts$qc.bool){
          export_tbls$qc_tbl <- as.character(tableHTML(qc_output$qc.stats, 
                                                       rownames = FALSE, 
                                                       border = 0))
        }
      }
        
#------------- PLOT SEPARATE TICS ---------------
      say('TIC plots saved as:')

      # Create variable for plot creation
      opts$plot.subtitle <- ifelse(opts$calibration.bool,
                                   paste(b, '| calibration:', opts$calibration.tag, '| level:', msn), 
                                   paste(b, '| levels:', msn))
      opts$plot.height <- (120*ceiling(length(unique(df$sample))/opts$plot.max.columns)) + 100
      opts$plot.width <- ifelse(length(samples) < opts$plot.max.columns,
                                   500*length(samples), 
                                   500*opts$plot.max.columns)
      
      img_tbl[nrow(img_tbl)+1,] <- plotTicSingle(df, 'RT', 'intensity', type = 'free', opts, start = opts$calibration.start,
                                                 end = opts$calibration.end, filename = paste0(opts$plots.dir, b, '_', msn, '_freeY.png'))
      img_tbl[nrow(img_tbl)+1,] <- plotTicSingle(df, 'RT', 'intensity', type = 'fixed', opts, start = opts$calibration.start,
                                                 end = opts$calibration.end, filename = paste0(opts$plots.dir, b, '_', msn, '_fixedY.png'))
      
#-------------- PLOT TICS OVERLAY ---------------
      say('Overlay TIC plots saved as:')
      img_tbl <- rbind(img_tbl, plotTicOverlay(df, 'RT', 'intensity', start = opts$calibration.start,
                                               end = opts$calibration.end, chunk.no = opts$plot.max.samples))
    
#------------------ QC METRICS ------------------
      # AUC
      export_tbls$sums_tbl <- rbind(export_tbls$sums_tbl, calculateAUC(df, remove.bg = FALSE, msn))
      say('AUC calculated and plot saved as:')
  
      img_tbl[nrow(img_tbl)+1,] <- export_tbls$sums_tbl |>
        dplyr::filter(Metric == 'AUC', Level == msn) |>
        plotAUC('Sample', 'Value')
    
      # TIC fluctuations
      if(tolower(msn) == 'ms') {
        data_fluc <- ticFluc(df, cutoff = 10)
        export_tbls$sums_tbl <- rbind(export_tbls$sums_tbl, data_fluc$export)
        
        say('TIC fluctuations calculated and plot saved as:')
        img_tbl[nrow(img_tbl)+1,] <- plotFluctuations(data_fluc$data, limits = 10)
      
      # Correlation matrix
        say('Sample correlation calculated and plot saved as:')
        img_tbl[nrow(img_tbl)+1,] <- plotCorrelationMatrix(df, 'RT', 'intensity')
      }
      
      
#--------------- INJECTION AMOUNT ---------------
      if(tolower(msn) == 'ms' & opts$calibration.bool) {
        export_tbls$inject_tbl_exp <- calculateAUC(subset(df, RT > opts$calibration.start & RT < opts$calibration.end), 
                                                   remove.bg = TRUE, msn) |>
          select(Sample, Value) |>
          rename('AUC' = 'Value') |>
          calculateInjection(coefs = opts$calibration.curve, max = opts$calibration.max)
      }
    }
    
#----------------- HTML OUTPUT ------------------
    # Rename colnames in sums
    export_tbls$sums_tbl <- export_tbls$sums_tbl |>
      spread(key = Level, value = Value)
    
    # Create HTML ouput
    createHTMLreport(header = paste('<b>Sample batch:</b>', b, '\n<br>'), 
                     title = 'QC sample analysis', img.table = img_tbl,
                     export.tables = export_tbls, output.dir = opts$sample.dir,
                     qc = opts$qc.bool, cal = opts$calibration.bool)
  }
},

error = function(err) {
  e(err)
  quit()
})

diff <- Sys.time() - time_start

say(paste('Done!', new_tic, 'new sample data and', new_qc, 
          'new QC peptide data processed in', round(as.numeric(diff), 2), units(diff), '\n'))