# ***********************************************
# Title       : Sample injection calculation
# Description : Calculates the injection volume
#               based on previously acquired
#               calibration curve
# Author      : Karolina Krystofova
# Date        : 2024/11/10
# Version     : 3.2.0
# ***********************************************

# File paths
args <- commandArgs(trailingOnly=TRUE)

logfile <- args[1]
file_config <- args[2]
fcn_src <- args[3]
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
  for(f in list.files(file_fcn, full.names = TRUE)) {
    source(f)
  }
  
  for(l in libs) {
    suppressWarnings(suppressMessages(library(l, character.only = TRUE, quietly = TRUE, 
                             warn.conflicts = FALSE, verbose = FALSE), 
                     classes = c('message', 'warning')))
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
    
    # Create table where info about images is stored
    img_tbl <- img_tbl_template
    errors <- list()
    
    # Variable to be changed during the first pass
    first_pass <- TRUE

#------------------ MS LEVEL LOOP ---------------
    for (msn in opts$MS.levels$name) {
      # Get pure sample names
      samples <- getSampleNames(b)
      
      say(paste0('Processing MS level ', msn, ' [',
                 which(opts$MS.levels == msn), '/', length(opts$MS.levels$name), ']'))
    
#----------------- CALCULATE TICS ---------------
      out <- getTICandQC(samples, msn, opts)
      
      samples.not.calc <- setdiff(samples, unique(out$tic.df$sample))
      samples <- samples[!samples %in% samples.not.calc]
      
      # Info for plot creation
      opts$plot.height <- (120*ceiling(length(samples)/opts$plot.max.columns)) + 100
      opts$plot.width <- ifelse(length(samples) < opts$plot.max.columns,
                                500*length(samples), 
                                500*opts$plot.max.columns)
      
      # Save info about which plots could not have been created
      errors[[length(errors)+1]] <- samples.not.calc
      names(errors)[length(errors)] <- msn
      
      if(length(samples) == 0) next
      
      opts$new.tic <- out$new.tic
      opts$new.bpc <- out$new.bpc

      if(opts$qc.bool & msn == 'MS') {
        opts$new.qc <- out$new.qc
        out$qc.plots$width <- opts$plot.width
        img_tbl <- rbind(img_tbl, out$qc.plots)
      }

      # Create a factor column out of the sample column
      # Order levels based on datetime
      out$tic.df <- arrange(out$tic.df, time, RT)
      out$tic.df$sample <- factor(out$tic.df$sample, levels = unique(out$tic.df$sample))
      
      if(msn == 'MS') {
        out$bpc.df <- arrange(out$bpc.df, time, RT)
        out$bpc.df$sample <- factor(out$bpc.df$sample, levels = unique(out$bpc.df$sample))
      }

#------------- GET CALIBRATION TAG --------------
      if(tolower(msn) == 'ms' & opts$calibration.bool) {
        say('Obtaining calibration curve information')
        opts <- modifyList(opts, getCalibrationInfo(samples, out$tic.df$time, opts))
        img_tbl[nrow(img_tbl)+1,] <- getCalibrationCurveInfo(opts)
        say(paste('Calibration curve', opts$calibration.tag, 'used for further processing'))
      }

#----------------- SAMPLE TABLE -----------------
      if(first_pass){
        export_tbls <- createExportTables(out$tic.df, opts$MS.levels, removed = samples.not.calc)
       
        if(opts$qc.bool & opts$qc.identify & length(nrow(out$qc.stats)) > 0) {
          export_tbls$qc_tbl <- as.character(tableHTML(out$qc.stats, 
                                                       rownames = FALSE, 
                                                       border = 0))
        }
      }
      
#--------------- INJECTION AMOUNT ---------------
      if(tolower(msn) == 'ms' & opts$calibration.bool) {
        say(paste('Calculating injection amount based on calibration curve', opts$calibration.tag))
        export_tbls <- modifyList(export_tbls, 
                                  calculateAUC(subset(out$tic.df, RT > opts$calibration.start & RT < opts$calibration.end), 
                                               remove.bg = TRUE, msn) |>
          select(Sample, Value) |>
          rename('AUC' = 'Value') |>
          calculateInjection(equation = opts$calibration.equation, max = opts$calibration.max))
        
        out$cal.df <- addCalibrationForPlotting(out$tic.df, export_tbls$inject_tbl, opts)
      }
        
#------------------ PLOT TICS -------------------
      say('Generating individual TIC plots')
      
      # Create variable for plot creation
      opts$plot.subtitle <- ifelse(opts$calibration.bool,
                                   paste(b, '| calibration:', opts$calibration.tag, '| level:', msn), 
                                   paste(b, '| levels:', msn))
      
      img_tbl[nrow(img_tbl)+1,] <- plotChromatogramSingle(out$tic.df, 'RT', 'intensity', cr.type = 'TIC', type = 'free', opts, 
                                                 start = opts$calibration.start, end = opts$calibration.end, 
                                                 filename = paste0(removeSpecialCharacters(b), '_', msn, '_tic_freeY.png'))
      img_tbl[nrow(img_tbl)+1,] <- plotChromatogramSingle(out$tic.df, 'RT', 'intensity', cr.type = 'TIC', type = 'fixed', opts,
                                                 start = opts$calibration.start, end = opts$calibration.end, 
                                                 filename = paste0(removeSpecialCharacters(b), '_', msn, '_tic_fixedY.png'))
      
      if(tolower(msn) == 'ms' & opts$calibration.bool){
        say('Generating TIC plots with calibration sample overlay')
        img_tbl[nrow(img_tbl)+1,] <- plotTicCalibration(out$cal.df, 'RT', 'intensity', opts, 
                                                   start = opts$calibration.start, end = opts$calibration.end, 
                                                   filename = paste0(removeSpecialCharacters(b), '_', msn, '_calibration_samples.png'))
      }
      
      say('Generating overlay TIC plots')
      img_tbl <- rbind(img_tbl, plotChromatogramOverlay(out$tic.df, 'RT', 'intensity', cr.type = 'TIC', start = opts$calibration.start,
                                                        end = opts$calibration.end, chunk.no = opts$plot.max.samples))
      say('TIC plots created')

#------------------- PLOT BPCS -------------------
      if(tolower(msn) == 'ms') {
        say('Generating individual BPC plots')
        img_tbl[nrow(img_tbl)+1,] <- plotChromatogramSingle(out$bpc.df, 'RT', 'intensity', cr.type = 'BPC', type = 'free', opts,
                                                            start = opts$calibration.start, end = opts$calibration.end, 
                                                            filename = paste0(removeSpecialCharacters(b), '_', msn, '_bpc_freeY.png'))
        img_tbl[nrow(img_tbl)+1,] <- plotChromatogramSingle(out$bpc.df, 'RT', 'intensity', cr.type = 'BPC', type = 'fixed', opts,
                                                            start = opts$calibration.start, end = opts$calibration.end, 
                                                            filename = paste0(removeSpecialCharacters(b), '_', msn, '_bpc_fixedY.png'))
        
        say('Generating overlay BPC plots')
        img_tbl <- rbind(img_tbl, plotChromatogramOverlay(out$tic.df, 'RT', 'intensity', cr.type = 'BPC', start = opts$calibration.start,
                                                          end = opts$calibration.end, chunk.no = opts$plot.max.samples))
        say('BPC plots created')
      }
      
#------------------ QC METRICS ------------------
      # AUC
      export_tbls$sums_tbl <- rbind(export_tbls$sums_tbl, calculateAUC(out$tic.df, remove.bg = FALSE, msn))
      
      say('Generating AUC plots')
      img_tbl[nrow(img_tbl)+1,] <- export_tbls$sums_tbl |>
        dplyr::filter(Metric == 'AUC', Level == msn) |>
        plotAUC('Sample', 'Value')
    
      # TIC fluctuations
      if(tolower(msn) == 'ms') {
        say('Generating TIC fluctuation plot')
        df_fluc <- ticFluc(out$tic.df, cutoff = opts$plot.fluctuation.threshold)
        export_tbls$sums_tbl <- rbind(export_tbls$sums_tbl, df_fluc$export)
        
        img_tbl[nrow(img_tbl)+1,] <- plotFluctuations(df_fluc$data, limits = opts$plot.fluctuation.threshold)

      # Correlation matrix
        say('Generating sample correlation plot')
        img_tbl[nrow(img_tbl)+1,] <- plotCorrelationMatrix(out$tic.df, 'RT', 'intensity')
      }
      
      say('QC metrics plots created')
      if(first_pass) first_pass <- FALSE
    }
    
#----------------- SAMPLE PLOTS -----------------
    sample_plots <- tibble::tibble(sample = character(),
                                   plot = list(),
                                   width = numeric(),
                                   height = numeric())
    samples <- getSampleNames(b)
    
    say('Generating individual sample plots for sample:')
    
    for (s in samples) {
      say(s, type = 'input')
      
      dfs <- getSampleData(s, opts)
      if(!is.null(dfs)) {
        tmp <- plotSampleData(s, dfs, x = 'RT', y = 'intensity', opts)
        sample_plots[nrow(sample_plots)+1,] <- tmp
        img_tbl[nrow(img_tbl)+1,] <- tibble(name = s, plot = tmp$plot, filename = paste0(s, '.png'),
                                            width = tmp$width, height = tmp$height, type = 'Individual samples',
                                            level = 'Individual samples', full.path = NA, include.in.report = FALSE)
      }
    }
    
    img_tbl <- rbind(img_tbl, createJoinedSamplePlots(sample_plots, opts))
    say('Plots for individual samples created')
    
#-------------- CREATE EXPORT DIRS --------------
    opts <- modifyList(opts, createExportDirectory(removeSpecialCharacters(b)))
    
#-------------- SAVE CREATED PLOTS --------------
    img_tbl$filename <- removeDoubleUnderscore(img_tbl$filename)
    say('Saving plots')
    img_tbl <- saveCreatedPlots(img_tbl)
    
#----------------- HTML OUTPUT ------------------
    # Rename colnames in sums
    export_tbls$sums_tbl <- export_tbls$sums_tbl |>
      spread(key = Level, value = Value)
    
    # Create HTML ouput
    cal_header <- ifelse(opts$calibration.bool,
                         paste0('<br><b>Calibration curve:</b> ', opts$calibration.tag, '\n'), '')
    say('Generating HTML output')
    createHTMLreport(header = paste0('<b>Sample batch:</b> ', b, '\n', cal_header,'<br>'), 
                     title = 'QC sample analysis', img.table = subset(img_tbl, include.in.report),
                     export.tables = export_tbls, output.dir = opts$sample.dir,
                     qc = (opts$qc.bool & opts$qc.identify), cal = opts$calibration.bool,
                     error.list = errors)
    
    # Remove files that have been created only as an intermediate step
    invisible(file.remove(list.files(opts$plots.dir, pattern = "^REMOVE_", full.names = TRUE)))

#------------------ TSV OUTPUT ------------------
    if(opts$calibration.bool) {
      say('Generating table output')
      saveExportTable(export_tbls$inject_tbl, cal.tag = opts$calibration.tag)
    }
  }
},

error = function(err) {
  e(err)
  quit()
})

diff <- Sys.time() - time_start

say(paste('Done!', opts$new.tic, 'new sample data and', opts$new.qc, 
          'new QC peptide data processed in', round(as.numeric(diff), 2), units(diff), '\n'))