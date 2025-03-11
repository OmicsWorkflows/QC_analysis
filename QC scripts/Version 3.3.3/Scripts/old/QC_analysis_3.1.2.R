# ***********************************************
# Title       : Sample injection calculation
# Description : Calculates the injection volume
#               based on previously acquired
#               calibration curve
# Author      : Karolina Krystofova
# Date        : 2024/10/31
# Version     : 3.1.2
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
    
#--------------- RAW FILE LOOKUP ----------------
    # Get pure sample names
    samples <- getSampleNames(b)

    # Info for plot creation
    opts$plot.height <- (120*ceiling(length(samples)/opts$plot.max.columns)) + 100
    opts$plot.width <- ifelse(length(samples) < opts$plot.max.columns,
                              500*length(samples), 
                              500*opts$plot.max.columns)
    
#------------------ MS LEVEL LOOP ---------------
    # Create table where info about images is stored
    img_tbl <- img_tbl_template
    
    # Variable to be changed during the first pass
    first_pass <- TRUE
    
    for (msn in opts$MS.levels$name) {
      say(paste0('Processing MS level ', msn, ' [',
                 which(opts$MS.levels == msn), '/', length(opts$MS.levels$name), ']'))
    
#----------------- CALCULATE TICS ---------------
      out <- getTICandQC(samples, msn, opts)
      
      samples.not.calc <- setdiff(samples, unique(out$tic.df$sample))
      samples <- samples[!samples %in% samples.not.calc]
      
      opts$plot.height <- (120*ceiling(length(samples)/opts$plot.max.columns)) + 100
      opts$plot.width <- ifelse(length(samples) < opts$plot.max.columns,
                                500*length(samples), 
                                500*opts$plot.max.columns)
      
      opts$new.tic <- out$new.tic
      opts$new.bpc <- out$new.bpc

      if(opts$qc.bool & msn == 'MS') {
        opts$new.qc <- out$new.qc
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
        opts <- modifyList(opts, getCalibrationInfo(samples, out$tic.df$time, opts))
        img_tbl[nrow(img_tbl)+1,] <- getCalibrationCurveInfo(opts)
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
        export_tbls <- modifyList(export_tbls, 
                                  calculateAUC(subset(out$tic.df, RT > opts$calibration.start & RT < opts$calibration.end), 
                                               remove.bg = TRUE, msn) |>
          select(Sample, Value) |>
          rename('AUC' = 'Value') |>
          calculateInjection(coefs = opts$calibration.curve, max = opts$calibration.max))
        
        say(paste('Calibration curve', opts$calibration.tag, 'used for further processing'))
        
        out$cal.df <- addCalibrationForPlotting(out$tic.df, export_tbls$inject_tbl, opts)
      }
        
#------------- PLOT SEPARATE TICS ---------------
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
        img_tbl[nrow(img_tbl)+1,] <- plotTicCalibration(out$cal.df, 'RT', 'intensity', opts, 
                                                   start = opts$calibration.start, end = opts$calibration.end, 
                                                   filename = paste0(removeSpecialCharacters(b), '_', msn, '_calibration_samples.png'))
      }
      
      say('Individual TIC plots created')
#------------- PLOT SEPARATE BPCS ---------------
      if(tolower(msn) == 'ms') {
        img_tbl[nrow(img_tbl)+1,] <- plotChromatogramSingle(out$bpc.df, 'RT', 'intensity', cr.type = 'BPC', type = 'free', opts,
                                                            start = opts$calibration.start, end = opts$calibration.end, 
                                                            filename = paste0(removeSpecialCharacters(b), '_', msn, '_bpc_freeY.png'))
        img_tbl[nrow(img_tbl)+1,] <- plotChromatogramSingle(out$bpc.df, 'RT', 'intensity', cr.type = 'BPC', type = 'fixed', opts,
                                                            start = opts$calibration.start, end = opts$calibration.end, 
                                                            filename = paste0(removeSpecialCharacters(b), '_', msn, '_bpc_fixedY.png'))
        say('Individual BPC plots created')
      }
      
#---------------- OVERLAY PLOTS -----------------
      img_tbl <- rbind(img_tbl, plotChromatogramOverlay(out$tic.df, 'RT', 'intensity', cr.type = 'TIC', start = opts$calibration.start,
                                                        end = opts$calibration.end, chunk.no = opts$plot.max.samples))
      say('Overlay TIC plots created')
      
      if(tolower(msn) == 'ms') {
        img_tbl <- rbind(img_tbl, plotChromatogramOverlay(out$tic.df, 'RT', 'intensity', cr.type = 'BPC', start = opts$calibration.start,
                                                          end = opts$calibration.end, chunk.no = opts$plot.max.samples))
        say('Overlay BPC plots created')
      }
      
#------------------ QC METRICS ------------------
      # AUC
      export_tbls$sums_tbl <- rbind(export_tbls$sums_tbl, calculateAUC(out$tic.df, remove.bg = FALSE, msn))

      img_tbl[nrow(img_tbl)+1,] <- export_tbls$sums_tbl |>
        dplyr::filter(Metric == 'AUC', Level == msn) |>
        plotAUC('Sample', 'Value')
      say('AUC calculated and plot created')
    
      # TIC fluctuations
      if(tolower(msn) == 'ms') {
        df_fluc <- ticFluc(out$tic.df, cutoff = opts$plot.fluctuation.threshold)
        export_tbls$sums_tbl <- rbind(export_tbls$sums_tbl, df_fluc$export)
        
        img_tbl[nrow(img_tbl)+1,] <- plotFluctuations(df_fluc$data, limits = opts$plot.fluctuation.threshold)
        say('TIC fluctuations calculated and plot created')

      # Correlation matrix
        img_tbl[nrow(img_tbl)+1,] <- plotCorrelationMatrix(out$tic.df, 'RT', 'intensity')
        say('Sample correlation calculated and plot created')
      }
      
      if(first_pass) first_pass <- FALSE
    }
    
#----------------- SAMPLE PLOTS -----------------
    sample_plots <- tibble::tibble(sample = character(),
                                   plot = list(),
                                   width = numeric(),
                                   height = numeric())
    samples <- getSampleNames(b)
    
    for (s in samples) {
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
    img_tbl <- saveCreatedPlots(img_tbl)
    
#----------------- HTML OUTPUT ------------------
    # Rename colnames in sums
    export_tbls$sums_tbl <- export_tbls$sums_tbl |>
      spread(key = Level, value = Value)
    
    # Create HTML ouput
    cal_header <- ifelse(opts$calibration.bool,
                         paste0('<br><b>Calibration curve:</b> ', opts$calibration.tag, '\n'), '')
    
    createHTMLreport(header = paste0('<b>Sample batch:</b> ', b, '\n', cal_header,'<br>'), 
                     title = 'QC sample analysis', img.table = subset(img_tbl, include.in.report),
                     export.tables = export_tbls, output.dir = opts$sample.dir,
                     qc = (opts$qc.bool & opts$qc.identify), cal = opts$calibration.bool)
    
    # Remove files that have been created only as an intermediate step
    invisible(file.remove(list.files(opts$plots.dir, pattern = "^REMOVE_", full.names = TRUE)))
  }
},

error = function(err) {
  e(err)
  quit()
})

diff <- Sys.time() - time_start

say(paste('Done!', new_tic, 'new sample data and', new_qc, 
          'new QC peptide data processed in', round(as.numeric(diff), 2), units(diff), '\n'))