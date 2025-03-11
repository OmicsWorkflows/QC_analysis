#------------------ Loading files ------------------
# Function for getting already calculated files and those that need recalculation
alreadyCalculatedFiles <- function(s, bool = FALSE, dir, suffix = '.csv'){
  if (bool) {
    load <- NULL
    calc <- s
  } else {
    searched <- sapply(s, function(x) {
      paste0(x, suffix) %in% list.files(dir)
    })
    
    load <- names(searched)[which(searched == TRUE)]
    calc <- names(searched)[which(searched == FALSE)]
  }
  
  return(list(load = load, calc = calc))
}

# Load files that have been calculated previously
loadCalculatedTable <- function(x, dir, sample.name = 'sample', suffix = '.csv') {
  out <- list()
  
  for (xi in x) {
    fn <- paste0(dir, xi, suffix)
    if(file.exists(fn)){
      tmp <- read.table(fn, header = TRUE, sep = ',', strip.white = TRUE) |>
        mutate(!!sample.name := xi)
      out[[length(out)+1]] <- tmp
    }
  }
  
  return(do.call(rbind, out))
  
}

#------------ Function for redirecting -------------
getTICandQC <- function(x, msn, opts) {
  # For .raw and .mzML files
  if (opts$file.format %in% c('.mzML', '.raw')) {
    out <- list()
    
    if (opts$file.format == '.mzML') {
      out <- modifyList(out, getTICmzML(x, type = 'MS'))
      if (msn == 'MS') out <- modifyList(out, getBPCmzML(x))
    } else {
      out <- modifyList(out, getTICraw(x, type = 'MS'))
      if (msn == 'MS') out <- modifyList(out, getBPCraw(x))
    }

    # Get QC peptides?
    if (msn == 'MS' & opts$qc.bool) {
      out <- modifyList(out, getQCs(x, opts))
      lvls <- out$tic.df |> dplyr::select(sample, time) |> distinct() |> arrange(time) |> pull(sample)
      out$qc.df$sample <- factor(out$qc.df$sample, levels = lvls)
      out$qc.df <- arrange(out$qc.df, sample)
      opts <- modifyList(opts, initiateQcPlotVariables(out$qc.df, by = 'QC'))
      
      qc.output <- plotQc(out$qc.df, opts, by = 'QC')
      
      out$qc.plots <- qc.output$plot
      out$qc.stats <- qc.output$qc.stats
    }
    
  } else if(opts$file.format == '.d') {
    out <- getTICtimsD(x, type = 'MS')
    if(msn == 'MS' & opts$qc.bool) {
      lvls <- out$tic.df |> dplyr::select(sample, time) |> distinct() |> arrange(time) |> pull(sample)
      out$qc.df$sample <- factor(out$qc.df$sample, levels = lvls)
      out$qc.df <- arrange(out$qc.df, sample)
      out$qc.plots <- plotQcTims(out$qc.df, opts, by = 'QC')
    }
    
  } else if(opts$file.format == '.csv') { # If the file format is .csv, only load the previously calculated data
    out <- list()
    tic.subdir <- addSlash(paste0(opts$tic.dir, msn))
    
    # Load TIC data
    say('Locating previously calculated TIC data')
    
    out$tic.df <- loadCalculatedTable(x, tic.subdir)
    out$new.tic <- 0
    
    say(paste('TIC data loaded from', tic.subdir))
    
    # Load BPC data
    if(msn == 'MS') {
      say('Locating previously calculated BPC data')
      
      out$bpc.df <- loadCalculatedTable(x, opts$bpc.dir)
      out$new.bpc <- 0
      
      say(paste('BPC data loaded from', opts$bpc.dir))
    }

    # Load QC data
    if(msn == 'MS' & opts$qc.bool) {
      say('Locating previously calculated QC data')
      
      out$new.qc <- 0
      
      out$qc.df <- loadCalculatedTable(x, opts$qc.dir)

      lvls <- out$tic.df |> dplyr::select(sample, time) |> distinct() |> arrange(time) |> pull(sample)
      out$qc.df$sample <- factor(out$qc.df$sample, levels = lvls)
      out$qc.df <- arrange(out$qc.df, sample)
      
      opts <- modifyList(opts, initiateQcPlotVariables(out$qc.df, by = 'QC'))
      qc.output <- plotQc(out$qc.df, opts, by = 'QC')
      
      out$qc.plots <- qc.output$plot
      out$qc.stats <- qc.output$qc.stats
    }
  }
  
  # Look for other types of traces
  if(msn == 'MS') {
    other.dirs <- c('mpc_data', 'eic_data', 'LC_data')
    other.dirs <- other.dirs[other.dirs %in% list.dirs(opts$file.output, full.names = FALSE, recursive = FALSE)]
    
    say('Locating data for additional traces')
    
    if(length(other.dirs) > 0) {

      other.traces <- list()
      other.dirs.full <- paste0(opts$file.output, other.dirs, '/')
      
      for (dir in other.dirs.full) {
        subdirs <- file.info(list.files(dir, full.names = TRUE, include.dirs = TRUE)) |>
          dplyr::filter(isdir) |>
          rownames_to_column('folder') |>
          pull(folder)
        
        type <- sub(opts$file.output, '', dir)
        type <- toupper(sub('_data/', '', type))
        
        # If there are subdirectories, access them and read data from them
        if(length(subdirs) > 0) {
          for(subdir in subdirs) {
            subtype <- sub(dir, '', subdir)
            subtype <- sub('_', ' ', subtype)
            
            subdir.full <- paste0(subdir, '/')
            d <- loadCalculatedTable(x, subdir.full, suffix = '.csv')
            if(!is.null(d)){
              say(paste('Reading data of type', paste0(type, ' (', subtype, ')'),
                        'from', subdir))
              
              d$variable <- colnames(d)[2]
              colnames(d)[2] <- 'y'
              d$type <- type
              d$subtype <- subtype
              
              other.traces[[length(other.traces)+1]] <- d
            } else {
              say(paste('No data of subtype', subtype, 'of type', type, 'found'))
            }
          }
        } else {
          d <- loadCalculatedTable(x, dir, suffix = '.csv')
          
          if(!is.null(d)){
            say(paste('Reading data of type', type,
                      'from', dir))
          
            d$type <- type
            d$subtype <- NA
            other.traces[[length(other.traces)+1]] <- d
          } else {
            say(paste('No data of type', type, 'found'))
          }
        }
      }
      
      if (length(other.traces) > 0) {
        out$other.traces <- do.call(rbind, other.traces)
      }

    } else {
      say('No additional data found')
    }
  }
  
  return(out)
}

#---------------- .mzML processing -----------------
getDateFromMzML <- function(x, date.format) {
  r <- runInfo(openMSfile(x))$startTimeStamp
  f <- format(as.POSIXct(r, format = date.format),
              '%Y-%m-%d %H:%M:%S')
  return(f)
}

processTICfileMzML <- function(x) {
  d <- data.frame(RT = x@rtime/60,
                  intensity = x@intensity,
                  row.names = NULL)
  return(d)
}

# Calculate chromatograms from data that haven't been previously calculated
calculateChromatogramMzML <- function(x, ms.level = 1, dir, output.dir, aggregation,
                                      config, sample.name = 'sample') {
  
  if(aggregation == 'sum') {
    type <- 'TIC'
  } else if(aggregation == 'max'){
    type <- 'BPC'
  }
  
  say(paste0('Reading .mzML files for ', type,' computation from:'))
  say(dir, type = 'input')
  
  # read .mzML files
  fn <- paste0(dir, x, '.mzML')
  readMS <- readMSData(fn, mode = "onDisk")
  
  # return success message
  say(paste0('.mzML files read succesfully'))
  
  # calculate chromatogram
  paste('Calculating', type, 'chromatogram data')
  t <- suppressWarnings(chromatogram(readMS, aggregationFun = aggregation, msLevel = ms.level,
                                     BPPARAM = BiocParallel::SerialParam()))
  
  say(paste(type, 'data calculated and saved as:'))
  
  # save as data frame
  l <- lapply(1:ncol(t), function(i) {
    tmp <- t[,i]
    date <- getDateFromMzML(paste0(dir, x[i], '.mzML'), config$file.date)
    d <- processTICfileMzML(tmp) |>
      mutate(time = date)
    
    d$intensity[which(is.na(d$intensity))] <- 0
    
    saveDataCsv(x[i], d, output.dir)
    
    d <- mutate(d, !!sample.name := x[i])
    return(d)
  })
  
  return(do.call(rbind, l))
}

# Wrapper function for mzMLs
getTICmzML <- function(x, type, ms.level = msn, suffix.csv = '.csv',
                       config = opts) {
  
  output.dir <- config$tic.dir
  counter <- config$new.tic
  
  if(type == 'MS') {
    ms.dict <- config$MS.levels
    ms.val <- ms.dict$value[ms.dict$name == ms.level]
    subdir <- addSlash(paste0(output.dir, ms.level))
  } else {
    subdir <- addSlash(paste0(output.dir, type))
    ms.val <- ms.level
  }
  
  # Check if files should be reanalyzed
  if(opts$file.reanalyze) {
    say('Initiating de novo TIC analysis of all samples')
    flist <- list(load = NULL,
                  calc = x)
  } else {
    say('Locating previously calculated TIC data')
    flist <- alreadyCalculatedFiles(x, dir = subdir, suffix = suffix.csv)
  }

  df <- NULL
  
  if(length(flist$load) > 0) {
    df <- loadCalculatedTable(flist$load, subdir, suffix = suffix.csv)
    say(paste('TIC data loaded from', subdir))
  }
  
  if(length(flist$calc) > 0) {
    calc <- calculateChromatogramMzML(flist$calc, 
                                      dir = config$file.input, 
                                      output.dir = subdir, 
                                      ms.level = ms.val, 
                                      aggregation = 'sum',
                                      config)
    df <- rbind(df, calc)
    
    counter <- counter + length(flist$calc)
  }
  
  return(list(tic.df = df,
              new.tic = counter))
}

#----------------- .raw processing -----------------
# calculate TICs from .raw files
calculateChromatogramRaw <- function(x, ms.level = msn, sample.name = 'sample', 
                                     cr.type, output.dir, config) {
  input.dir <- config$file.input
  
  say(paste('Reading .raw file for', toupper(cr.type),'chromatogram computation from', input.dir, 'for sample:'))
  say(x, type = 'input')
  
  fn <- paste0(input.dir, x, '.raw')
  cr <- readChromatogram(fn, filter = ms.level, type = cr.type)
  
  #success message
  say('.raw file read succesfully')
  
  date_read <- readFileHeader(fn)$`Creation date`
  date <- format(as.POSIXct(date_read, format = config$file.date),
                 '%Y-%m-%d %H:%M:%S')
  
  # Create output data frame
  d <- data.frame(RT = as.numeric(cr$times),
                  intensity = cr$intensities,
                  level = ms.level,
                  time = date,
                  row.names = NULL)
  
  d$intensity[which(is.na(d$intensity))] <- 0
  
  say(paste(toupper(cr.type), 'chromatogram data saved as:'))
  saveDataCsv(x, d, output.dir)
  d <- mutate(d, !!sample.name := x)
  
  return(d)
}

# Wrapper function for Thermo .raw
getTICraw <- function(x, type, ms.level = msn, suffix.csv = '.csv', 
                      config = opts) {
  
  output.dir <- config$tic.dir
  counter <- config$new.tic

  if(type == 'MS') {
    ms.dict <- config$MS.levels
    ms.val <- ms.dict$value[ms.dict$name == ms.level]
    subdir <- addSlash(paste0(output.dir, ms.level))
  } else {
    subdir <- addSlash(paste0(output.dir, type))
  }
  
  # Check if files should be reanalyzed
  if(opts$file.reanalyze) {
    say('Initiating de novo TIC analysis of all samples')
    flist <- list(load = NULL,
                  calc = x)
  } else {
    say('Locating previously calculated TIC data')
    flist <- alreadyCalculatedFiles(x, dir = subdir, suffix = suffix.csv)
  }

  df <- NULL

  if (length(flist$load) > 0) {
    df <- loadCalculatedTable(flist$load, subdir, suffix = suffix.csv)
    say(paste('TIC data loaded from', subdir))
  }
  
  if (length(flist$calc) > 0) {
    calc <- lapply(flist$calc, function(x) calculateChromatogramRaw(x, output.dir = subdir, cr.type = 'tic', config = config))
    calc <- do.call(rbind, calc)
    df <- rbind(df, calc)
    
    counter <- counter + length(flist$calc)
  }
  
  return(list(tic.df = df,
              new.tic = counter))
  
}

#------------------ .d processing ------------------
# Wrapper function for Bruker TIMS .d
getTICtimsD <- function(x, type, ms.level = msn,
                        suffix.csv = '.csv', config = opts) {
  
  if(type == 'MS') {
    qc.bool <- config$qc.bool & ms.level == 'MS'
      
    counter <- config$new.tic
    qc.counter <- config$new.qc
    bpc.counter <- config$new.bpc
    
    ms.val <- config$MS.levels$value[config$MS.levels$name == ms.level]
    subdir <- addSlash(paste0(config$tic.dir, ms.level))
    
    # Check if files should be reanalyzed
    if(opts$file.reanalyze) {
      say('Initiating de novo TIC analysis of all samples')
      flist <- list(load = NULL,
                    calc = x)
    } else {
      say('Locating previously calculated TIC data')
      flist <- alreadyCalculatedFiles(x, dir = subdir, suffix = suffix.csv)
    }
    
    if(qc.bool) {
      if(config$qc.reprocess) {
        flist.qc <- list(load = NULL,
                         calc = x)
        sapply(x, function(xi){
          fni <- paste0(config$qc.dir, xi, suffix.csv)
          if(file.exists(fni)) file.remove(fni)
        }, USE.NAMES = FALSE)
      } else {
        flist.qc <- alreadyCalculatedFiles(x, dir = config$qc.dir, suffix = suffix.csv)
      }
    }
    
    df <- NULL
    qc.df <- NULL
    bpc.df <- NULL
    
    # Reload previously calculated data
    if(length(flist$load) > 0) {
      df <- loadCalculatedTable(flist$load, subdir, suffix = suffix.csv)
      say(paste('TIC data loaded from', subdir))

      if(msn == 'MS'){
        bpc.df <- loadCalculatedTable(flist$load, config$bpc.dir, suffix = suffix.csv)
        say(paste('BPC data loaded from', config$bpc.dir))
      }
      
    }
    
    if(qc.bool) {
      if(length(flist.qc$load) > 0) {
        qc.df <- loadCalculatedTable(flist.qc$load, config$qc.dir, suffix = suffix.csv)
        say(paste('QC data loaded from', config$qc.dir))
      }
    }
    
    # If QC XICs should be calculated, but TIC has been calculated previously
    # it is necessary to open the same file again
    if(qc.bool) {
      if(length(setdiff(flist.qc$calc, flist$calc)) > 0){
        flist$calc <- c(flist$calc, setdiff(flist.qc$calc, flist$calc))
      }
    }
    
    # Calculate new data
    if(length(flist$calc) > 0) {
      py.script <- 'readRawBrukerTimsRustModified2.py'
      
      res <- sapply(flist$calc, function(y) {
        say(paste('Reading .d file from', config$file.input, 'for sample:'))
        say(y, type = 'input')
        
        qc.only <- ifelse(y %in% flist$load, 'True', 'False')
        y.full <- paste0(y, '.d')
        
        qc.path <- ifelse(qc.bool, config$qc.reference.matrix, '')
        qc.output <- ifelse(qc.bool, config$qc.dir, '')
        bpc.output <- ifelse(ms.level == 'MS', config$bpc.dir, '')
        
        command <- paste(config$python,
                         paste0('"', config$src, 'Python/', py.script, '"'), 
                         paste0('"',config$file.input,'"'), 
                         y.full,
                         ms.level, 
                         paste0('"', qc.path, '"'),
                         paste0('"', subdir, '"'),
                         paste0('"', bpc.output,'"'),
                         paste0('"', qc.output,'"'),
                         qc.only)
        result <- system(command, intern = TRUE, wait = TRUE)
        
        if (result == 'error') {
          say('Warning!')
          say(paste0('TIC curve could not be calculated for sample ', y,', ', py.script, ' encountered an error'), 
              type = 'warn')
          
          return(FALSE)
          
        } else {
          if (qc.only == 'False') {
            say('TIC data calculated and saved as:')
            say(paste0(subdir, y, '.csv'), type = 'output')
            
            if (msn == 'MS') {
              say('BPC data calculated and saved as:')
              say(paste0(config$bpc.dir, y, '.csv'), type = 'output')
            }
          }
          
          if(qc.bool) {
            say('QC data calculated and saved as:')
            say(paste0(config$qc.dir, y, '.csv'), type = 'output')
          }
          
          return(TRUE)
          
        }
      })
      
      success <- names(which(res == TRUE))
      flist$calc <- flist$calc[flist$calc %in% success]

      calc <- loadCalculatedTable(flist$calc, subdir, suffix = suffix.csv)
      df <- rbind(df, calc)
      counter <- counter + length(flist$calc)
      
      if(msn == 'MS') {
        bpc.calc <- loadCalculatedTable(flist$calc, config$bpc.dir, suffix = suffix.csv)
        bpc.df <- rbind(bpc.df, bpc.calc)
        bpc.counter <- bpc.counter + length(flist$calc)
      }
      
      if(qc.bool){
        flist.qc$calc <- flist.qc$calc[flist.qc$calc %in% success]
      
        qc.calc <- loadCalculatedTable(flist.qc$calc, config$qc.dir, suffix = suffix.csv)
        qc.df <- rbind(qc.df, qc.calc)
        qc.counter <- qc.counter + length(flist.qc$calc)
      }
      
      
    }
    
    return(list(tic.df = df,
                new.tic = counter,
                qc.df = qc.df,
                new.qc = qc.counter,
                bpc.df = bpc.df,
                new.bpc = bpc.counter))
    
  } else {
    say('Error!')
    say('Selected calculation type for TIMS .d files not available', type = 'error')
    quit()
  }
}

#---------- Create chromatogram plots -----------
# Just ggplot settings for a single chromatogram plot
plotChromatogramDefault <- function(d, subtitle, title, y.var = 'Intensity', start = 0, end = 0){
  p <- d |>
    ggplot(aes(x = x, y = y, color = sample)) +
    annotate('rect', xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
             color = 'transparent', alpha = 0.15) +
    geom_line() +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    expand_limits(y = 0)
  
  # add title to the plot?
  if(!is.na(title)) p <- p + ggtitle(title)
  
  # add subtitle to the plot? (while adding other labels)
  if(!is.na(subtitle)) {
    p <- p + labs(x = 'RT [min]', y = y.var, subtitle = subtitle) 
  } else {
    p <- p + labs(x = 'RT [min]', y = y.var)
  }
  
  return(p)
}

# Plot a fixed or free axis chromatogram plot
plotChromatogramSingle <- function(d, x, y, cr.type, type, subtitle, title = 'Auto', tic.opts = opts,
                          start = 0, end = 0, colors = 'black', strip.position = 'top', filename) {
  colnames(d)[which(colnames(d) == x)] <- 'x'
  colnames(d)[which(colnames(d) == y)] <- 'y'
  
  d$sample <- plotAdjustNames(d$sample, tic.opts$plot.prefix, tic.opts$plot.suffix)
  
  # Adjust color character vector
  if(length(colors) != length(unique(d$sample))) colors <- rep(colors, times = length(d$sample))
  
  if(type == 'free') {
    scales.type <- 'free_y'
  } else {
    scales.type <- 'fixed'
  }
  
  if(title == 'Auto') title <- paste('Individual', cr.type, 'curves with', type, 'Y axis')
  
  # change upper y limit so that curve doesn't reach the top of the plot
  if(end != 0) {
    if(type == 'fixed') {
      d.limit <- d |>
        dplyr::filter(x >= start & x <= end)
      d$max.y <- max(d.limit$y, na.rm = TRUE)*1.05
      d <- d |> mutate(y = ifelse(y > max.y, max.y, y))
    } else {
      d.limit <- d |>
        dplyr::filter(x >= start & x <= end) |>
        group_by(sample) |>
        summarise(max.y = max(y, na.rm = TRUE)*1.05)
      d <- merge(d, d.limit) |>
        mutate(y = ifelse(y > max.y, max.y, y))
    }
  } else {
    if(type == 'fixed') {
      d$max.y <- max(d$y, na.rm = TRUE)*1.05
    } else {
      d <- d |>
        group_by(sample) |>
        mutate(max.y = max(y, na.rm = TRUE)*1.05) |>
        ungroup()
    }
  }
  
  # create plot
  p <- plotChromatogramDefault(d, tic.opts$plot.subtitle, title, start = start, end = end) +
    scale_color_manual(values = colors) +
    theme(legend.position = 'none') +
    geom_rect(aes(xmin = min(x, na.rm = TRUE), xmax = max(x, na.rm = TRUE), 
                  ymin = 0, ymax = max.y, group = sample), alpha = 0) +
    facet_wrap(.~sample, ncol = tic.opts$plot.max.columns, scales = scales.type, 
               shrink = FALSE, strip.position = strip.position)
  
  tab <- tibble(name = paste0(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)),' Y axis'),
                link = paste(tolower(cr.type), type, tolower(msn), sep = '_'),
                plot = list(p),
                filename = filename, 
                width = tic.opts$ plot.width, 
                height = tic.opts$plot.height,
                type = paste('Individual', cr.type, 'curves'), 
                level = msn,
                full.path = NA,
                include.in.report = TRUE)
  return(tab)
}

plotChromatogramOverlaySingle <- function(d, type, subtitle, col, start, end, pos) {
  if(end != 0) {
    d.limit <- d |>
      dplyr::filter(x >= start & x <= end)
    d$max.y <- max(d.limit$y, na.rm = TRUE)*1.05
    d <- d |> mutate(y = ifelse(y > max.y, max.y, y))
  } else {
    d$max.y <- max(d$y, na.rm = TRUE)*1.05
  }
      
  if(length(unique(d$sample)) > length(col)){
    pal <- colorRampPalette(col)
    plot.cols <- pal(length(unique(d$sample)))
  } else {
    plot.cols <- col[1:length(unique(d$sample))]
  }
  
  plotChromatogramDefault(d, subtitle, paste('All', type, 'curves overlayed')) +
    guides(color = guide_legend(byrow = TRUE, ncol = ifelse(pos == 'right', 1, length(unique(d$sample))), 
                                title = 'Sample', override.aes = list(linewidth = 1))) +
    annotate('rect', xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
             color = 'transparent', alpha = 0.15) +
    geom_rect(aes(xmin = min(x, na.rm = TRUE), xmax = max(x, na.rm = TRUE), 
                  ymin = 0, ymax = max.y, group = sample), alpha = 0, show.legend = FALSE) +
    scale_color_manual(values = plot.cols) +
    theme(legend.position = pos,
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black'),
          legend.spacing.y = unit(0, "points"))
}

# Plot overlay of all TICs
plotChromatogramOverlay <- function(d, x, y, cr.type, start = 0, end = 0, chunk.no = opts$plot.max.samples,
                           colors = calibrationColors$overlayColor, tic.opts = opts,
                           calibration.bool = opts$calibration.bool, batch = opts$file.batch, tag = opts$calibration.tag, 
                           ms = msn, filename = 'auto') {

  d$sample <- plotAdjustNames(d$sample, tic.opts$plot.prefix, tic.opts$plot.suffix)
  
  l <- split(unique(d$sample), ceiling(seq_along(unique(d$sample))/chunk.no))
  
  if(filename == 'auto') {
    sub_suffix <- ifelse(calibration.bool, paste(' | calibration:', tag, '| level:', msn), paste(' | level:', msn))
    legend.pos <- 'right'
    width <- 800
  } else {
    sub_suffix <- tag
    legend.pos <- 'bottom'
    width <- 500
  }
  
  tbs <- lapply(1:length(l), function(i) {
    path <- ifelse(filename == 'auto', 
                   paste0(removeSpecialCharacters(batch), '_', ms, '_', tolower(cr.type), '_overlay', i, '.png'),
                   filename)
    tmp <- d |>
      dplyr::filter(sample %in% unlist(l[i]))
    colnames(tmp)[which(colnames(tmp) == x)] <- 'x'
    colnames(tmp)[which(colnames(tmp) == y)] <- 'y'
    
    sub <- ifelse(length(l) == 1, paste0(batch, sub_suffix), 
                  paste0(batch, ' [', i, '/', length(l), ']', sub_suffix))
    
    p <- plotChromatogramOverlaySingle(tmp, cr.type, sub, colors, start, end, legend.pos)
    
    # Note this plot in image table
    tab <- tibble(name = paste0(i, '/', length(l)),
                  link = paste(tolower(cr.type), i, tolower(ms), sep = '_'),
                  plot = list(p),
                  filename = path, 
                  width = width, 
                  height = 500, 
                  type = paste(cr.type, 'curves overlay'), 
                  level = ms,
                  full.path = NA,
                  include.in.report = TRUE)
    return(tab)
  })
  
  return(do.call(rbind, tbs))
}

# Plot TICs with calibration samples TIC
plotTicCalibration <- function(d, x, y, tic.opts = opts,
                               start = 0, end = 0, colors = 'black', 
                               strip.position = 'top', filename) {
  colnames(d)[which(colnames(d) == x)] <- 'x'
  colnames(d)[which(colnames(d) == y)] <- 'y'
  
  for (n in names(tic.opts)){
    assign(n, tic.opts[[n]])
  }
  
  if(end != 0) {
    d.limit <- d |>
      dplyr::filter(x >= start & x <= end) |>
      group_by(sample) |>
      summarise(max.y = max(y, na.rm = TRUE)*1.05)
    d <- merge(d, d.limit) |>
      mutate(y = ifelse(y > max.y, max.y, y))
  }
  
  d$sample.name <- plotAdjustNames(d$sample, plot.prefix, plot.suffix)
  d$sample <- d$type
  
  col.d <- d[,c('color', 'sample')] |> distinct() |>
    arrange(sample)
  
  p <- plotChromatogramDefault(d, plot.subtitle, 
                      paste('Individual TIC curves with calibration samples'), 
                      start = start, end = end) +
    scale_color_manual(values = col.d$color,
                       breaks = col.d$sample) +
    theme(legend.position = 'bottom',
          legend.title = element_blank(),
          #legend.key.size = unit(0.5, 'lines'),
          legend.background = element_rect(color = 'black'),
          legend.spacing.y = unit(0, "points")) + 
    guides(color = guide_legend(override.aes = list(linewidth = 1))) +
    facet_wrap(.~sample.name, ncol = plot.max.columns, 
               scales = 'free_y', shrink = FALSE, strip.position = strip.position)
  
  # Note this plot in image table
  tab <- tibble(name = 'With calibration',
                link = 'tic_with_calibration',
                plot = list(p),
                filename = filename, 
                width = plot.width, 
                height = plot.height,
                type = 'Individual TIC curves', 
                level = msn,
                full.path = NA,
                include.in.report = TRUE)
  return(tab)
}
