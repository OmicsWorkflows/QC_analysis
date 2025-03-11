# Retrieve all dataframes associated with given sample
getSampleData <- function(x, opts) {
  data <- list()
  
  for (msn in opts$MS.levels$name) {
    tmp <- loadCalculatedTable(x, paste0(opts$tic.dir, msn, '/'), suffix = '.csv')
    if(!is.null(tmp)) {
      tmp$type <- paste('TIC', msn)
      tmp$color <- ticColors$colors[ticColors$type == 'TIC' & ticColors$level == msn]
    }
    
    data[[length(data)+1]] <- tmp
    if(!is.null(tmp)) names(data)[length(data)] <- paste('TIC', msn)
    
    if(msn == 'MS') {
      tmp <- loadCalculatedTable(x, paste0(opts$bpc.dir, '/'), suffix = '.csv')
      if(!is.null(tmp)) {
        tmp$type <- 'BPC'
        tmp$color <- ticColors$colors[ticColors$type == 'BPC']
      }
      data[[length(data)+1]] <- tmp
      if(!is.null(tmp)) names(data)[length(data)] <- 'BPC'
    }
  }
  
  if (length(data) > 0) {
    other.dirs <- c('eic_data', 'mpc_data', 'LC_data')
    other.dirs <- other.dirs[other.dirs %in% list.dirs(opts$file.output, full.names = FALSE, recursive = FALSE)]
    
    if(length(other.dirs) > 0) {
      other.dirs.full <- paste0(opts$file.output, other.dirs, '/')
      
      for (dir in other.dirs.full) {
        subdirs <- file.info(list.files(dir, full.names = TRUE, include.dirs = TRUE)) |>
          dplyr::filter(isdir) |>
          rownames_to_column('folder') |>
          pull(folder)
        
        type <- sub(opts$file.output, '', dir)
        type <- toupper(sub('_data/', '', type))
        
        if(length(subdirs) > 0) {
          for(subdir in subdirs) {
            subtype <- sub(dir, '', subdir)
            subtype <- sub('_', ' ', subtype)
            
            subdir.full <- paste0(subdir, '/')
            tmp <- loadCalculatedTable(x, subdir.full, suffix = '.csv')
            
            if(!is.null(tmp)) {
              colors <- lcColors[lcColors$type == type,]
              tmp$type <- paste0(type, ' (', subtype, ')')
              if(nrow(colors) > 1) {
                tmp$color <- colors$colors[colors$subtype == subtype]
              } else {
                tmp$color <- colors$colors
              }
            }
            
            data[[length(data)+1]] <- tmp
            if(!is.null(tmp)) names(data)[length(data)] <- paste0(type, ' (', subtype, ')')
          }
        } else {
          tmp <- loadCalculatedTable(x, dir, suffix = '.csv')
          
          if(!is.null(tmp)) {
            tmp$color <- colors$colors[colors$type == type]
            tmp$type <- type
          }
          
          data[[length(data)+1]] <- tmp
          if(!is.null(tmp)) names(data)[length(data)] <- type
        }
      }
    
    }
    
    if (opts$qc.bool) {
      data$`QC peptides` <- loadCalculatedTable(x, opts$qc.dir, suffix = '.csv')
      data$`QC peptides` <- merge(data$`QC peptides`, qcColors)
    }
    
    if (exists('df_fluc')) data$`TIC fluctuation` <- df_fluc$data |> 
      dplyr::filter(sample == x) |>
      mutate(color = 'mediumseagreen')
    
    return(data)
  } else {
    return(NULL)
  }
}

# Create plot for a single sample
plotSampleData <- function(s, l, opts) {
  # Create a plot for TICs
  tic.names <- c(paste('TIC', opts$MS.levels$name), 'BPC')
  tic.names <- if(length(tic.names) > 2) tic.names[c(1, length(tic.names), 2:(length(tic.names)-1))]

  tic.tmp <- do.call(rbind, l[which(names(l) %in% tic.names)])
  colnames(tic.tmp)[which(colnames(tic.tmp) == 'RT')] <- 'x'
  colnames(tic.tmp)[which(colnames(tic.tmp) == 'intensity')] <- 'y'
  
  tic.tmp$type <- factor(tic.tmp$type, levels = tic.names)
  
  start <- opts$calibration.start
  end <- opts$calibration.end
  
  if(end != 0) {
    tmp.limit <- tic.tmp |>
      dplyr::filter(x >= start & x <= end) |>
      group_by(type) |>
      summarise(max.y = max(y, na.rm = TRUE)*1.05)
    tic.tmp <- merge(tic.tmp, tmp.limit) |>
      mutate(y = ifelse(y > max.y, max.y, y))
  } else {
    tic.tmp <- tic.tmp |>
        group_by(type) |>
        mutate(max.y = max(y, na.rm = TRUE)*1.05) |>
        ungroup()
  }
  
  tic.tmp <- tic.tmp |>
    mutate(sample = factor(type, levels = levels(type)))

  tmp.p <- plotChromatogramDefault(tic.tmp, title = NA, subtitle = NA, start = start, end = end) +
    theme(legend.position = 'none') +
    geom_rect(aes(xmin = min(x, na.rm = TRUE), xmax = max(x, na.rm = TRUE), 
                  ymin = 0, ymax = max.y, group = sample), alpha = 0) +
    facet_grid(sample~., scales = 'free_y') +
    scale_color_manual(values = unique(tic.tmp$color))
  
  plots <- tibble(row = 1, column = 1, plot = list(tmp.p))
  
  tic.p <- plots |>
    ggplot(aes(x = column, y = row)) +
    geom_plot_npc(aes(npcx = 0, npcy = 0, label = plot, vp.width = 1, vp.height = 1)) +
    theme_void() +
    labs(subtitle = 'Chromatograms', title = s) +
    theme(strip.text = element_blank(),
          plot.background = element_rect(fill = 'white', color = 'white'))
  
  # Create LC trace plots
  additional.traces <- names(l)[!names(l) %in% c('QC peptides', 'BPC',
                                                 paste('TIC', opts$MS.levels$name),
                                                 'TIC fluctuation')]
  if(length(additional.traces) > 0) {
    plot_list <- list()
    
    for (trace in additional.traces){
      tmp <- l[[trace]]
      y.var <- colnames(tmp)[2]
      colnames(tmp)[2] <- 'y'
      colnames(tmp)[1] <- 'x'
      
      i <- which(additional.traces == trace)
      
      tmp.p <- plotChromatogramDefault(tmp, title = s, y.var = y.var, subtitle = NA, start = 0, end = 0) +
        facet_grid(type~., scales = 'free_y') +
        scale_color_manual(values = unique(tmp$color)) +
        theme(legend.position = 'none',
              plot.title = element_blank())
      
      plot_list[[length(plot_list)+1]] <- tibble(row = i,
                                                 plot = list(tmp.p))
    }
    
    plots <- do.call(rbind, plot_list)
    plots$column <- 1
    
    lc.p <- plots |>
      ggplot(aes(x = column, y = row)) +
      geom_plot_npc(aes(npcx = 0, npcy = 0, label = plot, vp.width = 1, vp.height = 1)) +
      facet_grid(rows = vars(row), cols = vars(column), scales = 'free') +
      theme_void() +
      labs(subtitle = 'Additional traces') +
      theme(strip.text = element_blank(),
            plot.background = element_rect(fill = 'white', color = 'white'),
            panel.spacing = unit(0, "lines"))
  }
  
  # Create a plot for QC peptides
  if(opts$qc.bool) {
    tmp.p <- createSeparateQCplot(s, l$`QC peptides`, which(samples == s), opts$qc.identify,
                                 opts, by = 'sample', batch = b, inset.bool = opts$qc.identify)$plot[[1]] +
      theme(plot.subtitle = element_blank())
    plots <- tibble(row = 1, column = 1, plot = list(tmp.p))
    
    qc.p <- plots |>
      ggplot(aes(x = column, y = row)) +
      geom_plot_npc(aes(npcx = 0, npcy = 0, label = plot, vp.width = 1, vp.height = 1)) +
      theme_void() +
      labs(subtitle = 'QC peptides') +
      theme(strip.text = element_blank(),
            plot.background = element_rect(fill = 'white', color = 'white'))
  }
  
  # Composite plot
  # different included plots and heights based on available data for a given sample
  if(opts$qc.bool) {
    if(length(additional.traces) > 0) {
      all.p <- arrangeGrob(tic.p, lc.p, qc.p, ncol = 1, 
                                        heights = c(length(unique(tic.tmp$type))*1.25,
                                        length(additional.traces)*1.5,
                                        length(unique(l$`QC peptides`$QC))))
    } else {
      all.p <- arrangeGrob(tic.p, qc.p, ncol = 1, 
                           heights = c(length(unique(tic.tmp$type))*1.25,
                                       length(unique(l$`QC peptides`$QC))))
    }
    plot.height <- ((length(unique(tic.tmp$type))+
                       length(additional.traces)+
                       length(unique(l$`QC peptides`$QC)))*100)+100
  } else {
    if(length(additional.traces) > 0) {
      all.p <- arrangeGrob(tic.p, lc.p, ncol = 1, 
                           heights = c(length(unique(tic.tmp$type))*1.25,
                                       length(additional.traces)*1.5))
    } else {
      all.p <- arrangeGrob(tic.p, ncol = 1)
    }
    plot.height <- ((length(unique(tic.tmp$type))+
                       length(additional.traces))*100)+100
  }
  
  out <- tibble::tibble(sample = s,
                        plot = list(all.p),
                        width = 500,
                        height = plot.height)
  
  return(out)
}

createJoinedSamplePlots <- function(plots, opts) {
  plots$number <- as.integer(rownames(plots))
  plots$column <- ((plots$number+1) %% opts$plot.max.columns)+1
  plots$row <- rep(1:ceiling(nrow(plots)/opts$plot.max.columns), each = opts$plot.max.columns)[1:nrow(plots)]
  
  plot_list <- list()
  
  # Iterate per row
  for(r in unique(plots$row)){
    tmp <- plots |>
      dplyr::filter(row == r)
    tmp_plot <- do.call(arrangeGrob, c(tmp$plot, nrow = 1))
    tmp_out <- tibble(name = list(tmp$sample),
                      link = list(tmp$sample),
                      plot = list(tmp_plot),
                      filename = paste0('REMOVE_row', r, '.png'),
                      width = sum(tmp$width),
                      height = max(tmp$height),
                      type = 'Individual samples',
                      level = 'Individual samples',
                      full.path = NA,
                      include.in.report = TRUE)
    plot_list[[length(plot_list)+1]] <- tmp_out
  }
  
  plot_list <- do.call(rbind, plot_list)
  return(plot_list)
}
