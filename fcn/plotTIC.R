# Plots TIC
plotTicDefault <- function(d, subtitle, title, start = 0, end = 0){
  p <- d |>
    ggplot(aes(x = x, y = y, color = sample)) +
    annotate('rect', xmin = start, xmax = end, ymin = -Inf, ymax = Inf,
             color = 'transparent', alpha = 0.15) +
    geom_line() +
    theme_bw() +
    labs(x = 'RT [min]', y = 'Intensity', subtitle = subtitle, title = title) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    expand_limits(y = 0)
  return(p)
}

# Plot TIC with fixed axes
plotTicFixed <- function(d, x, y, tic.opts = tic_opts,
                         start = 0, end = 0, calibration.bool = cal_bool,
                         filename = paste0(plots_dir, b, '_', msn, '_fixedY.png')) {
  colnames(d)[which(colnames(d) == x)] <- 'x'
  colnames(d)[which(colnames(d) == y)] <- 'y'
  
  for (n in names(tic.opts)){
    assign(n, tic.opts[[n]])
  }
  
  p <- plotTicDefault(d, subtitle, 'Individual TIC curves with fixed Y axis', start = start, end = end) +
    scale_color_manual(values = rep('black', times = length(d$sample))) +
    theme(legend.position = 'none') + 
    facet_wrap(.~sample, ncol = n.col, shrink = FALSE, strip.position = "top")
  
  ggsave(filename, plot = p, scale = 4,
         width = width, height = height, units = "px", limitsize = FALSE)
  
  say(filename, type = 'output')
  
  tab <- img_tbl_template
  tab[1,] <- c('Fixed Y axis', filename, width, height, 
           'Individual TIC curves', msn)
  return(tab)
}

plotTicFree <- function(d, x, y, tic.opts = tic_opts,
                        start = 0, end = 0, calibration.bool = cal_bool,
                        filename = paste0(plots_dir, b, '_', msn, '_freeY.png')) {
  colnames(d)[which(colnames(d) == x)] <- 'x'
  colnames(d)[which(colnames(d) == y)] <- 'y'
  
  for (n in names(tic.opts)){
    assign(n, tic.opts[[n]])
  }
  
  p <- plotTicDefault(d, subtitle, 'Individual TIC curves with free Y axis', start = start, end = end) +
    scale_color_manual(values = rep('black', times = length(d$sample))) +
    theme(legend.position = 'none') + 
    facet_wrap(.~sample, ncol = n.col, scales = 'free_y', shrink = FALSE, strip.position = "top")
  
  ggsave(filename, plot = p, scale = 4,
         width = width, height = height, units = "px", limitsize = FALSE)
  
  say(filename, type = 'output')
  
  tab <- img_tbl_template
  tab[1,] <- c('Free Y axis', filename, width, height, 
               'Individual TIC curves', msn)
  return(tab)
}

plotTicOverlaySingle <- function(d, subtitle) {
  pal <- colorRampPalette(calibrationColors$overlayColor)
  
  plotTicDefault(d, subtitle, 'All TIC curves overlayed') +
    guides(color = guide_legend(byrow = TRUE, ncol = 1, title = 'Sample', override.aes = list(linewidth = 1))) +
    scale_color_manual(values = pal(length(unique(d$sample)))) +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black'),
          legend.spacing.y = unit(0, "points"))
}

# Plot overlay of all TICs
plotTicOverlay <- function(l, d, x, y,
                           start = 0, end = 0, calibration.bool = cal_bool, 
                           batch = b, tag = '', output.dir = plots_dir, ms = msn) {
  sub_suffix <- ifelse(cal_bool, paste('| calibration:', tag, '| level:', msn), paste('| level:', msn))
  
  tbs <- lapply(1:length(l), function(i) {
    tmp <- d |>
      dplyr::filter(sample %in% unlist(l[i]))
    colnames(tmp)[which(colnames(tmp) == x)] <- 'x'
    colnames(tmp)[which(colnames(tmp) == y)] <- 'y'
    
    sub <- ifelse(length(l) == 1, paste(b, sub_suffix), 
                  paste(b, 
                        paste0('[', i, '/', length(l), ']'),
                        sub_suffix))
    
    p <- plotTicOverlaySingle(tmp, sub)
    path <- paste0(output.dir, b, '_', ms, '_overlay', i, '.png')
    ggsave(path, plot = p, scale = 4, width = 800, height = 500, 
           units = "px", limitsize = FALSE)
    
    say(path, type = 'output')
    
    # Note this plot in image table
    tab <- img_tbl_template
    tab[1,] <- c(paste0(i, '/', length(l)), path, 
                 800, 500, 'TIC curves overlay', ms)
    return(tab)
  })
  
  return(do.call(rbind, tbs))
}

# Plot all fluctuations
plotFluctuations <- function(d, limits, tic.opts = tic_opts, ms.level = msn,
                             filename = paste0(plots_dir, b, '_', msn, '_fluctuation.png')) {
  fl_list <- lapply(unique(d$sample), function(s) {
    tmp <- dplyr::filter(d, sample == s)
    fl <- ticFluc(tmp$intensity)
    d <- data.frame(sample = s,
                    fluctuations = fl,
                    rt = tmp$RT[2:length(tmp$RT)])
    return(d)
  })
  
  fl_df <- do.call(rbind, fl_list) |>
    mutate(y = log2(fluctuations),
           x = rt)
  
  for (n in names(tic.opts)){
    assign(n, tic.opts[[n]])
  }
  
  p <- plotTicDefault(fl_df, subtitle, 'TIC fluctuation between scans') +
    geom_hline(yintercept = log2(limits), color = 'red', linetype = 'dashed') +
    geom_hline(yintercept = -log2(limits), color = 'red', linetype = 'dashed') +
    scale_color_manual(values = rep('black', times = length(unique(d$sample)))) +
    theme(legend.position = 'none') + 
    facet_wrap(.~sample, ncol = n.col, scales = 'fixed', strip.position = "top") +
    labs(y = 'Fluctuations (log2)')
  
  ggsave(filename, plot = p, scale = 4, 
         width = width, height = height, units = "px", limitsize = FALSE)
  
  say(filename, type = 'output')
  
  # Note this plot in image table
  tab <- img_tbl_template
  tab[1,] <- c('TIC fluctuation', filename, width, height, 'QC metrics', ms.level)
  return(tab)
  
}

# Plot AUC
plotAUC <- function(d, x, y, 
                    subtitle = ifelse(cal_bool, paste(b, '| calibration:', tag, '| level': msn), 
                                      paste(b, '| level:', msn)),
                    ms.level = msn, filename = paste0(plots_dir, b, '_', msn, '_AUC.png')) {
  colnames(d)[which(colnames(d) == x)] <- 'x'
  colnames(d)[which(colnames(d) == y)] <- 'y'
  
  p <- ggplot(d, aes(x = x, y = y)) +
    labs(x = 'Sample', y = 'AUC', title = 'AUC comparison', subtitle = subtitle) +
    geom_col(position = position_dodge()) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  height <- 350
  ggsave(filename, plot = p, scale = 4, 
         width = 500, height = height, units = "px", limitsize = FALSE)
  
  say(filename, type = 'output')
  
  # Note this plot in image table
  tab <- img_tbl_template
  tab[1,] <- c('AUC', filename, 500, height, 'QC metrics', ms.level)
  return(tab)
}
