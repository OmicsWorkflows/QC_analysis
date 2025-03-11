# Function for trapezoid area calculation
trapezoidArea <- function(i_start, i_end, x, y) {
  width <- x[i_end] - x[i_start]
  height <- (as.numeric(y[i_start]) + as.numeric(y[i_end]))/2
  area <- width*height
  return(area)
}

# AUC function utilizing the trapezoid area method
ticAUC <- function(rt, intensity) {
  sum(sapply(1:(length(rt)-1), function(x) trapezoidArea(x, x+1, rt, intensity)), na.rm = TRUE)
}

# Calculates fluctuations
fluc <- function(i, x) {
  diff <- x[i]/x[i+1]
  return(diff)
}

# Wrapper for fluctuation calculation
ticFluc <- function(x, cutoff) {
  l <- lapply(unique(x$sample), function(s) {
    tmp <- dplyr::filter(x, sample == s)
    fl <- sapply(1:(length(tmp$intensity)-1), function(i) fluc(i, tmp$intensity))
    d <- data.frame(sample = s,
                    fluctuations = fl,
                    rt = tmp$RT[2:length(tmp$RT)])
    return(d)
  })
  
  data <- do.call(rbind, l)
  
  f <- paste0('Signal fall (', cutoff, 'x) count')
  j <- paste0('Signal jump (', cutoff, 'x) count')
  
  export <- data |>
    group_by(sample) |>
    summarise(Fall = sum(fluctuations <= 1/cutoff),
              Jump = sum(fluctuations >= cutoff)) |>
    as.data.frame() |>
    gather(-sample, key = 'Level', value = 'Value') |>
    mutate(Metric = paste0('Signal fluctuation (', cutoff, 'x) in MS1 count')) |>
    rename(Sample = sample)
  
  return(list(data = data,
              export = export,
              cutoff = cutoff))
}

# Wrapper for calculating AUC and inserting them into an export table
calculateAUC <- function(x, remove.bg, ms) {
  if(remove.bg){
    bg <- sapply(unique(x$sample), function (s) {
      tmp <- x |>
        dplyr::filter(sample == s)
      val <- trapezoidArea(1, nrow(tmp), tmp$RT, tmp$intensity)
      return(val)
    })
  }
  
  # Calculate areas under curve
  d <- x |>
    group_by(sample) |>
    summarise(Metric = 'AUC',
              Level = ms,
              Value = ticAUC(RT, intensity)) |>
    rename(Sample = sample) |>
    as.data.frame() |>
    ungroup()
  
  if(remove.bg) {
    d$Value <- round(d$Value - bg)
  } else {
    d$Value <- round(d$Value)
  }
  
  return(d)
}

# Calculate TIC quartiles
calculateQuartiles <- function(d, x, y, cumsum = TRUE) {
  d$x.norm <- d[,x] - min(d[,x], na.rm = TRUE)
  if(cumsum) {
    d[,y] <- cumsum(d[,y])/sum(d[,y], na.rm = TRUE)
  } else {
    d[,y] <- d[,y]/max(d[,y], na.rm = TRUE)
  }
  
  x.diff <- max(d$x.norm, na.rm = TRUE)
    
  qs <- sapply(c(0.25,0.5,0.75,1), function(qv) d$x.norm[max(which(d[,y]<=qv))])/x.diff
  return(list(cumulative = d[,y],
              quartiles = qs))
}

# Calculate RT quartiles
calculateRTquartiles <- function(d, x, y) {
  d$x.norm <- d[,x] - min(d[,x], na.rm = TRUE)
  d$cumulative <- cumsum(d[,y])/sum(d[,y], na.rm = TRUE)
  x.diff <- max(d$x.norm, na.rm = TRUE)
  
  qs <- sapply(c(0.25,0.5,0.75,1), function(qv) d$x.norm[max(which(d$cumulative<=qv))])/x.diff
  return(list(cumulative_intensity = d$cumulative,
              quartiles = qs))
}

# Plot all fluctuations
plotFluctuations <- function(d, limits, tic.opts = opts, subtitle = opts$plot.subtitle,
                             ms.level = msn,
                             filename = paste0(removeSpecialCharacters(opts$file.batch), '_', msn, '_fluctuation.png')) {
  
  d <- mutate(d, y = log2(fluctuations), x = rt)
  
  for (n in names(tic.opts)){
    assign(n, tic.opts[[n]])
  }
  
  d$sample <- plotAdjustNames(d$sample, plot.prefix, plot.suffix)
  subtitle <- paste0(subtitle, '\nred dashed line denotes log2(', limits,') and -log2(', limits,') limits')
  
  p <- plotChromatogramDefault(d, subtitle, 'TIC fluctuation between scans') +
    geom_hline(yintercept = log2(limits), color = 'red', linetype = 'dashed') +
    geom_hline(yintercept = -log2(limits), color = 'red', linetype = 'dashed') +
    scale_color_manual(values = rep('mediumseagreen', times = length(unique(d$sample)))) +
    theme(legend.position = 'none') + 
    facet_wrap(.~sample, ncol = plot.max.columns, scales = 'fixed', strip.position = "top") +
    labs(y = 'Fluctuations (log2)')

  # Note this plot in image table
  tab <- tibble(name = 'TIC fluctuation',
                link = 'tic_fluctuation',
                plot = list(p),
                filename = filename, 
                width = plot.width, 
                height = plot.height, 
                type = 'QC metrics',
                level = ms.level,
                full.path = NA,
                include.in.report = TRUE)
  return(tab)
}

# Split barchart labels into rows
splitIntoRows <- function(x, length = 20){
  chunks <- character()
  
  while (nchar(x) > length) {
    current.length <- length
    current <- substr(x, current.length, current.length)
    step <- 1
    
    while (!current %in% spec_chars) {
      current.length.plus <- current.length + step
      current.length.minus <- current.length - step
      
      if(substr(x, current.length.plus, current.length.plus) %in% spec_chars) {
        current.length <- current.length.plus
        current <- substr(x, current.length,  current.length)
      } else if(substr(x, current.length.minus, current.length.minus) %in% spec_chars) {
        current.length <- current.length.minus
        current <- substr(x, current.length,  current.length)
      } else {
        step <- step + 1
      }
    }
    
    chunks[length(chunks)+1] <- substr(x, 1, current.length)
    x <- substr(x, current.length+1, nchar(x))
  }
  
  chunks[length(chunks)+1] <- x
  return(paste(chunks, collapse = '\n'))
}

# Plot AUC
plotAUC <- function(d, x, y, 
                    subtitle = opts$plot.subtitle, ms.level = msn, 
                    filename = paste0(removeSpecialCharacters(opts$file.batch), '_', msn, '_AUC.png'),
                    tic.opts = opts) {
  
  colnames(d)[which(colnames(d) == x)] <- 'x'
  colnames(d)[which(colnames(d) == y)] <- 'y'
  
  d$x <- plotAdjustNames(d$x, tic.opts$plot.prefix, tic.opts$plot.suffix)
  #d$x <- sapply(as.character(d$x), splitIntoRows, USE.NAMES = FALSE)
  
  if(length(unique(d$x)) > tic.opts$plot.bar.max.samples) {
    n.facet <- ceiling(length(unique(d$x)) / tic.opts$plot.bar.max.samples)
    d$facet <- unlist(sapply(1:n.facet, function(n){
      if(n*tic.opts$plot.bar.max.samples <= length(unique(d$x))) {
        return(rep(n, times = tic.opts$plot.bar.max.samples))
      } else {
        return(rep(n, times = length(unique(d$x))-((n-1)*tic.opts$plot.bar.max.samples)))
      }
    }))
  } else {
    d$facet <- 1
  }

  p <- ggplot(d, aes(x = x, y = y)) +
    labs(x = 'Sample', y = 'AUC', title = 'AUC comparison', subtitle = subtitle) +
    geom_col(position = position_dodge()) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.background = element_blank(),
          strip.text = element_blank()) +
    facet_wrap(facet~., ncol = 1, scales = 'free_x')
  
  height <- 400 * length(unique(d$facet))
  
  if(length(unique(d$x)) > 30 & length(unique(d$x)) < tic.opts$plot.bar.max.samples) {
    width <- 500 + 10*(max(table(d$facet))-30)
  } else {
    width <- 500
  }
  
  if(width > 800) width <- 800

  # Note this plot in image table
  tab <- tibble(name = 'AUC',
                link = paste('auc', tolower(ms.level), sep = '_'),
                plot = list(p),
                filename = filename, 
                width = width, 
                height = height, 
                type = 'QC metrics', 
                level = ms.level,
                full.path = NA,
                include.in.report = TRUE)
  return(tab)
}

# Plot corr matrix
plotCorrelationMatrix <- function(d, x, y, rt.decimals = 1,
                                  subtitle = opts$plot.subtitle, ms.level = msn, 
                                  filename = paste0(removeSpecialCharacters(opts$file.batch), '_', msn, '_corr.png'),
                                  tic.opts = opts) {

  colnames(d)[which(colnames(d) == x)] <- 'x'
  colnames(d)[which(colnames(d) == y)] <- 'y'
  
  d$sample <- plotAdjustNames(d$sample, tic.opts$plot.prefix, tic.opts$plot.suffix)
  
  d <- suppressMessages(d |> 
    mutate(x = round(x, digits = rt.decimals)) |>
    group_by(x, sample) |>
    summarise(y = mean(y, na.rm = TRUE)) |>
    as.data.frame() |>
    spread(key = 'sample', value = 'y'))
  rownames(d) <- d$x
  d$x <- NULL

  m <- cor(d) |> as.data.frame()
  m$x <- rownames(m)
  rownames(m) <- NULL
  m <- gather(m, -x, key = 'y', value = 'corr.coef') |>
    mutate(text.col = ifelse(corr.coef < 0.75 & corr.coef > -0.75, 'grey10', 'white'))
  
  p <- ggplot(m, aes(x = x, y = y, fill = corr.coef, label = round(corr.coef, digits = 2))) +
    geom_tile() +
    labs(title = 'Sample correlation matrix', subtitle = subtitle) +
    #geom_text(aes(color = text.col), size = 2.5) +
    scale_y_discrete(limits = rev) +
    scale_color_identity() +
    theme_bw() +
    scale_fill_gradientn(colours = c('#b2182b', '#ef8a62', '#fddbc7', '#f7f7f7', 
                                     '#d1e5f0', '#67a9cf', '#2166ac'), 
                         limit = c(min(m$corr.coef), max(m$corr.coef)), space = "Lab", 
                         name="Pearson\ncorrelation") +
    coord_fixed() +
    guides(fill = guide_colourbar(barwidth = 10)) +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.title = element_text(vjust = 1, hjust = 0.5),
          legend.position = 'bottom')
  
  if(length(colnames(d)) <= 20) {
    width <- 500
    height <- 500
  } else {
    width <- 500 + 10 * (length(colnames(d)) - 20)
    height <- 500 + 10 * (length(colnames(d)) - 20)
  }
  
  if(width > 800) width <- 800
  if(height > 800) height <- 800

  # Note this plot in image table
  tab <- tibble(name = 'Sample correlation matrix',
                link = 'corr_matrix',
                plot = list(p),
                filename = filename, 
                width = width, 
                height = height,
                type = 'QC metrics',
                level = ms.level,
                full.path = NA,
                include.in.report = TRUE)
  return(tab)
}

