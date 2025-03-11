plotLCtraces <- function(d, t, current.batch, opts) {
  tmp <- d |>
    dplyr::filter(type == t)
  
  tmp$sample <- plotAdjustNames(tmp$sample, opts$plot.prefix, opts$plot.suffix)
  
  subtypes <- unique(tmp$subtype)
  
  if(length(subtypes) > 1) { # create several plots and join them or just one?
    subplot.list <- list()
    
    # Create separate plots for each subplot
    for (i in 1:length(subtypes)) {
      subt <- subtypes[i]
      tmp2 <- tmp |>
        dplyr::filter(subtype == subt)
      y.var <- unique(tmp2$variable)
      y.var <- paste0(toupper(substr(y.var, 1, 1)), substr(y.var, 2, nchar(y.var)))
      
      colnames(tmp2)[which(colnames(tmp2) == 'RT')] <- 'x'
      
      # Get appropriate plot color from dataframe
      type.colors <- lcColors |>
        dplyr::filter(type == t)
      
      if(nrow(type.colors) == 1) {
        color <- type.colors$colors
      } else {
        color <- type.colors$colors[which(type.colors$subtype == subt)]
      }
      
      # Plot with or without title
      if(i == 1) {
        p <- plotChromatogramDefault(tmp2, subt, t, y.var = y.var)
      } else {
        p <- plotChromatogramDefault(tmp2, subt, NA, y.var = y.var)
      }
      
      p <- p +
        scale_color_manual(values = rep(color, times = length(unique(tmp2$sample)))) +
        theme(legend.position = 'none') +
        facet_wrap(.~sample, ncol = opts$plot.max.columns, 
                   scales = 'fixed', shrink = FALSE, strip.position = 'top')
      
      tab <- tibble(plot = list(p),
                    subtype = subt,
                    height = (120*ceiling(length(samples)/opts$plot.max.columns)) + 100)
      
      subplot.list[[i]] <- tab
    }
    
    subplots <- do.call(rbind, subplot.list)
    
    subplots$row <- as.integer(rownames(subplots))
    subplots$column <- 1

    # Get plot size
    height <- sum(subplots$height) + 100
    width <- ifelse(length(unique(tmp2$sample)) < opts$plot.max.columns,
              500*length(unique(tmp2$sample)), 
              500*opts$plot.max.columns)
    
    p_final <- subplots |>
      ggplot(aes(x = column, y = row)) +
      geom_plot_npc(aes(npcx = 0, npcy = 0, label = plot, vp.width = 1, vp.height = 1)) +
      facet_grid(rows = vars(row), cols = vars(column), scales = 'free') +
      theme_void() +
      theme(strip.text = element_blank(),
            plot.background = element_rect(fill = 'white', color = 'white'))
  } else {
    y.var <- unique(tmp$variable)
    y.var <- paste0(toupper(substr(y.var, 1, 1)), substr(y.var, 2, nchar(y.var)))
    
    colnames(tmp)[which(colnames(tmp) == 'RT')] <- 'x'
    
    # Get appropriate plot color from dataframe
    type.colors <- lcColors |>
      dplyr::filter(type == t)
    
    color <- type.colors$colors
    
    p_final <- plotChromatogramDefault(tmp, subtypes, t, y.var = y.var)  +
      scale_color_manual(values = rep(color, times = length(unique(tmp$sample)))) +
      theme(legend.position = 'none') +
      facet_wrap(.~sample, ncol = opts$plot.max.columns,
                 scales = 'fixed', shrink = FALSE, strip.position = 'top')
   
    # Get plot size
    height <- (120*ceiling(length(unique(tmp$sample))/opts$plot.max.columns)) + 100
    width <- ifelse(length(unique(tmp$sample)) < opts$plot.max.columns,
                    500*length(unique(tmp$sample)), 
                    500*opts$plot.max.columns)
  }
  
  tab <- tibble(name = t,
                link = tolower(sub(' ', '_', t, fixed = TRUE)),
                plot = list(p_final),
                filename = paste0(removeSpecialCharacters(current.batch), '_', 
                                  sub(' ', '_', t, fixed = TRUE),'.png'),
                width = width, 
                height = height, 
                type = 'Additional traces', 
                level = 'MS',
                full.path = NA,
                include.in.report = TRUE)
  return(tab)
}
