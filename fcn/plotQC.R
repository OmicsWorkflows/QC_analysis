# Create variables for plot outputs
initiateQcPlotVariables <- function(d) {
  tib <- tibble(plot_separate = list(), 
                plot_merge = list(), 
                sample = character(),
                height = numeric())
  
  d <- length(unique(d$QC))
  
  template <- data.frame(QC = numeric(),
                                matrix_correlation = numeric(),
                                idsl_auc = numeric(),
                                sym = numeric(),
                                SN_ratio = numeric(),
                                QC_peak_y = numeric(),
                                width = numeric(),
                                sample = character())
  
  assign('qc_plot_tibble', tib, envir = parent.frame())
  assign('qc_subplot_row', d, envir = parent.frame())
  assign('qc_stats', template, envir = parent.frame())
}

# Create a dataframe for creating an inset plot
createInputDfForInsetPlotting <- function(x) {
  df <- x |>
    mutate(rt = ifelse(rt < peak_start | rt > peak_end, NA, rt)) |>
    dplyr::filter(!is.na(rt)) |>
    group_by(QC) |>
    mutate(rt_at_max = ifelse(is.na(matrix_correlation), NA, rt[which.max(intensity)]),
           max_intensity = max(intensity, na.rm = TRUE),
           idsl_auc = ifelse(is.na(matrix_correlation), NA, peakAreaCalculator(rt, intensity)),
           sym = NA,
           SN_ratio = ifelse(is.na(matrix_correlation), NA, SNRbaseline(intensity, baseline)),
           width = ifelse(is.na(matrix_correlation), NA, peak_end - peak_start)) |>
    as.data.frame()
  return(df)
}


# Wrapper function for creating a separate QC plot
createSeparateQCplot <- function(x, df){
  tmp <- dplyr::filter(df, sample == x) |>
    group_by(QC) |>
    mutate(peak_y = max(intensity, na.rm = TRUE),
           peak_start = ifelse(is.na(peak_start), min(rt, na.rm = TRUE), peak_start),
           peak_end = ifelse(is.na(peak_end), min(rt, na.rm = TRUE), peak_end)) |>
    ungroup() |>
    as.data.frame()
  
  createInputDfForInsetPlotting(tmp)
  
}

# Create inset for iRT plots
getInset <- function(df){
  p <- ggplot(data = df |> 
                group_by(QC),
              aes(x = rt, y = intensity)) +
    geom_line() +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
    geom_vline(aes(xintercept = rt_at_max), color = 'red') +
    theme_bw() +  ## makes everything smaller
    expand_limits(y = 0) +
    theme(panel.background = element_rect(fill="white"),  ## white plot background 
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(size=rel(0.7)), ## tiny axis text
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank())
  return(p)
}

createQCsubplot <- function(df, mark.peak = TRUE) {
  p <- df |>
    ggplot(aes(x = rt, y = intensity, color = QC)) +
    geom_plot_npc(data = inset_coords, aes(npcx = x, npcy = y, label = plot, vp.width = 0.25, vp.height = 0.9)) +
    geom_line(na.rm = TRUE) +
    facet_grid(factor(QC)~sample, scales = 'free', space = 'free_x') +
    theme_bw() +
    theme(legend.position = 'none', 
          strip.text.y = element_blank()) +
    coord_cartesian(xlim = c(0, max(temp_QC$rt)), expand = FALSE) +
    geom_text(data = temp_QC_stats, aes(x = max(temp_QC$rt)*0.01, y = y_label, label = paste('peptide ID:', QC)),
              size = 3, hjust = 0, vjust = 1.5, color = 'black') +
    labs(x = 'RT', y = 'Intensity')
  
  if(mark.peak){
    p <- p +
      geom_point(data = temp_QC_inset, aes(x = QC_peak_x, y = max_intensity * 1.1), color = 'black', shape = "\u2193", size = 3) +
      geom_text(data = temp_QC_stats, aes(x = max(temp_QC$rt)*0.75, y = y_label, 
                                          label = label), 
                size = 3, lineheight = 0.75, color = 'black', hjust = 1, vjust = 1.1)
  }
  
  if(i == 1) QC_subplot <- QC_subplot + labs(title = 'QC peptides', subtitle = fn)
  if(i == 2) QC_subplot <- QC_subplot + labs(title = ' ', subtitle = ' ')
  
  # Create QC plots for individual export
  QC_subplot_separate <- QC_subplot + labs(title = 'QC peptides', subtitle = temp_sample_name) + theme(strip.text.x = element_blank())
  
  say(temp_sample_name, type = 'output')
}