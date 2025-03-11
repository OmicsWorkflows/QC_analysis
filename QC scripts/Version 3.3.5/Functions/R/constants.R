# Libraries to load
libs <- c('signal', 'MSnbase', 'mzR', 'ggplot2', 'dplyr', 'tidyr', 'tableHTML', 'grid',
          'RCurl', 'purrr', 'ggpmisc', 'IDSL.IPA', 'xcms', 'rawrr', 'tibble', 'gridExtra')

# New processed files
new_tic <- 0
new_qc <- 0

# Expected injected amounts
expected_amounts <- c('000', '002', '005', '010', '020', '050', '100', '200')

# Dataframe containing designated colors for calibration lines
calibrationColors <- data.frame(amounts = c('0 ng', '2 ng', '5 ng', '10 ng', '20 ng', '50 ng', '100 ng', '200 ng'),
                                colors = c('#ed1c24', '#fb7400', '#809100', '#008400', 
                                            '#00877d', '#0808e9', '#682e8e', '#665d99'),
                                overlayColor = c('#ed1c24', '#fb7400', '#809100', '#008400', 
                                                 '#00877d', '#0808e9', '#682e8e', '#E349A2'))

# Image table template
img_tbl_template <- tibble::tibble(name = character(),
                                   link = character(),
                                   plot = list(),
                                   filename = character(),
                                   width = numeric(),
                                   height = numeric(),
                                   type = character(),
                                   level = character(),
                                   full.path = character(),
                                   include.in.report = logical())

# Special characters for label splitting
spec_chars <- c('-',' ','_','/')

# Colors for curves
qcColors <- data.frame(QC = 1:11,
                       color = c('mediumvioletred', 'lightseagreen','lightslategrey','black',
                                  'firebrick', 'forestgreen', 'navy', 'gold3', 'mediumvioletred', 'lightseagreen','lightslategrey'))
lcColors <- data.frame(type = c('EIC', 'MPC', 'LC', 'LC'),
                       subtype = c(NA, NA, 'LP pump', 'NC pump'),
                       colors = c('darkorchid', 'darkorange3', 'navy', 'gold3'))
ticColors <- data.frame(type = c('TIC', 'TIC', 'BPC'),
                        level = c('MS', 'MS2', 'MS'),
                        colors = c('black', 'forestgreen', 'firebrick'))
