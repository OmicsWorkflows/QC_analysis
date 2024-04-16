# Libraries to load
libs <- c('signal', 'MSnbase', 'mzR', 'ggplot2', 'dplyr', 'tidyr', 'tableHTML',
          'RCurl', 'purrr', 'ggpmisc', 'IDSL.IPA', 'xcms', 'rawrr')

# New processed files
new_tic <- 0
new_qc <- 0

# Expected injected amounts
expected_amounts <- c('000', '002', '005', '010', '020', '050', '100', '200')

# Dataframe containing designated colors for calibration lines
calibrationColors <- data.frame(amounts = c(0, 2, 5, 10, 20, 50, 100, 200),
                                colors = c('#ed1c24', '#fb7400', '#809100', '#008400', 
                                            '#00877d', '#0808e9', '#682e8e', '#665d99'),
                                overlayColor = c('#ed1c24', '#fb7400', '#809100', '#008400', 
                                                 '#00877d', '#0808e9', '#682e8e', '#E349A2'))

# Image table template
img_tbl_template <- data.frame(name = character(),
                               address = character(),
                               width = numeric(),
                               height = numeric(),
                               type = character(),
                               level = character())
