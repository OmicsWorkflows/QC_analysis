# Adjust dataframe
adjustDfCalibration <- function(d, config = opts) {
  d$sample <- sapply(d$sample, function(x) names(samples)[which(samples == x)], 
                     USE.NAMES = FALSE)
  d$real_amount <- sapply(d$sample, function(x) {
    config$real.amounts$amount[which(config$real.amounts$file == x)]}, 
    USE.NAMES = FALSE)
  d$real_amount <- as.numeric(d$real_amount)
  d <- relocate(d, 'time', .after = last_col())
  d$sample <- paste(as.integer(d$sample), 'ng')
  d$sample <- factor(d$sample, levels = paste(as.numeric(config$amounts), "ng"))
  
  return(d)
}

# Summarize calibration samples
summaryCalibration <- function(d, config = processing.config) {
  tmp <- d |>
    filter(RT > config$tic.start & RT < config$tic.end) |>
    calculateAUC(TRUE, NA) |>
    select(Sample, Value) |>
    rename(sample = Sample,
           AUC = Value)
  
  rv <- d |>
    select(sample, real_amount) |>
    distinct()
  
  tmp$max <- d |>
    filter(RT > config$tic.start & RT < config$tic.end) |>
    group_by(sample) |>
    summarise(max = max(intensity, na.rm = TRUE)) |>
    pull(max)
  
  tmp <- merge(tmp, rv) |>
    arrange(sample)
  
  return(tmp)
}

# Calculate residuals
getCalibrationCurve <- function(d, form) {
  # Calculate regression based on formula
  form.string <- paste0('real_amount~', gsub('x', 'AUC', form))
  form.string <- gsub(' ', '', form.string)
  form.string <- gsub("AUC\\^(\\d+)", "I(AUC^\\1)", form.string)
  fit <- lm(as.formula(form.string), data = d, weights = 1/(real_amount+0.01))
  
  cf <- coef(fit)
  names(cf) <- gsub("I\\(AUC\\^(\\d+)\\)", "AUC^\\1", names(cf))
  names(cf) <- gsub("AUC", "x", names(cf))
  
  # Get parts of a formula
  form.trim <- gsub(' ', '', form)
  form.parts <- unlist(strsplit(form.trim, '+', fixed = TRUE))

  # Create final formula string
  form.final <- 'y ='
  
  for (p in form.parts) {
    current.cf <- cf[p]
    if(which(form.parts == p) > 1) {
      if(current.cf > 0) form.final <- paste(form.final, '+')
    }
    
    form.final <- paste(form.final, current.cf, '*', p)
  }
  
  if(cf['(Intercept)'] > 0) form.final <- paste(form.final, '+')
  form.final <- paste(form.final, cf['(Intercept)'])
  form.final <- gsub(" \\-(\\d+)", " \\- \\1", form.final)
  
  # Calculate residuals
  resids <- resid(fit)
  names(resids) <- d$sample
  d$resid.abs <- sapply(d$sample, function(x) resids[as.character(x)], USE.NAMES = FALSE)
  d <- d |>
    mutate(expected = real_amount - resid.abs,
           resid = ifelse(!is.na(resid.abs), 
                          paste(round(resid.abs/expected, 3)*100, '%'),
                          NA))|>
    mutate(resid = ifelse(real_amount == 0, NA, resid)) |>
    select(-max, -resid.abs, -expected) |>
    relocate(real_amount, .after = 'sample') |>
    mutate(real_amount = paste(real_amount, 'ng'))
  
  return(list(fit = fit,
              equation = form.final,
              d = d))
}

# Function for plotting calibration curve
plotCalibrationCurve <- function(curve.info, subtitle, filename){
  line.d <- data.frame(AUC = seq(0, max(curve.info$d$AUC), length.out = 1000)) |>
    mutate(log = log(AUC))
  line.d$amount <- predict(curve.info$fit, line.d)
  line.d <- dplyr::filter(line.d, amount >= 0)
  
  curve.info$d$amount <- as.numeric(trimws(sub('ng', '', curve.info$d$real_amount, fixed = TRUE)))

  # Round numbers in the equation for a prettier plot
  rounded.equation <- curve.info$equation
  big.numbers <- unlist(regmatches(curve.info$equation, gregexpr("\\d+\\.?\\d*", rounded.equation)))
  big.numbers <- big.numbers[grep('.', big.numbers, fixed = TRUE)]
  big.numbers <- unique(big.numbers)
  rounded.num <- sapply(big.numbers, function(x) as.character(round(as.numeric(x), digits = 5)))
  
  for (i in seq_along(rounded.num)) {
    patt <- names(rounded.num)[i]
    repl <- rounded.num[i]
    rounded.equation <- gsub(patt, repl, rounded.equation, fixed = TRUE)
  }
  
  # Create plot
  p <- ggplot(data = curve.info$d, aes(x = AUC, y = amount, color = factor(sample))) +
    geom_point(size = 1.5) +
    geom_line(inherit.aes = FALSE, data = line.d, aes(x = AUC, y = amount), linetype = 2) +
    scale_color_manual(values = opts$plot.colors) +
    guides(color = guide_legend(byrow = TRUE, nrow = 1)) +
    labs(x = 'Area under curve', y = 'Amount [ng]', title = 'Calibration curve', subtitle = subtitle) +
    annotate(geom = 'text', y = max(curve.info$d$amount), x = 0, hjust = 0, label = rounded.equation) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    scale_y_continuous(expand = c(0.02, 0.02)) +
    theme_bw() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.key.size = unit(0.5, 'lines'),
          legend.background = element_rect(color = 'black'),
          legend.spacing.y = unit(0, "points"))
  
  tab <- tibble(name = 'Calibration curve',
                link = 'calibration_curve',
                plot = list(p),
                filename = filename,
                width = 500, 
                height = 600, 
                type = 'Calibration curve', 
                level = '',
                full.path = NA,
                include.in.report = TRUE)
  
  return(tab)
}

# Create output table
createCalibrationTable <- function(processing.info, curve.info){
  d <- data.frame(Tag = processing.info$tag,
                  Method = processing.info$method,
                  Equation = curve.info$equation,
                  check.names = FALSE)
  
  add.d <- curve.info$d |>
    select(-real_amount) |>
    mutate(AUC = as.character(AUC)) |>
    pivot_longer(!sample) |>
    na.omit() |>
    mutate(name = ifelse(name == 'resid', 'residual', 'integral'),
           sample = sub(' ', '', sample),
           colname = paste(sample, name))

  d <- cbind(d, add.d |>
               dplyr::filter(name == 'residual') |>
               select(-sample, -name) |>
               spread(key = 'colname', value = 'value') |>
               as.data.frame() |>
               select(add.d$colname[add.d$name == 'residual']))
  
  d <- cbind(d, add.d |>
               dplyr::filter(name == 'integral') |>
               select(-sample, -name) |>
               spread(key = 'colname', value = 'value') |>
               as.data.frame() |>
               select(add.d$colname[add.d$name == 'integral']))
  
  d$'TIC curve start' <- processing.info$tic.start
  d$'TIC curve end' <- processing.info$tic.end
  d$'Real amounts' <- paste(curve.info$d$real_amount, collapse = ', ')
  d$Time <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  d$'Used from' <- ''
  
  return(d)
}

# Create a new calibration file
newCalibrationTable <- function(d, path){
  if(file.exists(path)) {
    current.table <- read.table(path, header = TRUE, sep = '\t', 
                                quote = '', check.names = FALSE)
    
    # If new table is missing any columns, fill them in
    missing.cols <- setdiff(colnames(current.table), colnames(d))
    
    for(col in missing.cols) {
      preceding <- colnames(current.table)[which(colnames(current.table) == col)-1]
      i <- which(colnames(d) == preceding)

      join.list <- list()
      if(i != 0) join.list[[length(join.list)+1]] <- d[,1:i]
      join.list[[length(join.list)+1]] <- data.frame(key = col, value = NA) |> spread('key', 'value')
      if(i != ncol(d)) join.list[[length(join.list)+1]] <- d[,(i+1):ncol(d)]
      
      d <- do.call(cbind, join.list)
    }
    
    # Check if columns in new table correspond to columns in an existing table
    diffs <- setdiff(colnames(d), colnames(current.table))

    # If there are new columns, also add them in the table
    if(length(diffs) > 0) {
      for(col in diffs){
        i <- which(colnames(d) == col)-1
        
        join.list <- list()
        if(i != 0) join.list[[length(join.list)+1]] <- current.table[,1:i]
        join.list[[length(join.list)+1]] <- data.frame(key = col, 
                                                       value = rep(NA, times = nrow(current.table)),
                                                       id = 1:nrow(current.table)) |> 
          spread('key', 'value') |> select(-id)
        if(i != ncol(current.table)) join.list[[length(join.list)+1]] <- current.table[,(i+1):ncol(current.table)]
        
        current.table <- do.call(cbind, join.list)
      }
    }
    
    # If a row with a given tag is already present, replace it
    if(d$Tag %in% current.table$Tag) {
      current.table <- current.table[-which(current.table$Tag == d$Tag),] 
      current.table <- rbind(current.table, d[1,])
    } else {
      current.table <- rbind(current.table, d)
    }
    
    write.table(current.table, file = path, sep = '\t', 
                quote = FALSE, append = FALSE, col.names = TRUE, row.names = FALSE)
  } else {
    current.table <- d
    write.table(current.table, file = path, 
                sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE)
  }

  say(path, type = 'output')
  return(current.table)
}

# Create plot for the comparison of calibration samples
plotCalibrationSamplesComparison <- function(d, method, filename){
  # Create a table containing current calibration and at most 10 last used calibrations
  most.current <- d[nrow(d),]
  filtered.df <- d |>
    dplyr::filter(`Used from` != '',
                  Method == method)
  
  final.df <- rbind(filtered.df, most.current) 
  
  if(nrow(final.df) > 10) {
    final.df <- final.df[(nrow(final.df)-10):nrow(final.df),]
  }
  
  # Table to long
  final.df <- final.df |>
    select(Tag, ends_with('integral'), 'Real amounts') |>
    gather(-Tag, -'Real amounts', key = 'Sample', value = 'Integral') |>
    dplyr::filter(!is.na(Integral)) |>
    mutate(Sample = sub(' integral', '', Sample),
           `Real amounts` = gsub('"', '', `Real amounts`, fixed = TRUE),
           Amount = Sample)
  
  # If real amounts differ from expected ones, note this info in table
  replace.amounts <- final.df |>
    select(Tag, 'Real amounts') |>
    distinct() |>
    dplyr::filter(!is.na(`Real amounts`))
  
  if(nrow(replace.amounts) > 0) {
    lapply(1:nrow(replace.amounts), function(i){
      tmp <- replace.amounts[i,]
      tag <- tmp$Tag
      into.cols <- final.df |>
        dplyr::filter(Tag == tag) |>
        pull(Amount)
      
      tmp.long <- tmp |>
        separate(`Real amounts`, into = into.cols, sep = ',') |>
        gather(-Tag, key = 'Sample', value = 'Real amounts')
      
      final.df$Amount[final.df$Tag == tag] <<- tmp.long$`Real amounts`
    })
  }

  # Normalize AUC to real amount and further change df for plotting
  final.df <- final.df |>
    select(-`Real amounts`) |>
    mutate(Amount = as.numeric(gsub('([a-z,A-Z]+)', '', Amount)),
           Expected = as.numeric(gsub('([a-z,A-Z]+)', '', Sample)),
           Integral = as.numeric(Integral),
           'AUC (norm)' = ifelse(Sample != '0ng', (Integral/Amount)*Expected, Integral))
  final.df$Sample <- factor(final.df$Sample, levels = unique(final.df$Sample))
  
  # Create plot and save
  p <- ggplot(final.df, aes(x = Tag, y = `AUC (norm)`, color = Sample, group = Sample)) + 
    geom_line() +
    scale_color_manual(values = calibrationColors$colors[sub(' ', '', calibrationColors$amounts) %in% final.df$Sample]) +
    geom_point() +
    theme_bw() +
    labs(title = 'Comparison of samples over different calibrations') +
    facet_grid(Sample~., scales = 'free_y') +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  tab <- tibble(name = 'Calibration comparison',
                link = 'calibration_comparison',
                plot = list(p),
                filename = filename, 
                width = 500, 
                height = 600, 
                type = 'Calibration comparison', 
                level = '',
                full.path = NA,
                include.in.report = TRUE)
  
  return(tab)
}

# What was maximum injection amount sample
getMaxInjectionSample <- function(d, colname.substr = 'integral'){
  tmp <- d |>
    select(contains(colname.substr)) |>
    t() |>
    na.omit() |>
    rownames()
  vals <- as.numeric(gsub('([^0-9])', '', tmp))
  return(max(vals))
}

# Get calibration tag for current batch of samples
getCalibrationInfo <- function(names, d, opts) {
  calibration.tbl <- read.table(opts$calibration.file, header = TRUE, sep = '\t',
                                 strip.white = TRUE, fill = TRUE, check.names = FALSE)
  calibration.tbl$Time <- format(as.POSIXct(calibration.tbl$Time, format = '%Y-%m-%d %H:%M'),
                                 '%Y-%m-%d %H:%M:%S')
  config <- list()
  
  if(opts$calibration.tag.read == 'default') {
    calibration.tbl <- calibration.tbl  |>
      dplyr::filter(`Used from` != '')
    
    if(opts$calibration.method == 'default') {
      default.m <- NULL
      
      mtds.all <- unique(calibration.tbl$Method) |>
        na.omit()
      
      for (m in mtds.all){
        m.lookup <- paste0('_', m, '_')
        mtd.match <- sum(grepl(m.lookup, names, fixed = TRUE))
        if (mtd.match == length(names)) {
          default.m <- m
          break
        }
      }
      
      if(is.null(default.m)) {
        say('Error!')
        say('No appropriate method found in the calibration curve table! Specify method name in configuration file or check file names',
            type = 'error')
        quit()
      } else {
        opts$calibration.method <- default.m
        config$calibration.method <- default.m
      }
    }
    
    time.vector <- d[,c('sample', 'time')] |> distinct() |>
      pull(time)
  
    mtds.all <- unique(calibration.tbl$Method) |>
      na.omit()
    calibration.tbl <- calibration.tbl |>
      dplyr::filter(Method == opts$calibration.method) |>
      arrange(`Used from`)
    
    if(nrow(calibration.tbl) == 0){
      say('Error!')
      say(paste('No calibration curve with selected method found! Possible method options:',
                paste(mtds.all, collapse = ', ')), type = 'error')
      quit()
    }
    
    default.tag <- NULL
    
    for(i in 1:nrow(calibration.tbl)){
      sum.time <- sum(time.vector > calibration.tbl$`Used from`[i])
      
      if(sum.time == length(time.vector)) {
        default.tag <- calibration.tbl$Tag[i]
      } else if (sum.time == 0) {
        break
      } else {
        say('Error!')
        say('Multiple calibration curves assigned to file batch! 
             Edit batches or specify calibration tag in configuration file',type = 'error')
        quit()
      }
    }
    
    if(is.null(default.tag)) {
      say('Error!')
      say('No appropriate calibration tag found! Specify calibration tag in configuration file', type = 'error')
      quit()
    } else {
      calibration.tbl.info <- dplyr::filter(calibration.tbl, Tag == default.tag)
      config$calibration.tag <- default.tag
      config$calibration.equation <- calibration.tbl.info$Equation
      config$calibration.start <- as.numeric(calibration.tbl.info$`TIC curve start`)
      config$calibration.end <- as.numeric(calibration.tbl.info$`TIC curve end`)
      config$calibration.max <- getMaxInjectionSample(calibration.tbl.info)
      return(config)
    }
    
  } else if (opts$calibration.tag.read == '') {
    say('Error!')
    say('Specify calibration tag in configuration file or leave at "default"', type = 'error')
    quit()
  } else if (!opts$calibration.tag.read %in% calibration.tbl$Tag) {
    say('Error!')
    say('No calibration with selected tag found', type = 'error')
    quit()
  } else {
    config$calibration.tag <- opts$calibration.tag.read
    calibration.tbl.info <- dplyr::filter(calibration.tbl, Tag == opts$calibration.tag.read)
    config$calibration.equation <- calibration.tbl.info$Equation
    config$calibration.start <- as.numeric(calibration.tbl.info$`TIC curve start`)
    config$calibration.end <- as.numeric(calibration.tbl.info$`TIC curve end`)
    config$calibration.max <- getMaxInjectionSample(calibration.tbl.info)
    return(config)
  }
}

# Retrieve calibration curve plot
getCalibrationCurveInfo <- function(opts) {
  where <- list.files(paste0(opts$calibration.path, opts$calibration.tag, '_calibration', '/plots'), full.names = TRUE)
  img.path <- where[grepl('curve', where)]
  
  tab <- tibble(name = 'Calibration curve',
                link = 'calibration_curve',
                plot = NA,
                filename = NA, 
                width = 500, 
                height = 600, 
                type = 'Calibration curve', 
                level = 'Calibration curve',
                full.path = img.path,
                include.in.report = TRUE)
  return(tab)
}

# Create injection amount table
calculateInjection <- function(d, equation, max) {
  eq.final <- gsub(' ', '', equation)
  eq.final <- sub('y=', '', eq.final)
  eq.final <- gsub('x', 'AUC', eq.final)
  p <- parse(text = eq.final)
  
  d$Amount <- suppressWarnings(sapply(d$AUC, function(AUC) eval(p)))

  d <- d |>
    mutate(Amount = round(Amount, digits = 3),
           Note = case_when(Amount > max ~ '*',
                            is.na(Amount) ~ '&dagger;',
                            TRUE ~ '')) |>
    rename('Integral' = AUC)  |>
    mutate(`Calibration tag` = opts$calibration.tag,
           `Integration time` = format(Sys.time(), '%Y-%m-%d %H:%M:%S')) |>
    relocate(`Calibration tag`, .after = Sample) |>
    relocate(Note, .after = last_col())

  d_txt <- as.character(tableHTML(d, rownames = FALSE, border = 0))
  
  # Add warning about extrapolated values
  if(any(d$Note == '*')) {
    d_txt <- gsub('<td id="tableHTML_column_6">*</td>', '<td bgcolor="#f59105" id="tableHTML_column_6">*</td>', d_txt, fixed = TRUE)
    d_txt <- paste0(d_txt,
                    '\n<p>*<span class="tab"></span>extrapolated amount value, integral lies outside calibration curve boundary</p>\n')
  }
  
  # Add warning about missing values
  if(any(d$Note == '&dagger;')) {
    d_txt <- gsub('<td id="tableHTML_column_6">&dagger;</td>', '<td bgcolor="#eb3f1c" id="tableHTML_column_6">&dagger;</td>', d_txt, fixed = TRUE)
    d_txt <- paste0(d_txt,
                    '\n<p>&dagger;<span class="tab"></span>no amount value, integral lies outside permitted calibration curve value range</p>\n')
  }
  
  d_txt <- sub('</p>\n\n<p>', '\n<br>', d_txt, fixed = TRUE)
  
  return(list(inject_tbl_exp = d_txt,
              inject_tbl = d))
}

# Modify df
addCalibrationForPlotting <- function(output.df, d, opts) {
  d <- d |>
    dplyr::select(Sample, Amount)
  
  output.df$color = 'black'
  output.df$type = 'sample'
  
  cal.path <- paste0(opts$calibration.path, 'tic_data/calibration/', opts$calibration.tag, '_full_calibration.csv')
  
  if(!file.exists(cal.path)) {
    cal.path <- paste0(opts$calibration.path, 'tic_data/calibration/', opts$calibration.tag, '_full_calibration.tsv')
  }
  
  if(!file.exists(cal.path)) {
    say('Error!')
    say(paste('No TIC data for calibration tagged', opts$calibration.path, 'found'))
    quit()
  }
  
  cal.d <- read.table(cal.path, header = TRUE, sep = ',', strip.white = TRUE) |>
    dplyr::select(RT, intensity, sample, time)

  cal.amounts <- unique(as.numeric(gsub('[^[:digit:]]', '', cal.d$sample)))
  
  d[, c('Lower', 'Higher')] <- t(apply(d, 1, function(row){
      a <- as.numeric(row[2])
      if(is.na(a)) {
        return(c(NA, NA)) }
      else {
        if(a < min(cal.amounts)) {
          low <- min(cal.amounts)
          high <- sort(cal.amounts)[2]
        } else if(a > max(cal.amounts)) {
          low <- sort(cal.amounts, decreasing = TRUE)[2]
          high <- max(cal.amounts)
        } else {
          low <- max(cal.amounts[cal.amounts < a])
          high <- min(cal.amounts[cal.amounts > a])
        }
        return(c(low, high))
      }
    }))
  
  d <- d |>
    gather(-Sample, -Amount, key = 'Type', value = 'Cal.amount') |>
    dplyr::filter(!is.na(Cal.amount))

  out <- lapply(1:nrow(d), function(i){
    current.amount <- d[i,'Cal.amount']

    tmp.cal <- cal.d |>
      mutate(sample = gsub('[^[:digit:]]', '', sample)) |>
      dplyr::filter(sample == current.amount)

    tmp <- data.frame(RT = tmp.cal$RT,
                      intensity = tmp.cal$intensity)
    tmp$sample <- d[i,'Sample']
    tmp$time <- tmp.cal$time
    tmp$color = calibrationColors$colors[gsub('[^[:digit:]]', '', calibrationColors$amounts) == current.amount]
    tmp$type = paste(current.amount, 'ng')

    return(tmp)
  })
  
  out <- do.call(rbind, out)
  out.final <- rbind(output.df, out)

  fct.type.lvls <- c('sample', paste(sort(as.numeric(gsub('[^[:digit:]]', '', calibrationColors$amounts))), 'ng'))
  out.final$type <- factor(out.final$type, levels = fct.type.lvls)
  out.final <- arrange(out.final, type)

  return(out.final)
}
