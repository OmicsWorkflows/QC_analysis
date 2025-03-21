# ***********************************************
# Title       : HTML output template
# Author      : Karolina Krystofova
# Date        : 2024/03/22
# Version     : 1.0.0
# ***********************************************

css <- 'body { 
  font-family: Open Sans, Arial, Helvetica, sans-serif;
  font-size: 14px;
}
table {
  table-layout: auto;
  width:850px;
  border: 0px;
}
.first {
  width: 75%
}
th {
  border-bottom: 1px solid #777777;
  background-color: #DCDCDC;
}
td, th {
  padding-right: 2px;
  padding-left: 2px;
  text-align: center;
}
td {
  white-space: pre;
}
.main{
  margin-top: 30px;
  margin-bottom: 100px;
  margin-right: 100px;
  margin-left: 300px;
  padding: 0px 10px;
}
.tab {
  display: inline-block;
  margin-left: 20px;
}
p {
  text-align: justify;
  text-justify: inter-word;
  max-width: 800px;
}
.sidenav {
  height: 100%;
  width: 250px;
  position: fixed; 
  z-index: 1;
  top: 0;
  left: 0;
  background-color: #DCDCDC;
    overflow-x: hidden;
  padding-top: 30px;
}
.sidenav a {
  text-decoration: none;
  color: black;
  display: block;
  padding: 3px 10px 3px 15px;
}
.sidenav a:hover {
  color: #B8B8B8;
}'

# Function for generating navigation links to plots for a given type
# tbl: image table
# x: type
navByLevel <- function(tbl, x, m = '') {
  output <- ''
  tmp <- tbl |> dplyr::filter(tbl$type == x) |>
    unnest(name) |>
    unnest(link)
  
  if(nrow(tmp) > 1) { # are there more plots of the same type?
    output <- paste0(output, '<li><a href = "#', 
                     tolower(gsub(' ', '_', x, fixed=TRUE)))
    if(m != '') output <- paste0(output, '_', tolower(gsub(' ', '_', m, fixed=TRUE)))
    output <- paste0(output, '">', x, '</a></li>')
    #tmp$link <- tolower(paste(gsub(' ', '_', tmp$name), tmp$link, sep = '_'))
    output <- paste0(output, '<ul>\n')
    
    for(n in tmp$name) {
      output <- paste0(output, '<li><a href = "#', tmp$link[tmp$name == n], '">',
                       n, '</a></li>\n')
    }
    
    output <- paste0(output, '</ul>\n')
  } else {
    output <- paste0(output, '<li><a href = "#', 
                     tolower(gsub(' ', '_', x, fixed=TRUE)))
    if(m != '') output <- paste0(output, '_', tolower(gsub(' ', '_', m, fixed=TRUE)))
    output <- paste0(output, '">', x, '</a></li>')
  }
  
  return(output)
}

# Function for generating part of the report body for a given type
# tbl: image table
# x: type
# m: MS level, for determining header level
bodyByLevel <- function(tbl, x, m) {
  output <- ''
  tmp <- tbl |> dplyr::filter(tbl$type == x)
  
  if(m == '' | m == 'Calibration curve') {
    output <- paste0(output, '<a id ="', tolower(gsub(' ', '_', x, fixed=TRUE)), '"><h2>', x, '</h2>\n')
  } else if(m != 'Individual samples') {
    output <- paste0(output, '<a id ="', tolower(gsub(' ', '_', x, fixed=TRUE)), 
                     '_', tolower(gsub(' ', '_', m, fixed=TRUE)), '"><h3>', x, '</h3>\n')
  } else {
    output <- paste0(output, '<a id ="', tolower(gsub(' ', '_', x, fixed=TRUE)), '"><h3>', x, '</h3>\n')
  }
  
  if(nrow(tmp) > 1) {
    tmp <- tmp |>
      unnest(name) |>
      unnest(link)
    #tmp$link <- tolower(paste(gsub(' ', '_', tmp$name), tmp$level, sep = '_'))
    
    for(n in tmp$name) {
      img_base64 <- base64Encode(readBin(tmp$full.path[tmp$name == n], "raw", 
                                         file.info(tmp$full.path[tmp$name == n])[1, "size"]), "txt")
      img_str <- paste0('<a id = "', tmp$link[tmp$name == n],'"><img src="data:image/png;base64,', 
                        img_base64, '" width="', tmp$width[tmp$name == n],
                        'px" height="', tmp$height[tmp$name == n],'px"></img><br>\n')
      output <- paste0(output, img_str)
    }
    
    output <- paste0(output, '</ul>\n')
    
  } else {
    img_base64 <- base64Encode(readBin(tmp$full.path, "raw", 
                                       file.info(tmp$full.path)[1, "size"]), "txt")
    img_str <- paste0('<img src="data:image/png;base64,', img_base64, '" width="', tmp$width, 'px" height="',
                      tmp$height,'px"></img><br>\n')
    output <- paste0(output, img_str)
  }
  
  return(output)
}

# Generatate html code for the sidebar navigation which includes links to images
# tbl: image table
# msn: iterates over all found MS levels
htmlNav <- function(tbl, msn){
  tmp <- tbl |> dplyr::filter(level == msn)
  
  output_str <- ''
  
  if(msn == 'Individual samples') { # sample graphs
    tmp$type <- sapply(1:nrow(tmp), function(i) unlist(tmp$name[i])[1])
    tmp <- unnest(tmp, name)
    output_str <- paste0(output_str, '<li><a href = "#', 
                         tolower(gsub(' ', '_', msn, fixed = TRUE)),
                         '">', msn,'</a></li>\n<ul>')
    tmp$link <- tolower(gsub(' ', '_', tmp$type, fixed = TRUE))
    output_str <- paste0(output_str, paste(unname(sapply(unique(tmp$name), 
                        function(x) paste0('<li><a href = "#',
                                           tmp$link[tmp$name == x], '">', x,'</a></li>'))), 
                        collapse = '\n'), '</ul>\n')
  } else if(all(tmp$type == tmp$level) | all(tmp$level == '')) { # calibration graphs etc.
    tmp <- tmp |>
      unnest(name) |>
      unnest(link)
    #tmp$link <- tolower(gsub(' ', '_', tmp$type, fixed = TRUE))
    output_str <- paste(output_str, paste(unname(sapply(unique(tmp$type), function(x) navByLevel(tmp, x))), collapse = '\n'))
  } else { # standard MS and MS/MS graphs
    tmp <- tmp |>
      unnest(name) |>
      unnest(link)
    #temp$link <- tolower(paste(gsub(' ', '_', temp$type, fixed = TRUE), temp$level, sep = '_'))
    output_str <- paste0(output_str, '<li><a href = "#', tolower(msn), '">', msn,'</a></li>\n<ul>', 
                         paste(unname(sapply(unique(tmp$type), function(x) navByLevel(tmp, x, msn))), collapse = '\n'),
                         '</ul>\n')
  }
  
  return(output_str)
}

htmlBody <- function(tbl, msn) {
  tmp <- tbl |> dplyr::filter(level == msn)
  
  output_str <- ''
  
  if (msn == 'Individual samples') {
    tmp$type <- sapply(1:nrow(tmp), function(i) unlist(tmp$name[i])[1])
    tmp$link <- tolower(gsub(' ', '_', tmp$type, fixed = TRUE))
    output_str <- paste0(output_str, '<a id ="', tolower(gsub(' ', '_', msn, fixed = TRUE)), '"><h2>', msn, '</h2>\n')
  } else if(all(tmp$type == tmp$level) | all(tmp$level == '')){
    #tmp$link <- tolower(gsub(' ', '_', tmp$type, fixed = TRUE))
  } else {
    #tmp$link <- tolower(paste(gsub(' ', '_', tmp$type, fixed = TRUE), tmp$level, sep = '_'))
    output_str <- paste0(output_str, '<a id ="', tolower(msn), '"><h2>', toupper(msn), ' level</h2>\n')
  }
  
  output_str <- paste(output_str, paste(unname(sapply(unique(tmp$type), function(x) bodyByLevel(tmp, x, msn))), collapse = '\n'))
  if (msn == 'Individual samples') output_str <- gsub('<h3>[a-zA-Z0-9_-]+</h3>', '', output_str)
  
  return(output_str)
}

# Strings from output table
outputTableNav <- function(table, cal.bool) {
  output_str <- '<ul>'
  
  if(cal.bool) output_str <- paste0(output_str, '<li><a href = "#injection_amount">Injection amount</a></li>\n')
  
  output_str <- paste0(output_str, '<li><a href = "#metrics_table">QC metrics</a></li>\n<ul>\n')
  output_str <- paste0(output_str, paste(sapply(unique(table$Metric), function(x) {
      paste0('<li><a href = "#', tolower(gsub(' ', '_', x)), '">', x, '</a></li>')
    }, USE.NAMES = FALSE), collapse = '\n'))
  
  output_str <- paste0(output_str, '</ul></ul>\n')
  
  return(output_str)
}

outputTableBody <- function(export.tables, cal.bool) {
  output_str <- ''
  
  if(cal.bool) output_str <- paste0(output_str, '<a id ="injection_amount"><h3>Injection amount</h3>\n',
                                    export.tables$inject_tbl_exp)
  
  output_str <- paste0(output_str, '<a id ="metrics_table"><h3>QC metrics</h3>\n')
  
  output_str <- paste0(output_str, paste(sapply(unique(export.tables$sums_tbl$Metric), function(x) {
    tmp <- export.tables$sums_tbl |> 
      dplyr::filter(Metric == x) |>
      dplyr::select(-Metric) |>
      dplyr::select(where(~ !all(is.na(.))))
    
    tmp_str <- paste0('<a id ="', tolower(gsub(' ', '_', x)),'">', x, '\n',
           as.character(tableHTML(tmp, rownames = FALSE, border = 0)))
    tmp_str <- gsub('td id="tableHTML_column_1"', 
                    'td class = "first" id="tableHTML_column_1"', 
                    tmp_str, fixed = TRUE)
    tmp_str <- gsub('th id="tableHTML_header_1"', 
                    'th class = "first" id="tableHTML_header_1"', 
                    tmp_str, fixed = TRUE)
    return(tmp_str)
  }, USE.NAMES = FALSE), collapse = '<br>'))
  
  return(output_str)
}

# Template for HTML file
htmlTemplate <- function(css, image_nav, table_nav, qc_nav,
                         title, analysis_info, current_date, 
                         sample_table, images, output_table, qc_table){
  paste0('<!DOCTYPE html>
  <html>
  <head>
    </head>
    <body>
       <style>',
       css,
       '</style>',
       # Side navigation
       '<div class="sidenav">
          <ul style="list-style: none;padding: 0px;">
            <li><a href="#home">Home</a></li>
            <li><a href="#samples">Processed samples</a></li>',
            image_nav,
            '<li><a href="#output">Output table</a></li>',
            table_nav,
            qc_nav,
          '</ul>
       </div>',
       # Main body
       '<div class = "main">
         <a id="home">
         <h1>', title,'</h1>
         <p>',
            analysis_info,
            '<b>Date: </b>', current_date,
         '</p>
         <a id="samples">
         <h2>Processed samples</h2>',
         sample_table,
         images,
         '<a id="output">
         <h2>Output table</h2>',
         output_table,
         qc_table,
       '</div>
    </body>
  </html>')
}

# Wrapper function for creating HTML output
createHTMLreport <- function(header, title, img.table, 
                             export.tables, output.dir, 
                             qc, cal, error.list = NULL, type = 'MS') {
  # Reorder types of plots in image table
  img.table$type <- factor(img.table$type,
                           levels = c('Calibration curve',
                                      'QC peptides',
                                      'Individual TIC curves',
                                      'Individual BPC curves',
                                      'TIC curves overlay',
                                      'BPC curves overlay',
                                      'Additional traces',
                                      'QC metrics',
                                      'Calibration comparison',
                                      'Individual samples'))
  
  img.table <- arrange(img.table, type)


  ### Prepare strings which will be inserted into the template
  # Date and time now
  current_date <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  
  ## Side navigation
  # Links to images
  image_nav <- paste(unname(sapply(unique(img.table$level), 
                                   function(x) htmlNav(img.table, x))), collapse = '\n')
  
  table_nav <- if(type == 'MS') {
    outputTableNav(export.tables$sums_tbl, cal)
  } else if(type == 'calibration') {
    ''
  }
  qc_nav <- ifelse(qc, '<li><a href = "#qc_peptide_table">QC peptide table</a></li>\n', '')
  
  # Body of report
  image_body <- paste(unname(sapply(unique(img.table$level), function(x) htmlBody(img.table, x))), collapse = '\n')
  
  table_body <- if(type == 'MS') {
    outputTableBody(export.tables, cal)
  } else if (type == 'calibration') {
    export.tables$sums_tbl_exp
  }
  qc_table <- ifelse(qc, paste('<a id="qc_peptide_table">', '<h2>QC peptide table</h2>', export.tables$qc_tbl, sep = '\n'), '')
  
  # Alter sample table
  sample_table <- export.tables$samples_tbl_exp
  error_string <- ''
  
  for(n in names(error.list)) {
    errors <- error.list[[n]]
    
    if(length(errors) > 0) {
      if(n == 'duplicated') {
        error_string <- paste0(error_string, '\n<p><b>Following sample merged during plotting due to identical names after prefix and suffix removal:</b>\n', paste0('<li>', errors, collapse = '\n'), '</p>')
      } else {
        error_string <- paste0(error_string, '\n<p><b>Issues during obtaining ', n, ' traces for samples:</b>\n', paste0('<li>', errors, collapse = '\n'), '</p>')
      }
    } 
  }
  
  if(error_string != '') {
    sample_table <- paste0(sample_table, '<h3>Errors during processing</h3>\n', error_string)
    
  }
  
  # Insert strings into HTML template
  html <- htmlTemplate(css, image_nav, table_nav, qc_nav, title,
                       header, current_date, sample_table,
                       image_body, table_body, qc_table)
  
  # Export HTML file
  if(type == 'calibration') {
    fn <- paste0(output.dir, opts$tag, '.html')
  } else {
    fn <- ifelse(cal, paste(opts$file.batch, opts$calibration.tag, sep = '_'), opts$file.batch)
    fn <- removeDoubleUnderscore(fn)
    fn <- gsub('?', 'x', fn, fixed = TRUE)
    fn <- gsub('x', 'x', fn, fixed = TRUE)
    fn <- paste0(output.dir, fn, '.html')
  }
  
  cat(html, file = fn, sep = '\n')
  
  say('Report saved as:')
  say(fn, type = 'output')
}
