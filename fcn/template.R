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

# Strings for separate MSn levels
navByType <- function(tbl, x) {
  output <- ''
  tmp <- tbl |> dplyr::filter(tbl$type == x)
  
  output <- paste0(output, '<li><a href = "#', unique(tmp$link), '">', x, '</a></li>')
  
  if(nrow(tmp) > 1) {
    tmp$link <- tolower(paste(gsub(' ', '_', tmp$name), tmp$level, sep = '_'))
    output <- paste0(output, '<ul>\n')
    
    for(n in tmp$name) {
      output <- paste0(output, '<li><a href = "#', tmp$link[tmp$name == n], '">',
                       n, '</a></li>\n')
    }
    
    output <- paste0(output, '</ul>\n')
  }
  
  return(output)
}

htmlByType <- function(tbl, x, m) {
  output <- ''
  tmp <- tbl |> dplyr::filter(tbl$type == x)
  
  if(m == '') {
    output <- paste0(output, '<a id ="', unique(tmp$link), '"><h2>', x, '</h2>\n')
  } else {
    output <- paste0(output, '<a id ="', unique(tmp$link), '"><h3>', x, '</h3>\n')
  }
  
  if(nrow(tmp) > 1) {
    tmp$link <- tolower(paste(gsub(' ', '_', tmp$name), tmp$level, sep = '_'))
    
    for(n in tmp$name) {
      img_base64 <- base64Encode(readBin(tmp$address[tmp$name == n], "raw", 
                                         file.info(tmp$address[tmp$name == n])[1, "size"]), "txt")
      img_str <- paste0('<a id = "', tmp$link[tmp$name == n],'"><img src="data:image/png;base64,', 
                        img_base64, '" width="', tmp$width[tmp$name == n],
                        'px" height="', tmp$height[tmp$name == n],'px"></img><br>\n')
      output <- paste0(output, img_str)
    }
    
    output <- paste0(output, '</ul>\n')
    
  } else {
    img_base64 <- base64Encode(readBin(tmp$address, "raw", 
                                       file.info(tmp$address)[1, "size"]), "txt")
    img_str <- paste0('<img src="data:image/png;base64,', img_base64, '" width="', tmp$width, 'px" height="',
                      tmp$height,'px"></img><br>\n')
    output <- paste0(output, img_str)
  }
  
  return(output)
}

htmlNav <- function(tbl, msn){
  temp <- tbl |> dplyr::filter(level == msn)
  
  output_str <- ''
  
  if(msn == '') {
    temp$link <- tolower(gsub(' ', '_', temp$type, fixed = TRUE))
  } else {
    temp$link <- tolower(paste(gsub(' ', '_', temp$type, fixed = TRUE), temp$level, sep = '_'))
    output_str <- paste0(output_str, '<li><a href = "#', tolower(msn), '">', toupper(msn),'</a></li>\n<ul>')
  }
  
  output_str <- paste(output_str, paste(unname(sapply(unique(temp$type), function(x) navByType(temp, x))), collapse = '\n'))
  
  if(msn != '') {
    output_str <- paste0(output_str, '</ul>\n')
  }
  
  return(output_str)
}

htmlBody <- function(tbl, msn){
  temp <- tbl |> dplyr::filter(level == msn)
  
  output_str <- ''
  
  if(msn == ''){
    temp$link <- tolower(gsub(' ', '_', temp$type, fixed = TRUE))
  } else {
    temp$link <- tolower(paste(gsub(' ', '_', temp$type, fixed = TRUE), temp$level, sep = '_'))
    output_str <- paste0(output_str, '<a id ="', tolower(msn), '"><h2> ', toupper(msn), ' level</h2>\n')
  }
  
  output_str <- paste(output_str, paste(unname(sapply(unique(temp$type), function(x) htmlByType(temp, x, msn))), collapse = '\n'))
  
  
  return(output_str)
  
}

# Template for HTML file
htmlTemplate <- function(css, image_nav, qc_nav,
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
createHTMLreport <- function(batch = b, img.table = img_tbl, export.tables = export_tbls, output.dir = output_dir, qc = plot_qc_bool) {
  if(qc){
    img.table$type <- factor(img.table$type,
                             levels = c('Calibration curve',
                                        'QC peptides',
                                        'Individual TIC curves',
                                        'TIC curves overlay',
                                        'QC metrics'))
  } else {
    img.table$type <- factor(img.table$type,
                             levels = c('Calibration curve',
                                        'Individual TIC curves',
                                        'TIC curves overlay',
                                        'QC metrics'))
  }
  
  img.table <- arrange(img.table, type)
  img.table$level <- toupper(img.table$level)
  
  # Prepare strings which will be inserted into the template
  analysis_info <- paste('<b>Sample batch:</b>', batch, '\n<br>')
  current_date <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
  image_nav <- paste(unname(sapply(unique(img.table$level), function(x) htmlNav(img.table, x))), collapse = '\n')
  qc_nav <- ifelse(qc, '<li><a href = "#qc_peptide_table">QC peptide table</a></li>\n', '')
  image_body <- paste(unname(sapply(unique(img.table$level), function(x) htmlBody(img.table, x))), collapse = '\n')
  qc_table <- ifelse(qc, paste('<a id="qc_peptide_table">', '<h2>QC peptide table</h2>', export.tables$qc_tbl_exp, sep = '\n'), '')
  
  # Insert strings into HTML template
  html <- htmlTemplate(css, image_nav, qc_nav, 'QC sample analysis',
                        analysis_info, current_date, export.tables$samples_tbl_exp,
                        image_body, export.tables$sums_tbl_exp, qc_table)
  
  # Export HTML file
  fn <- paste0(output.dir, 'report.html')
  cat(html, file = fn, sep='\n')
  
  say('Report saved as:')
  say(fn, type = 'output')
}
