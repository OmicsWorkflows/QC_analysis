# Return current time
timeNow <- function() paste0('\r', format(Sys.time(), '[%H:%M:%S]'))
timeLog <- function() format(Sys.time(), '[%Y-%m-%d %H:%M:%S]')

# Log current message
say <- function(message, log = logfile, type = 'time', add = '\n\r') {
  m <- paste0(' ', message)
  
  if(type == 'time'){
    prefCons <- timeNow()
    prefLog <- timeLog()
  } else {
    prefCons <- prefLog <- paste0('[', toupper(type), ']')
  }
  
  cat(prefCons, m, add, sep = '')
  
  if(file.exists(log)) {
    cat('\n', prefLog, m, file = log, append = TRUE, sep = '')
  } else {
    cat(prefLog, m, file = log, append = TRUE, sep = '')
  }
}

# Returns error if block of code fails
e <- function(err, logfile = logfile) {
  say('An error occured', log = logfile)
  message <- paste(c(err, recursive = TRUE), collapse = '\n[ERROR] ')
  say(message, log = logfile)
}