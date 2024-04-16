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
ticFluc <- function(x) {
  sapply(1:(length(x)-1), function(i) fluc(i, x))
}

# Wrapper for calculating AUC and inserting them into an export table
calculateAUC <- function(x, ms) {
  remove.bg <- tolower(ms) == 'ms'
  
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
    summarise(!!ms := ticAUC(RT, intensity)) |>
    mutate(Metric = 'AUC') |>
    rename(Sample = sample) |>
    as.data.frame() |>
    ungroup()
  
  if(remove.bg) {
    d[ms] <- round(d[ms] - bg)
  } else {
    d[ms] <- round(d[ms])
  }
  
  return(d)
}
