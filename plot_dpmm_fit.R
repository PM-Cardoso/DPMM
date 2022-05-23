library(tidyverse)
library(coda)

## plot method for dpmm_fit objects
plot.dpmm_fit <- function(x, trace = TRUE, density = TRUE, ...) {
  
  ## check inputs
  if(class(x) != "dpmm_fit") stop("'x' must be a dpmm_fit object")
  if(!is.logical(trace) | !is.logical(density)) stop("'trace' and 'density' must be logical")
  if(length(trace) > 1 | length(density) > 1) {
    print("'trace' or 'density' have more than one element, only the first will be used")
    trace <- trace[1]
    density <- density[1]
  }
  if(is.na(trace) | is.na(density)) stop("'trace' or density' can't be missing values")
  
  ## convert samples to mcmc object
  plot(as.mcmc(as.matrix(x$samples)), trace = trace, density = density)
}

