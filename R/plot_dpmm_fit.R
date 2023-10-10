#' Plot Parameter Trace plots
#'
#' This function plots trace plots of the DPMM estimated parameters.
#'
#' @param x object class 'dpmm_fit'
#' @param trace Logical parameter for plotting trace plots. (default = TRUE)
#' @param density Logical parameter for plotting density plots. (default = TRUE)
#' @param ... other arguments used by base plot.
#'
#' @return Plot
#'
#' @import coda
#'
#' @examples
#' \dontrun{
#' ## load dataset
#' data(dataset_1)
#' 
#' ## fit model
#' posteriors <- runModel(dataset_1, 
#'                        mcmc_iterations = 100,
#'                        L = 6, 
#'                        mcmc_chains = 2, 
#'                        standardise = TRUE)
#' ## plot trace plots
#' plot(posteriors)
#' }
#'
#'
#' @export
plot_dpmm_fit <- function(x, trace = TRUE, density = TRUE, ...) {
  
  ## check inputs
  if(!inherits(x, "dpmm_fit")) stop("'x' must be a dpmm_fit object")
  if(!is.logical(trace) | !is.logical(density)) stop("'trace' and 'density' must be logical")
  if(length(trace) > 1 | length(density) > 1) {
    print("'trace' or 'density' have more than one element, only the first will be used")
    trace <- trace[1]
    density <- density[1]
  }
  if(is.na(trace) | is.na(density)) stop("'trace' or density' must be specified")
  
  ## convert samples to mcmc object
  plot(as.mcmc(as.matrix(x$samples)), trace = trace, density = density)
}

