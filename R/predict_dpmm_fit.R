#' Predict from DPMM
#'
#' Predict missing values from the fitted DPMM.
#'
#' @param object object class 'dpmm_fit'
#' @param newdata dataframe with missingness
#' @param samples vector of iterations to be used in the posterior
#' @param ... other parameters used in 'posterior.dpmm'
#'
#' @return A list with n entries for n rows with missingness, each entry is a dataframe with the sampled missing values.
#'
#' @examples
#' ## load dataset
#' data(dataset_1)
#' 
#' ## fit model
#' posteriors <- runModel(dataset_1, 
#'                        mcmc_iterations = 100,
#'                        L = 6, 
#'                        mcmc_chains = 2, 
#'                        standardise = TRUE)
#'                        
#' ## introduce missing data
#' rows <- 501:550
#' dataset_missing <- dataset_1
#' dataset_missing_predict <- dataset_missing[rows,]
#' dataset_missing_predict[,1] <- as.numeric(NA)
#' 
#' # predict missing values
#' posteriors.dpmmfit <- predict(posteriors, 
#'                               dataset_missing_predict, 
#'                               samples = c(1:100))
#'                        
#'
#' @export
predict_dpmm_fit <- function(object, newdata, samples = seq(1000,2500, 100), ...) {
  
  ## check object is correct class
  if(!inherits(object, "dpmm_fit") && !inherits(object, "ggpairs.fit")) stop("'object' must be a 'dpmm_fit' or 'ggpairs.fit' object")
  
  ## if newdata is provided, check if the dataset is formatted the same was the dataset used during the DPMM fit.
  if(!missing(newdata)) {
    ## check newdata
    if(!is.data.frame(newdata)) stop("'newdata' must be a data.frame")
    if(!all(colnames(newdata) %in% colnames(object$dataset))) stop("Incorrect column names for 'newdata'")
    
    ## check all variables are present
    if(ncol(newdata) != ncol(object$dataset)) stop("'newdata' columns do not mathc 'object$dataset'")
    
    ## reorder newdata to match original data set
    newdata <- dplyr::select(newdata, !!colnames(object$dataset))
    
    ## check variables match
    if(!identical(
      summarise(newdata, across(everything(), class)),
      summarise(object$dataset, across(everything(), class))
    )) stop("'newdata' and 'object$dataset' have different column classes")
  }
  
  if (!missing(newdata)) {
    
    #:----------------------------------------------------------
    ## check which columns are continuous
    continuous <- dplyr::select(newdata, where(is.numeric))
    # need to standardise everything if object$standardise is present
    if (object$standardise == TRUE) {
      if (ncol(continuous) != 0) {
        for (i in 1:(ncol(continuous))) {
          continuous[,i] <- (continuous[,i]-object$mean_values[i])/object$sd_values[i]
        }
      }
    }
    
    #:----------------------------------------------------------
    ## check which columns are categorical
    discrete <- dplyr::select(newdata, where(Negate(is.numeric)))
    
    
    #:----------------------------------------------------------
    ## check the names of continuous and/or categorical variables
    if (ncol(continuous) != 0) {
      continuous <- continuous
      cont_vars <- colnames(continuous)
    } else {
      cont_vars <- NULL
    }
    if (ncol(discrete) != 0) {
      discrete <- discrete
      cat_vars <- colnames(discrete)
    } else {
      cat_vars <- NULL
    }
    
    # for the dataframe for prediction
    dataset_prediction <- cbind(continuous, discrete)
    
    ## is there any na? if not, stop
    if (min(rowSums(is.na(dataset_prediction))) == 0) {
      stop("'newdata' has rows with no missing values")
    }
    
    
  } else {
    #:---------------------------------------------------------------
    ## if newdata is missing, then just generate random samples from posterior
    
    #:----------------------------------------------------------
    ## check which columns are continuous
    continuous <- dplyr::select(object$dataset, where(is.numeric))
    
    #:----------------------------------------------------------
    ## check which columns are categorical
    discrete <- dplyr::select(object$dataset, where(Negate(is.numeric)))
    
    #:----------------------------------------------------------
    ## check the names of continuous and/or categorical variables
    if (ncol(continuous) != 0) {
      continuous[] <- as.numeric(NA)
      cont_vars <- colnames(continuous)
    } else {
      cont_vars <- NULL
    }
    if (ncol(discrete) != 0) {
      for (i in 1:ncol(discrete)) {
        discrete[,i] <- factor(NA, levels = levels(discrete[,i]))
      }
      cat_vars <- colnames(discrete)
    } else {
      cat_vars <- NULL
    }
    if (ncol(continuous) != 0) {
      continuous <- continuous
    }
    if (ncol(discrete) != 0) {
      discrete <- discrete
    }
    
    #:----------------------------------------------------------
    ## set up hypothetical row with all variables missing for random samples
    dataset_prediction <- cbind(continuous, discrete)[1,]
    
  }
  
  #:----------------------------------------------------------
  ## MCMC iterations to be used for sampling
  samples_object <- object$samples
  if (object$mcmc_chains == 1) {
    samples_object <- samples_object[samples,]
  } else {
    for (mcmc_chains_iter in 1:length(samples_object)) {
      samples_object[[mcmc_chains_iter]] <- samples_object[[mcmc_chains_iter]][samples,]
    }
  }
  
  #:----------------------------------------------------------
  ## sample posterior predictive distributions for missing values
  posterior <- posterior_dpmm(dataset_prediction, samples_object, seed = NULL, cont_vars = cont_vars, cat_vars = cat_vars, mcmc_chain = object$mcmc_chains)
  
  return(posterior)
  
}

