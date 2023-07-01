## load libraries
library(MASS)
library(tidyverse)
library(abind)
library(mvtnorm)



## prediction method for dpmm_fit object
predict.dpmm_fit <- function(object, newdata, samples = seq(1000,2500, 100), ...) {
  
  ## check object is correct class
  if(class(object) != "dpmm_fit" && class(object) != "ggpairs.fit") stop("'object' must be a 'dpmm_fit' or 'ggpairs.fit' object")
  
  ## if newdata is missing, then just generate random samples from posterior
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
    continuous <- dplyr::select(newdata, where(is.numeric))
    # need to standardise everything if object$standardise is present
    if (object$standardise == TRUE) {
      for (i in 1:(ncol(continuous))) {
        continuous[,i] <- (continuous[,i]-object$mean_values[i])/object$sd_values[i]
      }
    }
    discrete <- dplyr::select(newdata, where(Negate(is.numeric)))
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
    dataset_prediction <- cbind(continuous, discrete)
    
    ## is there any na? if not, stop
    if (min(rowSums(is.na(dataset_prediction))) == 0) {
      stop("'newdata' has entries with no missing values")
    }
    ## see if any of newdata has all data
    
    
  } else {
    ## extract continuous and discrete variables
    continuous <- dplyr::select(object$dataset, where(is.numeric))
    
    discrete <- dplyr::select(object$dataset, where(Negate(is.numeric)))
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
    dataset_prediction <- cbind(continuous, discrete)[1,]
    
  }
  
  
  source("posterior_dpmm.R")

  # model estimations
  samples_object <- object$samples
  if (object$mcmc_chains == 1) {
    samples_object <- samples_object[samples,]
  } else {
    for (mcmc_chains_iter in 1:length(samples_object)) {
      samples_object[[mcmc_chains_iter]] <- samples_object[[mcmc_chains_iter]][samples,]
    }
  }
  
  posterior <- posterior_dpmm(dataset_prediction, samples_object, seed = NULL, cont_vars = cont_vars, cat_vars = cat_vars, mcmc_chain = object$mcmc_chains)
  
  
  return(posterior)
  
}

