library(tidyverse)
library(coda)
library(GGally)

plot.ggpairs <- function(x, newdata, iterations, nburn, ggpairs_title = "",...) {

  if (missing(iterations) && missing(nburn)) stop("'nburn' or 'iterations' must be provided")
  if (!missing(iterations) && !missing(nburn)) stop("'nburn' or 'iterations' cannot be used at same time")
  if (!missing(iterations)) {
    if (!is.numeric(iterations)) stop("'iterations' must be 'numeric'")
  }
  if (!missing(nburn)) {
      if (!is.numeric(nburn)) stop("'nburn' must be 'numeric'")
  }
  if (!missing(newdata) && !missing(iterations)) stop("'newdata' can only be used with 'nburn'")
  if (!missing(newdata) && missing(nburn)) stop("'nburn' is required when supplying 'newdata'")
  if (class(x) != "dpmm_fit") stop("'x' must be class 'dpmm_fit'")
  
  if (missing(x))  {
    stop("'x' needs to be supplied")
  } else {
      
    discrete <- dplyr::select(x$dataset, where(Negate(is.numeric)))
    
    if(!missing(newdata)) {
      ## check newdata
      if(!is.data.frame(newdata)) stop("'newdata' must be a data.frame")
      if(!all(colnames(newdata) %in% colnames(x$dataset))) stop("Incorrect column names for 'x'")
      
      ## check all variables are present
      if(ncol(newdata) != ncol(x$dataset)) stop("'newdata' columns do not mathc 'x$dataset'")
      
      ## reorder newdata to match original data set
      newdata <- dplyr::select(newdata, !!colnames(x$dataset))
      
      ## check variables match
      if(!identical(
        summarise(newdata, across(everything(), class)),
        summarise(x$dataset, across(everything(), class))
      )) stop("'newdata' and 'x$dataset' have different column classes")
      
      if (x$mcmc_chains == 1) {
        number_samples <- nrow(x$samples)
      } else {
        number_samples <- nrow(x$samples[[1]])
      }
      
      if ((nrow(newdata) > length((nburn+1):number_samples))) stop("not enough samples to compare against the dataset")
      
      iterations <- seq(nburn, number_samples, 1)
      
      if (x$mcmc_chains == 1) {
        samples <- x$samples[ceiling(iterations),]
      } else {
        samples <- x$samples
        for (i in 1:x$mcmc_chains) {
          samples[[i]] <- samples[[i]][ceiling(iterations),]
        }
      }
      # class(samples) <- c("tbl_df", "tbl", "data.frame")
      
      object <- list(dataset = x$dataset, L = x$L, samples = samples, standardise = x$standardise, mcmc_chains = x$mcmc_chains)
      
      
      
    } else {
      
      #in case only one posterior samples is given
      if (x$mcmc_chains == 1) {
        number_samples <- nrow(x$samples)
      } else {
        number_samples <- nrow(x$samples[[1]])
      }
      
      if (x$mcmc_chains == 1) {
        if (min(ceiling(iterations)) < 1 || max(ceiling(iterations)) > max(nrow(x$samples))) stop("number of iterations does not match mcmc samples")
      } else {
        if (min(ceiling(iterations)) < 1 || max(ceiling(iterations)) > max(nrow(x$samples[[1]]))) stop("number of iterations does not match mcmc samples")
      }
      
      if (x$mcmc_chains == 1) {
        samples <- x$samples[ceiling(iterations),]
      } else {
        samples <- x$samples
        for (i in 1:x$mcmc_chains) {
          samples$samples[[i]] <- samples$samples[[i]][ceiling(iterations),]
        }
      }
      # class(samples) <- c("tbl_df", "tbl", "data.frame")
      
      object <- list(dataset = x$dataset, L = x$L, samples = samples, standardise = x$standardise, mcmc_chains = x$mcmc_chains)
       
    }
      
    
    
    class(object) <- "ggpairs.fit"
    
    
    
    source("predict_dpmm_fit.R")
    
    if (!missing(nburn)) {
      if (x$mcmc_chains == 1) {
        posteriors <- predict.dpmm_fit(object, samples = seq(1, nrow(object$samples), 1))
      } else {
        posteriors <- predict.dpmm_fit(object, samples = seq(1, nrow(object$samples[[1]]), 1))
      }
    } else {
      posteriors <-  predict.dpmm_fit(object, samples = seq(1, length(iterations), 1))
    }
    
    
    if (!missing(newdata)) {
      
      # sample the same number of values as the dataset we put into the function
      posteriors[[1]] <- posteriors[[1]] %>% sample_n(nrow(newdata))
      
      if (ncol(discrete) != 0) {
        for (i in 1:ncol(discrete)) {
          posteriors[[1]][,colnames(discrete)[i]] <- factor(posteriors[[1]][,colnames(discrete)[i]] , 
                                                            levels = 1:length(unique(posteriors[[1]][,colnames(discrete)[i]])),
                                                            labels = levels(discrete[,i]))
        }
      }
      
      newdata_plot <- cbind(newdata, Data = "Original")
      
      posteriors_plot <- cbind(posteriors[[1]], Data = "DPMM")
      
      #if values were standardised, bring them back to normal
      if (class(x) == "dpmm_fit") {
        if (x$standardise == TRUE) {
          continuous <- dplyr::select(x$dataset, where(is.numeric))
          for (i in 1:ncol(continuous)) {
            posteriors_plot[,colnames(continuous)[i]] <- (posteriors_plot[,colnames(continuous)[i]]*x$sd_values[i])+x$mean_values[i]
          }
        }
      } else {
        if (x[[1]]$standardise == TRUE) {
          continuous <- dplyr::select(x[[1]]$dataset, where(is.numeric))
          for (i in 1:ncol(continuous)) {
            posteriors_plot[,colnames(continuous)[i]] <- (posteriors_plot[,colnames(continuous)[i]]*x[[1]]$sd_values[i])+x[[1]]$mean_values[i]
          }
        }
      }
      
      
      
      combination <- rbind(newdata_plot, posteriors_plot) %>%
        mutate(Data = factor(Data))
      
      
      my_dens_lower <- function(data, mapping, ...) {
        ggplot(mapping=mapping) +
          geom_density2d(data = filter(data, Data == "Original"), linewidth = 1.2, alpha = 0.7, colour = "#00BFC4") +
          geom_density2d(data = filter(data, Data == "DPMM"), linewidth = 1.2, alpha = 0.7, colour = '#F8766D')
      }
      
      my_dens_diagonal <- function(data, mapping, ...) {
        ggplot(data = data, mapping=mapping) +
          geom_density(aes(fill = Data), alpha = 0.3)
      }
      
      
      ggpairs(combination, columns = 1:(ncol(combination)-1), 
              aes(color = Data),
              showStrips = TRUE,
              lower = list(continuous = my_dens_lower, discrete = wrap(ggally_facetbar, position = "dodge"), combo = wrap(ggally_facetdensity,alpha=0.7)),
              diag = list(continuous = my_dens_diagonal, discrete = wrap(ggally_barDiag, position = "dodge")),
              upper = NULL,
              legend = 1,
              title = ggpairs_title) +
        theme(legend.position = 'bottom',
              panel.border = element_rect(fill = NA),
              panel.grid.major = element_blank())
      
    } else {
      
      if (ncol(discrete) != 0) {
        for (i in 1:ncol(discrete)) {
          posteriors[[1]][,colnames(discrete)[i]] <- factor(posteriors[[1]][,colnames(discrete)[i]] , 
                                                            levels = 1:length(unique(posteriors[[1]][,colnames(discrete)[i]])),
                                                            labels = levels(discrete[,i]))
        }
      }
      
      
      
      posteriors_plot <- cbind(posteriors[[1]], Data = "DPMM")
      
      combination <- posteriors_plot %>%
        mutate(Data = factor(Data))
      
      # combination %>% ggpairs()
      
      
      my_dens_lower <- function(data, mapping, ...) {
        ggplot(mapping=mapping) +
          geom_density2d(data = filter(data, Data == "DPMM"), linewidth = 1.2, alpha = 0.7, colour = '#F8766D')
      }
      
      my_dens_diagonal <- function(data, mapping, ...) {
        ggplot(data = data, mapping=mapping) +
          geom_density(aes(fill = Data), alpha = 0.3)
      }
      
      
      ggpairs(combination, columns = 1:(ncol(combination)-1),
              aes(color = Data),
              showStrips = TRUE,
              lower = list(continuous = my_dens_lower, discrete = wrap(ggally_facetbar, position = "dodge"), combo = wrap(ggally_facetdensity,alpha=0.7)),
              diag = list(continuous = my_dens_diagonal, discrete = wrap(ggally_barDiag, position = "dodge")),
              upper = NULL,
              title = ggpairs_title) +
        theme(panel.border = element_rect(fill = NA),
              panel.grid.major = element_blank())
      
      
    }
    
    
  }
  
  
}
