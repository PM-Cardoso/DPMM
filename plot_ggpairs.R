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
  
  
  if (missing(x))  {
    stop("'x' needs to be supplied")
  } else {
    if (class(x) == "dpmm_fit") {
      
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
        
        
        number_samples <- nrow(x$samples)
        
        if ((nrow(newdata) > length((nburn+1):number_samples))) stop("not enough samples to compare against the dataset")
        
        
        samples <- x$samples[ceiling(seq(nburn+1, number_samples, length.out = nrow(newdata))),]
        class(samples) <- c("tbl_df", "tbl", "data.frame")
        
        object <- list(dataset = x$dataset, L = x$L, samples = samples, standardise = x$standardise)
        
        
        
      } else {
        
        #in case only one posterior samples is given
        number_samples <- nrow(x$samples)
        if (min(ceiling(iterations)) < 1 || max(ceiling(iterations)) > max(nrow(x$samples))) stop("number of iterations does not match mcmc samples")
        
        samples <- x$samples[ceiling(iterations),]
        class(samples) <- c("tbl_df", "tbl", "data.frame")
        
        object <- list(dataset = x$dataset, L = x$L, samples = samples, standardise = x$standardise)
         
      }
      
      
    } else {
      #in case there is more than one dataset given
      if (!is.list(x)) stop("'x' must be of class(dpmm_fit) or 'list'")
      samples <- NULL
      
      discrete <- dplyr::select(x[[1]]$dataset, where(Negate(is.numeric)))
      
      for (number_list in 1:length(x)) {
        
        
        ## if newdata is missing, then just generate random samples from posterior
        if(!missing(newdata)) {
          ## check newdata
          if(!is.data.frame(newdata)) stop("'newdata' must be a data.frame")
          if(!all(colnames(newdata) %in% colnames(x[[number_list]]$dataset))) {
            title <-paste0("Incorrect column names for 'x ", number_list,"'")
            stop(title)
          }
          
          ## check all variables are present
          if(ncol(newdata) != ncol(x[[number_list]]$dataset)) {
            title <- paste0("'newdata' columns do not mathc 'x[[",number_list,"]]$dataset'") 
            stop(title)
          }
          
          ## reorder newdata to match original data set
          newdata <- dplyr::select(newdata, !!colnames(x[[number_list]]$dataset))
          
          ## check variables match
          if(!identical(
            summarise(newdata, across(everything(), class)),
            summarise(x[[number_list]]$dataset, across(everything(), class))
          )) {
            title <- paste0("'newdata' and 'x[[",i,"]]$dataset' have different column classes")
            stop(title)
          }
          
          ## samples
          
          number_samples <- nrow(x[[number_list]]$samples)
          if (number_samples <= nburn) stop("number of samples must be bigger than nburn")
          samples <- rbind(samples, x[[number_list]]$samples[(nburn+1):number_samples,])
          
          
        } else {
          
          if (min(ceiling(iterations)) < 1 || max(ceiling(iterations)) > max(nrow(x[[number_list]]$samples))) stop("number of iterations does not match mcmc samples")
         
          samples <- rbind(samples, x[[number_list]]$samples[iterations,])
          
        }
        
        
        
        
      }
      
      if (!missing(nburn)) {
        if (nrow(newdata) > nrow(samples)) stop("not enough samples to compare against the dataset")
        
        samples <- samples[ceiling(seq(1, nrow(samples), length.out = nrow(newdata))),]
      }
      
      samples <- as.data.frame(samples)
      class(samples) <- c("tbl_df", "tbl", "data.frame")
      object <- list(dataset = x[[1]]$dataset, L = x[[1]]$L, samples = samples, standardise = x[[1]]$standardise)
      
      
      
    }
    
    
    class(object) <- "ggpairs.fit"
    
    
    
    source("predict_dpmm_fit.R")
    
    if (!missing(nburn)) {
      posteriors <- predict.dpmm_fit(object, samples = seq(1, nrow(object$samples), 1))
    } else {
      posteriors <-  predict.dpmm_fit(object, samples = seq(1, length(iterations), 1))
    }
    
    
    if (!missing(newdata)) {
      
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
      
      combination %>% ggpairs()
      
      
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
