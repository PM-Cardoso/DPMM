


log_sum_exp <- function(lx) {
  
  ## extract maximum of logged values
  mX <- max(lx)
  
  ## return answer
  out <- mX + log(sum(exp(lx - mX)))
  out
}





posterior_dpmm <- function(patient, samples, seed = NULL, cont_vars = NULL, cat_vars = NULL) {
  
  #libraries
  library(MASS)
  library(tidyverse)
  library(abind)
  
  # seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (!is.null(cat_vars)) {
    ndisc = length(cat_vars)
    # ndiscdim = max(unique(as.numeric(patient[,cat_vars])))
    
    ndiscdim <- NULL
    for (i in 1:length(cat_vars)) {
      ndiscdim <- c(ndiscdim, length(unique(levels(patient[,cat_vars[i]]))))
    }
    ndiscdim <- max(ndiscdim)
  }
  
  
  # extract weights
  postW <- samples %>%
    select(starts_with("v[")) %>%
    apply(1, function(v) {
      w <- numeric(length(v) + 1)
      w[1] <- v[1]
      for(i in 2:length(v)) {
        w[i] <- v[i] * prod(1 - v[1:(i - 1)])
      }
      w[length(w)] <- prod(1 - v)
      w
    }) %>%
    t() %>%
    as_tibble() %>%
    mutate(iter = 1:n()) %>%
    gather(var, value, -iter) %>%
    mutate(var = as.numeric(gsub("V", "", var))) %>%
    arrange(iter, var) %>%
    group_by(iter) %>%
    nest() %>%
    mutate(data = map(data, "value")) %>%
    rename(w = data)
  
  # this allows for dpmm with continuous/categorical/continuous + categorical
  # so we check if vars names are null or not
  if (!is.null(cont_vars)) {
    postMu <- samples %>%
      select(starts_with("muL")) %>%
      mutate(iter = 1:n()) %>%
      gather(var, value, -iter) %>%
      mutate(var = gsub(" |muL\\[|\\]", "", var)) %>%
      separate(var, c("component", "dim"), sep = ",") %>%
      mutate_at(vars(c("component", "dim")), as.numeric) %>%
      arrange(iter, component, dim) %>%
      select(-dim) %>%
      group_by(iter, component) %>%
      nest() %>%
      mutate(data = map(data, "value")) %>%
      group_by(iter) %>%
      nest() %>%
      mutate(data = map(data, ~{
        pluck(., "data") %>%
          abind(along = 2) %>%
          t()
      })) %>%
      rename(muL = data)
    postTau <- samples %>%
      select(starts_with("tauL")) %>%
      mutate(iter = 1:n()) %>%
      gather(var, value, -iter) %>%
      mutate(var = gsub(" |tauL\\[|\\]", "", var)) %>%
      separate(var, c("dim1", "dim2", "component"), sep = ",") %>%
      mutate_at(vars(c("dim1", "dim2", "component")), as.numeric) %>%
      arrange(iter, component, dim2, dim1) %>%
      select(-dim1, -dim2) %>%
      group_by(iter, component) %>%
      nest() %>%
      mutate(data = map(data, "value")) %>%
      mutate(data = map(data, ~{
        matrix(., sqrt(length(.)), sqrt(length(.)))
      })) %>%
      group_by(iter) %>%
      nest() %>%
      mutate(data = map(data, ~{
        pluck(., "data") %>%
          abind(along = 3)
      })) %>%
      rename(tauL = data) 
  }
  if (!is.null(cat_vars)) {
    postPhi <- samples %>%
      select(starts_with("phiL")) %>%
      mutate(iter = 1:n()) %>%
      gather(var, value, -iter) %>%
      mutate(var = gsub(" |phiL\\[|\\]", "", var)) %>%
      separate(var, c("dim1", "dim2", "component"), sep = ",") %>%
      mutate_at(vars(c("dim1", "dim2","component")), as.numeric) %>%
      arrange(iter, component, dim1, dim2) %>%
      select(-dim1, -dim2) %>%
      group_by(iter, component) %>%
      nest() %>%
      mutate(data = map(data, "value")) %>%
      mutate(data = map(data, ~{
        matrix(., nrow = ndisc, ncol = ndiscdim, byrow = TRUE)
      })) %>%
      group_by(iter) %>%
      nest() %>%
      mutate(data = map(data, ~{
        pluck(., "data") %>%
          abind(along = 3)
      })) %>%
      rename(phiL = data)
  }
  
  
  # output of the function
  ordered_patient <- patient[,c(cont_vars,cat_vars)]
  # I decided to do the following
  #    input: several patients with several variables
  #    output: a list of patients, with matrix of missing vars
  
  # number of patients in data
  number_vars <- length(c(cont_vars))
  number_patients <- nrow(ordered_patient)
  number_samples <- nrow(samples)
  
  posterior <- vector(mode = "list", length = number_patients)
  
  for (i in 1:number_patients) {
    missing_vars <- ordered_patient[i,] %>%
      is.na() %>%
      colSums()
    
    
    missing_vars_names <- colnames(ordered_patient)[missing_vars > 0]
    
    posterior[[i]] <- matrix(as.numeric(NA), nrow = number_samples, ncol = length(missing_vars_names))
    
  }
  
  for (iteration_patient in 1:number_patients) {
    
    current_patient <- ordered_patient[iteration_patient,]
    
    if (!is.null(cat_vars)) {
      current_patient[cat_vars] <- as.numeric(current_patient[cat_vars])
    }
    
    #replicate patient into number of samples
    current_data <- as.data.frame(lapply(current_patient, rep, number_samples)) %>%
      mutate(iter = 1:n()) %>%
      gather(var, value, -iter) %>%
      arrange(iter) %>%
      group_by(iter) %>%
      nest() %>%
      mutate(data = map(data, "value"))  %>%
      rename(patient = data)
    
    # put together all posterior samples
    post <- current_data %>%
      inner_join(postW, by = "iter")
    if (!is.null(cont_vars)) {
      post <- post %>%
        inner_join(postMu, by = "iter") %>%
        inner_join(postTau, by = "iter")
    }
    if (!is.null(cat_vars)) {
      post <- post %>%
        inner_join(postPhi, by = "iter")
    }
    
    # draw conditional for only cont, only cat, cont and cat
    if (!is.null(cont_vars) & is.null(cat_vars)) {
      #   cont
      preds <- post %>%
        mutate(w = pmap(list(w, muL, tauL, patient), function(w, mu, tau, patient) {
          # vars without missing values
          variables <- which(!is.na(patient))
          
          if (!is_empty(variables)) {
            
            # marg matrices
            marg_mean <- mu[, variables, drop = FALSE]
            marg_sigma <- tau[variables, variables, , drop = FALSE]
            
            lmarg <- map_dbl(1:length(w), function(i, x, mu, sigma) {
              mvtnorm::dmvnorm(x, mu[i, ], solve(sigma[, , i]), log = TRUE) + log(w[i])
            }, x = patient[variables], mu = marg_mean, sigma = marg_sigma)
            lmarg <- log_sum_exp(lmarg)
            
            ## conditional for x2 for all components
            z_given_x_l <- map_dbl(1:length(w), function(i, x, mu, sigma, denom) {
              mvtnorm::dmvnorm(x, mu[i,], solve(sigma[, , i]), log = TRUE) + log(w[i]) - denom
            }, x = patient[variables], mu = marg_mean, sigma = marg_sigma, denom = lmarg)
            z_given_x <- exp(z_given_x_l)
            
            # draw cluster from adjusted probs
            draw <- which.max(rmultinom(n = 1, size = 1, prob = z_given_x))
            
          } else {
            
            draw <- which.max(rmultinom(n = 1, size = 1, prob = w))
            
          }
          
          })) %>%
        mutate(muL = map2(w, muL, function(w, mu) {
          mu[w, ]
        })) %>%
        mutate(tauL = map2(w, tauL, function(w, tau) {
          tau[, , w]
        })) %>%
        mutate(tauL = map(tauL, function(tau) {
          tauL <- matrix(unlist(tauL), ncol = number_vars)
          
          for (i in 1:sqrt(lengths(tauL))) {
            tauL[-i,i] = tauL[i,-i]
          }
          tauL
        })) %>%
        mutate(value = pmap(list(patient, muL, tauL), function(patient, mu, tau) {
          # vars without missing values
          variables_given <- which(!is.na(patient))
          # vars with missing values
          variables_dep <- which(is.na(patient))
          
          # conditional distribution for missing vars
          density <- condMVNorm::condMVN(mu, solve(tau), dep = variables_dep, given = variables_given,
                             X.given = patient[variables_given], check.sigma = FALSE)
          # draw for missing var
          preds <- cbind(t(MASS::mvrnorm(1, density$condMean, density$condVar))) %>%
            as.data.frame()
          preds
        }))
      
    } else if (is.null(cont_vars) & !is.null(cat_vars)) {
      #   cat
      preds <- post %>%
        mutate(w = pmap(list(w, phiL, patient), function(w, phi, patient) {
          variables <- which(!is.na(patient))
          
          if (!is_empty(variables)) {
            lmarg <- map_dbl(1:length(w), function(i, x, phi) {
              value = log(w[i])
              for (l in 1:length(variables)) {
                value <- value + log(phi[l ,x[variables[l]] , i])
              }
              value
            }, x = patient, phi = phi)
            lmarg <- log_sum_exp(lmarg)
            
            z_given_x_l <- map_dbl(1:length(w), function(i, x, phi, denom) {
              value <- log(w[i]) - denom
              for (l in 1:length(variables)) {
                value <- value + log(phi[l ,x[variables[l]] , i])
              }
              value
            }, x = patient, phi = phi, denom = lmarg)
            
            z_given_x <- exp(z_given_x_l)
            
            draw <- which.max(rmultinom(n = 1, size = 1, prob = z_given_x))
            
          } else {
            draw <- which.max(rmultinom(n = 1, size = 1, prob = w))
          }
        })) %>%
        mutate(phiL = map2(w, phiL, function(w, psi) {
          psi[,,w]
        })) %>%
        mutate(value = pmap(list(patient, phiL), function(patient, phi) {
          variables <- which(is.na(patient))
          
          preds <- NULL
          
          for (i in 1:length(variables)) {
            preds <- preds %>%
              cbind(rcat(1,phi[variables[i],]))
          }
          preds <- preds %>% as.data.frame()
          preds
        }))
        
      
    } else {
      #   cont and cat
      preds <- post %>%
        mutate(w = pmap(list(w, muL, tauL, phiL, patient), function(w, mu, tau, phi, patient) {
          # which cont vars to marginalize 
          continuous_patient <- patient[1:length(cont_vars)]
          continuous_vars_given <- which(!is.na(continuous_patient))
          continuous_vars_dep <- which(is.na(continuous_patient))
          # which cat vars to iterate over
          categorical_patient <- patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))]
          categorical_vars <- which(!is.na(categorical_patient))
          
          if (!is_empty(continuous_vars_given) & !is_empty(categorical_vars)) {
            # if there is continuous vars and categorical vars for w-conditional
            
            
            # marg matrices (with var that is missing)
            marg_mean <- mu[, continuous_vars_given, drop = FALSE]
            marg_sigma <- tau[continuous_vars_given, continuous_vars_given, , drop = FALSE]
            
            lmarg <- map_dbl(1:length(w), function(i, xcont, xcat, mu, sigma, phi) {
              # calculate the density for continuous
              value <- mvtnorm::dmvnorm(xcont, mu[i, ], solve(sigma[, , i]), log = TRUE) + log(w[i])
              # calculate the density for categorical
              for (l in 1:length(categorical_vars)) {
                value <- value + log(phi[l ,xcat[categorical_vars[l]] , i])
              }
              value
              
            }, xcont = patient[continuous_vars_given], xcat = patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))], mu = marg_mean, sigma = marg_sigma, phi = phi)
            lmarg <- log_sum_exp(lmarg)
            
            z_given_x_l <- map_dbl(1:length(w), function(i, xcont, xcat, mu, sigma, phi, denom) {
              # calculate the density for continuous
              value <- mvtnorm::dmvnorm(xcont, mu[i, ], solve(sigma[, , i]), log = TRUE) + log(w[i]) - denom
              # calculate the density for categorical
              for (l in 1:length(categorical_vars)) {
                value <- value + log(phi[l ,xcat[categorical_vars[l]] , i])
              }
              value
            }, xcont = patient[continuous_vars_given], xcat = patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))], mu = marg_mean, sigma = marg_sigma, phi = phi, denom = lmarg)
            
            z_given_x <- exp(z_given_x_l)
            
            # draw cluster from adjusted probs
            draw <- which.max(rmultinom(n = 1, size = 1, prob = z_given_x))
            
          } else if (!is_empty(continuous_vars_given) & is_empty(categorical_vars)) {
            # if there is continuous vars  for w-conditional
            
            # marg matrices (with var that is missing)
            marg_mean <- mu[, continuous_vars_given, drop = FALSE]
            marg_sigma <- tau[continuous_vars_given, continuous_vars_given, , drop = FALSE]
            
            lmarg <- map_dbl(1:length(w), function(i, xcont, mu, sigma) {
              # calculate the density for continuous
              value <- mvtnorm::dmvnorm(xcont, mu[i, ], solve(sigma[, , i]), log = TRUE) + log(w[i])
            }, xcont = patient[continuous_vars_given], mu = marg_mean, sigma = marg_sigma)
            lmarg <- log_sum_exp(lmarg)
            
            z_given_x_l <- map_dbl(1:length(w), function(i, xcont, mu, sigma, denom) {
              # calculate the density for continuous
              value <- mvtnorm::dmvnorm(xcont, mu[i, ], solve(sigma[, , i]), log = TRUE) + log(w[i]) - denom
            }, xcont = patient[continuous_vars_given], mu = marg_mean, sigma = marg_sigma, denom = lmarg)
            
            z_given_x <- exp(z_given_x_l)
            
            # draw cluster from adjusted probs
            draw <- which.max(rmultinom(n = 1, size = 1, prob = z_given_x))
            
            
          } else if (is_empty(continuous_vars_given) & !is_empty(categorical_vars)) {
            # if there is categorical vars  for w-conditional
            
            lmarg <- map_dbl(1:length(w), function(i, xcat, psi_var1) {
              # calculate the density for continuous
              value = log(w[i])
              # calculate the density for categorical
              for (l in 1:length(categorical_vars)) {
                value <- value + log(phi[l ,xcat[categorical_vars[l]] , i])
              }
              value
            }, xcat = patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))], phi = phi)
            lmarg <- log_sum_exp(lmarg)
            
            z_given_x_l <- map_dbl(1:length(w), function(i, xcat, psi_var1, denom) {
              # calculate the density for continuous
              value = log(w[i]) - denom
              # calculate the density for categorical
              for (l in 1:length(categorical_vars)) {
                value <- value + log(phi[l ,xcat[categorical_vars[l]] , i])
              }
              value
            }, xcat = patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))], phi = phi, denom = lmarg)
            
            z_given_x <- exp(z_given_x_l)
            
            # draw cluster from adjusted probs
            draw <- which.max(rmultinom(n = 1, size = 1, prob = z_given_x))
            
            
            
          } else {
            
            draw <- which.max(rmultinom(n = 1, size = 1, prob = w))
            
          }
            
          })) %>%
        mutate(muL = map2(w, muL, function(w, mu) {
          mu[w, ]
        })) %>%
        mutate(tauL = map2(w, tauL, function(w, tau) {
          tau[, , w]
        })) %>%
        mutate(tauL = map(tauL, function(tau) {
          tauL <- matrix(unlist(tauL), ncol = number_vars)
          
          for (i in 1:sqrt(lengths(tauL))) {
            tauL[-i,i] = tauL[i,-i]
          }
          tauL
        })) %>%
        mutate(phiL = map2(w, phiL, function(w, psi) {
          psi[,,w]
        })) %>%
        mutate(value = pmap(list(patient, muL, tauL, phiL), function(patient, mu, tau, phi) {
          # which cont vars to marginalize 
          continuous_patient <- patient[1:length(cont_vars)]
          continuous_vars_given <- which(!is.na(continuous_patient))
          continuous_vars_dep <- which(is.na(continuous_patient))
          # which cat vars to iterate over
          categorical_patient <- patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))]
          categorical_vars <- which(is.na(categorical_patient))
          
          preds <- NULL
          # if there is variables to be estimated
          if (!is_empty(continuous_vars_dep)) {
            # conditional distribution for missing vars
            density <- condMVNorm::condMVN(mu, solve(tau), dep = continuous_vars_dep, given = continuous_vars_given,
                                           X.given = patient[continuous_vars_given], check.sigma = FALSE)
            # draw for missing var
            preds <- cbind(t(MASS::mvrnorm(1, density$condMean, density$condVar)))
          }
          
          if (!is_empty(categorical_vars)) {
            for (i in 1:length(categorical_vars)) {
              if (length(categorical_vars) == 1) {
                preds <- preds %>%
                  cbind(rcat(1,phi))
              } else {
                preds <- preds %>%
                  cbind(rcat(1,phi[i,]))
              }
              
            }
          }
          
          preds <- as.data.frame(preds)
          preds
        }))
    }
      
      
    preds <- preds %>%
        select(iter, value) %>%
        unnest(cols = value) %>%
        ungroup() %>%
        select(-iter) %>%
        set_names(colnames(current_patient)[is.na(current_patient)])
        
    
    for (posterior_columns in 1:ncol(posterior[[iteration_patient]])) {
      
      posterior[[iteration_patient]][,posterior_columns] <- unlist(preds[,posterior_columns])
    
    }
    posterior[[iteration_patient]] <- posterior[[iteration_patient]] %>%
      as.data.frame() %>%
      set_names(colnames(current_patient)[which(is.na(current_patient))])
    
  }
  
  return(posterior)
  
}

  
