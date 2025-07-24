#' Sample from DPMM
#'
#' Sample missing values from the fitted DPMM
#'
#' @param patient dataset with missing values.
#' @param samples vector of iterations to be used in the posterior.
#' @param seed specify seed to be used. (default = NULL)
#' @param cont_vars names of continuous variables
#' @param cat_vars names of categorical variables
#' @param mcmc_chain MCMC posterior samples of DPMM parameters
#'
#' @return A list with n entries of n rows with missingness
#' 
#'
#' @import condMVNorm
#' @import mvtnorm
#' @import abind
#' @importFrom tidyr gather
#' @importFrom tidyr nest
#' @importFrom tidyr separate
#' @importFrom tidyr unnest
#' @importFrom tidyr pivot_longer
#' @importFrom purrr pmap
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @importFrom purrr map2
#' @importFrom purrr pluck
#' @importFrom rlang set_names
#' @importFrom utils data
#' 
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
#'                        
#' ## introduce missing data
#' rows <- 501:550
#' dataset_missing <- dataset_1
#' dataset_missing_predict <- dataset_missing[rows,]
#' dataset_missing_predict[,1] <- as.numeric(NA)
#' 
#' # predict missing values
#' posteriors.dpmmfit <- predict_dpmm_fit(posteriors, 
#'                                        dataset_missing_predict, 
#'                                        samples = c(1:100))
#' }
#' 
#'
#' @export
posterior_dpmm <- function(patient, samples, seed = NULL, cont_vars = NULL, cat_vars = NULL, mcmc_chain = NULL) {


  #:----------------------------------------------------------
  # check whether seed is provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #:----------------------------------------------------------
  # check whether cat_vars is provided
  if (!is.null(cat_vars)) {
    ndisc = length(cat_vars)
    # ndiscdim = max(unique(as.numeric(patient[,cat_vars])))
    
    ndiscdim <- NULL
    for (i in 1:length(cat_vars)) {
      ndiscdim <- c(ndiscdim, length(unique(levels(patient[,cat_vars[i]]))))
    }
    ndiscdim <- max(ndiscdim)
  }
  
  #:----------------------------------------------------------
  ## Iterate through all chains fitted in the DPMM
  for (mcmc_chains in 1:mcmc_chain) {
    
    # Extract the parameter posterior samples for the current chain 
    if (mcmc_chain == 1) {
      samples_chain <- samples
    } else {
      samples_chain <- samples[[mcmc_chains]]
    }
    
    #:----------------------------------------------------------
    ## Extract components weights
    v_mat <- samples_chain %>%
      as_tibble(.name_repair = 'unique') %>%
      select(starts_with("v[")) %>%
      as.matrix()
    
    stick_breaking <- function(v) {
      K <- ncol(v)
      N <- nrow(v)
      w <- matrix(0, nrow = N, ncol = K + 1)
      cumprod_1mv <- matrix(1, nrow = N, ncol = K)
      
      for (k in 2:K) {
        cumprod_1mv[, k] <- cumprod_1mv[, k - 1] * (1 - v[, k - 1])
      }
      w[, 1:K] <- v * cumprod_1mv
      w[, K + 1] <- apply(1 - v, 1, prod)
      return(w)
    }
    
    postW <- stick_breaking(v_mat) %>%
      as_tibble(.name_repair = "unique") %>%
      mutate(iter = row_number()) %>%
      tidyr::pivot_longer(-iter, names_to = "var", values_to = "value") %>%
      mutate(var = as.integer(gsub("V", "", var))) %>%
      arrange(iter, var) %>%
      group_by(iter) %>%
      summarise(w = list(value), .groups = "drop")
    
    #:----------------------------------------------------------
    ## If cont_vars is provided (and hence the model is fitted to continuous variables)
    if (!is.null(cont_vars)) {
      ### Extract posterior samples of component means
      mu_cols <- samples_chain %>%
        as_tibble(.name_repair = "unique") %>%
        select(starts_with("muL[")) %>%
        as.matrix()
      
      n_iter <- nrow(mu_cols)
      n_param <- ncol(mu_cols)
      
      postMu <- tibble(value = as.vector(mu_cols),
                       iter = rep(1:n_iter, times = n_param),
                       var = rep(colnames(mu_cols), each = n_iter)) %>%
        mutate(var = gsub("muL\\[|\\]", "", var)) %>%
        separate(var, c("component", "dim"), sep = ",", convert = TRUE) %>%
        arrange(iter, component, dim) %>%
        group_by(iter, component) %>%
        summarise(mu = list(value), .groups = "drop") %>%
        group_by(iter) %>%
        summarise(muL = list(t(simplify2array(mu))), .groups = "drop")
      
      ### Extract posterior samples of component precision matrices
      tau_cols <- samples_chain %>%
        as_tibble(.name_repair = "unique") %>%
        select(starts_with("tauL[")) %>%
        as.matrix()
      
      n_iter <- nrow(tau_cols)
      n_param <- ncol(tau_cols)
      
      postTau <- tibble(value = as.vector(tau_cols),
                        iter = rep(1:n_iter, times = n_param),
                        var = rep(colnames(tau_cols), each = n_iter)) %>%
        mutate(var = gsub("tauL\\[|\\]", "", var)) %>%
        separate(var, c("dim1", "dim2", "component"), sep = ",", convert = TRUE) %>%
        arrange(iter, component, dim2, dim1) %>%
        group_by(iter, component) %>%
        summarise(tau = list(matrix(value, nrow = max(dim1), ncol = max(dim2), byrow = TRUE)), .groups = "drop") %>%
        group_by(iter) %>%
        summarise(tauL = list(abind(tau, along = 3)), .groups = "drop")
      
    }
    
    #:----------------------------------------------------------
    ## If cat_vars is provided (and hence the model is fitted to categorical variables)
    if (!is.null(cat_vars)) {
      ### Extract posterior samples of component level probabilities
      phi_cols <- samples_chain %>%
        as_tibble(.name_repair = "unique") %>%
        select(starts_with("phiL[")) %>%
        as.matrix()
      
      n_iter <- nrow(phi_cols)
      n_param <- ncol(phi_cols)
      
      postPhi <- tibble(value = as.vector(phi_cols),
                        iter = rep(1:n_iter, times = n_param),
                        var = rep(colnames(phi_cols), each = n_iter)) %>%
        mutate(var = gsub("phiL\\[|\\]", "", var)) %>%
        separate(var, c("dim1", "dim2", "component"), sep = ",", convert = TRUE) %>%
        arrange(iter, component, dim1, dim2) %>%
        group_by(iter, component) %>%
        summarise(phi = list(matrix(value, nrow = ndisc, ncol = ndiscdim, byrow = TRUE)), .groups = "drop") %>%
        group_by(iter) %>%
        summarise(phiL = list(abind(phi, along = 3)), .groups = "drop")
      
    }
    
    #:----------------------------------------------------------
    ## Organise data in format needed for the function
    
    ## rows being used
    ordered_patient <- patient[,c(cont_vars,cat_vars)]
    
    ## The output is in the following shape
    ###    input: several patients with several variables
    ###    output: a list of patients, with matrix of missing vars
    
    number_vars <- length(c(cont_vars)) ## number of continuous variables
    number_patients <- nrow(ordered_patient) ## number of patients
    number_samples <- nrow(samples_chain)  ## number of iterations
    
    ## output list
    posterior <- vector(mode = "list", length = number_patients)
    
    ## for each entry on the list (each row)
    for (i in 1:number_patients) {
      ## check which variables are missing
      missing_vars <- ordered_patient[i,] %>%
        is.na() %>%
        colSums()
      
      ## collect the name of the variables missing
      missing_vars_names <- colnames(ordered_patient)[missing_vars > 0]
      
      ## add matrix of sampled values with the same number of rows are the model iterations
      posterior[[i]] <- matrix(as.numeric(NA), nrow = number_samples, ncol = length(missing_vars_names))
      
    }
    
    
    #:----------------------------------------------------------
    ## Iterate through each of the rows with missing data provided
    for (iteration_patient in 1:number_patients) {
      
      ## Extract data from current row
      current_patient <- ordered_patient[iteration_patient,]
      
      ## if current row has categorical variables, turn them into numerical values
      if (!is.null(cat_vars)) {
        current_patient[cat_vars] <- as.numeric(current_patient[cat_vars])
      }
      
      ## replicate row into number of samples
      current_data <- as.data.frame(lapply(current_patient, rep, number_samples)) %>%
        mutate(iter = 1:n()) %>%
        gather(var, value, -iter) %>%
        arrange(iter) %>%
        group_by(iter) %>%
        nest() %>%
        mutate(data = map(data, "value"))  %>%
        rename(patient = data)
      
      ## put together all posterior samples necessary for the prediction
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
      
      
      #:----------------------------------------------------------
      ## Iterate through each of the rows with missing data provided
      
      # draw conditional for only cont, only cat, cont and cat
      if (!is.null(cont_vars) & is.null(cat_vars)) {
        
        #:----------------------------------------------------------
        ## Only continuous variables in the dataset
        preds <- post %>%
          ## Adjust component weights depending on available values and choose one option
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
          ## Keep only means for the chosen component
          mutate(muL = map2(w, muL, function(w, mu) {
            mu[w, ]
          })) %>%
          ## Keep only precision matrices for the chosen component
          mutate(tauL = map2(w, tauL, function(w, tau) {
            tau[, , w]
          })) %>%
          ## Make sure precision matrices have equal lower and upper triangles (rounding error)
          mutate(tauL = map(tauL, function(tau) {
            tauL <- matrix(unlist(tauL), ncol = number_vars)
            
            for (i in 1:sqrt(length(tauL))) {
              tauL[-i,i] = tauL[i,-i]
            }
            tauL
          })) %>%
          ## Sample missing values
          mutate(value = pmap(list(patient, muL, tauL), function(patient, mu, tau) {
            
            # vars without missing values
            variables_given <- which(!is.na(patient))
            # vars with missing values
            variables_dep <- which(is.na(patient))
            
            # conditional distribution for missing vars
            density <- condMVNorm::condMVN(mu, solve(tau), dep = variables_dep, given = variables_given,
                                           X.given = patient[variables_given], check.sigma = FALSE)
            # draw for missing var
            preds <- cbind(t(mvrnorm(1, density$condMean, density$condVar))) %>%
              as.data.frame()
            preds
          }))
        
      } else if (is.null(cont_vars) & !is.null(cat_vars)) {
        ## Only categorical variables in the dataset
        preds <- post %>%
          ## Adjust component weights depending on available values and choose one option
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
          ## Keep only probabilities for the chosen component
          mutate(phiL = map2(w, phiL, function(w, psi) {
            psi[,,w]
          })) %>%
          ## Sample missing values
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
        ## Both continuous and categorical variables in the dataset
        preds <- post %>%
          ## Adjust component weights depending on available values and choose one option
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
                if (!isSymmetric(solve(sigma[, , i]))) {
                  value = as.numeric(NA)
                } else {
                  # calculate the density for continuous
                  value <- mvtnorm::dmvnorm(xcont, mu[i, ], solve(sigma[, , i]), log = TRUE) + log(w[i])
                  # calculate the density for categorical
                  for (l in 1:length(categorical_vars)) {
                    value <- value + log(phi[l ,xcat[categorical_vars[l]] , i])
                  }
                }
                value
                
              }, xcont = patient[continuous_vars_given], xcat = patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))], mu = marg_mean, sigma = marg_sigma, phi = phi)
              lmarg <- log_sum_exp(lmarg)
              
              z_given_x_l <- map_dbl(1:length(w), function(i, xcont, xcat, mu, sigma, phi, denom) {
                if (!isSymmetric(solve(sigma[, , i]))) {
                  value = NA
                } else {
                  # calculate the density for continuous
                  value <- mvtnorm::dmvnorm(xcont, mu[i, ], solve(sigma[, , i]), log = TRUE) + log(w[i]) - denom
                  # calculate the density for categorical
                  for (l in 1:length(categorical_vars)) {
                    value <- value + log(phi[l ,xcat[categorical_vars[l]] , i])
                  }
                }
                
                value
              }, xcont = patient[continuous_vars_given], xcat = patient[(length(cont_vars)+1):(length(cont_vars)+length(cat_vars))], mu = marg_mean, sigma = marg_sigma, phi = phi, denom = lmarg)
              
              z_given_x <- exp(z_given_x_l)
              
              z_given_x[is.na(z_given_x)] <- 0
              
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
          ## Keep only means for the chosen component
          mutate(muL = map2(w, muL, function(w, mu) {
            mu[w, ]
          })) %>%
          ## Keep only precision matrices for the chosen component
          mutate(tauL = map2(w, tauL, function(w, tau) {
            tau[, , w]
          })) %>%
          ## Make sure precision matrices have equal lower and upper triangles (rounding error)
          mutate(tauL = map(tauL, function(tau) {
            tauL <- matrix(unlist(tauL), ncol = number_vars)
            
            for (i in 1:sqrt(length(tauL))) {
              tauL[-i,i] = tauL[i,-i]
            }
            tauL
          })) %>%
          ## Keep only probabilities for the chosen component
          mutate(phiL = map2(w, phiL, function(w, psi) {
            psi[,,w]
          })) %>%
          ## Sample missing values
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
              preds <- cbind(t(mvrnorm(1, density$condMean, density$condVar)))
            }
            
            if (!is_empty(categorical_vars)) {
              for (i in 1:length(categorical_vars)) {
                if (length(cat_vars) == 1) {
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
      
      #:----------------------------------------------------------
      ## Extract only the values imputed
      preds <- preds %>%
        select(iter, value) %>%
        unnest(cols = value) %>%
        ungroup() %>%
        select(-iter) %>%
        set_names(colnames(current_patient)[is.na(current_patient)])
      
      #:----------------------------------------------------------
      ## Add the values to the list of estimated values
      for (posterior_columns in 1:ncol(posterior[[iteration_patient]])) {

        posterior[[iteration_patient]][,posterior_columns] <- unlist(preds[,posterior_columns])

      }
      
      posterior[[iteration_patient]] <- posterior[[iteration_patient]] %>%
        as.data.frame() %>%
        set_names(colnames(current_patient)[which(is.na(current_patient))])

    }
  
    
    #:----------------------------------------------------------
    ## Add the values to the output list
    if (mcmc_chains == 1) {
      # if on the first chain, just use the actual list
      predictions_final <- posterior
      
    } else {
      # if on other chains, append other chain values to the bottom of the matrix
      for (i in 1:length(predictions_final)) {
        
        predictions_final[[i]] <- rbind(
          predictions_final[[i]],
          posterior[[i]]
        )
        
      }
      
    }
  
  }
  
  
  return(predictions_final)
  
}

  
