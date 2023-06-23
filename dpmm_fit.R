## libraries
library(nimble)
library(MASS)
library(tidyverse)
library(abind)


## function to create and run nimble model
runModel <- function(dataset, mcmc_iterations = 2500 , L = 10, standardise = TRUE) {
  
  ## check inputs
  if(!is.numeric(L)) stop("'L' must be 'numeric'")
  if(L < 3) stop("'L' has to be 3 or more")
  if(!is.data.frame(dataset)) stop("'dataset' must be 'data.frame'")
  if(!is.numeric(mcmc_iterations)) stop("'mcmc_iterations' must be 'numeric'")
  if(!is.numeric(L)) stop("'L' must be 'numeric'")
  if(!is.logical(standardise)) stop("'standardise' must be 'logical'")
  
  continuous <- dplyr::select(dataset, where(is.numeric))
  discrete <- dplyr::select(dataset, where(Negate(is.numeric)))
  dataset <- cbind(continuous, discrete)
  if (standardise == TRUE) {
    if (ncol(continuous) != 0) {
      library(matrixStats)
      mean_values <- as.vector(continuous %>% colMeans(na.rm = TRUE)) %>% 
        as.data.frame() %>%
        t()
      colnames(mean_values) <- colnames(continuous)
      sd_values <- as.vector(colSds(as.matrix(continuous), na.rm = TRUE)) %>% 
        as.data.frame() %>%
        t()
      colnames(sd_values) <- colnames(continuous)
      continuous <- dplyr::select(dataset, where(is.numeric)) %>%
        apply(2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    }
  }
  if (ncol(continuous) != 0) {
    continuous <- as.matrix(continuous)
  }
  if (ncol(discrete) != 0) {
    discrete <- mutate(discrete, across(everything(), as.character)) %>%
      mutate(across(everything(), factor)) %>%
      mutate(across(everything(), as.numeric))
    discrete <- as.matrix(discrete)
  }
  
  
  
  if (ncol(continuous) > 0 & ncol(discrete) >0) {
    # model for both discrete and continuous vars
    rows_complete <- which(complete.cases(dataset))
    rows_incomplete <- which(!complete.cases(dataset))
    
    if (!is_empty(rows_incomplete)) {
      # if there is missing values
      if (length(rows_incomplete) == 1) {
        if (sum(!complete.cases(continuous[rows_incomplete,])) == 0) {
          rows_continuous_complete <- rows_incomplete
          rows_continuous_incomplete <- c()
        } else {
          rows_continuous_incomplete <- rows_incomplete
          rows_continuous_complete <- c()
        }
        if (sum(!complete.cases(discrete[rows_incomplete,])) == 0) {
          rows_discrete_complete <- rows_incomplete
          rows_discrete_incomplete <- c()
        } else {
          rows_discrete_incomplete <- rows_incomplete
          rows_discrete_complete <- c()
        }
        
        
      } else {
        
        rows_continuous_complete <- which(complete.cases(continuous[rows_incomplete,]))
        rows_continuous_incomplete <- which(!complete.cases(continuous[rows_incomplete,]))
        rows_discrete_complete <- which(complete.cases(discrete[rows_incomplete,]))
        rows_discrete_incomplete <- which(!complete.cases(discrete[rows_incomplete,]))
        
      }
      
    
      
      if (!is_empty(rows_continuous_incomplete) & !is_empty(rows_discrete_incomplete)) {
        # if there is missing values for continuous and discrete
        
        if (sum(length(rows_complete)) == 1) {
          discrete_complete <- discrete[-rows_incomplete,] %>%
            as.data.frame() %>%
            t()
          continuous_complete <- continuous[-rows_incomplete,] %>%
            as.data.frame() %>%
            t()
        } else {
          discrete_complete <- discrete[-rows_incomplete,] %>%
            as.data.frame()
          continuous_complete <- continuous[-rows_incomplete,] %>%
            as.data.frame()
        }
        
        if (sum(length(rows_incomplete)) == 1) {
          discrete_incomplete <- discrete[rows_incomplete,] %>%
            as.data.frame() %>%
            t()
          continuous_incomplete <- continuous[rows_incomplete,] %>%
            as.data.frame() %>%
            t()
        } else {
          discrete_incomplete <- discrete[rows_incomplete,] %>%
            as.data.frame()
          continuous_incomplete <- continuous[rows_incomplete,] %>%
            as.data.frame()
        }
        
        
        consts <- list(
          N = length(rows_complete),
          Nmiss = length(rows_incomplete),
          L = L,
          ndisc = ncol(discrete),
          ndiscdim = as.vector(apply(as.matrix(discrete_complete), 2, function(x) length(unique(x))))
        )
        
        tau0 <- apply(as.matrix(continuous_complete), 2, function(x) range(x,na.rm = TRUE))
        
        data <- list(
          x_disc = as.matrix(discrete_complete),
          x_disc_miss = as.matrix(discrete_incomplete),
          x_cont = as.matrix(continuous_complete),
          x_cont_miss = as.matrix(continuous_incomplete),
          mu0 = apply(as.matrix(continuous_complete), 2, function(x) mean(x, na.rm = TRUE)),
          tau0 = base::solve(diag(apply(tau0, 2, diff)^2)),
          R0 = base::solve(cov(as.matrix(continuous_complete))) / ncol(as.matrix(continuous_complete)),
          kappa0 = ncol(as.matrix(continuous_complete)),
          delta = matrix(rep(1, consts$ndisc * max(consts$ndiscdim)), nrow = consts$ndisc)
        )
        
        code <- nimbleCode({
          
          ## likelihood terms
          for (i in 1:N) {
            ## DPMM for continuous
            z[i] ~ dcat(w[1:L])
            x_cont[i, ] ~ dmnorm(muL[z[i], ], prec = tauL[, , z[i]])
            
            for (j in 1:ndisc) {
              if (length(ndiscdim) == 1) {
                x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim , z[i]])
              } else {
                x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , z[i]])
              }
            }
          }
          
          for (i in 1:Nmiss) {
            zmiss[i] ~ dcat(w[1:L])
            x_cont_miss[i, ] ~ dmnorm(muL[zmiss[i], ], prec = tauL[, , zmiss[i]])
            
            for (j in 1:ndisc) {
              if (length(ndiscdim) == 1) {
                x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim , zmiss[i]])
              } else {
                x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , zmiss[i]])
              }
            }
          }
          
          # priors for DPMM
          alpha ~ dgamma(shape = 2, rate = 1)
          for(i in 1:(L - 1)) {
            v[i] ~ dbeta(1, alpha)
          }
          w[1:L] <- stick_breaking(v[1:(L - 1)])
          
          # hyperpriors for continuous
          R1[, ] ~ dwish(R0[, ], kappa0)
          kappa1 ~ T(dexp(0.1), kappa0, )
          for(i in 1:L) {
            muL[i, ] ~ dmnorm(mu0[], prec = tau0[, ])
            tauL[, , i] ~ dwish(R1[, ], kappa1)
          }
          for(j in 1:ndisc) {
            for(i in 1:L) {
              if (length(ndiscdim) == 1) {
                phiL[j, 1:ndiscdim, i] ~ ddirch(delta[j, 1:ndiscdim])
              } else {
                phiL[j, 1:ndiscdim[j], i] ~ ddirch(delta[j, 1:ndiscdim[j]])
              }
            }
          }
        })
        
        ## sample initial values
        initFn <- function(L, N, Nmiss, mu0, tau0, R0, kappa0, ndiscdim, x_cont_miss, x_disc_miss) {
          
          # DPMM clustering
          alpha <- rgamma(1, shape = 2, rate = 1)
          v <- rbeta(L - 1, 1, alpha)
          w <- v[1]
          for(i in 2:(L - 1)) {
            w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
          }
          w <- c(w, prod(1 - v))
          z <- rcat(N, w)
          zmiss <- rcat(Nmiss, w)
          
          # DPMM for continuous
          R1 <- rWishart(1, kappa0, R0)[, , 1]
          kappa1 <- kappa0 - 1
          while(kappa1 < kappa0) {
            kappa1 <- rexp(1, 0.1)
          }
          muL <- mvrnorm(L, mu0, tau0)
          tauL <- rWishart(L, kappa1, R1)
          
          # DPMM discrete
          phiL <- map(1:L, function(i, ndiscdim, m, L) {
            p <- map(ndiscdim, function(n, m, L) {
              out <- numeric(m)
              out[1:n] <- nimble::rdirch(1, rep(1, n))
              out
            }, m = m, L = L)
            p <- do.call("rbind", p)
            p
          }, ndiscdim = ndiscdim, m = max(ndiscdim), L = L)
          phiL <- abind(phiL, along = 3)
          phiL <- array(phiL, dim = c(length(ndiscdim), max(ndiscdim), L))
          
          
          x_cont_miss1 <- x_cont_miss
          x_cont_miss1[!is.na(x_cont_miss)] <- NA
          x_cont_miss1[is.na(x_cont_miss)] <- 0
          x_disc_miss1 <- x_disc_miss
          x_disc_miss1[!is.na(x_disc_miss)] <- NA
          x_disc_miss1[is.na(x_disc_miss)] <- 1
          
          inits <- list(
            alpha = alpha,
            v = v,
            w = w,
            z = z,
            zmiss = zmiss,
            R1 = R1,
            kappa1 = kappa1,
            muL = muL,
            tauL = tauL,
            phiL = phiL,
            x_cont_miss = x_cont_miss1,
            x_disc_miss = x_disc_miss1
          )
          inits
        }
        
        model <- nimbleModel(
          code = code,
          constants = consts,
          data = data,
          inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim, data$x_cont_miss, data$x_disc_miss)
        )
        
        #compile the model
        cmodel <- compileNimble(model)
        
        #set monitors
        config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL", "x_cont_miss", "x_disc_miss"), thin = 1, print = FALSE)
        
        
        source("conditional_RW.R")
        source("conditional_RW_block.R")
        
        
        ## add custom sampler
        for(i in 1:nrow(data$x_cont_miss)) {
          
          if(sum(is.na(data$x_cont_miss[i, ])) == 1) {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'conditional_RW',
              control = list(scale = 1, index = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          } else if (sum(is.na(data$x_cont_miss[i, ])) == ncol(continuous) || sum(is.na(data$x_cont_miss[i, ])) == 0) {
            #nothing happens because we just want draws from the component
          } else {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'conditional_RW_block',
              control = list(scale = 1, indices = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          }
        }
        
        
      } else {
        if (!is_empty(rows_continuous_incomplete)) {
          # if there is missing values for continuous

          if (sum(c(length(rows_complete),length(rows_continuous_complete))) == 1) {
            discrete_complete <- discrete[-rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame() %>%
              t()
            continuous_complete <- continuous[-rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame() %>%
              t()
          } else {
            discrete_complete <- discrete[-rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame()
            continuous_complete <- continuous[-rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame() 
          }
          if (length(rows_continuous_incomplete) == 1) {
            discrete_incomplete <- discrete[rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame() %>%
              t()
            continuous_incomplete <- continuous[rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame() %>%
              t()
          } else {
            discrete_incomplete <- discrete[rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame()
            continuous_incomplete <- continuous[rows_incomplete[rows_continuous_incomplete],] %>%
              as.data.frame() 
          }
          
          consts <- list(
            N = length(rows_complete),
            Nmiss = length(rows_incomplete),
            L = L,
            ndisc = ncol(discrete),
            ndiscdim = as.vector(apply(as.matrix(discrete_complete), 2, function(x) length(unique(x))))
          )
          
          tau0 <- apply(as.matrix(continuous_complete), 2, range)
          
          data <- list(
            x_disc = as.matrix(discrete_complete),
            x_disc_miss = as.matrix(discrete_incomplete),
            x_cont = as.matrix(continuous_complete),
            x_cont_miss = as.matrix(continuous_incomplete),
            mu0 = apply(as.matrix(continuous_complete), 2, mean),
            tau0 = base::solve(diag(apply(tau0, 2, diff)^2)),
            R0 = base::solve(cov(as.matrix(continuous_complete))) / ncol(as.matrix(continuous_complete)),
            kappa0 = ncol(as.matrix(continuous_complete)),
            delta = matrix(rep(1, consts$ndisc * max(consts$ndiscdim)), nrow = consts$ndisc)
          )
          
          
          code <- nimbleCode({
            
            ## likelihood terms
            for (i in 1:N) {
              ## DPMM for continuous
              z[i] ~ dcat(w[1:L])
              x_cont[i, ] ~ dmnorm(muL[z[i], ], prec = tauL[, , z[i]])
              
              for (j in 1:ndisc) {
                if (length(ndiscdim) == 1) {
                  x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim , z[i]])
                } else {
                  x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , z[i]])
                }
              }
            }
            
            for (i in 1:Nmiss) {
              zmiss[i] ~ dcat(w[1:L])
              x_cont_miss[i, ] ~ dmnorm(muL[zmiss[i], ], prec = tauL[, , zmiss[i]])
              
              for (j in 1:ndisc) {
                if (length(ndiscdim) == 1) {
                  x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim , zmiss[i]])
                } else {
                  x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , zmiss[i]])
                }
              }
            }
            
            # priors for DPMM
            alpha ~ dgamma(shape = 2, rate = 1)
            for(i in 1:(L - 1)) {
              v[i] ~ dbeta(1, alpha)
            }
            w[1:L] <- stick_breaking(v[1:(L - 1)])
            
            # hyperpriors for continuous
            R1[, ] ~ dwish(R0[, ], kappa0)
            kappa1 ~ T(dexp(0.1), kappa0, )
            for(i in 1:L) {
              muL[i, ] ~ dmnorm(mu0[], prec = tau0[, ])
              tauL[, , i] ~ dwish(R1[, ], kappa1)
            }
            for(j in 1:ndisc) {
              for(i in 1:L) {
                if (length(ndiscdim) == 1) {
                  phiL[j, 1:ndiscdim, i] ~ ddirch(delta[j, 1:ndiscdim])
                } else {
                  phiL[j, 1:ndiscdim[j], i] ~ ddirch(delta[j, 1:ndiscdim[j]])
                }
              }
            }
          })
          
          
          ## sample initial values
          initFn <- function(L, N, Nmiss, mu0, tau0, R0, kappa0, ndiscdim, x_cont_miss) {
            
            # DPMM clustering
            alpha <- rgamma(1, shape = 2, rate = 1)
            v <- rbeta(L - 1, 1, alpha)
            w <- v[1]
            for(i in 2:(L - 1)) {
              w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
            }
            w <- c(w, prod(1 - v))
            z <- rcat(N, w)
            zmiss <- rcat(Nmiss, w)
            
            # DPMM for continuous
            R1 <- rWishart(1, kappa0, R0)[, , 1]
            kappa1 <- kappa0 - 1
            while(kappa1 < kappa0) {
              kappa1 <- rexp(1, 0.1)
            }
            muL <- mvrnorm(L, mu0, tau0)
            tauL <- rWishart(L, kappa1, R1)
            
            # DPMM discrete
            phiL <- map(1:L, function(i, ndiscdim, m, L) {
              p <- map(ndiscdim, function(n, m, L) {
                out <- numeric(m)
                out[1:n] <- nimble::rdirch(1, rep(1, n))
                out
              }, m = m, L = L)
              p <- do.call("rbind", p)
              p
            }, ndiscdim = ndiscdim, m = max(ndiscdim), L = L)
            phiL <- abind(phiL, along = 3)
            phiL <- array(phiL, dim = c(length(ndiscdim), max(ndiscdim), L))
            
            
            x_cont_miss1 <- x_cont_miss
            x_cont_miss1[!is.na(x_cont_miss)] <- NA
            x_cont_miss1[is.na(x_cont_miss)] <- 0
            
            inits <- list(
              alpha = alpha,
              v = v,
              w = w,
              z = z,
              zmiss = zmiss,
              R1 = R1,
              kappa1 = kappa1,
              muL = muL,
              tauL = tauL,
              phiL = phiL,
              x_cont_miss = x_cont_miss1
            )
            inits
          }
          
          model <- nimbleModel(
            code = code,
            constants = consts,
            data = data,
            inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim, data$x_cont_miss)
          )
          
          #compile the model
          cmodel <- compileNimble(model)
          
          #set monitors
          config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL", "x_cont_miss"), thin = 1, print = FALSE)
          
          
          source("conditional_RW.R")
          source("conditional_RW_block.R")
          
          
          ## add custom sampler
          ## add custom sampler
          for(i in 1:nrow(data$x_cont_miss)) {
            
            if(sum(is.na(data$x_cont_miss[i, ])) == 1) {
              target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
              config$removeSampler(target)
              config$addSampler(
                target = target,
                type = 'conditional_RW',
                control = list(scale = 1, index = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
              )
            } else if (sum(is.na(data$x_cont_miss[i, ])) == ncol(continuous) || sum(is.na(data$x_cont_miss[i, ])) == 0) {
              #nothing happens because we just want draws from the component
            } else {
              target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
              config$removeSampler(target)
              config$addSampler(
                target = target,
                type = 'conditional_RW_block',
                control = list(scale = 1, indices = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
              )
            }
          }
          
          
        } else {
          # if there is missing values for discrete
          tau0 <- apply(as.matrix(continuous), 2, range)
          
          
          if (sum(c(length(rows_complete),length(rows_discrete_complete))) == 1) {
            discrete_complete <- discrete[-rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame() %>%
              t()
            continuous_complete <- continuous[-rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame() %>%
              t()
          } else {
            discrete_complete <- discrete[-rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame()
            continuous_complete <- continuous[-rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame() 
          }
          if (length(rows_discrete_incomplete) == 1) {
            discrete_incomplete <- discrete[rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame() %>%
              t()
            continuous_incomplete <- continuous[rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame() %>%
              t()
          } else {
            discrete_incomplete <- discrete[rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame()
            continuous_incomplete <- continuous[rows_incomplete[rows_discrete_incomplete],] %>%
              as.data.frame() 
          }
          
          
          consts <- list(
            N = length(rows_complete),
            Nmiss = length(rows_incomplete),
            L = L,
            ndisc = ncol(discrete),
            ndiscdim = as.vector(apply(as.matrix(discrete_complete), 2, function(x) length(unique(x))))
          )
          
          
          data <- list(
            x_disc = as.matrix(discrete_complete),
            x_disc_miss = as.matrix(discrete_incomplete),
            x_cont = as.matrix(continuous_complete),
            x_cont_miss = as.matrix(continuous_incomplete),
            mu0 = apply(as.matrix(continuous), 2, mean),
            tau0 = base::solve(diag(apply(tau0, 2, diff)^2)),
            R0 = base::solve(cov(as.matrix(continuous))) / ncol(as.matrix(continuous)),
            kappa0 = ncol(as.matrix(continuous)),
            delta = matrix(rep(1, consts$ndisc * max(consts$ndiscdim)), nrow = consts$ndisc)
          )
          
          code <- nimbleCode({
            
            ## likelihood terms
            for(i in 1:N) {
              ## DPMM for continuous
              z[i] ~ dcat(w[1:L])
              x_cont[i, ] ~ dmnorm(muL[z[i], ], prec = tauL[, , z[i]])
              
              for (j in 1:ndisc) {
                if (length(ndiscdim) == 1) {
                  x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim , z[i]])
                  } else {
                  x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , z[i]])
                  }
                }
              }
            
            for(i in 1:Nmiss) {
              ## DPMM for continuous
              zmiss[i] ~ dcat(w[1:L])
              x_cont_miss[i, ] ~ dmnorm(muL[zmiss[i], ], prec = tauL[, , zmiss[i]])
              
              for (j in 1:ndisc) {
                if (length(ndiscdim) == 1) {
                  x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim , zmiss[i]])
                } else {
                  x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , zmiss[i]])
                }
              }
            }
            
            # priors for DPMM
            alpha ~ dgamma(shape = 2, rate = 1)
            for(i in 1:(L - 1)) {
              v[i] ~ dbeta(1, alpha)
            }
            w[1:L] <- stick_breaking(v[1:(L - 1)])
            
            # hyperpriors for continuous
            R1[, ] ~ dwish(R0[, ], kappa0)
            kappa1 ~ T(dexp(0.1), kappa0, )
            for(i in 1:L) {
              muL[i, ] ~ dmnorm(mu0[], prec = tau0[, ])
              tauL[, , i] ~ dwish(R1[, ], kappa1)
              
            }
            for(j in 1:ndisc) {
              for(i in 1:L) {
                if (length(ndiscdim) == 1) {
                  phiL[j, 1:ndiscdim, i] ~ ddirch(delta[j, 1:ndiscdim])
                } else {
                  phiL[j, 1:ndiscdim[j], i] ~ ddirch(delta[j, 1:ndiscdim[j]])
                }
              }
            }
            
          })
          
          ## sample initial values
          initFn <- function(L, N, Nmiss, mu0, tau0, R0, kappa0, ndiscdim, x_disc_miss) {
            
            # DPMM clustering
            alpha <- rgamma(1, shape = 2, rate = 1)
            v <- rbeta(L - 1, 1, alpha)
            w <- v[1]
            for(i in 2:(L - 1)) {
              w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
            }
            w <- c(w, prod(1 - v))
            z <- rcat(N, w)
            zmiss <- rcat(Nmiss, w)
            
            # DPMM for continuous
            R1 <- rWishart(1, kappa0, R0)[, , 1]
            kappa1 <- kappa0 - 1
            while(kappa1 < kappa0) {
              kappa1 <- rexp(1, 0.1)
            }
            muL <- mvrnorm(L, mu0, tau0)
            tauL <- rWishart(L, kappa1, R1)
            
            # DPMM discrete
            phiL <- map(1:L, function(i, ndiscdim, m, L) {
              p <- map(ndiscdim, function(n, m, L) {
                out <- numeric(m)
                out[1:n] <- nimble::rdirch(1, rep(1, n))
                out
              }, m = m, L = L)
              p <- do.call("rbind", p)
              p
            }, ndiscdim = ndiscdim, m = max(ndiscdim), L = L)
            phiL <- abind(phiL, along = 3)
            phiL <- array(phiL, dim = c(length(ndiscdim), max(ndiscdim), L))
            
            x_disc_miss1 <- x_disc_miss
            x_disc_miss1[!is.na(x_disc_miss1)] <- NA
            x_disc_miss1[is.na(x_disc_miss1)] <- 1
            
            inits <- list(
              alpha = alpha,
              v = v,
              w = w,
              z = z,
              zmiss = zmiss,
              R1 = R1,
              kappa1 = kappa1,
              muL = muL,
              tauL = tauL,
              phiL = phiL,
              x_disc_miss = x_disc_miss1
            )
            inits
          }
          
          model <- nimbleModel(
            code = code,
            constants = consts,
            data = data,
            inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim, data$x_disc_miss)
          )
          
          #compile the model
          cmodel <- compileNimble(model)
          
          #set monitors
          config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL", "x_disc_miss"), thin = 1, print = FALSE)
          
          
        }
      }
      
      
    } else {
      
      consts <- list(
        N = nrow(dataset),
        L = L,
        ndisc = ncol(discrete),
        ndiscdim = as.vector(apply(as.matrix(discrete), 2, function(x) length(unique(x))))
      )
      
      tau0 <- apply(as.matrix(continuous), 2, range)
      
      data <- list(
        x_disc = as.matrix(discrete),
        x_cont = as.matrix(continuous),
        mu0 = apply(as.matrix(continuous), 2, mean),
        tau0 = base::solve(diag(apply(tau0, 2, diff)^2)),
        R0 = base::solve(cov(as.matrix(continuous))) / ncol(as.matrix(continuous)),
        kappa0 = ncol(as.matrix(continuous)),
        delta = matrix(rep(1, consts$ndisc * max(consts$ndiscdim)), nrow = consts$ndisc)
      )
      
      code <- nimbleCode({
        
        ## likelihood terms
        for(i in 1:N) {
          ## DPMM for continuous
          z[i] ~ dcat(w[1:L])
          x_cont[i, ] ~ dmnorm(muL[z[i], ], prec = tauL[, , z[i]])
          
          for (j in 1:ndisc) {
            if (length(ndiscdim) == 1) {
              x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim , z[i]])
            } else {
              x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , z[i]])
            }
          }
        }
        
        # priors for DPMM
        alpha ~ dgamma(shape = 2, rate = 1)
        for(i in 1:(L - 1)) {
          v[i] ~ dbeta(1, alpha)
        }
        w[1:L] <- stick_breaking(v[1:(L - 1)])
        
        # hyperpriors for continuous
        R1[, ] ~ dwish(R0[, ], kappa0)
        kappa1 ~ T(dexp(0.1), kappa0, )
        for(i in 1:L) {
          muL[i, ] ~ dmnorm(mu0[], prec = tau0[, ])
          tauL[, , i] ~ dwish(R1[, ], kappa1)
          
        }
        for(j in 1:ndisc) {
          for(i in 1:L) {
            if (length(ndiscdim) == 1) {
              phiL[j, 1:ndiscdim, i] ~ ddirch(delta[j, 1:ndiscdim])
            } else {
              phiL[j, 1:ndiscdim[j], i] ~ ddirch(delta[j, 1:ndiscdim[j]])
            }
          }
        }
        
      })
      
      
      ## sample initial values
      initFn <- function(L, N, mu0, tau0, R0, kappa0, ndiscdim) {
        
        # DPMM clustering
        alpha <- rgamma(1, shape = 2, rate = 1)
        v <- rbeta(L - 1, 1, alpha)
        w <- v[1]
        for(i in 2:(L - 1)) {
          w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
        }
        w <- c(w, prod(1 - v))
        z <- rcat(N, w)
        
        # DPMM for continuous
        R1 <- rWishart(1, kappa0, R0)[, , 1]
        kappa1 <- kappa0 - 1
        while(kappa1 < kappa0) {
          kappa1 <- rexp(1, 0.1)
        }
        muL <- mvrnorm(L, mu0, tau0)
        tauL <- rWishart(L, kappa1, R1)
        
        # DPMM discrete
        phiL <- map(1:L, function(i, ndiscdim, m, L) {
          p <- map(ndiscdim, function(n, m, L) {
            out <- numeric(m)
            out[1:n] <- nimble::rdirch(1, rep(1, n))
            out
          }, m = m, L = L)
          p <- do.call("rbind", p)
          p
        }, ndiscdim = ndiscdim, m = max(ndiscdim), L = L)
        phiL <- abind(phiL, along = 3)
        phiL <- array(phiL, dim = c(length(ndiscdim), max(ndiscdim), L))
        
        inits <- list(
          alpha = alpha,
          v = v,
          w = w,
          z = z,
          R1 = R1,
          kappa1 = kappa1,
          muL = muL,
          tauL = tauL,
          phiL = phiL
        )
        inits
      }
      
      model <- nimbleModel(
        code = code,
        constants = consts,
        data = data,
        inits = initFn(consts$L, consts$N, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim)
      )
      
      #compile the model
      cmodel <- compileNimble(model)
      
      #set monitors
      config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL"), thin = 1, print = FALSE)
      
          
    }
    
    
    
  } else {
    if (ncol(continuous) > 0) {
      # model for only continuous
      continuous_complete <- continuous[complete.cases(continuous),] %>%
        as.data.frame()
      if (is.null(nrow(continuous[!complete.cases(continuous),]))) {
        continuous_incomplete <- continuous[!complete.cases(continuous),] %>%
          as.data.frame() %>%
          t()
      } else {
        continuous_incomplete <- continuous[!complete.cases(continuous),] %>%
          as.data.frame()
      }
      
      if (nrow(continuous_incomplete) != 0) {
        # if some missing values
        
        # load library
        source("conditional_RW.R")
        source("conditional_RW_block.R")
        
        
        consts <- list(
          N = nrow(continuous_complete),
          Nmiss = nrow(continuous_incomplete),
          L = L
        )
        
        tau0 <- apply(as.matrix(continuous_complete), 2, range)
        
        data <- list(
          x_cont = as.matrix(continuous_complete),
          
          x_cont_miss = as.matrix(continuous_incomplete),
          
          mu0 = apply(as.matrix(continuous_complete), 2, mean),
          tau0 = base::solve(diag(apply(tau0, 2, diff)^2)),
          R0 = base::solve(cov(as.matrix(continuous_complete))) / ncol(as.matrix(continuous_complete)),
          kappa0 = ncol(as.matrix(continuous_complete))
        )
        
        code <- nimbleCode({
          
          ## likelihood terms
          for (i in 1:N) {
            ## DPMM for continuous
            z[i] ~ dcat(w[1:L])
            x_cont[i, ] ~ dmnorm(muL[z[i], ], prec = tauL[, , z[i]])
            
          }
          
          for (i in 1:Nmiss) {
            zmiss[i] ~ dcat(w[1:L])
            x_cont_miss[i, ] ~ dmnorm(muL[zmiss[i], ], prec = tauL[, , zmiss[i]])
          }
          
          # priors for DPMM
          alpha ~ dgamma(shape = 2, rate = 1)
          for(i in 1:(L - 1)) {
            v[i] ~ dbeta(1, alpha)
          }
          w[1:L] <- stick_breaking(v[1:(L - 1)])
          
          # hyperpriors for continuous
          R1[, ] ~ dwish(R0[, ], kappa0)
          kappa1 ~ T(dexp(0.1), kappa0, )
          for(i in 1:L) {
            muL[i, ] ~ dmnorm(mu0[], prec = tau0[, ])
            tauL[, , i] ~ dwish(R1[, ], kappa1)
          }
        })
        
        ## sample initial values
        initFn <- function(L, N, Nmiss, mu0, tau0, R0, kappa0, x_cont_miss) {
          
          # DPMM clustering
          alpha <- rgamma(1, shape = 2, rate = 1)
          v <- rbeta(L - 1, 1, alpha)
          w <- v[1]
          for(i in 2:(L - 1)) {
            w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
          }
          w <- c(w, prod(1 - v))
          z <- rcat(N, w)
          zmiss <- rcat(Nmiss, w)
          
          # DPMM for continuous
          R1 <- rWishart(1, kappa0, R0)[, , 1]
          kappa1 <- kappa0 - 1
          while(kappa1 < kappa0) {
            kappa1 <- rexp(1, 0.1)
          }
          muL <- mvrnorm(L, mu0, tau0)
          tauL <- rWishart(L, kappa1, R1)
          
          x_cont_miss1 <- x_cont_miss
          x_cont_miss1[!is.na(x_cont_miss)] <- NA
          x_cont_miss1[is.na(x_cont_miss)] <- 0
          
          inits <- list(
            alpha = alpha,
            v = v,
            w = w,
            z = z,
            zmiss = zmiss,
            R1 = R1,
            kappa1 = kappa1,
            muL = muL,
            tauL = tauL,
            x_cont_miss = x_cont_miss1
          )
          inits
        }
        
        
        model <- nimbleModel(
          code = code,
          constants = consts,
          data = data,
          inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, data$x_cont_miss)
        )
        
        
        #compile the model
        cmodel <- compileNimble(model)
        
        #set monitors
        config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "x_cont_miss"), thin = 1, print = FALSE)
        
        
        ## add custom sampler
        for(i in 1:nrow(data$x_cont_miss)) {
          
          if(sum(is.na(data$x_cont_miss[i, ])) == 1) {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'conditional_RW',
              control = list(scale = 1, index = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          } else if (sum(is.na(data$x_cont_miss[i, ])) == ncol(continuous) || sum(is.na(data$x_cont_miss[i, ])) == 0) {
            #nothing happens because we just want draws from the component
          } else {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'conditional_RW_block',
              control = list(scale = 1, indices = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          }
        }
        
        
      } else {
        # if no missing values
        
        consts <- list(
          N = nrow(continuous_complete),
          L = L
        )
        
        tau0 <- apply(as.matrix(continuous_complete), 2, range)
        
        data <- list(
          x_cont = as.matrix(continuous_complete),
          mu0 = apply(as.matrix(continuous_complete), 2, mean),
          tau0 = base::solve(diag(apply(tau0, 2, diff)^2)),
          R0 = base::solve(cov(as.matrix(continuous_complete))) / ncol(as.matrix(continuous_complete)),
          kappa0 = ncol(as.matrix(continuous_complete))
        )
        
        
        code <- nimbleCode({
          
          ## likelihood terms
          for(i in 1:N) {
            ## DPMM for continuous
            z[i] ~ dcat(w[1:L])
            x_cont[i, ] ~ dmnorm(muL[z[i], ], prec = tauL[, , z[i]])
            
          }
          
          # priors for DPMM
          alpha ~ dgamma(shape = 2, rate = 1)
          for(i in 1:(L - 1)) {
            v[i] ~ dbeta(1, alpha)
          }
          w[1:L] <- stick_breaking(v[1:(L - 1)])
          
          # hyperpriors for continuous
          R1[, ] ~ dwish(R0[, ], kappa0)
          kappa1 ~ T(dexp(0.1), kappa0, )
          for(i in 1:L) {
            muL[i, ] ~ dmnorm(mu0[], prec = tau0[, ])
            tauL[, , i] ~ dwish(R1[, ], kappa1)
          }
        })
        
        ## sample initial values
        initFn <- function(L, N, mu0, tau0, R0, kappa0) {
          
          # DPMM clustering
          alpha <- rgamma(1, shape = 2, rate = 1)
          v <- rbeta(L - 1, 1, alpha)
          w <- v[1]
          for(i in 2:(L - 1)) {
            w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
          }
          w <- c(w, prod(1 - v))
          z <- rcat(N, w)
          
          # DPMM for continuous
          R1 <- rWishart(1, kappa0, R0)[, , 1]
          kappa1 <- kappa0 - 1
          while(kappa1 < kappa0) {
            kappa1 <- rexp(1, 0.1)
          }
          muL <- mvrnorm(L, mu0, tau0)
          tauL <- rWishart(L, kappa1, R1)
          
          
          inits <- list(
            alpha = alpha,
            v = v,
            w = w,
            z = z,
            R1 = R1,
            kappa1 = kappa1,
            muL = muL,
            tauL = tauL
          )
          inits
        }
        
        
        model <- nimbleModel(
          code = code,
          constants = consts,
          data = data,
          inits = initFn(consts$L, consts$N, data$mu0, data$tau0, data$R0, data$kappa0)
        )
        
        
        #compile the model
        cmodel <- compileNimble(model)
        
        #set monitors
        config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha"), thin = 1, print = FALSE)
        
      }
      
    } else {
      # model for only categorical
      discrete_complete <- discrete[complete.cases(discrete),] %>%
        as.data.frame()
      if (is.null(nrow(discrete[!complete.cases(discrete),]))) {
        discrete_incomplete <- discrete[!complete.cases(discrete),] %>%
          as.data.frame() %>%
          t()
      } else {
        discrete_incomplete <- discrete[!complete.cases(discrete),] %>%
          as.data.frame()
      }
      
      
      if (nrow(discrete_incomplete) != 0) {
        # if some missing values
        
        consts <- list(
          N = nrow(discrete_complete),
          Nmiss = nrow(discrete_incomplete),
          L = L,
          ndisc = ncol(discrete_complete),
          ndiscdim = as.vector(apply(as.matrix(discrete_complete), 2, function(x) length(unique(x))))
        )
        
        data <- list(
          x_disc = as.matrix(discrete_complete),
          x_disc_miss = as.matrix(discrete_incomplete),
          delta = matrix(rep(1, consts$ndisc * max(consts$ndiscdim)), nrow = consts$ndisc)
        )
        
        code <- nimbleCode({
          
          ## likelihood terms
          for(i in 1:N) {
            ## DPMM for continuous
            z[i] ~ dcat(w[1:L])
            
            for (j in 1:ndisc) {
              if (length(ndiscdim) == 1) {
                x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim , z[i]])
              } else {
                x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , z[i]])
              }
            }
          }
          
          for(i in 1:Nmiss) {
            ## DPMM for continuous
            zmiss[i] ~ dcat(w[1:L])
            
            for (j in 1:ndisc) {
              if (length(ndiscdim) == 1) {
                x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim , zmiss[i]])
              } else {
                x_disc_miss[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , zmiss[i]])
              }
            }
          }
          
          # priors for DPMM
          alpha ~ dgamma(shape = 2, rate = 1)
          for(i in 1:(L - 1)) {
            v[i] ~ dbeta(1, alpha)
          }
          w[1:L] <- stick_breaking(v[1:(L - 1)])
          
          for(j in 1:ndisc) {
            for(i in 1:L) {
              if (length(ndiscdim) == 1) {
                phiL[j, 1:ndiscdim, i] ~ ddirch(delta[j, 1:ndiscdim])
              } else {
                phiL[j, 1:ndiscdim[j], i] ~ ddirch(delta[j, 1:ndiscdim[j]])
              }
            }
          }
        })
        
        ## sample initial values
        initFn <- function(L, N, Nmiss, ndiscdim, x_disc_miss) {
          
          # DPMM clustering
          alpha <- rgamma(1, shape = 2, rate = 1)
          v <- rbeta(L - 1, 1, alpha)
          w <- v[1]
          for(i in 2:(L - 1)) {
            w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
          }
          w <- c(w, prod(1 - v))
          z <- rcat(N, w)
          zmiss <- rcat(Nmiss, w)
          
          
          # DPMM discrete
          phiL <- map(1:L, function(i, ndiscdim, m, L) {
            p <- map(ndiscdim, function(n, m, L) {
              out <- numeric(m)
              out[1:n] <- nimble::rdirch(1, rep(1, n))
              out
            }, m = m, L = L)
            p <- do.call("rbind", p)
            p
          }, ndiscdim = ndiscdim, m = max(ndiscdim), L = L)
          phiL <- abind(phiL, along = 3)
          phiL <- array(phiL, dim = c(length(ndiscdim), max(ndiscdim), L))
          
          
          x_disc_miss1 <- x_disc_miss
          x_disc_miss1[!is.na(x_disc_miss)] <- NA
          x_disc_miss1[is.na(x_disc_miss)] <- 1
          
          
          inits <- list(
            alpha = alpha,
            v = v,
            w = w,
            z = z,
            zmiss = zmiss,
            phiL = phiL,
            x_disc_miss = x_disc_miss1
          )
          inits
        }
        
        model <- nimbleModel(
          code = code,
          constants = consts,
          data = data,
          inits = initFn(consts$L, consts$N, consts$Nmiss, consts$ndiscdim, data$x_disc_miss)
        )
        
        #compile the model
        cmodel <- compileNimble(model)
        
        #set monitors
        config <- configureMCMC(cmodel, monitors = c("v","z", "alpha", "phiL", "x_disc_miss"), thin = 1, print = FALSE)
        

      } else {
      
        consts <- list(
          N = nrow(discrete),
          L = L,
          ndisc = ncol(discrete_complete),
          ndiscdim = as.vector(apply(as.matrix(discrete), 2, function(x) length(unique(x))))
        )
        
        data <- list(
          x_disc = as.matrix(discrete),
          delta = matrix(rep(1, consts$ndisc * max(consts$ndiscdim)), nrow = consts$ndisc)
        )
        
        code <- nimbleCode({
          
          ## likelihood terms
          for(i in 1:N) {
            ## DPMM for continuous
            z[i] ~ dcat(w[1:L])
            
            for (j in 1:ndisc) {
              if (length(ndiscdim) == 1) {
                x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim , z[i]])
              } else {
                x_disc[i, j] ~ dcat(phiL[j , 1:ndiscdim[j] , z[i]])
              }
            }
          }
          
          # priors for DPMM
          alpha ~ dgamma(shape = 2, rate = 1)
          for(i in 1:(L - 1)) {
            v[i] ~ dbeta(1, alpha)
          }
          w[1:L] <- stick_breaking(v[1:(L - 1)])
          
          for(j in 1:ndisc) {
            for(i in 1:L) {
              if (length(ndiscdim) == 1) {
                phiL[j, 1:ndiscdim, i] ~ ddirch(delta[j, 1:ndiscdim])
              } else {
                phiL[j, 1:ndiscdim[j], i] ~ ddirch(delta[j, 1:ndiscdim[j]])
              }
            }
          }
        })
        
        ## sample initial values
        initFn <- function(L, N, ndiscdim) {
          
          # DPMM clustering
          alpha <- rgamma(1, shape = 2, rate = 1)
          v <- rbeta(L - 1, 1, alpha)
          w <- v[1]
          for(i in 2:(L - 1)) {
            w <- c(w, v[i] * prod(1 - v[1:(i - 1)]))
          }
          w <- c(w, prod(1 - v))
          z <- rcat(N, w)
          
          
          # DPMM discrete
          phiL <- map(1:L, function(i, ndiscdim, m, L) {
            p <- map(ndiscdim, function(n, m, L) {
              out <- numeric(m)
              out[1:n] <- nimble::rdirch(1, rep(1, n))
              out
            }, m = m, L = L)
            p <- do.call("rbind", p)
            p
          }, ndiscdim = ndiscdim, m = max(ndiscdim), L = L)
          phiL <- abind(phiL, along = 3)
          phiL <- array(phiL, dim = c(length(ndiscdim), max(ndiscdim), L))
          
          inits <- list(
            alpha = alpha,
            v = v,
            w = w,
            z = z,
            phiL = phiL
          )
          inits
        }
        
        model <- nimbleModel(
          code = code,
          constants = consts,
          data = data,
          inits = initFn(consts$L, consts$N, consts$ndiscdim)
        )
        
        
        #compile the model
        cmodel <- compileNimble(model)
        
        #set monitors
        config <- configureMCMC(cmodel, monitors = c("v","z", "alpha", "phiL"), thin = 1, print = FALSE)
        
      }
    }
    
  }
  ## print config
  print(config)
    
  # build the model
  built <- buildMCMC(config)
  # compile model
  cbuilt <- compileNimble(built)
  
  # run model
  cbuilt$run(niter = mcmc_iterations, reset = TRUE)
  
  # collect samples
  samples <- as.matrix(cbuilt$mvSamples) %>%
    as_tibble()
  
  if (standardise == TRUE) {
    if (ncol(continuous) != 0) {
      output <- list(dataset = dataset, L = L, samples = samples, standardise = standardise, mean_values = mean_values, sd_values = sd_values)
    } else {
      output <- list(dataset = dataset, L = L, samples = samples, standardise = standardise)
    }
  } else {  
    output <- list(dataset = dataset, L = L, samples = samples, standardise = standardise)
  }
  
  class(output) <- "dpmm_fit"
  
  return(output)
}















