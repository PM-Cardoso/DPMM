#' Fit DPMM
#'
#' This function fits the Dirichlet Process Mixture Model (DPMM).
#'
#' @param dataset dataframe with continuous and/or categorical variables.
#' @param mcmc_iterations Number of MCMC iterations.
#' @param thinning Interval for collecting MCMC samples.
#' @param L Number of DPMM components fitted.
#' @param mcmc_chains Number of MCMC chains fitted.
#' @param standardise Should continuous variables be standardised. default = TRUE
#'
#' @return
#' \subsection{Output: List of class 'dpmm_fit'}{
#'    \describe{
#'      \item{dataset}{dataframe with 1 row but defined the same way as the dataset fitted.}
#'      \item{L}{Number of DPMM components fitted.}
#'      \item{mcmc_chains}{Number of MCMC chains fitted.}
#'      \item{samples}{MCMC samples.}
#'      \item{standardise}{TRUE or FALSE whether values were standardised.}
#'      \item{mean_values}{Mean values for covariates standardised.}
#'      \item{sd_values}{SD values for covariates standardised.}
#'    }
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom stats complete.cases
#' @importFrom stats cov
#' @importFrom stats rWishart
#' @importFrom stats rbeta
#' @importFrom stats rexp
#' @importFrom stats rgamma
#' @importFrom stats rmultinom
#' @importFrom stats sd
#' @importFrom stats setNames
#' @importFrom stats var
#' @import nimble
#' @import abind
#' @import synthpop
#' @importFrom tidyselect vars_select_helpers
#' @importFrom rlang is_empty
#' @importFrom purrr map
#' @importFrom MASS mvrnorm
#' @importFrom matrixStats colSds
#' @import dplyr
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
#' }
#' 
#' 
#' @export
runModel <- function(dataset, mcmc_iterations = 2500, thinning = 1, L = 10, mcmc_chains = 2, standardise = TRUE) {
  
  #:-------------------------------------------------------------------
  ## check inputs
  if(!is.numeric(L)) stop("'L' must be 'numeric'")
  if(L < 3) stop("'L' has to be 3 or more")
  if(!is.data.frame(dataset)) stop("'dataset' must be 'data.frame'")
  if(!is.numeric(mcmc_iterations)) stop("'mcmc_iterations' must be 'numeric'")
  if(!is.numeric(thinning)) stop("'thinning' must be 'numeric'")
  if(!is.logical(standardise)) stop("'standardise' must be 'logical'")

  #:-------------------------------------------------------------------
  ## check which columns are continuous or categorical
  continuous <- dplyr::select(dataset, tidyselect::vars_select_helpers$where(is.numeric))
  discrete <- dplyr::select(dataset, tidyselect::vars_select_helpers$where(Negate(is.numeric)))
  dataset <- cbind(continuous, discrete)
  
  #:-------------------------------------------------------------------
  ## check if columns should be standardised
  if (standardise == TRUE) {
    if (ncol(continuous) != 0) {
      mean_values <- as.vector(continuous %>% colMeans(na.rm = TRUE)) %>%
        as.data.frame() %>%
        t()
      colnames(mean_values) <- colnames(continuous)
      sd_values <- as.vector(colSds(as.matrix(continuous), na.rm = TRUE)) %>%
        as.data.frame() %>%
        t()
      colnames(sd_values) <- colnames(continuous)
      continuous <- dplyr::select(dataset, tidyselect::vars_select_helpers$where(is.numeric)) %>%
        apply(2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    }
  }
  
  #:-------------------------------------------------------------------
  ## set up continuous and categorical variables
  if (ncol(continuous) != 0) {
    continuous <- as.matrix(continuous)
  }
  if (ncol(discrete) != 0) {
    discrete <- mutate(discrete, across(everything(), as.character)) %>%
      mutate(across(everything(), factor)) %>%
      mutate(across(everything(), as.numeric))
    discrete <- as.matrix(discrete)
  }

  
  #:-------------------------------------------------------------------
  ## if the dataset includes both continuous and categorical variables
  if (ncol(continuous) > 0 & ncol(discrete) >0) {
    
    
    #:-------------------------------------------------------------------
    ## check which rows are complete and which have missingness
    rows_complete <- which(complete.cases(dataset))
    rows_incomplete <- which(!complete.cases(dataset))

    if (!is_empty(rows_incomplete)) {
      
      #:-------------------------------------------------------------------
      ## If there are rows with incomplete data
      
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


      
      #:-------------------------------------------------------------------
      ## If there are rows with incomplete continuous and categorical predictors
      
      if (!is_empty(rows_continuous_incomplete) & !is_empty(rows_discrete_incomplete)) {

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
          
          # calculate log-density
          logDens ~ dnorm(0, 1)    ## this distribution does not matter
          
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
            x_disc_miss = x_disc_miss1,
            logDens = 0
          )
          inits
        }

        # adding this so that the logProb isn't -Inf
        logProb = "-Inf"

        while (logProb == "-Inf") {

          model <- nimbleModel(
            code = code,
            constants = consts,
            data = data,
            inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim, data$x_cont_miss, data$x_disc_miss)
          )

          logProb = model$calculate()

        }

        #compile the model
        cmodel <- compileNimble(model)

        #set monitors
        config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL", "x_cont_miss", "x_disc_miss", "logDens"), thin = 1, print = FALSE)

        ## add custom sampler
        for(i in 1:nrow(data$x_cont_miss)) {

          if(sum(is.na(data$x_cont_miss[i, ])) == 1) {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'sampler_conditional_RW',
              control = list(scale = 1, index = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          } else if (sum(is.na(data$x_cont_miss[i, ])) == ncol(continuous) || sum(is.na(data$x_cont_miss[i, ])) == 0) {
            #nothing happens because we just want draws from the component
          } else {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'sampler_conditional_RW_block',
              control = list(scale = 1, indices = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          }
        }
        
        #:-------------------------------------------------------------------
        ## If there are rows with incomplete continuous
      } else {
        if (!is_empty(rows_continuous_incomplete)) {

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
            
            # calculate log-density
            logDens ~ dnorm(0, 1)    ## this distribution does not matter
            
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
              x_cont_miss = x_cont_miss1,
              logDens = 0
            )
            inits
          }


          # adding this so that the logProb isn't -Inf
          logProb = "-Inf"

          while (logProb == "-Inf") {

            model <- nimbleModel(
              code = code,
              constants = consts,
              data = data,
              inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim, data$x_cont_miss)
            )

            logProb = model$calculate()

          }

          #compile the model
          cmodel <- compileNimble(model)

          #set monitors
          config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL", "x_cont_miss", "logDens"), thin = 1, print = FALSE)


          ## add custom sampler
          for(i in 1:nrow(data$x_cont_miss)) {

            if(sum(is.na(data$x_cont_miss[i, ])) == 1) {
              target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
              config$removeSampler(target)
              config$addSampler(
                target = target,
                type = 'sampler_conditional_RW',
                control = list(scale = 1, index = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
              )
            } else if (sum(is.na(data$x_cont_miss[i, ])) == ncol(continuous) || sum(is.na(data$x_cont_miss[i, ])) == 0) {
              #nothing happens because we just want draws from the component
            } else {
              target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
              config$removeSampler(target)
              config$addSampler(
                target = target,
                type = 'sampler_conditional_RW_block',
                control = list(scale = 1, indices = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
              )
            }
          }
          
          #:-------------------------------------------------------------------
          ## If there are rows with incomplete categorical predictors
          
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
            
            # calculate log-density
            logDens ~ dnorm(0, 1)    ## this distribution does not matter
            
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
              x_disc_miss = x_disc_miss1,
              logDens = 0
            )
            inits
          }


          # adding this so that the logProb isn't -Inf
          logProb = "-Inf"

          while (logProb == "-Inf") {

            model <- nimbleModel(
              code = code,
              constants = consts,
              data = data,
              inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim, data$x_disc_miss)
            )

            logProb = model$calculate()

          }

          #compile the model
          cmodel <- compileNimble(model)

          #set monitors
          config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL", "x_disc_miss", "logDens"), thin = 1, print = FALSE)

          
        }
      }
      
      #:-------------------------------------------------------------------
      ## If there are only complete rows with continuous and categorical predictors
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
        
        # calculate log-density
        logDens ~ dnorm(0, 1)    ## this distribution does not matter
        
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
          phiL = phiL,
          logDens = 0
        )
        inits
      }

      # adding this so that the logProb isn't -Inf
      logProb = "-Inf"

      while (logProb == "-Inf") {

        model <- nimbleModel(
          code = code,
          constants = consts,
          data = data,
          inits = initFn(consts$L, consts$N, data$mu0, data$tau0, data$R0, data$kappa0, consts$ndiscdim)
        )

        logProb = model$calculate()

      }

      #compile the model
      cmodel <- compileNimble(model)

      #set monitors
      config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "phiL", "logDens"), thin = 1, print = FALSE)
      
    }

  } else {
    
    #:-------------------------------------------------------------------
    ## if the dataset includes only continuous
    
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

      
      #:-------------------------------------------------------------------
      ## If there are rows with incomplete continuous predictors
      
      if (nrow(continuous_incomplete) != 0) {
        # if some missing values

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
          
          # calculate log-density
          logDens ~ dnorm(0, 1)    ## this distribution does not matter
          
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
            x_cont_miss = x_cont_miss1,
            logDens = 0
          )
          inits
        }


        # adding this so that the logProb isn't -Inf
        logProb = "-Inf"

        while (logProb == "-Inf") {

          model <- nimbleModel(
            code = code,
            constants = consts,
            data = data,
            inits = initFn(consts$L, consts$N, consts$Nmiss, data$mu0, data$tau0, data$R0, data$kappa0, data$x_cont_miss)
          )

          logProb = model$calculate()

        }



        #compile the model
        cmodel <- compileNimble(model)

        #set monitors
        config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "x_cont_miss", "logDens"), thin = 1, print = FALSE)


        ## add custom sampler
        for(i in 1:nrow(data$x_cont_miss)) {

          if(sum(is.na(data$x_cont_miss[i, ])) == 1) {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'sampler_conditional_RW',
              control = list(scale = 1, index = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          } else if (sum(is.na(data$x_cont_miss[i, ])) == ncol(continuous) || sum(is.na(data$x_cont_miss[i, ])) == 0) {
            #nothing happens because we just want draws from the component
          } else {
            target = paste0("x_cont_miss[", i, ", 1:", ncol(data$x_cont_miss), "]")
            config$removeSampler(target)
            config$addSampler(
              target = target,
              type = 'sampler_conditional_RW_block',
              control = list(scale = 1, indices = which(is.na(data$x_cont_miss[i, ])), adapt = TRUE)
            )
          }
        }
        
        #:-------------------------------------------------------------------
        ## If there are only complete rows of continuous predictors
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
          
          # calculate log-density
          logDens ~ dnorm(0, 1)    ## this distribution does not matter
          
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
            tauL = tauL,
            logDens = 0
          )
          inits
        }


        # adding this so that the logProb isn't -Inf
        logProb = "-Inf"

        while (logProb == "-Inf") {

          model <- nimbleModel(
            code = code,
            constants = consts,
            data = data,
            inits = initFn(consts$L, consts$N, data$mu0, data$tau0, data$R0, data$kappa0)
          )

          logProb = model$calculate()

        }


        #compile the model
        cmodel <- compileNimble(model)

        #set monitors
        config <- configureMCMC(cmodel, monitors = c("muL", "tauL", "v","z", "alpha", "logDens"), thin = 1, print = FALSE)
        
      }

      
      #:-------------------------------------------------------------------
      ## if the dataset includes only categorical
      
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

      
      #:-------------------------------------------------------------------
      ## If there are rows with incomplete categorical predictors
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
          
          # calculate log-density
          logDens ~ dnorm(0, 1)    ## this distribution does not matter
          
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
            x_disc_miss = x_disc_miss1,
            logDens = 0
          )
          inits
        }

        # adding this so that the logProb isn't -Inf
        logProb = "-Inf"

        while (logProb == "-Inf") {

          model <- nimbleModel(
            code = code,
            constants = consts,
            data = data,
            inits = initFn(consts$L, consts$N, consts$Nmiss, consts$ndiscdim, data$x_disc_miss)
          )

          logProb = model$calculate()

        }

        #compile the model
        cmodel <- compileNimble(model)

        #set monitors
        config <- configureMCMC(cmodel, monitors = c("v","z", "alpha", "phiL", "x_disc_miss", "logDens"), thin = 1, print = FALSE)
        
        
        #:-------------------------------------------------------------------
        ## If there are only rows with complete categorical predictors
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
          
          # calculate log-density
          logDens ~ dnorm(0, 1)    ## this distribution does not matter
          
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
            phiL = phiL,
            logDens = 0
          )
          inits
        }

        # adding this so that the logProb isn't -Inf
        logProb = "-Inf"

        while (logProb == "-Inf") {

          model <- nimbleModel(
            code = code,
            constants = consts,
            data = data,
            inits = initFn(consts$L, consts$N, consts$ndiscdim)
          )

          logProb = model$calculate()

        }


        #compile the model
        cmodel <- compileNimble(model)

        #set monitors
        config <- configureMCMC(cmodel, monitors = c("v","z", "alpha", "phiL", "logDens"), thin = 1, print = FALSE)
        
      }
    }

  }
  
  ## add custom sampler for log density
  config$removeSamplers('logDens')   ## remove sampler assigned to 'logDens'
  config$addSampler(target = 'logDens', type = 'sumLogPostDens')   ## add our custom sampler
  
  ## print config
  print(config)

  # build the model
  built <- buildMCMC(config)
  # compile model
  cbuilt <- compileNimble(built)

  # run model
  run <- runMCMC(cbuilt,
                 niter = mcmc_iterations,
                 nburnin = 0,
                 thin = thinning,
                 nchains = mcmc_chains,
                 progressBar = TRUE,
                 summary = FALSE,
                 samplesAsCodaMCMC = TRUE)

  # collect samples
  samples <- run

  # create output
  if (standardise == TRUE) {
    if (ncol(continuous) != 0) {
      output <- list(dataset = synthpop::syn(dataset, k = 1, print.flag = FALSE)$syn, L = L, mcmc_chains = mcmc_chains, thinning = thinning, samples = samples, standardise = standardise, mean_values = mean_values, sd_values = sd_values)
    } else {
      output <- list(dataset = synthpop::syn(dataset, k = 1, print.flag = FALSE)$syn, L = L, mcmc_chains = mcmc_chains, thinning = thinning, samples = samples, standardise = standardise)
    }
  } else {
    output <- list(dataset = synthpop::syn(dataset, k = 1, print.flag = FALSE)$syn, L = L, mcmc_chains = mcmc_chains, thinning = thinning, samples = samples, standardise = standardise)
  }

  class(output) <- "dpmm_fit"

  return(output)
}

