#' Vectorised draw of a mixture component per row
#'
#' Draws one component index per row of a log-likelihood matrix. Fully
#' vectorised (no per-row R-level loop): normalises each row with
#' 'rowLogSumExps', then locates the drawn component via a row-wise
#' cumulative-probability threshold instead of calling 'sample()' once per row.
#'
#' @param loglik a matrix (patients x components) of unnormalised log weights
#'   (log-likelihood + log prior weight per component). 'NA' entries (e.g. from
#'   a non-invertible covariance) are treated as zero probability.
#'
#' @return an integer vector, one drawn component index per row of 'loglik'.
#'
#' @keywords internal
#' @noRd
draw_component <- function(loglik) {
  loglik[is.na(loglik)] <- -Inf
  lmarg <- rowLogSumExps(loglik)
  probs <- exp(loglik - lmarg)
  cum_probs <- rowCumsums(probs)
  u <- runif(nrow(loglik))
  draw <- rowSums(cum_probs < u) + 1L
  pmin(draw, ncol(loglik))
}

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
#' @return A list with n entries of n rows with missingness, in the same row order as 'patient'.
#'
#' @import mvtnorm
#' @importFrom matrixStats rowLogSumExps
#' @importFrom matrixStats rowCumsums
#' @importFrom stats rnorm
#' @importFrom stats runif
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
#'                        standardise = TRUE,
#'                        seed = 1)
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
#'                                        samples = c(1:100),
#'                                        seed = 1)
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
    ## number of real categories per categorical variable. 'phiL' is stored
    ## in a rectangular array padded out to max(ndiscdim_vec) columns (nimble
    ## requires a fixed-size array across all L components/variables), so
    ## variables with fewer levels than the widest one have unused padding
    ## columns beyond their real level count -- these must never be read
    ## when drawing a new categorical value below.
    ndiscdim_vec <- vapply(cat_vars, function(v) length(unique(levels(patient[, v]))), integer(1))
  }

  n_cont <- length(cont_vars)
  n_cat <- length(cat_vars)

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

    ## w's columns are already in the correct 1..K+1 order by construction
    ## (stick_breaking() builds them in that order directly), so no column-name
    ## parsing is needed to recover component order -- just split rows to a list
    postW <- tibble(iter = seq_len(nrow(v_mat)), w = asplit(stick_breaking(v_mat), 1))

    #:----------------------------------------------------------
    ## If cont_vars is provided (and hence the model is fitted to continuous variables)
    if (!is.null(cont_vars)) {
      ### Extract posterior samples of component means
      mu_cols <- samples_chain %>%
        as_tibble(.name_repair = "unique") %>%
        select(starts_with("muL[")) %>%
        as.matrix()

      ## parse each column's (component, dim) indices once, instead of
      ## re-deriving them via a pivot_longer/group_by round trip per element
      mu_idx <- do.call(rbind, strsplit(gsub("muL\\[|\\]", "", colnames(mu_cols)), ",\\s*"))
      mode(mu_idx) <- "integer"
      n_iter <- nrow(mu_cols)
      L_mu <- max(mu_idx[, 1])
      d_mu <- max(mu_idx[, 2])

      muL_list <- vector("list", n_iter)
      for (it in seq_len(n_iter)) {
        m <- matrix(NA_real_, nrow = L_mu, ncol = d_mu)
        m[cbind(mu_idx[, 1], mu_idx[, 2])] <- mu_cols[it, ]
        muL_list[[it]] <- m
      }
      postMu <- tibble(iter = seq_len(n_iter), muL = muL_list)

      ### Extract posterior samples of component precision matrices, and
      ### immediately invert them to covariance matrices here -- once per
      ### (iteration, component) -- so that no solve() is ever repeated
      ### per patient/group further down.
      tau_cols <- samples_chain %>%
        as_tibble(.name_repair = "unique") %>%
        select(starts_with("tauL[")) %>%
        as.matrix()

      tau_idx <- do.call(rbind, strsplit(gsub("tauL\\[|\\]", "", colnames(tau_cols)), ",\\s*"))
      mode(tau_idx) <- "integer"
      n_iter <- nrow(tau_cols)
      d1_tau <- max(tau_idx[, 1])
      d2_tau <- max(tau_idx[, 2])
      L_tau <- max(tau_idx[, 3])

      covL_list <- vector("list", n_iter)
      for (it in seq_len(n_iter)) {
        a <- array(NA_real_, dim = c(d1_tau, d2_tau, L_tau))
        a[cbind(tau_idx[, 1], tau_idx[, 2], tau_idx[, 3])] <- tau_cols[it, ]
        cov_a <- array(NA_real_, dim = c(d1_tau, d2_tau, L_tau))
        for (l in seq_len(L_tau)) {
          ## symmetrise defensively (sampled precision matrices can carry tiny
          ## floating-point asymmetry) before inverting
          tau_l <- matrix(a[, , l], nrow = d1_tau, ncol = d2_tau)
          tau_l <- (tau_l + t(tau_l)) / 2
          cov_a[, , l] <- solve(tau_l)
        }
        covL_list[[it]] <- cov_a
      }
      postCov <- tibble(iter = seq_len(n_iter), covL = covL_list)

    }

    #:----------------------------------------------------------
    ## If cat_vars is provided (and hence the model is fitted to categorical variables)
    if (!is.null(cat_vars)) {
      ### Extract posterior samples of component level probabilities
      phi_cols <- samples_chain %>%
        as_tibble(.name_repair = "unique") %>%
        select(starts_with("phiL[")) %>%
        as.matrix()

      phi_idx <- do.call(rbind, strsplit(gsub("phiL\\[|\\]", "", colnames(phi_cols)), ",\\s*"))
      mode(phi_idx) <- "integer"
      n_iter <- nrow(phi_cols)
      d1_phi <- max(phi_idx[, 1])
      d2_phi <- max(phi_idx[, 2])
      L_phi <- max(phi_idx[, 3])

      phiL_list <- vector("list", n_iter)
      for (it in seq_len(n_iter)) {
        a <- array(NA_real_, dim = c(d1_phi, d2_phi, L_phi))
        a[cbind(phi_idx[, 1], phi_idx[, 2], phi_idx[, 3])] <- phi_cols[it, ]
        phiL_list[[it]] <- a
      }
      postPhi <- tibble(iter = seq_len(n_iter), phiL = phiL_list)

    }

    #:----------------------------------------------------------
    ## Organise data in format needed for the function

    ## rows being used, kept in the exact order supplied by the caller
    ordered_patient <- patient[, c(cont_vars, cat_vars), drop = FALSE]

    number_patients <- nrow(ordered_patient) ## number of patients
    number_samples <- nrow(samples_chain)  ## number of iterations

    ## numeric matrix version of the patient data (categorical columns
    ## recoded to their underlying integer factor levels, exactly as the
    ## model was fit on) -- used for fast vectorised indexing below
    patient_num <- do.call(cbind, lapply(ordered_patient, function(col) {
      if (is.factor(col)) as.numeric(col) else as.numeric(col)
    }))
    colnames(patient_num) <- colnames(ordered_patient)
    na_mask <- is.na(patient_num)

    #:----------------------------------------------------------
    ## Group patients by identical missingness pattern: every patient in a
    ## group shares the same given/missing variable indices, so the
    ## per-component marginal/conditional setup (the expensive part) is
    ## computed once per group per iteration and reused across every patient
    ## in that group, instead of once per patient.
    pattern_key <- apply(na_mask, 1, function(r) paste(which(r), collapse = ","))
    groups <- split(seq_len(number_patients), pattern_key)

    group_info <- lapply(groups, function(idx) {
      missing_idx <- which(na_mask[idx[1], ])
      given_idx <- which(!na_mask[idx[1], ])
      list(
        patients      = idx,
        missing_idx   = missing_idx,
        cont_given    = given_idx[given_idx <= n_cont],
        cont_dep      = missing_idx[missing_idx <= n_cont],
        cat_given_rel = given_idx[given_idx > n_cont] - n_cont,
        cat_dep_rel   = missing_idx[missing_idx > n_cont] - n_cont
      )
    })

    ## preallocate one result array per group: [iteration, missing variable, patient];
    ## every draw below fills into these pre-sized arrays instead of growing
    ## objects (rbind/cbind) inside the sampling loop
    group_results <- lapply(group_info, function(gi) {
      array(NA_real_, dim = c(number_samples, length(gi$missing_idx), length(gi$patients)))
    })

    #:----------------------------------------------------------
    ## Loop over posterior iterations once; everything that only depends on
    ## the iteration (component weights, means, covariances, category
    ## probabilities) is extracted a single time per iteration and shared
    ## across every patient/group processed within that iteration, rather
    ## than being re-derived per patient as before.
    for (it in seq_len(number_samples)) {

      w <- postW$w[[it]]
      L <- length(w)
      log_w <- log(w)

      if (n_cont > 0) {
        mu <- postMu$muL[[it]]
        covL <- postCov$covL[[it]]
        ## materialise each component's covariance matrix once per iteration
        ## (not once per group) -- reused for both the weight calculation and
        ## the conditional draw below
        Sigma_by_l <- vector("list", L)
        for (l in seq_len(L)) Sigma_by_l[[l]] <- matrix(covL[, , l], nrow = d1_tau, ncol = d2_tau)
      }
      if (n_cat > 0) {
        phi <- postPhi$phiL[[it]]
      }

      for (g in seq_along(group_info)) {

        gi <- group_info[[g]]
        pts <- gi$patients
        n_group <- length(pts)

        #:----------------------------------------------------------
        ## per-component log-likelihood for every patient in the group,
        ## vectorised across patients: one dmvnorm()/phi lookup per
        ## component covers the whole group, instead of one per patient.
        ll <- matrix(log_w, nrow = n_group, ncol = L, byrow = TRUE)

        if (length(gi$cont_given) > 0) {
          Xg <- patient_num[pts, gi$cont_given, drop = FALSE]
          for (l in seq_len(L)) {
            sigma_gg <- Sigma_by_l[[l]][gi$cont_given, gi$cont_given, drop = FALSE]
            ll[, l] <- ll[, l] + mvtnorm::dmvnorm(Xg, mean = mu[l, gi$cont_given], sigma = sigma_gg, log = TRUE)
          }
        }

        if (length(gi$cat_given_rel) > 0) {
          Xg_cat <- patient_num[pts, n_cont + gi$cat_given_rel, drop = FALSE]
          for (l in seq_len(L)) {
            for (j in seq_along(gi$cat_given_rel)) {
              ll[, l] <- ll[, l] + log(phi[gi$cat_given_rel[j], Xg_cat[, j], l])
            }
          }
        }

        draws <- draw_component(ll)

        #:----------------------------------------------------------
        ## draw the missing values, batched by chosen component: patients in
        ## the same group that land on the same component this iteration
        ## share the conditional covariance/Cholesky factor and are drawn in
        ## a single vectorised step (rather than one mvrnorm()/rcat() call
        ## per patient).
        for (l in unique(draws)) {

          bucket <- which(draws == l)  ## positions within the group
          n_bucket <- length(bucket)
          bucket_patients <- pts[bucket]

          if (length(gi$cont_dep) > 0) {
            mu_l <- mu[l, ]
            if (length(gi$cont_given) > 0) {
              sigma_dd <- Sigma_by_l[[l]][gi$cont_dep, gi$cont_dep, drop = FALSE]
              sigma_dg <- Sigma_by_l[[l]][gi$cont_dep, gi$cont_given, drop = FALSE]
              sigma_gg <- Sigma_by_l[[l]][gi$cont_given, gi$cont_given, drop = FALSE]
              B <- sigma_dg %*% solve(sigma_gg)
              cond_var <- sigma_dd - B %*% t(sigma_dg)

              Xg_bucket <- patient_num[bucket_patients, gi$cont_given, drop = FALSE]
              centred <- sweep(Xg_bucket, 2, mu_l[gi$cont_given], "-")
              cond_mean <- sweep(centred %*% t(B), 2, mu_l[gi$cont_dep], "+")
            } else {
              cond_var <- Sigma_by_l[[l]][gi$cont_dep, gi$cont_dep, drop = FALSE]
              cond_mean <- matrix(mu_l[gi$cont_dep], nrow = n_bucket, ncol = length(gi$cont_dep), byrow = TRUE)
            }
            cond_var <- (cond_var + t(cond_var)) / 2
            chol_var <- chol(cond_var)
            z <- matrix(rnorm(n_bucket * length(gi$cont_dep)), nrow = n_bucket)
            draws_cont <- cond_mean + z %*% chol_var

            var_pos <- match(gi$cont_dep, gi$missing_idx)
            group_results[[g]][it, var_pos, bucket] <- t(draws_cont)
          }

          if (length(gi$cat_dep_rel) > 0) {
            for (j in seq_along(gi$cat_dep_rel)) {
              ## only read the variable's real categories -- 'phi' is padded
              ## out to the widest categorical variable's level count, and
              ## those padding columns must never be sampled from (see note
              ## on 'ndiscdim_vec' above)
              n_levels <- ndiscdim_vec[gi$cat_dep_rel[j]]
              probs <- phi[gi$cat_dep_rel[j], seq_len(n_levels), l]
              vals <- sample.int(n_levels, size = n_bucket, replace = TRUE, prob = probs)
              var_pos <- match(n_cont + gi$cat_dep_rel[j], gi$missing_idx)
              group_results[[g]][it, var_pos, bucket] <- vals
            }
          }
        }
      }
    }

    #:----------------------------------------------------------
    ## Unpack each group's result array back into the per-patient output
    ## list, indexed by the patient's *original* row position -- so the
    ## returned list is always in the same order 'patient' was supplied in,
    ## regardless of how patients were grouped/reordered internally above.
    posterior <- vector(mode = "list", length = number_patients)
    for (g in seq_along(group_info)) {
      gi <- group_info[[g]]
      var_names <- colnames(ordered_patient)[gi$missing_idx]
      n_vars <- length(gi$missing_idx)
      for (p_local in seq_along(gi$patients)) {
        p <- gi$patients[p_local]
        mat <- matrix(group_results[[g]][, , p_local, drop = FALSE], nrow = number_samples, ncol = n_vars)
        colnames(mat) <- var_names
        posterior[[p]] <- as.data.frame(mat)
      }
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
