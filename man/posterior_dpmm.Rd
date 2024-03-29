% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_dpmm.R
\name{posterior_dpmm}
\alias{posterior_dpmm}
\title{Sample from DPMM}
\usage{
posterior_dpmm(
  patient,
  samples,
  seed = NULL,
  cont_vars = NULL,
  cat_vars = NULL,
  mcmc_chain = NULL
)
}
\arguments{
\item{patient}{dataset with missing values.}

\item{samples}{vector of iterations to be used in the posterior.}

\item{seed}{specify seed to be used. (default = NULL)}

\item{cont_vars}{names of continuous variables}

\item{cat_vars}{names of categorical variables}

\item{mcmc_chain}{MCMC posterior samples of DPMM parameters}
}
\value{
A list with n entries of n rows with missingness
}
\description{
Sample missing values from the fitted DPMM
}
\examples{
\dontrun{
## load dataset
data(dataset_1)

## fit model
posteriors <- runModel(dataset_1, 
                       mcmc_iterations = 100,
                       L = 6, 
                       mcmc_chains = 2, 
                       standardise = TRUE)
                       
## introduce missing data
rows <- 501:550
dataset_missing <- dataset_1
dataset_missing_predict <- dataset_missing[rows,]
dataset_missing_predict[,1] <- as.numeric(NA)

# predict missing values
posteriors.dpmmfit <- predict_dpmm_fit(posteriors, 
                                       dataset_missing_predict, 
                                       samples = c(1:100))
}


}
