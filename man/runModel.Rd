% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runModel.R
\name{runModel}
\alias{runModel}
\title{Fit DPMM}
\usage{
runModel(
  dataset,
  mcmc_iterations = 2500,
  thinning = 1,
  L = 10,
  mcmc_chains = 2,
  standardise = TRUE
)
}
\arguments{
\item{dataset}{dataframe with continuous and/or categorical variables.}

\item{mcmc_iterations}{Number of MCMC iterations.}

\item{thinning}{Interval for collecting MCMC samples.}

\item{L}{Number of DPMM components fitted.}

\item{mcmc_chains}{Number of MCMC chains fitted.}

\item{standardise}{Should continuous variables be standardised. default = TRUE}
}
\value{
\subsection{Output: List of class 'dpmm_fit'}{
   \describe{
     \item{dataset}{dataframe with 1 row but defined the same way as the dataset fitted.}
     \item{L}{Number of DPMM components fitted.}
     \item{mcmc_chains}{Number of MCMC chains fitted.}
     \item{samples}{MCMC samples.}
     \item{standardise}{TRUE or FALSE whether values were standardised.}
     \item{mean_values}{Mean values for covariates standardised.}
     \item{sd_values}{SD values for covariates standardised.}
   }
}
}
\description{
This function fits the Dirichlet Process Mixture Model (DPMM).
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
}


}
