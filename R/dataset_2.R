#' Synthetic dataset 2
#'
#' A synthetic dataset comprised of two clusters (1:500, 501:1000), where X1 values overlap for both clusters, X2 values do not overlap for both clusters and X3 (categorical) is distinct between both clusters. Both X2 and X3 are informative of the cluster (not X1).
#'
#' @format ## `dataset_2`
#' A data frame with 1,000 rows and 3 columns:
#' \describe{
#'   \item{X1}{Continuous variable, mean = 8, sd = 1}
#'   \item{X2}{Continuous variable, mean = 2, sd = 1; mean = 9, sd = 1}
#'   \item{X3}{Categorical variable, levels = A or B, distinct for each cluster}
#' }
#' @source Synthetic
"dataset_2"
