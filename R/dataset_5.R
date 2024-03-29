#' Synthetic dataset 5
#'
#' A synthetic dataset comprised of two clusters (1:500, 501:1000), where X1 values overlap for both clusters, X2 values do not overlap for both clusters, X3 (categorical) is distinct between both clusters, and X4 (categorical) has 4 levels and each cluster has two levels. Both X2/X3/X4 are informative of the cluster (not X1).
#'
#' @format `dataset_5`
#' A data frame with 1,000 rows and 3 columns:
#' \describe{
#'   \item{X1}{Continuous variable, mean = 8, sd = 1}
#'   \item{X2}{Continuous variable, mean = 2, sd = 1; mean = 9, sd = 1}
#'   \item{X3}{Categorical variable, levels = A or B, distinct for each cluster}
#'   \item{X4}{Categorical variable, levels = 1/2 or 3/4, distinct for each cluster}
#' }
#' @source Synthetic
"dataset_5"
