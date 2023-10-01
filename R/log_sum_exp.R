log_sum_exp <- function(lx) {
  
  ## extract maximum of logged values
  mX <- max(lx, na.rm = TRUE)
  
  ## return answer
  out <- mX + log(sum(exp(lx[!is.na(lx)] - mX)))
  out
}