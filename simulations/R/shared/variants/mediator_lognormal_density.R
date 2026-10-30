# Shared verbatim implementation of compute_f_M() used by 4 simulation scenarios.
# Extracted without rewriting the function body.

compute_f_M <- function(data, beta_m) {
  # Calculate linear predictors for the mean of the normal distribution
  meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
  # Total predicted mean
  means <- beta_m["A"] * data$A + meanL + beta_m["b0"] + beta_m["AL1"]*data$A*data$L1
  
  # Compute the log-normal density for each M
  # Convert the mean to a log scale for log-normal parameters
  f_M <- dlnorm(data$M, meanlog = means, sdlog = beta_m["sd"])
  return(f_M)
}
