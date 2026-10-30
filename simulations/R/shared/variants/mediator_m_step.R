# Shared verbatim implementation of M_stepf_M() used by 2 simulation scenarios.
# Extracted without rewriting the function body.

M_stepf_M <- function(data, beta_m) {
  # Define an internal function for optimization
  neg_log_likelihood <- function(beta_m, data) {
    # Extract weights
    weights <- data$w_norm
    # Calculate the mean for each sample using provided beta_m
    meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
    means <- beta_m["A"] * data$A + meanL + beta_m["b0"] + beta_m["AL1"]*data$A*data$L1
    
    # Compute log-likelihood for each data point using the truncated normal distribution
    # Log transformation is done using the log() function
    log_likelihoods <- log(dlnorm(data$M, meanlog = means, sdlog = beta_m["sd"]))
    
    # Incorporate weights into the likelihood calculation
    weighted_log_likelihood <- weights * log_likelihoods
    
    # Return the negative of the weighted log-likelihood sum
    return(-sum(weighted_log_likelihood))
  }
  
  optim_results <- optim(par = beta_m, fn = neg_log_likelihood, data = data)
  
  # Results
  return(beta_m_new = optim_results$par)
}
