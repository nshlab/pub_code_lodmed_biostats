# Shared verbatim implementation of sample_density_proposal() used by 4 simulation scenarios.
# Extracted without rewriting the function body.

sample_density_proposal <- function(n, data, LOD, beta_m) {
  # Linear predictor for L variables
  meanL <- as.matrix(data[c("L1", "L2", "L3")]) %*% beta_m[c("L1", "L2", "L3")]
  # Mean calculation for the normal distribution
  means <- beta_m["A"] * data$A + meanL + beta_m["b0"] + beta_m["AL1"]*data$A*data$L1
  # Initialize an empty vector to store valid samples
  valid_samples <- numeric(0)
  
  # Log-normal sampling loop until we have enough valid samples
  while (length(valid_samples) < n) {
    # Draw more samples to fill the request
    remaining <- n - length(valid_samples)
    # Draw samples from a log-normal distribution
    samples <- rlnorm(remaining, meanlog = means, sdlog = beta_m["sd"])
    # Append valid samples (less than LOD)
    valid_samples <- c(valid_samples, samples[samples <= LOD])
  }
  
  # Return only the first n valid samples
  return(valid_samples[1:n])
}
