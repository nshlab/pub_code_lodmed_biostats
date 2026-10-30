# Shared verbatim implementation of choose_ci() used by 2 simulation scenarios.
# Extracted without rewriting the function body.

choose_ci <- function(data, p_cen = 0.5, k = 5, gamma_seq = seq(0, 30, length.out = 200), B1 = 150, B2 = 100, nominal_coverage = 0.95,
                      type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  n <- nrow(data)  # Original sample size
  original_estimate <- compute_eff(data, 1:n, type = type, measure = measure)  # Estimate from original data
  
  # Initialize variables
  coverage_rates <- numeric(length(gamma_seq))
  optimal_c <- NA
  opt_m <- NA
  early_stop <- FALSE  # Flag to indicate early stopping
  
  # Loop over gamma_seq manually to allow early stopping
  for (i in seq_along(gamma_seq)) {
    if (early_stop) break  # Exit the loop if early stop is triggered
    
    gamma <- gamma_seq[i]
    
    # Compute c_i
    c_i <- (1 + gamma * exp(- k*p_cen)) / (1 + gamma)
    m_i <- floor(n^c_i)  # Compute m_i
    
    # Parallelize the first-stage bootstrap (B1 iterations)
    coverage_rate <- future_sapply(1:B1, function(b1) {
      first_stage_indices <- sample(1:n, n, replace = TRUE)  # Resample data
      #first_stage_estimate <- compute_eff(data, first_stage_indices, type = type, measure = measure)
      
      # Second-stage (nested) bootstrap
      nested_res <- gcom_boots_ci_m_out_n(data[first_stage_indices, ], m_i, n_bootstrap = B2, type = type, measure = measure)
      nested_ci <- nested_res$ci
      
      # Check if the original estimate falls within the CI
      if (original_estimate >= nested_ci[1] && original_estimate <= nested_ci[2]) 1 else 0
    }, future.seed = TRUE)  # Ensure reproducibility
    
    # Calculate the coverage rate for this gamma
    coverage_rates[i] <- mean(coverage_rate)
    
    # Check if the coverage rate exceeds the nominal coverage
    if (coverage_rates[i] >= nominal_coverage) {
      optimal_c <- c_i
      opt_m <- m_i
      early_stop <- TRUE  # Trigger early stopping
      cat("Stopping early at gamma =", gamma, "with coverage rate =", coverage_rates[i], "\n")
    }
    
    cat("Iteration", i, "gamma =", gamma, "coverage rate =", coverage_rates[i], "\n")
  }
  
  # If no gamma achieved the nominal coverage, select the one with the highest coverage rate
  if (is.na(optimal_c)) {
    opt_gamma <- gamma_seq[which.max(coverage_rates)]
    optimal_c <- (1 + opt_gamma * exp(- k*p_cen)) / (1 + opt_gamma)
    opt_m <- floor(n^optimal_c)
    cat("No gamma achieved the nominal coverage. Selecting gamma with the highest coverage rate:", optimal_c, "\n")
  }
  
  return(list(
    optimal_c = optimal_c,
    coverage_rates = data.frame(c = gamma_seq, coverage_rate = coverage_rates),
    opt_m = opt_m
  ))
}
