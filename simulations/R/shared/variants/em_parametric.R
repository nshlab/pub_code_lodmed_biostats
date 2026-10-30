# Shared verbatim implementation of EM_algorithm() used by 2 simulation scenarios.
# Extracted without rewriting the function body.

EM_algorithm <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 50, tol = 1e-20, verbose = TRUE) {
  iter <- 0; diff <- Inf; rel_diff <- Inf; beta_m <- beta_m_0 
  
  data_imputed <- Impute_step(data, beta_m_0, S, LOD)
  mod_pred <- mod_update_glm(data_imputed)
  
  start_time <- Sys.time()  # Track runtime
  
  while (iter < max_iter && diff > tol && rel_diff > tol) {
    iter_time <- Sys.time()
    
    # W-Step
    data_weted <- Wet_step(data_imputed, mod_pred, beta_m_new = beta_m, beta_m_0 = beta_m_0)
    
    # Maximization Step
    beta_m_new <- M_stepf_M(data = data_weted, beta_m)
    
    # Compute convergence criteria
    diff <- sum((beta_m_new - beta_m)^2)
    
    cen_coef_old <- mod_pred$cen_coef
    out_coef_old <- mod_pred$out_coef
    
    # Parameter update
    beta_m <- beta_m_new
    mod_pred <- mod_update_glm(data_weted)
    
    # Update diff with coefficient changes
    diff <- diff + sum((mod_pred$out_coef - out_coef_old)^2)
    
    # Fixed denominator issue in rel_diff calculation
    rel_diff <- diff / (sum(beta_m^2) + sum(out_coef_old^2) + 1e-8)  # Prevent division by zero
    # rel_diff <- diff / (sum(beta_m^2) + 1e-8)  # Prevent division by zero
    iter_time_end <- Sys.time()
    iter_duration <- round(difftime(iter_time_end, iter_time, units = "secs"), 4)
    
    if (verbose) {
      cat("Iteration:", iter, 
          "| Absolute Diff:", diff, 
          "| Relative Diff:", rel_diff,
          "| Time:", iter_duration, "s\n")
    }
    
    iter <- iter + 1
  }
  
  total_runtime <- round(difftime(Sys.time(), start_time, units = "secs"), 4)
  
  # Stopping condition logging
  if (iter >= max_iter) {
    message("EM Algorithm stopped: Maximum iterations reached (", max_iter, ").")
  } else if (diff <= tol || rel_diff <= tol) {
    message("EM Algorithm converged successfully in ", iter, " iterations.")
  }
  
  cat("Total Runtime:", total_runtime, "seconds\n")
  
  return(list(beta_m = beta_m, data_weted = data_weted))
}
