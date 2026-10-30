# Shared verbatim implementation of EM_algorithm() used by 2 simulation scenarios.
# Extracted without rewriting the function body.

EM_algorithm <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 100, tol = 1e-5, verbose = TRUE) {
  iter <- 0; diff <- Inf; rel_diff <- Inf
  
  # I-Step: Initial Imputation
  data_imputed <- Impute_step(data, beta_m_0, S, LOD)
  halmod_pred <- mod_update_hal(data = data_imputed)
  halbasis <- halmod_pred$basis_list
  wet_step <- Wet_step(data_imputed = data_imputed, mod_pred = halmod_pred, beta_m_0 = beta_m_0)
  data_weted <- wet_step$wet_dat
  
  start_time <- Sys.time()  # Track total runtime
  
  while (iter < max_iter && diff > tol && rel_diff > tol) {
    iter_time <- Sys.time()
    
    # W-Step
    wet_step <- Wet_step(data_imputed = data_weted, mod_pred = halmod_pred, beta_m_0 = beta_m_0)
    
    data_weted <- wet_step$wet_dat
    errors <- wet_step$dens_errors
    
    # M-Step
    errors_new <- compute_f_M_HAL(data = data_weted)$loss
    
    # Compute Convergence Criteria
    diff <- sum((errors_new - errors)^2)
    
    out_coef_old <- halmod_pred$out_coef
    
    # Parameter Update
    halmod_pred <- mod_update_hal(data_weted, basis_list = halbasis)
    
    # Update diff with coefficient changes
    #diff <- diff + sum((halmod_pred$out_coef - out_coef_old)^2)
    
    # Fixed denominator issue in rel_diff calculation
    # rel_diff <- diff / (sum(errors^2) + sum(out_coef_old^2) + 1e-8)  # Prevent division by zero
    rel_diff <- diff / (sum(errors^2) + 1e-8)  # Prevent division by zero
    
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
  
  return(list(data_weted = data_weted))
}
