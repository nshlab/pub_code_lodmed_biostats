# AUTO-REFACTORED FROM ORIGINAL APPLICATION CODE.
# Function bodies are preserved unless explicitly noted.

suppressPackageStartupMessages(library(future.apply))

gcom_boots_ci <- function(data, n_bootstrap = 1000, alpha = 0.05, 
                          type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  # Parallel bootstrap estimation using future_lapply()
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = nrow(data), replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    # Extract the corresponding estimate
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    
    return(NA)  # Handle unexpected measure values safely
  })
  
  # Convert to numeric vector and remove NAs
  boots_est <- na.omit(unlist(boots_est))
  
  # Check for empty bootstrap results
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA. Check data or computation methods.")
    return(NULL)
  }
  
  # Compute confidence intervals
  ci_lower <- quantile(boots_est, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boots_est, probs = 1 - alpha / 2, na.rm = TRUE)
  
  list(
    point_estimate = mean(boots_est, na.rm = TRUE),
    variance = var(boots_est, na.rm = TRUE),
    ci = c(lower = ci_lower, upper = ci_upper)
  )
}

gcom_boots_ci_m_out_n <- function(data, m, n_bootstrap = 1000, alpha = 0.05, 
                                  type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  boots_est <- numeric(n_bootstrap)
  
  # Parallel bootstrap estimation using future_lapply
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = m, replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    # Extract the corresponding estimate
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    
    return(NA)  # Handle unexpected measure values safely
  }, future.seed = TRUE)  # Ensure reproducibility
  
  boots_est <- na.omit(unlist(boots_est))  # Flatten and remove NAs
  
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA. Check data or computation methods.")
    return(NULL)
  }
  
  # Scale variance by (n/m) for m-out-of-n bootstrap correction
  point_estimate <- mean(boots_est, na.rm = TRUE)
  scaled_variance <- var(boots_est, na.rm = TRUE) * (m / nrow(data))
  ci <- quantile(boots_est, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  list(
    point_estimate = point_estimate,
    variance = scaled_variance,
    ci = ci
  )
}

compute_eff <- function(data, indices, type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  sample_data <- data[indices, ]  # Resample data
  bres <- if (type == "NDE") gcomp.nde(sample_data) else gcomp.nie(sample_data)
  
  # Extract the corresponding estimate
  if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
  if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
  if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
}

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

compute_se <- function(boots_est) {
  sd(boots_est, na.rm = TRUE)
}

compute_pvalue <- function(estimate, se, null_value = 0) {
  if (is.na(estimate) || is.na(se) || se == 0) {
    return(NA)
  }
  z_stat <- (estimate - null_value) / se
  p_value <- 2 * pnorm(-abs(z_stat))  # Two-tailed p-value
  return(p_value)
}

compute_ci_from_se <- function(estimate, se, alpha = 0.05) {
  z_crit <- qnorm(1 - alpha / 2)
  ci_lower <- estimate - z_crit * se
  ci_upper <- estimate + z_crit * se
  return(c(lower = ci_lower, upper = ci_upper))
}

gcom_boots_ci_with_se <- function(data, n_bootstrap = 1000, alpha = 0.05, 
                                  type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  # Parallel bootstrap estimation
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = nrow(data), replace = TRUE), ]
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    return(NA)
  })
  
  boots_est <- na.omit(unlist(boots_est))
  
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA.")
    return(NULL)
  }
  
  point_estimate <- mean(boots_est, na.rm = TRUE)
  se <- compute_se(boots_est)
  p_value <- compute_pvalue(point_estimate, se)
  ci_quantile <- quantile(boots_est, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  ci_normal <- compute_ci_from_se(point_estimate, se, alpha = alpha)
  
  list(
    point_estimate = point_estimate,
    se = se,
    variance = var(boots_est, na.rm = TRUE),
    p_value = p_value,
    ci_quantile = ci_quantile,
    ci_normal = ci_normal,
    boots_est = boots_est
  )
}

gcom_boots_ci_m_out_n_with_se <- function(data, m, n_bootstrap = 1000, alpha = 0.05, 
                                          type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = m, replace = TRUE), ]
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    return(NA)
  }, future.seed = TRUE)
  
  boots_est <- na.omit(unlist(boots_est))
  
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA.")
    return(NULL)
  }
  
  # Compute statistics with m-out-of-n correction
  point_estimate <- mean(boots_est, na.rm = TRUE)
  se <- compute_se(boots_est)
  scaled_variance <- var(boots_est, na.rm = TRUE) * (m / nrow(data))
  scaled_se <- sqrt(scaled_variance)
  p_value <- compute_pvalue(point_estimate, scaled_se)
  ci_quantile <- quantile(boots_est, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  ci_normal <- compute_ci_from_se(point_estimate, scaled_se, alpha = alpha)
  
  list(
    point_estimate = point_estimate,
    se = scaled_se,
    variance = scaled_variance,
    p_value = p_value,
    ci_quantile = ci_quantile,
    ci_normal = ci_normal,
    boots_est = boots_est
  )
}

# medoutcon-specific m-out-of-n helpers
moutn_multiplier_bootstrap <- function(eif, theta, B = 200, m = 100,  multipliers = "rademacher", ci_level = 0.95) {
  n <- length(eif)
  m <- min(m, n)  # Ensure m does not exceed n
  
  # Remove NAs from eif
  eif <- eif[!is.na(eif)]
  
  # Handle edge case where eif is empty or has only one value
  if (length(eif) < 2) {
    return(list(
      m_process = NA,
      lwr_ci = NA,
      upr_ci = NA,
      ci = c(NA, NA),
      c_alpha = NA
    ))
  }
  
  samples <- sample(1:length(eif), size = min(m, length(eif)), replace = FALSE)
  eif <- eif[samples]
  
  m_process <- numeric(B)
  se_eif <- sqrt(stats::var(eif) / n)
  
  # Handle edge case where se_eif is 0, NA, or very small
  if (is.na(se_eif) || se_eif < 1e-10) {
    se_eif <- 1e-10
  }
  
  for (b in 1:B) {
    # Step 1: Draw multipliers
    if (multipliers == "rademacher") {
      xi <- (2 * stats::rbinom(length(eif), 1, 0.5)) - 1
    } else if (multipliers == "gaussian") {
      xi <- rnorm(length(eif))
    }
    
    # Step 2: compute process M using multipliers
    m_process[b] <- mean((xi * eif) / se_eif)
  }
  
  # Remove any NAs before computing quantile
  m_process <- m_process[!is.na(m_process)]
  
  c_alpha <- unname(stats::quantile(m_process, ci_level, na.rm = TRUE))
  
  lwr_ci <- theta - c_alpha * se_eif
  upr_ci <- theta + c_alpha * se_eif
  return(list(m_process = m_process,
              lwr_ci = lwr_ci,
              upr_ci = upr_ci,
              se = se_eif,
              ci = c(lwr_ci, upr_ci),
              c_alpha = c_alpha))
}

os_compute_eff <- function(data, indices, type = c("NDE", "NIE"), fast = FALSE) {
  type <- match.arg(type)
  sample_data <- data[indices, ]  # Resample data
  
  # Choose function based on fast argument
  if (fast) {
    bres <- Onestep_med_effs_fast(dat_aug = sample_data, effect = type)
  } else {
    bres <- Onestep_med_effs(dat_aug = sample_data, effect = type)
  }
  
  # Extract estimate and EIF
  est_name <- if (type == "NDE") "nde_onestep_est" else "nie_onestep_est"
  
  return(list(
    eif = bres$eif,          # Assumes eif is always returned regardless of effect
    est = bres[[est_name]]
  ))
}

os_choose_ci <- function(data, original_estimate, p_cen = 0.5, k = 5, gamma_seq = seq(0.001, 30, length.out = 100), B1 = 25, B2 = 50000, 
                         nominal_coverage = 0.95, type = c("NDE", "NIE")) {
  n <- nrow(data)  # Original sample size
  
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
    m_i <- min(m_i, n)   # Ensure m_i does not exceed n
    
    # Parallelize the first-stage bootstrap (B1 iterations)
    coverage_rate <- future_sapply(1:B1, function(b1) {
      first_stage_indices <- sample(1:n, n, replace = TRUE)  # Resample data
      #first_stage_estimate <- compute_eff(data, first_stage_indices, type = type, measure = measure)
      
      # Second-stage (nested) bootstrap
      os_eff <- os_compute_eff(data = data, indices = first_stage_indices, type = type, fast = TRUE)
      nested_res <- moutn_multiplier_bootstrap(eif = os_eff$eif, theta = os_eff$est, B = B2, m = m_i)
      nested_ci <- nested_res$ci
      
      # Check if CI is valid (not NA) and if the original estimate falls within the CI
      if (any(is.na(nested_ci))) {
        return(0)  # Return 0 if CI is invalid
      } else if (original_estimate >= nested_ci[1] && original_estimate <= nested_ci[2]) {
        return(1)
      } else {
        return(0)
      }
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
    opt_m <- min(opt_m, n)  # Ensure opt_m does not exceed n
    cat("No gamma achieved the nominal coverage. Selecting gamma with the highest coverage rate:", optimal_c, "\n")
  }
  
  return(list(
    optimal_c = optimal_c,
    coverage_rates = data.frame(c = gamma_seq, coverage_rate = coverage_rates),
    opt_m = opt_m
  ))
}
