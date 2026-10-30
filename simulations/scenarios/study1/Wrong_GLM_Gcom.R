# Simulation scenario definition; settings preserved from the supplied repository.
# Scenario-specific function bodies are retained; only byte-identical shared helpers were centralized.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
source(here::here("simulations", "R", "shared", "dgp.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "simulation_helpers.R"), local = TRUE)
SCENARIO_ID <- "study1_Wrong_GLM_Gcom"
source(here::here("simulations", "R", "shared", "variants", "impute_numeric_indicator.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "onestep_mediation.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "weights_parametric.R"), local = TRUE)

library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(xgboost)
library(future.apply)


# NaiveImpute() moved to simulations/R/shared/simulation_helpers.R


#NaiveImpute(dataset, LOD = LOD_value)

###############################
## Density of M for samples below LOD
##
# Sampling function using Poisson distribution
sample_density_proposal <- function(n, data, LOD, beta_m) {
  # Linear predictor for L variables
  meanL <- as.matrix(data[c("L3")]) %*% beta_m[c("L3")]
  # Mean (rate) calculation for Poisson distribution
  lambda <- beta_m["A"] * data$A + meanL + beta_m["b0"]
  lambda <- pmax(lambda, 0)  # Ensure non-negative lambda
  
  # Initialize an empty vector to store valid samples
  valid_samples <- numeric(0)
  
  # Poisson sampling loop until we have enough valid samples
  while (length(valid_samples) < n) {
    # Draw more samples to fill the request
    remaining <- n - length(valid_samples)
    # Replicate lambda vector for multiple draws
    sampled_lambda <- rep(lambda, length.out = remaining)
    # Draw samples from Poisson
    samples <- rpois(remaining, lambda = sampled_lambda)
    # Append valid samples (less than or equal to LOD)
    valid_samples <- c(valid_samples, samples[samples <= LOD])
  }
  
  # Return only the first n valid samples
  return(valid_samples[1:n])
}

# Density function using Poisson distribution
compute_f_M <- function(data, beta_m) {
  # Calculate lambda (rate parameter) for the Poisson distribution
  meanL <- as.matrix(data[, .(L3)]) %*% beta_m[c("L3")]
  lambda <- beta_m["A"] * data$A + meanL + beta_m["b0"]
  lambda <- pmax(lambda, 0)  # Ensure non-negative lambda
  
  # Compute the Poisson density for each M
  f_M <- dpois(data$M, lambda = lambda)
  f_M[f_M < 1e-5] <- 1e-5
  return(f_M)
}

########################################
# Define the E-step
# Impute_step() loaded from a shared simulation implementation below.


#(data_impu <- Impute_step(data = dataset, beta_m, S = 10, LOD = LOD_value))


mod_update_glm <- function(data) {
  # Check and assign w_norm if NULL
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  # Censoring Modeling (parametric model, does not include Y)
  #cen_mod <- glm(C ~ L1 + L2 + L3 + A + M + A*M, data = data,family = binomial, weights = w_norm)
  
  #cen_prob <- predict(cen_mod, type = "response")
  # Outcome Modeling
  out_mod <- glm(Y ~ L3 + A + M, data = data, family = binomial, weights = w_norm)
  out_probY1 <- predict(out_mod, type = "response")
  out_probY0 <- 1 - out_probY1
  
  # Create a new column for predicted probabilities based on observed Y values
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)
  
  # Return predictions and coefficients
  return(list(
    #cen_prob = cen_prob, 
    out_prob = out_prob,
    #cen_coef = cen_mod$coefficients,
    out_coef = out_mod$coefficients
  ))
}


# mod_pred are from mod_update_glm funciton, which inlucdes both out_prob and cen_prob
# Wet_step() loaded from a shared simulation implementation below.


########################################
# Define the M-step, here is for M density
M_stepf_M <- function(data, beta_m) {

  neg_log_likelihood <- function(par, data) {
    weights <- data$w_norm
    eta <- par["b0"] + par["A"] * data$A + par["L3"] * data$L3
    lambda <- pmax(exp(eta), 1e-10)  # Lower bound to avoid 0
    
    log_lik <- dpois(data$M, lambda = lambda, log = TRUE)
    
    # Replace -Inf with large negative value for stability
    log_lik[!is.finite(log_lik)] <- -1e3
    
    return(-sum(weights * log_lik))
  }
  
  # Ensure parameters are finite
  beta_m_optim <- beta_m[c("b0", "A", "L3")]
  if (any(!is.finite(beta_m_optim))) stop("Initial parameters contain NA/Inf/NaN")
  
  # Optimization
  optim_results <- optim(par = beta_m_optim, fn = neg_log_likelihood, data = data, method = "BFGS")
  
  return(beta_m_new = optim_results$par)
}


# Example usage:
# QM <- -3 + 1.5 * A + 1.75 * L1 + 1.5 * L2 - 0.25 * L3,       0.25
# beta_m_start <- c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.5)
#M_stepf_M(data = data_weted, beta_m = beta_m_start)

# EM Algorithm
EM_algorithm <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 50, tol = 1e-5, verbose = TRUE) {
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


# gcomp.nde() moved to simulations/R/shared/simulation_helpers.R


# Perform bootstrap
# gcom_boots_ci() moved to simulations/R/shared/simulation_helpers.R


#######
#m-out-of-n bootstrap confidence intervals
#######
# Set up the parallel plan (e.g., multisession for all systems or multicore for Unix-like systems)

# Function for m-out-of-n bootstrap confidence intervals
# gcom_boots_ci_m_out_n() moved to simulations/R/shared/simulation_helpers.R


# Function to compute effect estimates
# compute_eff() moved to simulations/R/shared/simulation_helpers.R


# Function to choose c_i using double bootstrap ideas (parallelized with early stopping)
choose_ci <- function(data, p_cen = 0.5, k = 5, gamma_seq = seq(0, 30, length.out = 50), B1 = 150, B2 = 100, nominal_coverage = 0.85,
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


# gcomp.nie() moved to simulations/R/shared/simulation_helpers.R


# Onestep_med_effs() loaded from a shared simulation implementation below.


sim_iteration <- function(i, n_obs, k = 5, type = c("Trad", "EM-MLE", "Oracle"), censor_rate, S, beta_m) {
  type <- match.arg(type)
  beta_mres <- NULL  # No fitted EM parameters for Trad/Oracle.
  output_dir <- here::here("simulations", "results", "raw", SCENARIO_ID)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, paste0("Gcom_W2", type, "_n_", n_obs, "rate", censor_rate, ".txt"))
  data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
  simdata <- data$study_dat
  LOD_value <- as.numeric(data$lod)
  
  if (type == "Trad") {
    data_aug <- NaiveImpute(simdata, LOD_value)
  } else if (type == "EM-MLE") {
    em_result <- EM_algorithm(data = simdata, beta_m_0 = beta_m, S = S, LOD = LOD_value)
    data_aug <- em_result$data_weted
    beta_mres <- em_result$beta_m  
  } else if (type == "Oracle") {
    data_aug <- data$study_dat_full
  }
  
  # G-computation
  nde_gcom_glm <- gcomp.nde(data_aug)
  # m-out-of-n #p_cen = 0.5, k = 5, 
  choose_ciNDE <- choose_ci(data = data_aug, p_cen = censor_rate, k = k, type = "NDE", measure = "RD")
  CI.nderes <- gcom_boots_ci_m_out_n(data_aug, m = choose_ciNDE$opt_m, n_bootstrap = 1000, type = "NDE", measure = "RD")
  CI.nde <- CI.nderes$ci
  var.nde <- CI.nderes$variance
  # regular
  CI.nderes.reg <- gcom_boots_ci(data_aug, n_bootstrap = 1000, type = "NDE", measure = "RD")
  CI.nde.reg <- CI.nderes.reg$ci
  var.nde.reg <- CI.nderes.reg$variance
  
  #NIE
  nie_gcom_glm <- gcomp.nie(data_aug)
  # m-out-of-n 
  choose_ciNIE <- choose_ci(data_aug, p_cen = censor_rate, k = k, type = "NIE", measure = "RD")
  CI.nieres <- gcom_boots_ci_m_out_n(data_aug, m = choose_ciNIE$opt_m, n_bootstrap = 1000, type = "NIE", measure = "RD")
  CI.nie <- CI.nieres$ci
  var.nie <- CI.nieres$variance
  # regular
  CI.nieres.reg <- gcom_boots_ci(data_aug, n_bootstrap = 1000, type = "NIE", measure = "RD")
  CI.nie.reg <- CI.nieres.reg$ci
  var.nie.reg <- CI.nieres.reg$variance
  
  
  # Save results
  results <- c(i,nde_gcom_glm$nde.rd, CI.nde, var.nde, CI.nde.reg, var.nde.reg,
               nie_gcom_glm$nie.rd, CI.nie, var.nie, CI.nie.reg, var.nie.reg, beta_mres)
  
  write(results, file = output_file, ncolumns = length(results), append = TRUE)
  return(NULL)  # Returning NULL keeps memory usage low
}
