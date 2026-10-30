# Simulation scenario definition; settings preserved from the supplied repository.
# Scenario-specific function bodies are retained; only byte-identical shared helpers were centralized.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
source(here::here("simulations", "R", "shared", "dgp.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "simulation_helpers.R"), local = TRUE)
SCENARIO_ID <- "study1_LLoQ_Gcom"
source(here::here("simulations", "R", "shared", "variants", "em_parametric.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "impute_numeric_indicator.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "onestep_mediation.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "weights_parametric.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "update_outcome_glm.R"), local = TRUE)

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
sample_density_proposal <- function(n, data, LOD, beta_m) {
  # Linear predictor for L variables
  meanL <- as.matrix(data[c("L1", "L2", "L3")]) %*% beta_m[c("L1", "L2", "L3")]
  # Mean calculation for the normal distribution
  means <- beta_m["A"] * data$A + meanL + beta_m["b0"]
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
# Density Function to compute f(M)
compute_f_M <- function(data, beta_m) {
  # Calculate linear predictors for the mean of the normal distribution
  meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
  # Total predicted mean
  means <- beta_m["A"] * data$A + meanL + beta_m["b0"]
  
  # Compute the log-normal density for each M
  # Convert the mean to a log scale for log-normal parameters
  f_M <- dlnorm(data$M, meanlog = means, sdlog = beta_m["sd"])
  return(f_M)
}
########################################
# Define the E-step
# Impute_step() loaded from a shared simulation implementation below.


#(data_impu <- Impute_step(data = dataset, beta_m, S = 10, LOD = LOD_value))


# mod_update_glm() loaded from a shared simulation implementation below.


# mod_pred are from mod_update_glm funciton, which inlucdes both out_prob and cen_prob
# Wet_step() loaded from a shared simulation implementation below.


########################################
# Define the M-step, here is for M density
M_stepf_M <- function(data, beta_m) {
  # Define an internal function for optimization
  neg_log_likelihood <- function(beta_m, data) {
    # Extract weights
    weights <- data$w_norm
    # Calculate the mean for each sample using provided beta_m
    meanL <- as.matrix(data[, .(L1, L2, L3)]) %*% beta_m[c("L1", "L2", "L3")]
    means <- beta_m["A"] * data$A + meanL + beta_m["b0"]
    
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

# Example usage:
# QM <- -3 + 1.5 * A + 1.75 * L1 + 1.5 * L2 - 0.25 * L3,       0.25
# beta_m_start <- c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.5)
#M_stepf_M(data = data_weted, beta_m = beta_m_start)

# EM Algorithm
# EM_algorithm() loaded from a shared simulation implementation below.


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
choose_ci <- function(data, p_cen = 0.5, k = 5, gamma_seq = seq(0, 30, length.out = 100), B1 = 100, B2 = 50, nominal_coverage = 0.95,
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
  output_dir <- here::here("simulations", "results", "raw", SCENARIO_ID)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, paste0("Gcom_Trad", type, "_n_", n_obs, "rate", censor_rate, ".txt"))
  data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
  simdata <- data$study_dat
  LOD_value <- as.numeric(data$lod)
  
  if (type == "Trad") {
    data_aug <- NaiveImpute(simdata, LOD_value)
  } else if (type == "EM-MLE") {
    em_result <- EM_algorithm(simdata, beta_m, S = S, LOD = LOD_value)
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
               nie_gcom_glm$nie.rd, CI.nie, var.nie, CI.nie.reg, var.nie.reg)
  
  write(results, file = output_file, ncolumns = length(results), append = TRUE)
  return(NULL)  # Returning NULL keeps memory usage low
}
