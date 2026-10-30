# Simulation scenario definition; settings preserved from the supplied repository.
# Scenario-specific function bodies are retained; only byte-identical shared helpers were centralized.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
source(here::here("simulations", "R", "shared", "dgp.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "simulation_helpers.R"), local = TRUE)
SCENARIO_ID <- "study1_Hetero_CDE_HAL"
source(here::here("simulations", "R", "shared", "variants", "em_hal.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "impute_integer_indicator.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "onestep_mediation.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "weights_hal.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "select_gcomp_ci.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "mediator_lognormal_density.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "sample_truncated_lognormal.R"), local = TRUE)

library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(Rcpp)
library(medoutcon)
library(truncnorm)
library(xgboost)
library(zoo) 
# main dependencies 

###############################
## Density of M for samples below LOD
##
# sample_density_proposal() loaded from a shared simulation implementation below.

# Density Function to compute f(M)
# compute_f_M() loaded from a shared simulation implementation below.


compute_f_M_HAL <- function(data, K = 3, epsilon = 1e-8) {
  # Check for required variables
  required_vars <- c("M", "L1", "L2", "L3", "A")
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
  }
  
  # Ensure we have weights
  if (is.null(data$w_norm)) {
    data$w_norm <- 1
  }
  
  # Create fold assignments
  n <- nrow(data)
  fold_ids <- sample(rep(1:K, length.out = n))
  f_M_all <- numeric(n)
  
  for (k in 1:K) {
    # Split data
    train_data <- data[fold_ids != k, ]
    valid_data <- data[fold_ids == k, ]
    
    # Fit model on training data
    fit_k <- nadir::lnr_heteroskedastic_density(
      data = train_data,
      formula = M ~ L1 + L2 + L3 + A,
      mean_lnr = lnr_hal,
      var_lnr  = lnr_hal,
      mean_lnr_args = list(
        weights = train_data$w_norm,
        smoothness_orders = 1,
        max_degree        = 2,
        num_knots         = c(100, 75, 50),
        fit_control = list(
          cv_select       = TRUE,
          use_min         = TRUE,
          lambda.min.ratio= 1e-5,
          prediction_bounds = "default"
        )
      ),
      var_lnr_args = list(
        weights = train_data$w_norm,
        smoothness_orders = 1,
        max_degree        = 2,
        num_knots         = c(100, 75, 50),
        fit_control = list(
          cv_select       = TRUE,
          use_min         = TRUE,
          lambda.min.ratio= 1e-5,
          prediction_bounds = "default"
        )
      )
    )
    
    # Predict density on validation data
    f_M_valid <- fit_k(valid_data)
    f_M_valid <- pmax(f_M_valid, epsilon)  # stability
    
    # Store results
    f_M_all[fold_ids == k] <- f_M_valid
  }
  
  # Compute cross-fitted loss
  loss <- mean(-log(f_M_all))
  # loss <- sum(-log(f_M_all))
  
  
  return(list(
    f_M = f_M_all,
    loss = loss,
    fold_ids = fold_ids
  ))
}


########################################
# Define the E-step
# Impute_step() loaded from a shared simulation implementation below.


#############
# METHOD 2
# mod_update_hal() moved to simulations/R/shared/simulation_helpers.R


#system.time({
#  halmod_pred <- mod_update_hal(data = data_impu)
#})
# Access results
#halmod_pred$out_prob # Predicted probabilities for the outcome model
#halmod_pred$basis_list # Shared basis functions


# Wet_step() loaded from a shared simulation implementation below.

# (data_weted <- Wet_step(data_imputed = data_impu, mod_pred = halmod_pred, beta_m_0 = beta_m))


# EM Algorithm
# EM_algorithm() loaded from a shared simulation implementation below.


# Define initial guesses for beta_m parameters
beta_m_0 = c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.35)
# Run the EM algorithm
#dat_res <- EM_algorithm(data = dataset, beta_m_0 = beta_m_0)

#dat_aug <- dat_res$data_weted


# gcomp.nde() moved to simulations/R/shared/simulation_helpers.R


# Perform bootstrap
library(future.apply)
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
# choose_ci() loaded from a shared simulation implementation below.


# gcomp.nie() moved to simulations/R/shared/simulation_helpers.R


# Onestep_med_effs() loaded from a shared simulation implementation below.


sim_iteration <- function(i, n_obs, k = 5, type = c("Trad", "EM-MLE", "Oracle"), censor_rate, S, beta_m) {
  type <- match.arg(type)
  output_dir <- here::here("simulations", "results", "raw", SCENARIO_ID)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, paste0("GcomHAL_hetero4", type, "_n_", n_obs, "rate", censor_rate, ".txt"))
  data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
  simdata <- data$study_dat
  LOD_value <- as.numeric(data$lod)
  
  if (type == "Trad") {
    data_aug <- NaiveImpute(simdata, LOD_value)
  } else if (type == "EM-MLE") {
    em_result <- EM_algorithm(data = simdata,  beta_m_0 = beta_m, S=S, LOD = LOD_value)
    data_aug <- em_result$data_weted
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
  
  # One-step estimation
  #onestep.medeff <- Onestep_med_effs(dat_aug = data_aug)
  
  # Save results
  results <- c(i,nde_gcom_glm$nde.rd, CI.nde, var.nde, CI.nde.reg, var.nde.reg,
               nie_gcom_glm$nie.rd, CI.nie, var.nie, CI.nie.reg, var.nie.reg)
  
  write(results, file = output_file, ncolumns = length(results), append = TRUE)
  return(NULL)  # Returning NULL keeps memory usage low
}
