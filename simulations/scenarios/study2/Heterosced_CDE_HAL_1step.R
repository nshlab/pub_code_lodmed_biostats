# Simulation scenario definition; settings preserved from the supplied repository.
# Scenario-specific function bodies are retained; only byte-identical shared helpers were centralized.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
source(here::here("simulations", "R", "shared", "dgp.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "simulation_helpers.R"), local = TRUE)
SCENARIO_ID <- "study2_Heterosced_CDE_HAL_1step"
source(here::here("simulations", "R", "shared", "variants", "em_hal.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "impute_integer_indicator.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "weights_hal.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "mediator_lognormal_density.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "sample_truncated_lognormal.R"), local = TRUE)

library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(xgboost)
library(zoo) 
# main dependencies 
library(nadir)


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
  if (!"w_norm" %in% names(data) || is.null(data$w_norm)) {
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


#######
#m-out-of-n bootstrap confidence intervals
#######
library(future.apply)
# Set up the parallel plan (e.g., multisession for all systems or multicore for Unix-like systems)

moutn_multiplier_bootstrap <- function(eif, theta, B = 200, m = 100,  multipliers = "rademacher", ci_level = 0.95) {
  n <- length(eif)
  samples <- sample(1:n, size = m, replace = F)
  eif <- eif[samples]
  
  m_process <- numeric(B)
  se_eif <- sqrt(stats::var(eif) / n)
  
  for (b in 1:B) {
    # Step 1: Draw multipliers
    if (multipliers == "rademacher") {
      xi <- (2 * stats::rbinom(length(samples), 1, 0.5)) - 1
    } else if (multipliers == "gaussian") {
      xi <- rnorm(length(samples))
    }
    
    # Step 2: compute process M using multipliers
    m_process[b] <- mean((xi * eif) / se_eif)
  }
  
  c_alpha <- unname(stats::quantile(m_process, ci_level, na.rm = TRUE))
  
  
  lwr_ci <- theta - c_alpha * se_eif
  upr_ci <- theta + c_alpha * se_eif
  return(list(m_process = m_process,
              lwr_ci = lwr_ci,
              upr_ci = upr_ci,
              ci = c(lwr_ci, upr_ci),
              c_alpha = c_alpha))
}


# os_compute_eff() moved to simulations/R/shared/simulation_helpers.R


# Function to choose c_i using double bootstrap ideas (parallelized with early stopping)
os_choose_ci <- function(data, original_estimate, p_cen = 0.5, k = 5, gamma_seq = seq(0.001, 30, length.out = 100), B1 = 100, B2 = 50000, 
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
    
    # Parallelize the first-stage bootstrap (B1 iterations)
    coverage_rate <- future_sapply(1:B1, function(b1) {
      first_stage_indices <- sample(1:n, n, replace = TRUE)  # Resample data
      #first_stage_estimate <- compute_eff(data, first_stage_indices, type = type, measure = measure)
      
      # Second-stage (nested) bootstrap
      os_eff <- os_compute_eff(data = data, indices = first_stage_indices, type = type, fast = T)
      nested_res <- moutn_multiplier_bootstrap(eif = os_eff$eif, theta = os_eff$est, B = B2, m = m_i)
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


Onestep_med_effs <- function(dat_aug, effect = c("both", "NDE", "NIE")) {
  effect <- match.arg(effect)
  
  # Learners
  glm_lrnr <- Lrnr_glm_fast$new()
  lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
  lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
  lrn_earth <- Lrnr_earth$new()
  lrn_ranger <- Lrnr_ranger$new(
    probability = TRUE, num.trees = 200,
    mtry = 2, min.node.size = 5
  )
  lrn_xgb <- Lrnr_xgboost$new(
    objective = "binary:logistic",
    max_depth = 3, eta = 0.15, nrounds = 200,
    subsample = 0.9, colsample_bytree = 0.9,
    min_child_weight = 5, gamma = 0.1, verbosity = 0
  )
  lrn_gam <- Lrnr_gam$new()
  lrn_hal <- Lrnr_hal9001$new(
    family = "binomial", max_degree = 2, smoothness_orders = 1,
    num_knots = c(100, 75, 50),
    fit_control = list(
      cv_select = TRUE, use_min = TRUE,
      lambda.min.ratio = 1e-5, prediction_bounds = "default"
    )
  )
  
  lrnr_lib <- Stack$new(lrn_hal, lrn_earth, lrn_xgb)
  sl_lrnr <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_solnp$new())
  glm_lrnr <- Lrnr_glm_fast$new(family = binomial())
  
  # Define interactions
  lrnr_interact1 <- Lrnr_define_interactions$new(list(c("W_1", "W_3"), c("W_1", "W_2")))
  lrnr_interact2 <- Lrnr_define_interactions$new(list(c("A", "M_1"), c("A", "W_1")))
  lrnr_interact3 <- Lrnr_define_interactions$new(list(c("W_1", "W_3"), c("W_1", "W_2"), c("W_2", "W_3")))
  lrnr_interact4 <- Lrnr_define_interactions$new(list(c("A", "M_1"), c("A", "W_1"), c("M_1", "W_3")))
  
  interaction_pipeline_glmPS <- sl3::make_learner(sl3::Pipeline, lrnr_interact1, glm_lrnr)
  interaction_pipeline_glmPS2 <- sl3::make_learner(sl3::Pipeline, lrnr_interact3, glm_lrnr)
  interaction_pipeline_glmOUT <- sl3::make_learner(sl3::Pipeline, lrnr_interact2, glm_lrnr)
  interaction_pipeline_glmOUT2 <- sl3::make_learner(sl3::Pipeline, lrnr_interact4, glm_lrnr)
  
  result_list <- list()
  
  # ---- NDE ----
  if (effect %in% c("both", "NDE")) {
    nde_onestep <- medoutcon(
      W = dat_aug[, c("L1", "L2", "L3")],
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = dat_aug$M,
      Y = dat_aug$Y,
      obs_weights = dat_aug$w_norm,
      g_learners = interaction_pipeline_glmPS2,
      h_learners = lrn_xgb,
      b_learners = interaction_pipeline_glmOUT2,
      v_learners = lrn_xgb,
      effect = "direct",
      estimator = "onestep",
      estimator_args = list(cv_folds = 15L, max_iter = 5L, tiltmod_tol = 5)
    )
    
    nde_os <- summary(nde_onestep)
    result_list$eif <- nde_onestep$eif
    result_list$nde_onestep_est <- nde_os$param_est
    result_list$nde_onestep_var <- nde_os$var_est
    result_list$nde_onestep_ci_sd <- c(nde_os$lwr_ci, nde_os$upr_ci)
    result_list$nde_onestep_ci <- c(
      nde_os$param_est - 2.12 * sqrt(nde_os$var_est),
      nde_os$param_est + 2.12 * sqrt(nde_os$var_est)
    )
  }
  
  # ---- NIE ----
  if (effect %in% c("both", "NIE")) {
    nie_onestep <- medoutcon(
      W = dat_aug[, c("L1", "L2", "L3")],
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = dat_aug$M,
      Y = dat_aug$Y,
      obs_weights = dat_aug$w_norm,
      g_learners = interaction_pipeline_glmPS2,
      h_learners = lrn_xgb,
      b_learners = interaction_pipeline_glmOUT2,
      v_learners = lrn_xgb,
      effect = "indirect",
      estimator = "onestep",
      estimator_args = list(cv_folds = 15L, max_iter = 5L, tiltmod_tol = 5)
    )
    
    nie_os <- summary(nie_onestep)
    result_list$eif <- nie_onestep$eif
    result_list$nie_onestep_est <- nie_os$param_est
    result_list$nie_onestep_var <- nie_os$var_est
    result_list$nie_onestep_ci_sd <- c(nie_os$lwr_ci, nie_os$upr_ci)
    result_list$nie_onestep_ci <- c(
      nie_os$param_est - 2.12 * sqrt(nie_os$var_est),
      nie_os$param_est + 2.12 * sqrt(nie_os$var_est)
    )
  }
  
  return(result_list)
}


Onestep_med_effs_fast <- function(dat_aug, effect = c("both", "NDE", "NIE")) {
  effect <- match.arg(effect)
  
  # Learners
  glm_lrnr <- Lrnr_glm_fast$new()
  lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
  lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
  lrn_earth <- Lrnr_earth$new()
  lrn_ranger <- Lrnr_ranger$new(
    probability = TRUE, num.trees = 200,
    mtry = 2, min.node.size = 5
  )
  lrn_xgb <- Lrnr_xgboost$new(
    objective = "binary:logistic",
    max_depth = 3, eta = 0.15, nrounds = 200,
    subsample = 0.9, colsample_bytree = 0.9,
    min_child_weight = 5, gamma = 0.1, verbosity = 0
  )
  lrn_gam <- Lrnr_gam$new()
  lrn_hal <- Lrnr_hal9001$new(
    family = "binomial", max_degree = 2, smoothness_orders = 1,
    num_knots = c(100, 75, 50),
    fit_control = list(
      cv_select = TRUE, use_min = TRUE,
      lambda.min.ratio = 1e-5, prediction_bounds = "default"
    )
  )
  
  lrnr_lib <- Stack$new(lrn_hal, lrn_earth, lrn_xgb)
  sl_lrnr <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_solnp$new())
  glm_lrnr <- Lrnr_glm_fast$new(family = binomial())
  
  # Define interactions
  lrnr_interact1 <- Lrnr_define_interactions$new(list(c("W_1", "W_3"), c("W_1", "W_2")))
  lrnr_interact2 <- Lrnr_define_interactions$new(list(c("A", "M_1"), c("A", "W_1")))
  lrnr_interact3 <- Lrnr_define_interactions$new(list(c("W_1", "W_3"), c("W_1", "W_2"), c("W_2", "W_3")))
  lrnr_interact4 <- Lrnr_define_interactions$new(list(c("A", "M_1"), c("A", "W_1"), c("M_1", "W_3")))
  
  interaction_pipeline_glmPS <- sl3::make_learner(sl3::Pipeline, lrnr_interact1, glm_lrnr)
  interaction_pipeline_glmPS2 <- sl3::make_learner(sl3::Pipeline, lrnr_interact3, glm_lrnr)
  interaction_pipeline_glmOUT <- sl3::make_learner(sl3::Pipeline, lrnr_interact2, glm_lrnr)
  interaction_pipeline_glmOUT2 <- sl3::make_learner(sl3::Pipeline, lrnr_interact4, glm_lrnr)
  
  result_list <- list()
  
  # ---- NDE ----
  if (effect %in% c("both", "NDE")) {
    nde_onestep <- medoutcon(
      W = dat_aug[, c("L1", "L2", "L3")],
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = dat_aug$M,
      Y = dat_aug$Y,
      obs_weights = dat_aug$w_norm,
      g_learners = interaction_pipeline_glmPS2,
      h_learners = lrn_xgb,
      b_learners = interaction_pipeline_glmOUT2,
      effect = "direct",
      estimator = "onestep",
      estimator_args = list(cv_folds = 2L, max_iter = 5L, tiltmod_tol = 5)
    )
    
    nde_os <- summary(nde_onestep)
    result_list$eif <- nde_onestep$eif
    result_list$nde_onestep_est <- nde_os$param_est
    result_list$nde_onestep_var <- nde_os$var_est
    result_list$nde_onestep_ci_sd <- c(nde_os$lwr_ci, nde_os$upr_ci)
    result_list$nde_onestep_ci <- c(
      nde_os$param_est - 2.12 * sqrt(nde_os$var_est),
      nde_os$param_est + 2.12 * sqrt(nde_os$var_est)
    )
  }
  
  # ---- NIE ----
  if (effect %in% c("both", "NIE")) {
    nie_onestep <- medoutcon(
      W = dat_aug[, c("L1", "L2", "L3")],
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = dat_aug$M,
      Y = dat_aug$Y,
      obs_weights = dat_aug$w_norm,
      g_learners = interaction_pipeline_glmPS2,
      h_learners = lrn_xgb,
      b_learners = interaction_pipeline_glmOUT2,
      v_learners = lrn_xgb,
      effect = "indirect",
      estimator = "onestep",
      estimator_args = list(cv_folds = 2L, max_iter = 5L, tiltmod_tol = 5)
    )
    
    nie_os <- summary(nie_onestep)
    result_list$eif <- nie_onestep$eif
    result_list$nie_onestep_est <- nie_os$param_est
    result_list$nie_onestep_var <- nie_os$var_est
    result_list$nie_onestep_ci_sd <- c(nie_os$lwr_ci, nie_os$upr_ci)
    result_list$nie_onestep_ci <- c(
      nie_os$param_est - 2.12 * sqrt(nie_os$var_est),
      nie_os$param_est + 2.12 * sqrt(nie_os$var_est)
    )
  }
  
  return(result_list)
}


sim_iteration <- function(i, n_obs, k = 5, type = c("Trad", "EM-MLE", "Oracle"), censor_rate, S, beta_m) {
  type <- match.arg(type)
  output_dir <- here::here("simulations", "results", "raw", SCENARIO_ID)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, paste0("hetero1stepJul23_v9", type, "_n_", n_obs, "rate", censor_rate, ".txt"))
  data <- Study_dgp_Mcensor(n_obs = n_obs, censor_rate = censor_rate)
  simdata <- data$study_dat
  LOD_value <- as.numeric(data$lod)
  
  if (type == "Trad") {
    data_aug <- NaiveImpute(simdata, LOD_value)
  } else if (type == "EM-MLE") {
    em_result <- EM_algorithm(simdata, beta_m, S = S, LOD = LOD_value)
    data_aug <- em_result$data_weted
    beta_m <- em_result$beta_m  
  } else if (type == "Oracle") {
    data_aug <- data$study_dat_full
  }
  
  # One-step estimation
  onestep.nde <- Onestep_med_effs(dat_aug = data_aug, effect = "NDE")
  # m-out-of-n #p_cen = 0.5, k = 5, 
  choose_ciNDE <- os_choose_ci(data = data_aug, original_estimate = onestep.nde$nde_onestep_est, p_cen = censor_rate, k = k, type = "NDE")
  CI.nderes <- moutn_multiplier_bootstrap(eif = onestep.nde$eif, theta = onestep.nde$nde_onestep_est, 
                                          B = 2000, m = choose_ciNDE$opt_m)
  CI.nde <- CI.nderes$ci
  
  onestep.nie <- Onestep_med_effs(dat_aug = data_aug, effect = "NIE")
  #choose_ciNIE <- os_choose_ci(data = data_aug, original_estimate = onestep.nie$nie_onestep_est, p_cen = censor_rate, k = k, type = "NIE")
  CI.nieres <- moutn_multiplier_bootstrap(eif = onestep.nie$eif, theta = onestep.nie$nie_onestep_est,
                                          B = 2000, m = choose_ciNDE$opt_m)
  CI.nie <- CI.nieres$ci
  # Save results
  results <- c(i,
               onestep.nde$nde_onestep_est, CI.nde, onestep.nde$nde_onestep_ci_sd, onestep.nde$nde_onestep_var, 
               onestep.nie$nie_onestep_est, CI.nie, onestep.nie$nie_onestep_ci_sd, onestep.nie$nie_onestep_var, choose_ciNDE$optimal_c)
  
  
  write(results, file = output_file, ncolumns = length(results), append = TRUE)
  return(NULL)  # Returning NULL keeps memory usage low
}
