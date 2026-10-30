# Simulation scenario definition; settings preserved from the supplied repository.
# Scenario-specific function bodies are retained; only byte-identical shared helpers were centralized.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
source(here::here("simulations", "R", "shared", "dgp.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "simulation_helpers.R"), local = TRUE)
SCENARIO_ID <- "study1_EM1_ture_nuis_Gcom"
source(here::here("simulations", "R", "shared", "variants", "em_parametric.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "impute_numeric_indicator.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "mediator_m_step.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "onestep_mediation.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "weights_parametric.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "select_gcomp_ci.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "mediator_lognormal_density.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "update_outcome_glm.R"), local = TRUE)
source(here::here("simulations", "R", "shared", "variants", "sample_truncated_lognormal.R"), local = TRUE)

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
# sample_density_proposal() loaded from a shared simulation implementation below.

# Density Function to compute f(M)
# compute_f_M() loaded from a shared simulation implementation below.

########################################
# Define the E-step
# Impute_step() loaded from a shared simulation implementation below.


#(data_impu <- Impute_step(data = dataset, beta_m, S = 10, LOD = LOD_value))


# mod_update_glm() loaded from a shared simulation implementation below.


# mod_pred are from mod_update_glm funciton, which inlucdes both out_prob and cen_prob
# Wet_step() loaded from a shared simulation implementation below.


########################################
# Define the M-step, here is for M density
# M_stepf_M() loaded from a shared simulation implementation below.


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
# choose_ci() loaded from a shared simulation implementation below.


# gcomp.nie() moved to simulations/R/shared/simulation_helpers.R


# Onestep_med_effs() loaded from a shared simulation implementation below.


sim_iteration <- function(i, n_obs, k = 5, type = c("Trad", "EM-MLE", "Oracle"), censor_rate, S, beta_m) {
  type <- match.arg(type)
  beta_mres <- NULL  # No fitted EM parameters for Trad/Oracle.
  output_dir <- here::here("simulations", "results", "raw", SCENARIO_ID)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- file.path(output_dir, paste0("Gcom", type, "_n_", n_obs, "rate", censor_rate, ".txt"))
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
               nie_gcom_glm$nie.rd, CI.nie, var.nie, CI.nie.reg, var.nie.reg, beta_mres)
  
  write(results, file = output_file, ncolumns = length(results), append = TRUE)
  return(NULL)  # Returning NULL keeps memory usage low
}
