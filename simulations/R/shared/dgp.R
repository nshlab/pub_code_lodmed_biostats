# Clear the entire environment
# Load necessary package
library(data.table)
#' Simulate Study Data with Censoring Due to Limit of Detection (LOD)
#'
#' This function generates synthetic data for a study involving:
#' - Baseline covariates (L1–L3)
#' - Treatment assignment (A), either randomized or confounded
#' - A continuous post-treatment biomarker (M), subject to left censoring at a limit of detection (LOD)
#' - A binary outcome (Y) depending on exposure and biomarker
#'
#' Censoring is induced either by specifying a fixed LOD value or by targeting a desired censoring rate.
#'
#' @param n_obs Number of observations (default = 5000)
#' @param random Logical; whether treatment is randomized (default = FALSE)
#' @param lod Optional numeric; fixed limit of detection. If NULL, a quantile-based LOD is derived.
#' @param censor_rate If `lod` is NULL, sets LOD as the quantile corresponding to this proportion (default = 0.3)
#' @param return_dgp Logical; whether to return the underlying data-generating functions (default = TRUE)
#'
#' @return A list containing:
#' - `study_dat`: Observed data with censored M values ("LOD" where M < LOD)
#' - `study_dat_full`: Complete data without censoring
#' - `dgp_funs`: List of functions used to generate A, M, and Y (only if `return_dgp = TRUE`)
#' - `lod`: The value of the limit of detection used
#'
#' @examples
#' # Simulate randomized trial data with 1000 participants
#' sim_data <- Study_dgp_Mcensor(n_obs = 1000, random = TRUE)
#' head(sim_data$study_dat)
#'
#' # Simulate observational data with a fixed LOD of 5
#' sim_data2 <- Study_dgp_Mcensor(n_obs = 500, random = FALSE, lod = 5)
#'
#' # Simulate data with a target censoring rate of 40%
#' sim_data3 <- Study_dgp_Mcensor(n_obs = 2000, censor_rate = 0.4)
#'
#' # Access the full (uncensored) biomarker values
#' head(sim_data$study_dat_full$M)
#'
#' # Plot observed vs full M values to assess censoring
#' with(sim_data, plot(study_dat_full$M, col = ifelse(study_dat$C == 1, "blue", "red"),
#'      main = "Biomarker M: Observed (blue) vs Censored (red)"))
#'
#' @export
Study_dgp_Mcensor <- function(n_obs = 5000, random = FALSE, lod = NULL, censor_rate = 0.3, return_dgp = TRUE) {
  # 1. Generate baseline covariates
  L1 <- rbinom(n_obs, 1, 0.6)
  L2 <- rbinom(n_obs, 1, 0.5)
  L3 <- rbinom(n_obs, 1, 0.25)
  
  # 2. Treatment assignment: random or confounded
  gA <- function(L1, L2, L3, random) {
    if (random) {
      return(rep(0.5, length(L1)))
    } else {
      return(plogis(-1 + 0.5 * L1 + 1.25 * L2 + 0.75 * L3 - 1.25 * L1 * L3))
    }
  }
  gA_obs <- gA(L1, L2, L3, random)
  A <- rbinom(n_obs, 1, gA_obs)
  
  # 3. Generate biomarker values (M_full)
  QM <- function(A, L1, L2, L3) {
    linear_predictor <- -3 + 1.5 * A + 1.75 * L1 + 0.25 * A * L1 + 1.5 * L2 - 0.25 * L3
    return(list(QM = linear_predictor))
  }
  QM_obs <- QM(A, L1, L2, L3)
  M_full <- rlnorm(n_obs, QM_obs$QM, 0.25)
  
  # 4. Apply LOD censoring
  if (is.null(lod)) {
    lod <- quantile(M_full, probs = censor_rate)
  }
  C <- M_full >= lod   # Censoring indicator: 1 = observed, 0 = censored
  M <- M_full
  M[!C] <- "LOD"       # Replace censored values with "LOD" string
  
  # 5. Outcome generation based on full M values
  QY <- function(A, M, L1, L2, L3) {
    return(-1 + 2.5 * A + 1.75 * M + 0.5 * A * M - 2.25 * L1 - 1.75 * L2 - 1.5 * L3)
  }
  QY_obs <- QY(A, M_full, L1, L2, L3)
  Yprob <- plogis(QY_obs)
  Y <- rbinom(n_obs, 1, Yprob)
  
  # 6. Assemble datasets
  study_dat <- data.table::setDT(
    data.frame(id = seq_len(n_obs), L1, L2, L3, A, M, Y, C = as.numeric(C))
  )
  study_dat_full <- data.table::setDT(
    data.frame(id = seq_len(n_obs), L1, L2, L3, A, M = M_full, Y)
  )
  
  # 7. Return list with optional DGP functions
  dgp_funs <- list(gA = gA, QM = QM, QY = QY)
  if (return_dgp) {
    return(list(
      study_dat = study_dat,
      study_dat_full = study_dat_full,
      dgp_funs = dgp_funs,
      lod = lod
    ))
  } else {
    return(list(
      study_dat = study_dat,
      study_dat_full = study_dat_full
    ))
  }
}


# Executable DGP demo removed from the active shared helper.
# Exploratory demonstrations are excluded from this module.
#' Approximate Average Treatment Effect (ATE), Natural Direct and Indirect Effects (NDE/NIE)
#'
#' @param gen_data A function that generates data from a data-generating process.
#' @param n_truth Number of observations to generate for an "asymptotic" sample
#'   in which the plug-in estimator of the ATE will be used.
#' @param random Logical; indicates whether the study is a randomized trial.
#'
#' @return A list containing:
#' \describe{
#'   \item{psi_rd}{Risk difference estimate of the ATE.}
#'   \item{psi_rr}{Risk ratio estimate of the ATE.}
#'   \item{psi_or}{Odds ratio estimate of the ATE.}
#'   \item{var_eif}{Variance bound based on the efficient influence function.}
#'   \item{nde_rd}{Natural Direct Effect (risk difference).}
#'   \item{nde_rr}{Natural Direct Effect (risk ratio).}
#'   \item{nde_or}{Natural Direct Effect (odds ratio).}
#'   \item{nie_rd}{Natural Indirect Effect (risk difference).}
#'   \item{nie_rr}{Natural Indirect Effect (risk ratio).}
#'   \item{nie_or}{Natural Indirect Effect (odds ratio).}
#' }
#'
#' @examples
#'
#' # Run get_truth with Study_dgp_Mcensor DGP (using a smaller n for illustration)
#' truth <- get_truth(gen_data = Study_dgp_Mcensor,  n_truth = 1e4)
#' print(truth$psi_rd)
#' print(truth$nde_rd)
#'
#' @export

get_truth <- function(gen_data, random = FALSE, n_truth = 1e8) {
  # Compute large data set from data-generating mechanism
  dgp <- gen_data(n_obs = n_truth, random =random, return_dgp = TRUE)
  wow_so_much_data <- dgp$study_dat_full
  
  # Extract DGP helper functions
  gA <- dgp$dgp_funs$gA
  QM <- dgp$dgp_funs$QM
  QY <- dgp$dgp_funs$QY
  
  # Compute exposure (infection) mechanism in "asymptotic" sample
  gA_mech <- with(wow_so_much_data, gA(L1, L2, L3, random))
  
  # compute post-infection labs and biomarkers in "asymptotic" sample
  QM_Anat_mech <- with(wow_so_much_data, QM(A, L1, L2, L3))
  
  # compute post-infection labs and biomarkers in "asymptotic" sample for
  # counterfactual data where all units are exposed (A = 1)
  QM_Ais1_mech <- with(wow_so_much_data, QM(1, L1, L2, L3))
  M1_Ais1 <- rlnorm(n_truth, QM_Ais1_mech$QM, 0.25)
  
  # compute post-infection labs and biomarkers in "asymptotic" sample for
  # counterfactual data where all units are unexposed (A = 0)
  QM_Ais0_mech <- with(wow_so_much_data, QM(0, L1, L2, L3))
  M1_Ais0 <- rlnorm(n_truth, QM_Ais0_mech$QM, 0.25)
  
  # Compute outcome of interest in "asymptotic" sample
  Y_mech <- with(wow_so_much_data, QY(A, M, L1, L2, L3))
  QY_mech <- rbinom(n_truth, 1, plogis(Y_mech) )
  
  Y_Ais1_mech <- with(wow_so_much_data, QY(1, M1_Ais1, L1, L2, L3))
  Y_Ais0_mech <- with(wow_so_much_data, QY(0, M1_Ais0, L1, L2, L3))
  
  QY_Ais1_mech <- rbinom(n_truth, 1, plogis(Y_Ais1_mech) )
  QY_Ais0_mech <- rbinom(n_truth, 1, plogis(Y_Ais0_mech) )
  
  
  # Compute average treatment effect (ATE) via plug-in in "asymptotic" sample
  ate_rd <- mean(QY_Ais1_mech) - mean(QY_Ais0_mech)
  ate_rr <- mean(QY_Ais1_mech)/mean(QY_Ais0_mech)
  ate_or <- (mean(QY_Ais1_mech)/(1 - mean(QY_Ais1_mech)))/(mean(QY_Ais0_mech)/(1 - mean(QY_Ais0_mech)))
  
  # Compute direct and indirect effects
  Y_corss_mech <- with(wow_so_much_data, QY(1, M1_Ais0, L1, L2, L3))
  QY_corss_mech <- rbinom(n_truth, 1, plogis(Y_corss_mech) )
  
  nie_rd <- mean(QY_Ais1_mech) - mean(QY_corss_mech)  # Natural Indirect Effect RD
  nie_rr <- mean(QY_Ais1_mech)/mean(QY_corss_mech)  # Natural Indirect Effect RR
  nie_or <- (mean(QY_Ais1_mech)/(1 - mean(QY_Ais1_mech)))/(mean(QY_corss_mech)/(1-mean(QY_corss_mech)) ) # Natural Indirect Effect OR
  
  nde_rd <- mean(QY_corss_mech) - mean(QY_Ais0_mech)  # Natural Direct Effect RD
  nde_rr <- mean(QY_corss_mech)/mean(QY_Ais0_mech)  # Natural Direct Effect RD
  nde_or <- (mean(QY_corss_mech)/(1-mean(QY_corss_mech)) )/(mean(QY_Ais0_mech)/(1 - mean(QY_Ais0_mech)))  # Natural Direct Effect RD
  
  
  # Compute the variance bound via the EIF
  # EIF for ATE
  eif_Ais1 <- (wow_so_much_data$A == 1) / gA_mech * (wow_so_much_data$Y - QY_mech) + QY_Ais1_mech
  eif_Ais0 <- (wow_so_much_data$A == 0) / (1 - gA_mech) * (wow_so_much_data$Y - QY_mech) + QY_Ais0_mech
  eif_ate <- (eif_Ais1 - eif_Ais0) - ate_rd
  
  # EIF for NIE and NDE
  return(list(
    psi_rd = ate_rd,
    psi_rr = ate_rr,
    psi_or = ate_rr,
    var_eif = var(eif_ate),
    nde_rd = nde_rd,
    nde_rr = nde_rr,
    nde_or = nde_or,
    nie_rd = nie_rd,
    nie_rr = nie_rr,
    nie_or = nie_or
  ))
}

#(truth <- get_truth(gen_data = Study_dgp_Mcensor))



