simulation_scenarios <- list(
  study1_EM1_ture_nuis_Gcom = list(
    file = "simulations/scenarios/study1/EM1_ture_nuis_Gcom.R",
    replicates = 500L, n_obs = 2000L, k = 10L,
    type = "EM-MLE", censor_rate = 0.75, S = 5L,
    beta_m = c(b0 = -2.75, A = 1.75, L1 = 2, AL1 = 0.35, L2 = 1.25, L3 = -0.35, sd = 0.35),
    packages = c("data.table", "future.apply", "dplyr", "sl3", "hal9001", "medoutcon", "truncnorm", "xgboost")
  ),
  study1_Hetero_CDE_HAL = list(
    file = "simulations/scenarios/study1/Hetero_CDE_HAL.R",
    replicates = 1000L, n_obs = 2000L, k = 10L,
    type = "EM-MLE", censor_rate = 0.75, S = 5L,
    beta_m = c(b0 = -2.85, A = 1.65, L1 = 1.65, AL1 = 0.3, L2 = 1.35, L3 = -0.35, sd = 0.25),
    packages = c("data.table", "future.apply", "dplyr", "sl3", "hal9001", "Rcpp", "medoutcon", "truncnorm", "xgboost", "zoo", "nadir")
  ),
  study1_LLoQ_Gcom = list(
    file = "simulations/scenarios/study1/LLoQ_Gcom.R",
    replicates = 300L, n_obs = 1200L, k = 10L,
    type = "Trad", censor_rate = 0.75, S = 5L,
    beta_m = c(b0 = -2.75, A = 1.75, L1 = 2, L2 = 1.25, L3 = -0.35, sd = 0.35),
    packages = c("data.table", "future.apply", "dplyr", "sl3", "hal9001", "medoutcon", "truncnorm", "xgboost")
  ),
  study1_Wrong_GLM_Gcom = list(
    file = "simulations/scenarios/study1/Wrong_GLM_Gcom.R",
    replicates = 500L, n_obs = 120L, k = 2L,
    type = "EM-MLE", censor_rate = 0.25, S = 3L,
    beta_m = c(b0 = log(5), A = 1, L3 = 0.2),
    packages = c("data.table", "future.apply", "dplyr", "sl3", "hal9001", "medoutcon", "truncnorm", "xgboost")
  ),
  study2_Heterosced_CDE_HAL_1step = list(
    file = "simulations/scenarios/study2/Heterosced_CDE_HAL_1step.R",
    replicates = 1000L, n_obs = 300L, k = 8L,
    type = "EM-MLE", censor_rate = 0.5, S = 4L,
    beta_m = c(b0 = -2.95, A = 1.55, L1 = 1.85, AL1 = 0.25, L2 = 1.45, L3 = -0.3, sd = 0.24),
    packages = c("data.table", "future.apply", "dplyr", "sl3", "hal9001", "medoutcon", "truncnorm", "xgboost", "zoo", "nadir")
  ),
  study2_LLoQ_1step = list(
    file = "simulations/scenarios/study2/LLoQ_1step.R",
    replicates = 300L, n_obs = 300L, k = 6L,
    type = "Trad", censor_rate = 0.25, S = 3L,
    beta_m = c(b0 = -2.95, A = 1.55, L1 = 1.85, AL1 = 0.25, L2 = 1.45, L3 = -0.3, sd = 0.24),
    packages = c("data.table", "future.apply", "dplyr", "sl3", "hal9001", "medoutcon", "truncnorm", "xgboost", "ranger", "Rsolnp", "nadir")
  )
)
