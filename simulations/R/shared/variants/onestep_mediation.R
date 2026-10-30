# Shared verbatim implementation of Onestep_med_effs() used by 4 simulation scenarios.
# Extracted without rewriting the function body.

Onestep_med_effs <- function(dat_aug) {
  # Load necessary libraries
  # library(sl3)
  # library(tlverse)
  
  # Instantiate learners
  fglm_lrnr <- Lrnr_glm_fast$new()
  # penalized regressions:
  lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
  lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
  # spline regressions:
  #lrn_polspline <- Lrnr_polspline$new()
  lrn_earth <- Lrnr_earth$new()
  # tree-based methods
  lrn_ranger <- Lrnr_ranger$new()
  lrn_xgb <- Lrnr_xgboost$new()
  # generalized additive model (GAM) and Bayesian GLM
  lrn_gam <- Lrnr_gam$new()
  #lrn_bayesglm <- Lrnr_bayesglm$new()
  # fast highly adaptive lasso (HAL) implementation
  lrn_hal <- Lrnr_hal9001$new(max_degree = 2, num_knots = c(3,2), nfolds = 5)
  
  
  # Create learner library and instantiate super learner ensemble
  lrnr_lib <- Stack$new(fglm_lrnr, lrn_lasso, lrn_earth, lrn_xgb)
  sl_lrnr <- Lrnr_sl$new(learners = lrnr_lib, metalearner = Lrnr_nnls$new())
  
  
  nde_onestep <- medoutcon(
    W = dat_aug[, c("L1", "L2", "L3")],
    A = as.numeric(dat_aug$A),
    Z = NULL,
    M = dat_aug$M,
    Y = dat_aug$Y,
    obs_weights = dat_aug$w_norm,
    g_learners = sl_lrnr, # PS  # interaction here lasso of PS is not good
    h_learners = sl_lrnr, # modefied PS with Mediator.  hard to model P(A | M, L)
    b_learners = sl_lrnr, # Outcome regression
    effect = "direct",
    estimator = "onestep"
  )
  
  
  nie_onestep <- medoutcon(
    W = dat_aug[, c("L1", "L2", "L3")],
    A = as.numeric(dat_aug$A),
    Z = NULL,
    M = dat_aug$M,
    Y = dat_aug$Y,
    obs_weights = dat_aug$w_norm,
    g_learners = sl_lrnr,
    h_learners = sl_lrnr,
    b_learners = sl_lrnr,
    effect = "indirect",
    estimator = "onestep"
  )
  
  
  nde_os <- summary(nde_onestep)
  nie_os <- summary(nie_onestep)
  
  # Return the summaries as a named list
  return(list(
    nde_onestep_est = nde_os$param_est,
    nde_onestep_var = nde_os$var_est,
    nde_onestep_ci = c(nde_os$lwr_ci, nde_os$upr_ci),
    nie_onestep_est = nie_os$param_est,
    nie_onestep_var = nie_os$var_est,
    nie_onestep_ci = c(nie_os$lwr_ci, nie_os$upr_ci)
  ))
}
