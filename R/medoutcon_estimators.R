# Consolidated from onestep.R; estimator choice is now explicit.
# Wrapping nadir::super_learner() in an sl3 Lrnr template

# In this script, we define a new sl3 Lrnr as a wrapper around
# nadir::super_learner() so that software that uses sl3 may easily
# extend to using nadir::super_learner()

# dependencies
library(nadir)
library(sl3)
library(R6)
library(data.table) # Assuming 'data' is a data.table
library(dplyr)
library(sl3)
library(hal9001)
library(medoutcon)
library(truncnorm)
library(xgboost)
library(future.apply)
library(dplyr) 
library(ranger)
library(nadir)
library(earth)

# new class : nadir_super_learner_sl3_Lrnr --------------------------------
#
#   * Inherits from the sl3::Lrnr_base
#   * Initialize a member of this class via calling
#     `nadir_sl_Lrnr <- make_learner(nadir_super_learner_sl3_Lrnr,
#        ...)`
#   * At initialization, any extra parameters given will be passed
#     on to nadir::super_learner,
#     so the following should be passed:
#       * formulas (optional: if not specified, it is constructed from
#         task$nodes$outcome and task$nodes$covariates)
#       * learners = a list of nadir learners -- see ?nadir::learners
#       * extra_learner_args (optional) = a list of extra arguments for
#         individual learners.
#       * any other arguments nadir::super_learner() needs like
#         outcome_type, etc.
#   * Fit via running on a task:
#     `nadir_sl_Lrnr_fit <- nadir_sl_Lrnr$train(task)` given an sl3 task
#   * Predict after fitting:
#     `nadir_sl_Lrnr_fit$predict(task)`
#
nadir_super_learner_sl3_Lrnr <- R6Class(
  classname = "nadir_sl",
  inherit = sl3::Lrnr_base,
  portable = TRUE,
  class = TRUE,
  
  public = list(
    
    initialize = function(...) {
      # this captures all parameters to initialize and saves them as self$params
      params <- sl3::args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    # list properties your learner supports here.
    # Use sl3_list_properties() for a list of options
    .properties = c("binomial", 'categorical', 'density', 'ids', 'weights'),
    
    # .train takes task data and returns a fit object that can be used to generate predictions
    .train = function(task) {
      # generate an argument list from the parameters that were
      # captured when your learner was initialized.
      # this allows users to pass arguments directly to your ml function
      args <- self$params
      
      # get outcome variable type
      # preferring learner$params$outcome_type first, then task$outcome_type
      # args$outcome_type <- self$get_outcome_type(task)
      # should pass something on to your learner indicating outcome_type
      # e.g. family or objective
      
      # specific to nadir super learning:
      args$data <- task$get_data()
      if (is.null(args$formulas)) {
        args$formulas <- as.formula(paste0(task$nodes$outcome, " ~ ", paste0(task$nodes$covariates, collapse = " + ")))
      }
      
      # only add arguments on weights and offset
      # if those were specified when the task was generated
      # if (task$has_node("weights")) {
      #   args$weights <- task$weights
      # }
      
      # if (task$has_node("id")) {
      #   args$cluster_ids <- task$id
      # }
      
      # store the fit super learner
      fit_nadir_super_learner <- sl3:::call_with_args(nadir::super_learner, args)
      
      # return the fit object, which will be stored
      # in a learner object and returned from the call
      # to learner$predict
      return(fit_nadir_super_learner)
    },
    
    # .predict takes a task and returns predictions from that task
    .predict = function(task = NULL) {
      self$training_task
      self$training_outcome_type
      self$fit_object
      
      predictions <- predict(self$fit_object, task$X)
      return(predictions)
    },
    # list any packages required for your learner here.
    .required_packages = c("nadir")
  )
)



# nadir learner wrapper ---------------------------------------------

nadir_learner_sl3_Lrnr <- R6Class(
  classname = "nadir_sl",
  inherit = sl3::Lrnr_base,
  portable = TRUE,
  class = TRUE,
  
  public = list(
    
    initialize = function(...) {
      # this captures all parameters to initialize and saves them as self$params
      
      params <- sl3::args_to_list()
      super$initialize(params = params, ...)
    }
  ),
  private = list(
    # list properties your learner supports here.
    # Use sl3_list_properties() for a list of options
    .properties = c("binomial", 'categorical', 'density', 'ids', 'weights'),
    
    # .train takes task data and returns a fit object that can be used to generate predictions
    .train = function(task) {
      
      # generate an argument list from the parameters that were
      # captured when your learner was initialized.
      # this allows users to pass arguments directly to your ml function
      learner <- self$params$learner
      stopifnot(!is.null(learner))
      # self$params$learner <- NULL
      args <- self$params
      args[['learner']] <- NULL
      
      
      # get outcome variable type
      # preferring learner$params$outcome_type first, then task$outcome_type
      # args$outcome_type <- self$get_outcome_type(task)
      # should pass something on to your learner indicating outcome_type
      # e.g. family or objective
      
      # specific to nadir super learning:
      args$data <- task$get_data()
      if (is.null(args$formula)) {
        args$formula <- as.formula(paste0(task$nodes$outcome, " ~ ", paste0(task$nodes$covariates, collapse = " + ")))
      }
      
      
      
      # only add arguments on weights and offset
      # if those were specified when the task was generated
      # if (task$has_node("weights")) {
      #   args$weights <- task$weights
      # }
      
      # if (task$has_node("id")) {
      #   args$cluster_ids <- task$id
      # }
      
      # store the fit super learner
      fit_learner <- do.call(learner, args)
      
      
      # return the fit object, which will be stored
      # in a learner object and returned from the call
      # to learner$predict
      return(fit_learner)
    },
    
    # .predict takes a task and returns predictions from that task
    .predict = function(task = NULL) {
      self$training_task
      self$training_outcome_type
      self$fit_object
      
      predictions <- self$fit_object(task$X)
      return(predictions)
    },
    # list any packages required for your learner here.
    .required_packages = c("nadir")
  )
)




#######
#m-out-of-n bootstrap confidence intervals
#######
library(future.apply)



Onestep_med_effs <- function(dat_aug, effect = c("both", "NDE", "NIE"), L_vars = NULL, M_vars = NULL, estimator = c("onestep", "tmle")) {
  estimator <- match.arg(estimator)
  effect <- match.arg(effect)
  
  # Auto-detect covariates if not provided
  if (is.null(L_vars)) L_vars <- get_L_vars(dat_aug)
  if (is.null(M_vars)) M_vars <- get_M_vars(dat_aug)
  
  # Extract W and M as data frames/matrices
  if (is.data.table(dat_aug)) {
    W_data <- dat_aug[, ..L_vars]
    M_data <- dat_aug[, ..M_vars]
  } else {
    W_data <- dat_aug[, L_vars, drop = FALSE]
    M_data <- dat_aug[, M_vars, drop = FALSE]
  }
  
  # Ensure M_data is a data frame
  if (!is.data.frame(M_data)) {
    M_data <- as.data.frame(M_data)
  }
  
  # Drop zero-variance and linear-combo columns from W_data automatically
  W_df <- as.data.frame(W_data)
  keep_nonconst <- sapply(W_df, function(x) length(unique(na.omit(x))) > 1)
  if (any(!keep_nonconst)) {
    W_df <- W_df[, keep_nonconst, drop = FALSE]
  }
  if (requireNamespace("caret", quietly = TRUE)) {
    lc <- tryCatch(caret::findLinearCombos(W_df), error = function(e) NULL)
    if (!is.null(lc$remove) && length(lc$remove) > 0) {
      dropped <- names(W_df)[lc$remove]
      W_df <- W_df[, -lc$remove, drop = FALSE]
      warning("Dropped linear-combo columns from W: ", paste(dropped, collapse = ", "))
    }
  }
  W_data <- as.matrix(W_df)
  
  # arguments for the nadir::super_learner()
  learners <- list(
    logistic = nadir::lnr_glm,
    gam = nadir::lnr_gam,
    earth = nadir::lnr_earth,
    ranger = nadir::lnr_ranger,
    xgb = nadir::lnr_xgboost,
    ridge = nadir::lnr_glmnet,
    lasso = nadir::lnr_glmnet
  )
  
  
  extra_learner_args <- list(
    ridge = list(alpha = 0),
    lasso = list(alpha = 1))
  
  # instantiate a nadir_sl_Lrnr
  nadir_sl_Lrnr_binary <- make_learner(learner_class = nadir_super_learner_sl3_Lrnr,
                                       learners = learners,
                                       extra_learner_args = extra_learner_args,
                                       outcome_type = 'binary')
  lrnr_nadir_gam <- make_learner(learner_class = nadir_learner_sl3_Lrnr, 
                                 learner = nadir::lnr_gam, 
                                 family = binomial(link = 'logit')
  )
  
  lrnr_nadir_earth <- make_learner(learner_class = nadir_learner_sl3_Lrnr, 
                                   learner = nadir::lnr_earth
  )
  
  result_list <- list()
  
  # ---- NDE ----
  if (effect %in% c("both", "NDE")) {
    nde_onestep <- medoutcon(
      W = as.matrix(W_data),
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = as.matrix(M_data),
      Y = as.numeric(dat_aug$Y_scale),
      g_learners = lrnr_nadir_gam,
      h_learners = lrnr_nadir_gam,
      b_learners = lrnr_nadir_gam,
      obs_weights = dat_aug$w_norm,
      effect = "direct",
      estimator = estimator
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
      W = as.matrix(W_data),
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = as.matrix(M_data),
      Y = as.numeric(dat_aug$Y_scale),
      obs_weights = dat_aug$w_norm,
      g_learners = lrnr_nadir_gam,
      h_learners = lrnr_nadir_gam,
      b_learners = lrnr_nadir_gam,
      effect = "indirect",
      estimator = estimator
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

Onestep_med_effs_fast <- function(dat_aug, effect = c("both", "NDE", "NIE"), L_vars = NULL, M_vars = NULL, estimator = c("onestep", "tmle")) {
  estimator <- match.arg(estimator)
  effect <- match.arg(effect)
  
  # Auto-detect covariates if not provided
  if (is.null(L_vars)) L_vars <- get_L_vars(dat_aug)
  if (is.null(M_vars)) M_vars <- get_M_vars(dat_aug)
  
  # Extract W and M as data frames/matrices
  if (is.data.table(dat_aug)) {
    W_data <- dat_aug[, ..L_vars]
    M_data <- dat_aug[, ..M_vars]
  } else {
    W_data <- dat_aug[, L_vars, drop = FALSE]
    M_data <- dat_aug[, M_vars, drop = FALSE]
  }
  
  # Ensure M_data is a data frame
  if (!is.data.frame(M_data)) {
    M_data <- as.data.frame(M_data)
  }
  
  # Drop zero-variance and linear-combo columns from W_data automatically
  W_df <- as.data.frame(W_data)
  keep_nonconst <- sapply(W_df, function(x) length(unique(na.omit(x))) > 1)
  if (any(!keep_nonconst)) {
    W_df <- W_df[, keep_nonconst, drop = FALSE]
  }
  if (requireNamespace("caret", quietly = TRUE)) {
    lc <- tryCatch(caret::findLinearCombos(W_df), error = function(e) NULL)
    if (!is.null(lc$remove) && length(lc$remove) > 0) {
      dropped <- names(W_df)[lc$remove]
      W_df <- W_df[, -lc$remove, drop = FALSE]
      warning("Dropped linear-combo columns from W: ", paste(dropped, collapse = ", "))
    }
  }
  W_data <- as.matrix(W_df)
  
  # Learners - removed nfolds parameter from HAL
  lrn_earth <- Lrnr_earth$new()
  lrn_hal <- Lrnr_hal9001$new(max_degree = 2, num_knots = c(3, 2))  # Removed nfolds
  lrn_ranger <- Lrnr_ranger$new()
  
  lrnr_lib <- Stack$new(lrn_hal, lrn_ranger, lrn_earth)
  
  result_list <- list()
  
  # ---- NDE ----
  if (effect %in% c("both", "NDE")) {
    nde_onestep <- medoutcon(
      W = as.matrix(W_data),
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = as.matrix(M_data),
      Y = as.numeric(dat_aug$Y),
      obs_weights = dat_aug$w_norm,
      effect = "direct",
      estimator = estimator,
      estimator_args = list(cv_folds = 5L, max_iter = 2L, tiltmod_tol = 2)
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
      W = as.matrix(W_data),
      A = as.numeric(dat_aug$A),
      Z = NULL,
      M = as.matrix(M_data),
      Y = as.numeric(dat_aug$Y_scale),
      obs_weights = dat_aug$w_norm,
      effect = "indirect",
      estimator = estimator,
      estimator_args = list(cv_folds = 5L, max_iter = 2L, tiltmod_tol = 2)
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

onestep_analysis <- function(data_aug, k = 5, censor_rate, B1 = 25, B2 = 50000, B_bootstrap = 50000, estimator = c("onestep", "tmle"), L_vars = NULL) {
  estimator <- match.arg(estimator)
  
  # ---- NDE Analysis ----
  onestep.nde <- Onestep_med_effs(dat_aug = data_aug, effect = "NDE", L_vars = L_vars, estimator = estimator)
  
  # Skip if NDE estimation failed
  if (is.na(onestep.nde$nde_onestep_est)) {
    nde_results <- list(
      point_est = NA_real_,
      se_onestep = NA_real_,
      se_bootstrap = NA_real_,
      ci_onestep = c(NA_real_, NA_real_),
      ci_bootstrap = c(NA_real_, NA_real_),
      pvalue_onestep = NA_real_,
      pvalue_bootstrap = NA_real_
    )
  } else {
    # Choose optimal m using double bootstrap
    choose_ciNDE <- os_choose_ci(
      data = data_aug, 
      original_estimate = onestep.nde$nde_onestep_est, 
      p_cen = censor_rate, 
      k = k, 
      B1 = B1,
      B2 = B2,
      type = "NDE"
    )
    
    # Bootstrap CI using m-out-of-n
    CI.nderes <- moutn_multiplier_bootstrap(
      eif = onestep.nde$eif, 
      theta = onestep.nde$nde_onestep_est, 
      B = B_bootstrap, 
      m = choose_ciNDE$opt_m
    )
    
    # Extract SEs from onestep summary and bootstrap
    se_onestep <- sqrt(onestep.nde$nde_onestep_var)
    se_bootstrap <- CI.nderes$se
    
    # Compute p-values
    pvalue_onestep <- compute_pvalue(onestep.nde$nde_onestep_est, se_onestep)
    pvalue_bootstrap <- compute_pvalue(onestep.nde$nde_onestep_est, se_bootstrap)
    
    nde_results <- list(
      point_est = onestep.nde$nde_onestep_est,
      se_onestep = se_onestep,
      se_bootstrap = se_bootstrap,
      ci_onestep = onestep.nde$nde_onestep_ci,
      ci_bootstrap = CI.nderes$ci,
      pvalue_onestep = pvalue_onestep,
      pvalue_bootstrap = pvalue_bootstrap,
      optimal_m = choose_ciNDE$opt_m,
      optimal_c = choose_ciNDE$optimal_c
    )
  }
  
  # ---- NIE Analysis ----
  onestep.nie <- Onestep_med_effs(dat_aug = data_aug, effect = "NIE", L_vars = L_vars, estimator = estimator)
  
  # Skip if NIE estimation failed
  if (is.na(onestep.nie$nie_onestep_est)) {
    nie_results <- list(
      point_est = NA_real_,
      se_onestep = NA_real_,
      se_bootstrap = NA_real_,
      ci_onestep = c(NA_real_, NA_real_),
      ci_bootstrap = c(NA_real_, NA_real_),
      pvalue_onestep = NA_real_,
      pvalue_bootstrap = NA_real_
    )
  } else {
    # Use same m as NDE for consistency
    opt_m <- if (!is.na(nde_results$optimal_m)) nde_results$optimal_m else floor(nrow(data_aug)^0.5)
    
    # Bootstrap CI using m-out-of-n
    CI.nieres <- moutn_multiplier_bootstrap(
      eif = onestep.nie$eif, 
      theta = onestep.nie$nie_onestep_est, 
      B = B_bootstrap, 
      m = opt_m
    )
    
    # Extract SEs from onestep summary and bootstrap
    se_onestep <- sqrt(onestep.nie$nie_onestep_var)
    se_bootstrap <- CI.nieres$se
    
    # Compute p-values
    pvalue_onestep <- compute_pvalue(onestep.nie$nie_onestep_est, se_onestep)
    pvalue_bootstrap <- compute_pvalue(onestep.nie$nie_onestep_est, se_bootstrap)
    
    nie_results <- list(
      point_est = onestep.nie$nie_onestep_est,
      se_onestep = se_onestep,
      se_bootstrap = se_bootstrap,
      ci_onestep = onestep.nie$nie_onestep_ci,
      ci_bootstrap = CI.nieres$ci,
      pvalue_onestep = pvalue_onestep,
      pvalue_bootstrap = pvalue_bootstrap,
      optimal_m = opt_m
    )
  }
  
  # Return comprehensive results
  return(list(
    NDE = nde_results,
    NIE = nie_results,
    summary_table = data.frame(
      Effect = c("NDE", "NIE"),
      Point_Est = c(nde_results$point_est, nie_results$point_est),
      SE_Onestep = c(nde_results$se_onestep, nie_results$se_onestep),
      SE_Bootstrap = c(nde_results$se_bootstrap, nie_results$se_bootstrap),
      CI_Lower_Onestep = c(nde_results$ci_onestep[1], nie_results$ci_onestep[1]),
      CI_Upper_Onestep = c(nde_results$ci_onestep[2], nie_results$ci_onestep[2]),
      CI_Lower_Bootstrap = c(nde_results$ci_bootstrap[1], nie_results$ci_bootstrap[1]),
      CI_Upper_Bootstrap = c(nde_results$ci_bootstrap[2], nie_results$ci_bootstrap[2]),
      Pvalue_Onestep = c(nde_results$pvalue_onestep, nie_results$pvalue_onestep),
      Pvalue_Bootstrap = c(nde_results$pvalue_bootstrap, nie_results$pvalue_bootstrap)
    )
  ))
}
