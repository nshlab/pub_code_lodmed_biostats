# V6 fast numerical smoke tests for refactored sim_iteration() wrappers.
#
# Purpose: verify that legacy and refactored orchestration produce the same
# numerical result vector under the same DGP draw, while avoiding thousands of
# bootstrap/SL fits. Core statistical functions and the DGP are checked
# separately and strictly by tests/simulation/run_all.R.

source(here::here("tests", "simulation", "helpers.R"))

install_capture_write <- function(env) {
  env$.captured_result <- NULL
  env$write <- function(x, file = "", ncolumns = if (is.character(x)) 1 else 5,
                        append = FALSE, sep = " ") {
    env$.captured_result <- x
    invisible(NULL)
  }
}

install_fast_stubs <- function(env, study) {
  # Keep the real DGP. Stub only expensive downstream fitting/bootstrap pieces.
  env$EM_algorithm <- function(data, beta_m = NULL, beta_m_0 = NULL, S, LOD, ...) {
    # Legacy scenarios are not uniform: some call beta_m= and HAL scenarios
    # call beta_m_0=. Accept both while returning the same deterministic stub.
    beta_used <- if (!is.null(beta_m)) beta_m else beta_m_0
    list(data_weted = data, beta_m = beta_used)
  }
  env$NaiveImpute <- function(dataset, LOD, ...) dataset

  if (identical(study, "study1")) {
    env$gcomp.nde <- function(data, ...) list(nde.rd = 0.111)
    env$gcomp.nie <- function(data, ...) list(nie.rd = -0.022)
    env$choose_ci <- function(data, p_cen, k, type, measure, ...) {
      list(opt_m = max(2L, floor(nrow(data) * 0.8)))
    }
    env$gcom_boots_ci_m_out_n <- function(data, m, n_bootstrap, type, measure, ...) {
      if (identical(type, "NDE")) list(ci = c(0.10, 0.12), variance = 0.0011)
      else list(ci = c(-0.03, -0.01), variance = 0.0022)
    }
    env$gcom_boots_ci <- function(data, n_bootstrap, type, measure, ...) {
      if (identical(type, "NDE")) list(ci = c(0.09, 0.13), variance = 0.0033)
      else list(ci = c(-0.04, 0.00), variance = 0.0044)
    }
  } else if (identical(study, "study2")) {
    env$Onestep_med_effs <- function(dat_aug, effect, ...) {
      n <- nrow(dat_aug)
      if (identical(effect, "NDE")) {
        list(nde_onestep_est = 0.111, eif = rep(0, n),
             nde_onestep_ci_sd = c(0.10, 0.12), nde_onestep_var = 0.0011)
      } else {
        list(nie_onestep_est = -0.022, eif = rep(0, n),
             nie_onestep_ci_sd = c(-0.03, -0.01), nie_onestep_var = 0.0022)
      }
    }
    env$os_choose_ci <- function(data, original_estimate, p_cen, k, type, ...) {
      list(opt_m = max(2L, floor(nrow(data) * 0.8)), optimal_c = 0.75)
    }
    env$moutn_multiplier_bootstrap <- function(eif, theta, B, m, ...) {
      list(ci = theta + c(-0.01, 0.01))
    }
  } else stop("Unknown study: ", study)

  install_capture_write(env)
  invisible(env)
}

load_smoke_envs <- function(legacy_path, refactored_path) {
  old_env <- new.env(parent = globalenv())
  new_env <- new.env(parent = globalenv())
  source_legacy_definitions(legacy_path, old_env)
  source_refactored_scenario(refactored_path, new_env)
  source_dgp_definition(here::here("tests", "fixtures", "simulation", "DGP", "DGP.R.txt"), old_env)
  source_dgp_definition(here::here("simulations", "R", "shared", "dgp.R"), new_env)
  list(old = old_env, new = new_env)
}

run_sim_iteration_smoke <- function(legacy_path, refactored_path, study,
                                    seed = 20261005L, n_obs = 80L,
                                    censor_rate = 0.50, k = 2L, S = 2L,
                                    beta_m = c(b0=-2.75, A=1.75, L1=2,
                                               AL1=0.25, L2=1.25, L3=-0.35,
                                               sd=0.35)) {
  envs <- load_smoke_envs(legacy_path, refactored_path)
  install_fast_stubs(envs$old, study)
  install_fast_stubs(envs$new, study)

  # Refactored wrappers intentionally use SCENARIO_ID only to construct a
  # portable output path. Supply it explicitly in both environments so that
  # path refactoring is not confused with a numerical change.
  scenario_id <- tools::file_path_sans_ext(basename(refactored_path))
  envs$old$SCENARIO_ID <- scenario_id
  envs$new$SCENARIO_ID <- scenario_id

  # EM-MLE is used because several legacy Study 1 wrappers append beta_mres,
  # which is only created on that branch.
  set.seed(seed)
  envs$old$sim_iteration(1L, n_obs=n_obs, k=k, type="EM-MLE",
                         censor_rate=censor_rate, S=S, beta_m=beta_m)
  old_result <- envs$old$.captured_result

  set.seed(seed)
  envs$new$sim_iteration(1L, n_obs=n_obs, k=k, type="EM-MLE",
                         censor_rate=censor_rate, S=S, beta_m=beta_m)
  new_result <- envs$new$.captured_result

  if (is.null(old_result) || is.null(new_result)) stop("Smoke test did not capture a result vector.")
  same <- isTRUE(all.equal(old_result, new_result, tolerance=0, check.attributes=TRUE))
  list(passed=same, legacy_result=old_result, refactored_result=new_result,
       max_abs_diff=if (length(old_result)==length(new_result)) max(abs(as.numeric(old_result)-as.numeric(new_result)), na.rm=TRUE) else Inf)
}
