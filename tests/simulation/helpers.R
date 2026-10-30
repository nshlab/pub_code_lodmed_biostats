# Utilities for legacy-vs-refactored simulation regression tests (v5).
# Static checks load function definitions only: top-level executable code, package
# attachment, parallel plans, and Monte Carlo tails are deliberately not run.

strip_install_lines <- function(lines) lines[!grepl("install_github\\(", lines)]

read_definition_block <- function(path) {
  lines <- strip_install_lines(readLines(path, warn = FALSE))
  marker <- grep("^parallelly::availableCores\\(\\)\\s*$", lines)
  if (length(marker) == 1L) lines <- lines[seq_len(marker - 1L)]
  lines
}

# Evaluate only top-level assignments whose RHS is a function. This avoids
# historical top-level statements such as data$lod, plan(multisession), and
# library(Rsolnp), which are irrelevant to source-code equivalence and caused
# misleading test failures in v4.
source_function_definitions <- function(path, envir, legacy = FALSE) {
  lines <- if (legacy) read_definition_block(path) else readLines(path, warn = FALSE)
  exprs <- parse(text = lines, keep.source = FALSE)
  for (expr in exprs) {
    if (!is.call(expr) || !as.character(expr[[1L]]) %in% c("<-", "=")) next
    rhs <- expr[[3L]]
    if (!is.call(rhs) || !identical(as.character(rhs[[1L]]), "function")) next
    eval(expr, envir = envir)
  }
  invisible(envir)
}

source_legacy_definitions <- function(path, envir) {
  source_function_definitions(path, envir, legacy = TRUE)
}

# Refactored scenarios delegate shared/variant functions to these files. Load
# function definitions from the scenario and all repository-local source()
# dependencies without executing package loads or top-level scenario code.
source_refactored_scenario <- function(path, envir) {
  lines <- readLines(path, warn = FALSE)
  src <- grep('^source\\(here::here\\(', lines, value = TRUE)
  for (s in src) {
    parts <- regmatches(s, gregexpr('"[^"]+"', s))[[1L]]
    if (!length(parts)) next
    parts <- gsub('^"|"$', '', parts)
    dep <- do.call(here::here, as.list(parts))
    if (file.exists(dep)) source_function_definitions(dep, envir, legacy = FALSE)
  }
  source_function_definitions(path, envir, legacy = FALSE)
  invisible(envir)
}

function_fingerprint <- function(fun) paste(deparse(body(fun), width.cutoff = 500L), collapse = "\n")

compare_named_functions <- function(old_env, new_env, ignore = character()) {
  common <- setdiff(intersect(ls(old_env, all.names = TRUE), ls(new_env, all.names = TRUE)), ignore)
  common_fun <- common[vapply(common, function(x) is.function(get(x, old_env)) && is.function(get(x, new_env)), logical(1))]
  data.frame(function_name = common_fun,
             identical_body = vapply(common_fun, function(x) identical(function_fingerprint(get(x, old_env)), function_fingerprint(get(x, new_env))), logical(1)),
             stringsAsFactors = FALSE)
}

source_dgp_definition <- function(path, envir) source_function_definitions(path, envir, legacy = FALSE)

compare_dgp_once <- function(old_env, new_env, n_obs, censor_rate, seed = 20261005L) {
  old_exists <- exists("Study_dgp_Mcensor", envir = old_env, inherits = FALSE)
  new_exists <- exists("Study_dgp_Mcensor", envir = new_env, inherits = FALSE)
  old_fun_ok <- old_exists && is.function(get("Study_dgp_Mcensor", old_env, inherits = FALSE))
  new_fun_ok <- new_exists && is.function(get("Study_dgp_Mcensor", new_env, inherits = FALSE))
  if (!old_fun_ok || !new_fun_ok) stop("DGP loading failure")
  old_fun <- get("Study_dgp_Mcensor", old_env, inherits = FALSE)
  new_fun <- get("Study_dgp_Mcensor", new_env, inherits = FALSE)
  set.seed(seed); old <- old_fun(n_obs = n_obs, censor_rate = censor_rate)
  set.seed(seed); new <- new_fun(n_obs = n_obs, censor_rate = censor_rate)
  list(legacy_dgp_found=old_exists, refactored_dgp_found=new_exists,
       legacy_dgp_is_function=old_fun_ok, refactored_dgp_is_function=new_fun_ok,
       dgp_function_body_identical=identical(function_fingerprint(old_fun), function_fingerprint(new_fun)),
       lod_identical=isTRUE(all.equal(old$lod,new$lod,tolerance=0)),
       observed_data_identical=isTRUE(all.equal(old$study_dat,new$study_dat,tolerance=0)),
       full_data_identical=isTRUE(all.equal(old$study_dat_full,new$study_dat_full,tolerance=0)))
}

run_static_scenario_regression <- function(legacy_path, refactored_path, n_obs, censor_rate,
                                           seed=20261005L,
                                           ignore_functions=c("sim_iteration")) {
  old_env <- new.env(parent=globalenv()); new_env <- new.env(parent=globalenv())
  source_legacy_definitions(legacy_path, old_env)
  source_refactored_scenario(refactored_path, new_env)
  source_dgp_definition(here::here("tests","fixtures","simulation","DGP","DGP.R.txt"), old_env)
  source_dgp_definition(here::here("simulations","R","shared","dgp.R"), new_env)
  fn_cmp <- compare_named_functions(old_env,new_env,ignore=ignore_functions)
  dgp_cmp <- compare_dgp_once(old_env,new_env,n_obs,censor_rate,seed)
  list(function_comparison=fn_cmp, dgp_comparison=dgp_cmp,
       intentionally_excluded_functions=ignore_functions,
       all_common_function_bodies_identical=if(nrow(fn_cmp)) all(fn_cmp$identical_body) else NA,
       dgp_identical=all(unlist(dgp_cmp)))
}

check_optional_dependencies <- function(packages) {
  data.frame(package=packages, installed=vapply(packages, requireNamespace, logical(1), quietly=TRUE), row.names=NULL)
}
