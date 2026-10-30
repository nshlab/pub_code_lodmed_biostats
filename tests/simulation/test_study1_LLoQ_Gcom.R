if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
source(here::here("tests", "simulation", "helpers.R"), local = TRUE)

res <- run_static_scenario_regression(
  legacy_path = here::here("tests", "fixtures", "simulation", "Study1/LLoQ_Gcom.R.txt"),
  refactored_path = here::here("simulations", "scenarios", "study1/LLoQ_Gcom.R"),
  n_obs = 200L,
  censor_rate = 0.25,
  seed = 20261005L
)

print(res$function_comparison)
cat("Intentionally excluded from static body equality:", paste(res$intentionally_excluded_functions, collapse = ", "), "\n")
print(res$dgp_comparison)

if (!isTRUE(res$all_common_function_bodies_identical)) {
  stop("At least one common function body differs between legacy and refactored definitions.")
}
if (!isTRUE(res$dgp_identical)) {
  stop("DGP output differs between legacy and refactored scenario under the same seed.")
}

message("PASS: study1_LLoQ_Gcom")
