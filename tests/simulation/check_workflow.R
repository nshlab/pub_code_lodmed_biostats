# Orchestration integration test with a stochastic stand-in, not an estimator test.
# Requires only here, future, future.apply, parallelly.
if (!requireNamespace("here", quietly = TRUE)) stop("here is required")
here::i_am("tests/simulation/check_workflow.R")
source(here::here("simulations", "R", "workflow.R"))
cfg <- load_simulation_config()[[1L]]
cfg$replicates <- 3L; cfg$packages <- character()
# Avoid statistical dependencies while exercising the real workflow/checkpoints.
simulation_preflight <- function(scenario) cfg
load_simulation_definitions <- function(path, envir) {
  envir$sim_iteration <- eval(quote(function(i, ...) write(c(i, runif(2)), file = "unused")), envir)
  invisible(envir)
}
id <- paste0("workflow_test_", Sys.getpid())
scenario <- names(load_simulation_config())[[1L]]
out <- here::here("simulations", "results", "raw", scenario, id)
tryCatch({
  a <- run_simulation(scenario, seed = 12L, run_id = id)
  b <- run_simulation(scenario, seed = 12L, run_id = id)
  stopifnot(identical(a, b), all(vapply(a, function(x) x$ok, logical(1))))
  unlink(file.path(out, "replicate_000002.rds"))
  c <- run_simulation(scenario, seed = 12L, run_id = id)
  stopifnot(identical(a, c))
  rejected <- tryCatch({run_simulation(scenario, seed = 13L, run_id = id); FALSE}, error = function(e) TRUE)
  stopifnot(rejected)
  if (parallelly::availableCores() >= 2L) {
    unlink(file.path(out, sprintf("replicate_%06d.rds", 1:3)))
    d <- run_simulation(scenario, seed = 12L, workers = 2L, run_id = id)
    stopifnot(identical(a, d))
  }
  message("Workflow checkpoint/resume/seed checks passed (stand-in estimator).")
}, finally = unlink(out, recursive = TRUE))
