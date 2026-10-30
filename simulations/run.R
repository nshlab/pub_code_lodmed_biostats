# Run a simulation scenario.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
here::i_am("simulations/run.R")
source(here::here("simulations", "R", "workflow.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || identical(args[[1L]], "--list")) {
  cat(paste(names(load_simulation_config()), collapse = "\n"), "\n")
} else {
  scenario <- args[[1L]]
  flags <- args[-1L]
  if (any(!grepl("^--(replicates|workers|seed|run-id|preflight)=", flags))) stop("Unknown flag")
  keys <- sub("=.*$", "", flags)
  if (anyDuplicated(keys)) stop("Duplicate flag")
  opt <- setNames(sub("^[^=]*=", "", flags), sub("^--", "", keys))
  if ("preflight" %in% names(opt) && opt[["preflight"]] != "true") stop("Use --preflight=true")
  if ("preflight" %in% names(opt)) simulation_preflight(scenario)
  else {
    number <- function(key, default) if (key %in% names(opt)) as.numeric(opt[[key]]) else default
    run_simulation(scenario, replicates = number("replicates", NULL),
      workers = number("workers", 1L), seed = number("seed", 20261013L),
      run_id = if ("run-id" %in% names(opt)) opt[["run-id"]] else "default")
  }
}
