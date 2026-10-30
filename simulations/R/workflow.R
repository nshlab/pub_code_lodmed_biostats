# Simulation orchestration only. Estimators remain in shared/ and scenarios/.
load_simulation_config <- function() {
  e <- new.env(parent = baseenv())
  sys.source(here::here("simulations", "config", "scenarios.R"), e)
  e$simulation_scenarios
}

# Load definitions, not exploratory top-level expressions or parallel plans.
load_simulation_definitions <- function(path, envir) {
  for (expr in parse(path, keep.source = FALSE)) {
    if (!is.call(expr)) next
    op <- as.character(expr[[1L]])
    if (identical(op, "source")) {
      dep <- eval(expr[[2L]], envir)
      load_simulation_definitions(dep, envir)
    } else if (length(op) == 1L && op %in% c("<-", "=")) {
      rhs <- expr[[3L]]
      if (is.call(rhs) && identical(rhs[[1L]], as.name("function"))) eval(expr, envir)
    }
  }
  invisible(envir)
}

positive_integer <- function(x, label) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x < 1 || x != floor(x) || x > .Machine$integer.max) stop(label, " must be a positive integer")
  as.integer(x)
}

simulation_preflight <- function(scenario) {
  cfg <- load_simulation_config()[[scenario]]
  if (is.null(cfg)) stop("Unknown scenario: ", scenario)
  packages <- unique(c("here", "future", "future.apply", "parallelly", cfg$packages))
  status <- data.frame(package = packages,
    installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE))
  print(status, row.names = FALSE)
  if (!all(status$installed)) stop("Missing packages: ", paste(status$package[!status$installed], collapse = ", "))
  invisible(cfg)
}

atomic_save_rds <- function(object, path) {
  temp <- tempfile(pattern = ".partial-", tmpdir = dirname(path))
  on.exit(unlink(temp), add = TRUE)
  saveRDS(object, temp)
  if (!file.rename(temp, path)) stop("Could not commit result: ", path)
  invisible(path)
}

run_simulation <- function(scenario, replicates = NULL, workers = 1L,
                           seed = 20261013L, run_id = "default", resume = TRUE) {
  cfg <- simulation_preflight(scenario)
  if (!is.null(replicates)) cfg$replicates <- positive_integer(replicates, "replicates")
  workers <- positive_integer(workers, "workers")
  seed <- positive_integer(seed, "seed")
  if (seed > .Machine$integer.max - cfg$replicates) stop("seed + replicates exceeds integer range")
  if (length(run_id) != 1L || !grepl("^[A-Za-z0-9][A-Za-z0-9_-]*$", run_id)) stop("Invalid run_id")
  if (workers > parallelly::availableCores()) stop("workers exceeds available cores")
  for (pkg in cfg$packages) library(pkg, character.only = TRUE)
  e <- new.env(parent = globalenv())
  e$SCENARIO_ID <- scenario
  load_simulation_definitions(here::here(cfg$file), e)
  # Intercept the historical append-write; each worker owns one checkpoint.
  e$write <- function(x, ...) { e$.result <- x; invisible(NULL) }
  out <- here::here("simulations", "results", "raw", scenario, run_id)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  lock <- file.path(out, ".running")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Run is already active (or has a stale .running lock): ", out)
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  src <- sort(c(here::here(cfg$file),
    list.files(here::here("simulations", "R"), "\\.R$", recursive = TRUE, full.names = TRUE),
    here::here("simulations", "config", "scenarios.R")))
  hashes <- tools::md5sum(src)
  names(hashes) <- substring(src, nchar(here::here()) + 2L)
  meta <- list(scenario = scenario, config = cfg, seed = seed,
               R_version = as.character(getRversion()),
               package_versions = vapply(unique(c("here", "future", "future.apply", "parallelly", cfg$packages)),
                 function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
               rng = "L'Ecuyer-CMRG, seed + replicate ID", source_md5 = hashes)
  manifest <- file.path(out, "manifest.rds")
  if (file.exists(manifest)) {
    if (!resume) stop("run_id already exists; choose a new run_id")
    if (!identical(meta, readRDS(manifest))) stop("Configuration/source changed; choose a new run_id")
  } else atomic_save_rds(meta, manifest)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  # Outer parallelism only: bootstrap tasks remain sequential inside a replicate.
  if (workers == 1L) future::plan(future::sequential)
  else future::plan(list(future::tweak(future::multisession, workers = workers), future::sequential))
  ids <- seq_len(cfg$replicates)
  checkpoint <- function(i) file.path(out, sprintf("replicate_%06d.rds", i))
  done <- vapply(ids, function(i) {
    if (!file.exists(checkpoint(i))) return(FALSE)
    x <- tryCatch(readRDS(checkpoint(i)), error = function(err) NULL)
    is.list(x) && identical(x$ok, TRUE) && identical(x$id, i)
  }, logical(1))
  message(scenario, ": ", sum(done), " cached; ", sum(!done), " remaining")
  future.apply::future_lapply(ids[!done], function(i) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed + i)
    warnings <- character()
    e$.result <- NULL
    ans <- tryCatch(withCallingHandlers({
      do.call(e$sim_iteration, c(list(i = i), cfg[c("n_obs", "k", "type", "censor_rate", "S", "beta_m")]))
      if (!is.numeric(e$.result) || !length(e$.result) || e$.result[[1L]] != i) stop("No valid result vector captured")
      list(id = i, ok = TRUE, values = e$.result)
    }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(err) list(id = i, ok = FALSE, error = conditionMessage(err)))
    ans$warnings <- warnings
    atomic_save_rds(ans, checkpoint(i))
    invisible(NULL)
  }, future.seed = TRUE)
  records <- lapply(ids, function(i) readRDS(checkpoint(i)))
  status <- data.frame(replicate = ids,
    ok = vapply(records, function(x) isTRUE(x$ok), logical(1)),
    message = vapply(records, function(x) if (isTRUE(x$ok)) "" else x$error, character(1)))
  utils::write.csv(status, file.path(out, "status.csv"), row.names = FALSE)
  atomic_save_rds(records, file.path(out, "results.rds"))
  if (!all(status$ok)) stop(sum(!status$ok), " replicates failed; see status.csv. Resume retries failed replicates.")
  message("Completed: ", out)
  invisible(records)
}
