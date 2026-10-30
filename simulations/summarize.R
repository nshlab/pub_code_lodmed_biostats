# Summarize simulation results.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript simulations/summarize.R RUN_DIR TRUE_NDE TRUE_NIE")
truth <- suppressWarnings(as.numeric(args[2:3]))
if (any(!is.finite(truth))) stop("Supply finite numerical truths on the risk-difference scale")
out <- normalizePath(args[[1L]], mustWork = TRUE)
meta <- readRDS(file.path(out, "manifest.rds"))
records <- readRDS(file.path(out, "results.rds"))
if (length(records) != meta$config$replicates || !all(vapply(records, function(x) isTRUE(x$ok), logical(1))))
  stop("Incomplete/failed experiment: inspect status.csv and resume first")
gcomp <- startsWith(meta$scenario, "study1_")
fields <- if (gcomp) c("replicate", "nde", "nde_adjusted_lower", "nde_adjusted_upper", "nde_adjusted_var",
  "nde_standard_lower", "nde_standard_upper", "nde_standard_var",
  "nie", "nie_adjusted_lower", "nie_adjusted_upper", "nie_adjusted_var",
  "nie_standard_lower", "nie_standard_upper", "nie_standard_var") else
  c("replicate", "nde", "nde_adjusted_lower", "nde_adjusted_upper", "nde_standard_lower", "nde_standard_upper", "nde_standard_var",
    "nie", "nie_adjusted_lower", "nie_adjusted_upper", "nie_standard_lower", "nie_standard_upper", "nie_standard_var", "selected_c")
if (any(vapply(records, function(x) length(x$values) < length(fields), logical(1)))) stop("Unexpected result schema")
dat <- as.data.frame(do.call(rbind, lapply(records, function(x) unname(x$values[seq_along(fields)]))))
names(dat) <- fields
if (any(!is.finite(as.matrix(dat)))) stop("Nonfinite results: resolve before summarizing")
utils::write.csv(dat, file.path(out, "estimates.csv"), row.names = FALSE)
summaries <- lapply(seq_along(truth), function(j) {
  effect <- c("nde", "nie")[[j]]
  estimate <- dat[[effect]]
  adj <- dat[[paste0(effect, "_adjusted_lower")]] <= truth[[j]] & dat[[paste0(effect, "_adjusted_upper")]] >= truth[[j]]
  standard <- dat[[paste0(effect, "_standard_lower")]] <= truth[[j]] & dat[[paste0(effect, "_standard_upper")]] >= truth[[j]]
  R <- length(estimate)
  data.frame(scenario = meta$scenario, effect = toupper(effect), truth = truth[[j]], replicates = R,
    mean_estimate = mean(estimate), bias = mean(estimate - truth[[j]]),
    empirical_variance = var(estimate), mse = mean((estimate - truth[[j]])^2),
    bias_mcse = sd(estimate)/sqrt(R), adjusted_coverage = mean(adj), standard_coverage = mean(standard),
    adjusted_coverage_mcse = sqrt(mean(adj)*(1-mean(adj))/R),
    standard_coverage_mcse = sqrt(mean(standard)*(1-mean(standard))/R),
    mean_standard_variance = mean(dat[[paste0(effect, "_standard_var")]]))
})
utils::write.csv(do.call(rbind, summaries), file.path(out, "summary.csv"), row.names = FALSE)
message("Saved estimates.csv and summary.csv. Original named vectors remain in results.rds.")
