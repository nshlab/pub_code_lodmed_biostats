# Explicit application tasks, run from the repository root.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")
here::i_am("analysis/run.R")
root <- here::here()
if (!identical(normalizePath(getwd()), normalizePath(root)))
  stop("Run this command from the repository root.")
args <- commandArgs(trailingOnly = TRUE)
task <- if (length(args)) args[[1L]] else "--help"
if (length(args) > 1L) stop("Supply one task only")
if (task %in% c("--help", "--list")) {
  cat("Application tasks:\n",
      "  table2       Prepared CSVs -> NDE/NIE point estimates -> comparison/export\n",
      "  descriptive  Authorized observed-data CSV -> Table S1 draft and QC\n",
      "  figure2      Same authorized observed data -> descriptive RNA figures and count QC\n",
      "  qqplot       Conditional Gaussian RNA diagnostics (authorized observed data)\n",
      "Usage: Rscript analysis/run.R TASK\n", sep = "")
} else {
  scripts <- switch(task,
    table2 = c("05_record_session.R", "03_table2_point_estimates.R", "04_export_table2.R"),
    descriptive = "01_descriptive_analysis.R",
    figure2 = "02_figure2.R",
    qqplot = "02a_qplot_conditional.R",
    stop("Unknown task: ", task))
  for (script in scripts) source(file.path(root, "analysis", script), local = TRUE)
}
