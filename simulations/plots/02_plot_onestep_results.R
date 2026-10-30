# Reproduce the one-step simulation summary plots from precomputed results.

required <- c("here", "dplyr", "purrr", "readr", "stringr", "tibble", "tidyr",
              "ggplot2", "latex2exp", "MetBrewer")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

source(here::here("simulations", "plots", "plot_helpers.R"))

input_dir <- here::here("simulations", "results", "precomputed", "onestep")
output_dir <- here::here("simulations", "results", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
  Trad1step = "LLoQ/2 imputation",
  `EM-heterHAL1step` = "EM (4): heterosced. CDE, HAL",
  `EM-correctGLM1step` = "EM (1): true nuisances"
)

# Preserve the original plotting decision to omit n = 3000.
full_data <- load_precomputed_results(
  input_dir, method_labels, exclude_pattern = "n_3000"
) |>
  dplyr::mutate(Mean_Var = pmax(MSE - Mean_Bias^2, 0))

custom_color <- MetBrewer::met.brewer("Troy", n = 8, type = "discrete")
custom_colors <- c(custom_color[1], custom_color[3], custom_color[8])
custom_pch <- c(8, 17, 5)

make_plot <- function(effect) {
  summary_plot_fun(
    dplyr::filter(full_data, effect_type == effect),
    coverage_adjusted = "Cov_adj",
    coverage_standard = "Cov_stand",
    coverage_adjusted_label = "m-BS \\\\ Cov.",
    coverage_standard_label = "Wald \\\\ Cov.",
    method_colors = custom_colors,
    method_shapes = custom_pch
  )
}

p_nde <- make_plot("NDE")
p_nie <- make_plot("NIE")

ggplot2::ggsave(file.path(output_dir, "simulation_onestep_NDE.pdf"), p_nde,
                width = 12, height = 10)
ggplot2::ggsave(file.path(output_dir, "simulation_onestep_NIE.pdf"), p_nie,
                width = 12, height = 10)

saveRDS(full_data, here::here("simulations", "results", "precomputed", "onestep_summary_combined.rds"))
