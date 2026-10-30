# Reproduce the g-computation simulation summary plots from precomputed results.

required <- c("here", "dplyr", "purrr", "readr", "stringr", "tibble", "tidyr",
              "ggplot2", "latex2exp", "MetBrewer")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

source(here::here("simulations", "plots", "plot_helpers.R"))

input_dir <- here::here("simulations", "results", "precomputed", "gcomp")
output_dir <- here::here("simulations", "results", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

method_labels <- c(
  Traditional = "LLoQ/2 imputation",
  `EM-HAL-tetero` = "EM (4): heterosced. CDE, HAL",
  `EM-HAL` = "EM (5): homosced. CDE, HAL",
  `EM-correctGLM` = "EM (1): true nuisances",
  `EM-misspecGLM` = "EM (2): oracle dens., misspec. GLMs",
  `EM-misspecDens` = "EM (3): misspec. dens., misspec. GLMs"
)

full_data <- load_precomputed_results(input_dir, method_labels)

custom_color <- MetBrewer::met.brewer("Troy", n = 8, type = "discrete")
custom_colors <- c(custom_color[1], custom_color[3], custom_color[5],
                   custom_color[3], custom_color[5], custom_color[8])
custom_pch <- c(8, 15, 15, 17, 17, 5)

make_plot <- function(effect) {
  summary_plot_fun(
    dplyr::filter(full_data, effect_type == effect),
    coverage_adjusted = "Cov_moutn",
    coverage_standard = "Cov_stand",
    coverage_adjusted_label = "m-BS \\\\ Cov.",
    coverage_standard_label = "Std BS \\\\ Cov.",
    method_colors = custom_colors,
    method_shapes = custom_pch
  )
}

p_nde <- make_plot("NDE")
p_nie <- make_plot("NIE")

ggplot2::ggsave(file.path(output_dir, "simulation_gcomp_NDE.pdf"), p_nde,
                width = 12, height = 10)
ggplot2::ggsave(file.path(output_dir, "simulation_gcomp_NIE.pdf"), p_nie,
                width = 12, height = 10)

saveRDS(full_data, here::here("simulations", "results", "precomputed", "gcomp_summary_combined.rds"))
