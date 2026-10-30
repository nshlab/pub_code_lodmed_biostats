# Shared helpers for plotting precomputed simulation summaries.

extract_result_metadata <- function(filename) {
  parts <- stringr::str_match(
    basename(filename),
    "(.*)_n_(\\d+)rate(0\\.\\d+)\\.csv"
  )
  if (anyNA(parts[1, 2:4])) {
    stop("Could not parse simulation result filename: ", basename(filename))
  }
  tibble::tibble(
    method_name = parts[, 2],
    sample_size = as.integer(parts[, 3]),
    censoring_rate = as.numeric(parts[, 4]),
    file = filename
  )
}

read_two_effect_summary <- function(file, method_name, sample_size, censoring_rate) {
  dat <- readr::read_csv(file, show_col_types = FALSE)
  if (nrow(dat) != 2L) {
    stop("Expected exactly two rows (NDE, NIE) in: ", basename(file))
  }
  dplyr::mutate(
    dat,
    effect_type = c("NDE", "NIE"),
    method_name = method_name,
    sample_size = sample_size,
    censoring_rate = censoring_rate
  )
}

load_precomputed_results <- function(folder_path, method_labels, exclude_pattern = NULL) {
  files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  if (!is.null(exclude_pattern)) {
    files <- files[!grepl(exclude_pattern, basename(files))]
  }
  if (!length(files)) stop("No CSV result files found in: ", folder_path)

  metadata <- purrr::map_dfr(files, extract_result_metadata) |>
    dplyr::mutate(
      method_name = dplyr::recode(method_name, !!!method_labels)
    )

  purrr::pmap_dfr(
    metadata,
    ~ read_two_effect_summary(..4, ..1, ..2, ..3)
  )
}

summary_plot_fun <- function(full_data, coverage_adjusted, coverage_standard,
                             coverage_adjusted_label, coverage_standard_label,
                             method_colors, method_shapes, type = "b") {
  measures <- c("Mean_Bias", "Mean_Var", "MSE", coverage_standard, coverage_adjusted)

  comb_tbl <- full_data |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(measures),
      names_to = "measure",
      values_to = "value"
    ) |>
    dplyr::mutate(
      measure = factor(measure, levels = measures),
      sample_size = factor(sample_size),
      censoring_rate = factor(censoring_rate),
      y_intercept = dplyr::if_else(
        measure %in% c(coverage_standard, coverage_adjusted), 0.95, 0
      )
    )

  levels(comb_tbl$measure) <- c(
    latex2exp::TeX("Bias"),
    latex2exp::TeX("Variance"),
    latex2exp::TeX("MSE"),
    latex2exp::TeX(coverage_standard_label),
    latex2exp::TeX(coverage_adjusted_label)
  )

  ggplot2::ggplot(
    comb_tbl,
    ggplot2::aes(
      x = sample_size, y = value,
      colour = method_name, shape = method_name, group = method_name
    )
  ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(measure), cols = ggplot2::vars(censoring_rate),
      scales = "free_y", labeller = ggplot2::label_parsed
    ) +
    ggplot2::geom_jitter(
      alpha = 1, size = ifelse(type == "o", 1.5, 2),
      width = 0.1, height = 0
    ) +
    ggplot2::geom_line(alpha = 0.25) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = y_intercept), linetype = 2, alpha = 0.95
    ) +
    ggplot2::labs(
      x = "Sample Size", y = "Value (300 Replicates)",
      colour = "Method", shape = "Method"
    ) +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::scale_shape_manual(values = method_shapes) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}
