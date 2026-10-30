if (!requireNamespace("here", quietly = TRUE)) stop("here is required")
source(here::here("analysis", "helpers", "day3_window.R"))
# Distinguish Day 3 priority, lower-value backup, source-paired flags, ties and missing.
control <- data.frame(day2_value = c(1, 5, 2, NA), day3_value = c(9, NA, NA, NA),
  day4_value = c(0.5, 1, 2, NA), day2_cat2 = c("< LLOQ", ">= LLOQ", ">= LLOQ", NA),
  day3_cat2 = c(">= LLOQ", NA, NA, NA), day4_cat2 = c("< LLOQ", "< LLOQ", ">= LLOQ", NA))
b <- select_day3_window(control)
stopifnot(identical(b$selected_day, c(3L, 4L, 2L, NA_integer_)),
  identical(b$status, c("Observed", "<LLoQ", "<LLoQ", "Missing")))
message("Synthetic Day 3 selection controls passed.")
