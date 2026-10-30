if (!requireNamespace("here", quietly = TRUE)) stop("here is required")
source(here::here("analysis", "helpers", "rna_classification.R"))
x <- data.frame(day0_value = c(2.01, 2, 1.7, NA, NA, NA, NA, NA),
 day0_text = c(NA, NA, NA, "< LLOQ", "Detected", NA, ">= LLOQ", "NOT DONE"),
 day0_cat1 = c(NA, NA, NA, NA, NA, "NOT DETECTED", NA, NA))
stopifnot(identical(rna_status_at(x, 0), c("Observed", "<LLoQ", "<LLoQ",
 "<LLoQ", "Unclassified", "<LLoQ", "Unclassified", "Missing")))
source(here::here("analysis", "helpers", "day3_window.R"))
x <- data.frame(day2_value = c(5, NA), day3_value = c(NA, NA),
 day4_value = c(1.7, NA), day3_text = c("< LLOQ", NA),
 day2_cat1 = c(NA, "NOT DETECTED"))
a <- select_day3_window(x)
stopifnot(identical(a$selected_day, c(3L, 2L)), all(is.na(a$selected_value)),
 all(a$status == "<LLoQ"), identical(a$eligible_count, c(3L, 1L)))
message("PASS: threshold boundary, below-LLoQ numeric NA, unavailable and unknown-status assays; paired visit selection")
