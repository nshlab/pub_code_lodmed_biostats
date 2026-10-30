if (!requireNamespace("here", quietly = TRUE)) stop("here is required")
source(here::here("analysis", "helpers", "frozen_input.R"))
source(here::here("analysis", "helpers", "day3_window.R"))
d <- read_frozen_longitudinal()
a <- select_day3_window(d)
m <- suppressWarnings(as.numeric(d$M))
ix <- !is.na(m) & is.na(d$day3_value)
stopifnot(sum(ix) == 73L,
  sum(a$available_day2[ix]) == 41L, sum(a$available_day3[ix]) == 0L,
  sum(a$available_day4[ix]) == 41L, sum(a$eligible_count[ix] > 1L) == 9L,
  sum(a$eligible_count[ix] == 0L) == 0L,
  sum(a$selected_day[ix] == 2L) == 34L, sum(a$selected_day[ix] == 4L) == 39L,
  all(a$selected_value[ix] == m[ix]), sum(a$status == "Observed") == 355L,
  sum(a$status == "<LLoQ") == 356L, sum(a$status == "Missing") == 38L,
  sum(a$legacy_flag_disagreement) == 6L)
message("Frozen cohort Day 3 regression passed.")
