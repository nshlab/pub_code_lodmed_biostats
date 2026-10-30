select_day3_window <- function(d) {
  source(here::here("analysis", "helpers", "rna_classification.R"), local = TRUE)
  required <- unlist(lapply(c(2L, 3L, 4L), function(day) paste0("day", day, "_value")))
  if (!all(required %in% names(d))) stop("Missing assigned-visit window fields")
  statuses <- lapply(c(2L, 3L, 4L), function(day) rna_status_at(d, day))
  available2 <- statuses[[1]] != "Missing"
  available3 <- statuses[[2]] != "Missing"
  available4 <- statuses[[3]] != "Missing"
  # Mirror original AVAL_filled: Day 3 first; otherwise min(Day 2, Day 4).
  # When equal, Day 2 represents the same numeric value and receives the tie.
  selected_day <- rep(NA_integer_, nrow(d))
  selected_day[available3] <- 3L
  use2 <- !available3 & available2 & (!available4 | d$day2_value <= d$day4_value)
  use4 <- !available3 & available4 & (!available2 | d$day4_value < d$day2_value)
  selected_day[which(use2)] <- 2L
  selected_day[which(use4)] <- 4L
  # A single categorical-only backup is usable without inventing a value.
  selected_day[which(!available3 & available2 & !available4)] <- 2L
  selected_day[which(!available3 & available4 & !available2)] <- 4L
  if (any(!available3 & available2 & available4 &
          (is.na(d$day2_value) | is.na(d$day4_value))))
    stop("Cannot rank two backup measurements with missing numeric RNA; original selection needs review")
  selected_value <- rep(NA_real_, nrow(d))
  status <- rep("Missing", nrow(d))
  selected_flag <- rep(NA_character_, nrow(d))
  for (day in c(2L, 3L, 4L)) {
    ix <- which(selected_day == day)
    selected_value[ix] <- d[[paste0("day", day, "_value")]][ix]
    status[ix] <- statuses[[match(day, c(2L, 3L, 4L))]][ix]
    selected_flag[ix] <- ifelse(status[ix] == "<LLoQ", "< LLOQ",
      ifelse(status[ix] == "Observed", ">= LLOQ", NA_character_))
  }
  # Audit the historical independent Day3->Day2->Day4 flag fallback using
  # derived categories, without requiring any stored CAT2 column.
  flags <- lapply(statuses, function(x) ifelse(x == "<LLoQ", "< LLOQ",
    ifelse(x == "Observed", ">= LLOQ", NA_character_)))
  legacy_flag <- flags[[2]]
  ix <- which(is.na(legacy_flag)); legacy_flag[ix] <- flags[[1]][ix]
  ix <- which(is.na(legacy_flag)); legacy_flag[ix] <- flags[[3]][ix]
  mismatch <- !is.na(selected_flag) & !is.na(legacy_flag) & selected_flag != legacy_flag
  data.frame(selected_day = selected_day, selected_value = selected_value,
    selected_flag = selected_flag, status = status,
    available_day2 = available2, available_day3 = available3, available_day4 = available4,
    eligible_count = as.integer(available2) + as.integer(available3) + as.integer(available4),
    legacy_flag_disagreement = mismatch,
    eligibility_basis = "Original assigned-visit labels; actual collection dates unavailable",
    stringsAsFactors = FALSE)
}
