# Shared input validation only. Does not change or save participant-level records.
read_frozen_longitudinal <- function() {
  path <- here::here("data", "private", "activ2_749_TRUE_VISITS_INTERNAL.csv")
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  dat <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("", "NA"))
  required <- c("id", "A", "M", "Y")
  if (!all(required %in% names(dat))) stop("Missing core columns")
  # The supplied frozen file is wide, one participant per row. Do not deduplicate.
  if (nrow(dat) != 749L || anyNA(dat$id) || anyDuplicated(dat$id))
    stop("Expected 749 unique nonmissing participant IDs; source is not modified")
  if (anyNA(dat$A) || !all(dat$A %in% c(0, 1)) ||
      sum(dat$A == 0) != 371L || sum(dat$A == 1) != 378L)
    stop("Frozen cohort treatment counts differ from 371 placebo / 378 A/R")
  dat
}

read_application_longitudinal <- function() {
  mode <- Sys.getenv("LODMED_DATASET", unset = "candidate")
  if (!mode %in% c("candidate", "internal")) stop("LODMED_DATASET must be candidate or internal")
  if (mode == "internal") {
    dat <- read_frozen_longitudinal()
    dat$age <- dat$age_cens80
  } else {
    path <- here::here("data", "release_candidate", "activ2_749_public_candidate.csv")
    if (!file.exists(path)) stop("Missing private release candidate: ", path)
    dat <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                          na.strings = c("", "NA"))
    required <- c("id", "A", "Y", "age", "sex", "bmi_cat", "symptom_duration",
      "race", "ethnicity", "country", "vaccine", "DIAB", "CVD", "HPTN", "CURSMOKR",
      unlist(lapply(c(0, 2, 3, 4, 7, 14, 28), function(day) paste0("day", day, "_value"))))
    if (!setequal(names(dat), required)) stop("Candidate must contain exactly the audited 22 columns")
    if (nrow(dat) != 749L || anyNA(dat$id) || anyDuplicated(dat$id)) stop("Candidate cohort differs from 749 unique participants")
    if (anyNA(dat$A) || !all(dat$A %in% c(0, 1)) ||
        sum(dat$A == 0) != 371L || sum(dat$A == 1) != 378L) stop("Candidate treatment counts differ")
    source(here::here("analysis", "helpers", "day3_window.R"), local = TRUE)
    assessment <- select_day3_window(dat)
    if (any(assessment$status == "Unclassified")) stop("Cannot recover legacy M for unclassified source status")
    # M is excluded as redundant from the candidate. Exact reconstruction was
    # verified independently against all original 749 M strings before release.
    # This is deterministic source selection/categorical encoding, not imputation.
    dat$M <- as.character(assessment$selected_value)
    dat$M[assessment$status == "<LLoQ"] <- "<= LLOQ"
    dat$M[assessment$status == "Missing"] <- NA_character_
  }
  attr(dat, "input_mode") <- mode
  dat
}
