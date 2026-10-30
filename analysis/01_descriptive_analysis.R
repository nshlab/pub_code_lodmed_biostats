# Descriptive analysis of the application longitudinal data.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required")
source(here::here("analysis", "helpers", "frozen_input.R"), local = TRUE)
d <- read_application_longitudinal()
required <- c("age", "sex", "bmi_cat", "race", "ethnicity", "country",
              "vaccine", "symptom_duration", "DIAB", "CVD", "HPTN", "CURSMOKR", "Y")
if (!all(required %in% names(d))) stop("Unavailable Table S1 fields: ", paste(setdiff(required, names(d)), collapse = ", "))
groups <- list(Overall = rep(TRUE, nrow(d)), Placebo = d$A == 0, `A/R` = d$A == 1)
rows <- list(); limitations <- character()
fmt_binary <- function(x) {
  if (all(is.na(x))) return("NA")
  sprintf("%d (%.1f%%)", sum(x == 1, na.rm = TRUE), 100 * sum(x == 1, na.rm = TRUE)/sum(!is.na(x)))
}
fmt_mean <- function(x) if (all(is.na(x))) "NA" else sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
add_row <- function(label, fun) {
  vals <- vapply(groups, function(ix) fun(d[ix, , drop = FALSE]), character(1))
  rows[[length(rows) + 1L]] <<- data.frame(Characteristics = label,
    Overall = vals[1], Placebo = vals[2], `A/R` = vals[3], check.names = FALSE)
}
check_levels <- function(name, allowed) {
  unexpected <- setdiff(unique(d[[name]][!is.na(d[[name]])]), allowed)
  if (length(unexpected)) stop("Unexpected ", name, " codes: ", paste(unexpected, collapse = ", "))
}
for (name in c("vaccine", "DIAB", "CVD", "HPTN", "CURSMOKR")) check_levels(name, c("Yes", "No"))
check_levels("sex", c("F", "M"))
check_levels("bmi_cat", c("Normal", "Overweight", "Obese", "Underweight"))
check_levels("race", c("WHITE", "BLACK OR AFRICAN AMERICAN", "OTHER"))
check_levels("ethnicity", c("HISPANIC OR LATINO", "NOT HISPANIC OR LATINO"))
check_levels("country", c("USA", "Non-USA"))
if (!all(is.na(d$Y) | d$Y %in% c(0, 1))) stop("Unexpected outcome coding")
age <- suppressWarnings(as.numeric(d$age))
if (any(!is.na(d$age) & is.na(age) & d$age != ">80"))
  stop("Unrecognized censored age code")
# Local in-memory derivation only; source file is never written.
d$.age_ge60 <- ifelse(!is.na(d$age) & d$age == ">80", TRUE,
                     ifelse(is.na(age), NA, age >= 60))
add_row("N", function(x) as.character(nrow(x)))
add_row("Age >= 60 years", function(x) fmt_binary(x$.age_ge60))
add_row("Female", function(x) fmt_binary(x$sex == "F"))
for (value in c("Normal", "Overweight", "Obese")) local({
  level <- value
  add_row(paste0("BMI: ", level), function(x) fmt_binary(x$bmi_cat == level))
})
# Preserve the legacy combined row explicitly; raw missingness is separately audited.
add_row("BMI: Underweight / unclassified (NOT Table S1 underweight)", function(x) {
  value <- is.na(x$bmi_cat) | x$bmi_cat == "Underweight"
  sprintf("%d (%.1f%%)", sum(value), 100 * mean(value))
})
add_row("Race: White", function(x) fmt_binary(x$race == "WHITE"))
add_row("Race: Non-White (Black + Other combined)", function(x) fmt_binary(x$race != "WHITE"))
add_row("Hispanic or Latino", function(x) fmt_binary(x$ethnicity == "HISPANIC OR LATINO"))
add_row("USA", function(x) fmt_binary(x$country == "USA"))
add_row("Vaccinated", function(x) fmt_binary(x$vaccine == "Yes"))
add_row("Days from symptom onset to enrollment", function(x) fmt_mean(x$symptom_duration))
for (spec in list(c("Diabetes", "DIAB"), c("Cardiovascular disease", "CVD"),
                  c("Hypertension", "HPTN"), c("Current smoker", "CURSMOKR"))) local({
  label <- spec[1]; field <- spec[2]
  add_row(label, function(x) fmt_binary(x[[field]] == "Yes"))
})
# Baseline visit is original assigned Day 0 (no Day 3 window applied here).
# day0_value == legacy Day-0 AVAL.y; values are already on the source log10 scale.
# Summarize the supplied numeric values as recorded, including source below-LLoQ
# codes. Do not recreate historical substitution-to-1 or BASE fallback.
baseline_required <- "day0_value"
if (!all(baseline_required %in% names(d))) stop("Missing raw Day 0 fields")
baseline_raw <- as.character(d$day0_value)
d$.baseline_rna <- suppressWarnings(as.numeric(baseline_raw))
source(here::here("analysis", "helpers", "rna_classification.R"), local = TRUE)
status <- rna_status_at(d, 0L)
d$.baseline_status <- ifelse(status == "Observed", "Quantifiable",
  ifelse(status == "<LLoQ", "Below LLoQ", status))
# Validate the historical field mapping without changing either source field.
if ("AVAL.y" %in% names(d) && !isTRUE(all.equal(d$.baseline_rna,
    suppressWarnings(as.numeric(d$AVAL.y)), tolerance = 0, check.attributes = FALSE)))
  stop("Day 0 source mapping differs from legacy AVAL.y")
add_row("AN RNA day 0 (log10 copies/mL)", function(x) fmt_mean(x$.baseline_rna))
limitations <- c(limitations,
  "Baseline mean (SD) uses numeric day0_value exactly as stored, on the existing log10 scale, including source-coded below-LLoQ values 0.7/1.7. Those codes are not exact latent concentrations.",
  "Baseline variable mapping: original NSLRNA assigned Day 0 AVAL -> frozen day0_value (identical to AVAL.y); status derived at 2.0 log10 copies/mL; optional text/CAT1 availability is respected. No baseline visit window is introduced.",
  "Historical M_base used coalesced baseline/Day-3 LLoQ flags to substitute 1, then AVAL.y/BASE fallback. Those substitutions are intentionally NOT reproduced; this is the requested unreplaced source-value description, not the old substituted M_base summary.",
  "Missing numeric baseline values are excluded from mean/SD but categorical below-LLoQ results remain separately counted; no imputation, recoding of source values, re-log transformation or BASE fallback is performed.")
add_row("Hospitalization or death after Day 3", function(x) fmt_binary(x$Y))
limitations <- c(limitations,
  "Age >=60 is derived from age; the unavailable legacy age_get60 column cannot be checked directly.",
  "The legacy combined Underweight/unclassified row includes true underweight and missing BMI; inspect the separate category QC.",
  "The original dictionary lists several derived indicators not present in this 72-column CSV. Original dictionary bytes are retained, not rewritten.",
  "This is a descriptive draft, not a claim of exact agreement with every manuscript Table S1 cell.")
out <- here::here("results", "tables"); qcdir <- here::here("results", "qc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(qcdir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(do.call(rbind, rows), file.path(out, "tableS1.csv"), row.names = FALSE, na = "")
qc <- data.frame(variable = required, missing = vapply(d[required], function(x) sum(is.na(x)), integer(1)))
utils::write.csv(qc, file.path(qcdir, "tableS1_missingness.csv"), row.names = FALSE)
category_qc <- do.call(rbind, lapply(c("A", "sex", "bmi_cat", "race", "ethnicity", "country", "vaccine"), function(v) {
  x <- table(d[[v]], useNA = "always")
  data.frame(variable = v, category = names(x), count = as.integer(x))
}))
utils::write.csv(category_qc, file.path(qcdir, "tableS1_category_counts.csv"), row.names = FALSE)
writeLines(limitations, file.path(qcdir, "tableS1_limitations.txt"))
# Aggregate baseline availability, numeric statistics and assay status by group.
baseline_qc <- do.call(rbind, lapply(names(groups), function(group) {
  x <- d[groups[[group]], , drop = FALSE]
  n <- sum(!is.na(x$.baseline_rna))
  quant <- x$.baseline_rna[x$.baseline_status == "Quantifiable"]
  data.frame(group = group, n_participants = nrow(x),
    n_available_assay = sum(x$.baseline_status %in% c("Quantifiable", "Below LLoQ")),
    n_numeric = n, n_quantifiable = sum(x$.baseline_status == "Quantifiable"),
    n_below_lloq = sum(x$.baseline_status == "Below LLoQ"),
    n_below_lloq_without_numeric = sum(x$.baseline_status == "Below LLoQ" & is.na(x$.baseline_rna)),
    n_missing = sum(x$.baseline_status == "Missing"),
    n_unclassified = sum(x$.baseline_status == "Unclassified"),
    mean_source_numeric = if (n) mean(x$.baseline_rna, na.rm = TRUE) else NA_real_,
    sd_source_numeric = if (n > 1L) stats::sd(x$.baseline_rna, na.rm = TRUE) else NA_real_,
    mean_quantifiable_only = if (length(quant)) mean(quant) else NA_real_,
    sd_quantifiable_only = if (length(quant) > 1L) stats::sd(quant) else NA_real_)
}))
stopifnot(all(baseline_qc$n_quantifiable + baseline_qc$n_below_lloq +
  baseline_qc$n_missing + baseline_qc$n_unclassified == baseline_qc$n_participants))
utils::write.csv(baseline_qc, file.path(qcdir, "tableS1_baseline_rna.csv"), row.names = FALSE)
writeLines(c(
  "Baseline AN RNA: mean (SD) of available source numeric Day 0 values, already on the source log10 scale.",
  "Includes the unchanged numeric codes attached to below-LLoQ measurements; these are not exact uncensored concentrations.",
  "No replacement-to-1, BASE fallback, missing-value imputation or new log transformation was applied.",
  "See results/qc/tableS1_baseline_rna.csv for available/numeric/quantifiable/below-LLoQ/missing counts by treatment group."),
  file.path(out, "tableS1_notes.txt"))
message("Saved descriptive Table S1 with baseline source numeric RNA and group-specific availability QC.")
