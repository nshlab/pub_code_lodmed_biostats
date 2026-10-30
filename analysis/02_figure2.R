# Frozen longitudinal CSV -> descriptive RNA plots; no imputation/estimation.
if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
source(here::here("analysis", "helpers", "frozen_input.R"), local = TRUE)
d <- read_application_longitudinal()
source(here::here("analysis", "helpers", "day3_window.R"), local = TRUE)
assessment <- select_day3_window(d)
days <- c(0, 2, 3, 4, 7, 14, 28)
required <- unlist(lapply(days, function(day) paste0("day", day, "_value")))
if (!all(required %in% names(d))) stop("Required raw visit fields unavailable")
out <- here::here("results", "figures")
qc <- here::here("results", "qc")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(qc, recursive = TRUE, showWarnings = FALSE)
# Numeric >2.0: observed; <=2.0: at/below LLoQ. Numeric NA is resolved
# using available text/CAT1 evidence; unknown assay status is never silently missing.
source(here::here("analysis", "helpers", "rna_classification.R"), local = TRUE)
long <- do.call(rbind, lapply(days, function(day) data.frame(
  id = d$id, A = d$A, day = day,
  value = d[[paste0("day", day, "_value")]],
  status = rna_status_at(d, day))))
# Preserve raw visit summaries separately. Only the Day 3 assessment uses the window.
raw_long <- long
ix3 <- which(long$day == 3)
long$value[ix3] <- assessment$selected_value
long$status[ix3] <- assessment$status
levels_status <- c("Observed", "<LLoQ", "Missing", "Unclassified")
# Determine missingness after assigned-visit window selection, then check the
# frozen historical cohort. M never supplies measurement values.
missing_ids <- d$id[assessment$eligible_count == 0L]
if (!setequal(missing_ids, d$id[is.na(d$M)]))
  stop("Window-missing cohort differs from frozen historical M-missing cohort")
if (length(missing_ids) != 38L) stop("Historical Figure 2(b) cohort differs from 38")
pattern <- long[long$id %in% missing_ids, , drop = FALSE]
make_counts <- function(dat, cohort, denominator) {
  counts <- as.data.frame(table(day = factor(dat$day, levels = days),
    status = factor(dat$status, levels = levels_status)), stringsAsFactors = FALSE)
  names(counts)[3L] <- "count"
  counts$cohort <- cohort
  counts$denominator <- denominator
  counts$percent <- 100 * counts$count/denominator
  stopifnot(all(tapply(counts$count, counts$day, sum) == denominator))
  counts
}
counts <- rbind(make_counts(long, "All 749 participants", 749L),
                make_counts(pattern, "Historical 38 M-missing participants", 38L))
utils::write.csv(counts, file.path(qc, "figure2_measurement_counts.csv"), row.names = FALSE)
raw_counts <- rbind(make_counts(raw_long, "All 749 participants", 749L),
  make_counts(raw_long[raw_long$id %in% missing_ids, , drop = FALSE], "Historical 38 M-missing participants", 38L))
utils::write.csv(raw_counts, file.path(qc, "figure2_raw_visit_counts.csv"), row.names = FALSE)
# Observed Day 3 ASSESSMENT: source measurements from assigned Days 2, 3 or 4.
# The frozen M column is not used to supply any plotted RNA value.
observed <- assessment$selected_value[assessment$status == "Observed"]
if (!length(observed)) stop("No above-threshold observed Day 3 assessments")
p_a <- ggplot2::ggplot(data.frame(rna = observed), ggplot2::aes(x = rna)) +
  ggplot2::geom_histogram(binwidth = 0.5, fill = "#7fcdbb", colour = "white", alpha = 0.9) +
  ggplot2::geom_density(ggplot2::aes(y = ggplot2::after_stat(count) * 0.5),
    colour = "#253494", linewidth = 0.7) +
  ggplot2::labs(title = "Distribution of observed Day 3 AN RNA assessments (> 2.0 log10 copies/mL)",
    subtitle = paste0("Assigned Day 3 +/-1-day window; n = ", length(observed)),
    x = "RNA (source log10 copies/mL)", y = "Count") + ggplot2::theme_minimal()
p_b <- ggplot2::ggplot(pattern, ggplot2::aes(x = day, y = factor(id))) +
  ggplot2::geom_point(ggplot2::aes(colour = status), size = 2) +
  ggplot2::scale_colour_manual(values = c(Observed = "#0072B2", `<LLoQ` = "#E69F00",
    Missing = "#999999", Unclassified = "#CC79A7"),
    labels = c(Observed = "Observed (>2.0)", `<LLoQ` = "At/below LLoQ (<=2.0)",
      Missing = "Missing", Unclassified = "Unclassified")) +
  ggplot2::scale_x_continuous(breaks = days, limits = c(0, 28)) +
  ggplot2::labs(x = "Study Day", y = "Participant", colour = "Measurement Status",
    title = "Measurement Pattern of Viral RNA Among 38 Missing Participants") +
  ggplot2::theme_bw() + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
for (ext in c("pdf", "png")) {
  ggplot2::ggsave(file.path(out, paste0("figure2a.", ext)), p_a, width = 7, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(out, paste0("figure2b.", ext)), p_b, width = 8, height = 7, dpi = 300)
}
historical_numeric <- suppressWarnings(as.numeric(d$M))
subset73 <- !is.na(historical_numeric) & is.na(d$day3_value)
raw_eligibility <- function(ix, label) data.frame(cohort = label,
  metric = c("assigned_Day2_available", "assigned_Day3_available", "assigned_Day4_available",
    "no_available_measurement_in_window", "multiple_available_measurements",
    "historical_numeric_M_without_matching_selected_measurement",
    "historical_numeric_M_with_matching_selected_measurement"),
  count = c(sum(assessment$available_day2[ix]), sum(assessment$available_day3[ix]),
    sum(assessment$available_day4[ix]), sum(assessment$eligible_count[ix] == 0L),
    sum(assessment$eligible_count[ix] > 1L),
    sum(!is.na(historical_numeric[ix]) & (is.na(assessment$selected_value[ix]) |
      historical_numeric[ix] != assessment$selected_value[ix]), na.rm = TRUE),
    sum(!is.na(historical_numeric[ix]) & historical_numeric[ix] == assessment$selected_value[ix], na.rm = TRUE)))
raw_audit <- rbind(raw_eligibility(rep(TRUE, nrow(d)), "All 749 participants"),
                   raw_eligibility(subset73, "73 previously questioned participants"))
utils::write.csv(raw_audit, file.path(qc, "figure2_day3_window_eligibility.csv"), row.names = FALSE)
final_counts <- function(ix, label) {
  day <- ifelse(is.na(assessment$selected_day[ix]), "No measurement", paste0("Day ", assessment$selected_day[ix]))
  tab <- table(factor(day, levels = c("Day 2", "Day 3", "Day 4", "No measurement")))
  data.frame(cohort = label, classification = names(tab), count = as.integer(tab), denominator = sum(ix))
}
utils::write.csv(rbind(final_counts(rep(TRUE, nrow(d)), "All 749 participants"),
  final_counts(subset73, "73 previously questioned participants")),
  file.path(qc, "figure2_day3_window_final.csv"), row.names = FALSE)
legacy_status <- ifelse(is.na(d$M), "Missing", ifelse(d$M == "<= LLOQ", "<LLoQ", "Observed"))
metrics <- data.frame(metric = c("historical_numeric_M", "historical_numeric_M_without_nominal_Day3",
  "observed_window_assessments", "legacy_M_status_disagreements", "legacy_independent_flag_disagreements",
  "historical_M_missing", "CSV_visit_count", "actual_collection_dates_available"),
  value = c(sum(!is.na(historical_numeric)), sum(subset73), length(observed),
    sum(assessment$status != legacy_status), sum(assessment$legacy_flag_disagreement),
    length(missing_ids), length(days), 0L))
utils::write.csv(metrics, file.path(qc, "figure2_definition_audit.csv"), row.names = FALSE)
writeLines(c(
  "RNA status is derived in memory: >2.0 observed; <=2.0 at/below LLoQ. Optional text/CAT1 resolves numeric-NA assay availability; no stored CAT2 is read.",
  "Figure 2 is descriptive only; no imputation or statistical correction is performed.",
  "Day 3 assessment uses actual recorded RNA from original assigned Days 2/3/4: prefer Day 3, otherwise the lower Day 2/4 value (tie Day 2).",
  "The original implementation used AVISIT labels. Actual calendar collection dates are absent from the CSV; calendar-date eligibility cannot be independently verified.",
  "The 73 questioned numeric M values all match selected Day 2/4 source measurements. No model-filled values are needed to account for them; model-provenance flags are not stored in this CSV.",
  "Quantification status is paired with the selected source measurement; the original independent flag fallback differs in six participants, but final observed/below-LLoQ/missing categories agree with frozen M for all 749.",
  "Only Day 3 uses the grace window. Other visits retain their own source values; raw-visit counts are saved separately.",
  "Panel b preserves the historical 38 M-missing participants and the seven available visits. Additional original visits are unavailable and are not invented."),
  file.path(qc, "figure2_limitations.txt"))
message("Day 3 window applied using original assigned-visit logic. See results/qc/ for eligibility, final classifications and limitations.")
