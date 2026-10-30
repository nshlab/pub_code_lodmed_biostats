# Export point estimates only for Table 2.
input <- "results/tables/table2_v9_comparison.csv"
if (!file.exists(input)) stop("Run analysis/03_table2_point_estimates.R first")
d <- read.csv(input, stringsAsFactors=FALSE)
if (nrow(d) != 18L || anyNA(d$estimate) || any(nzchar(d$error))) stop("Incomplete point estimates; inspect comparison and error fields")
d$percent_mediated <- NA_real_
for (est in unique(d$estimator)) for (v in unique(d$variant)) {
  ii <- which(d$estimator == est & d$variant == v)
  nde <- d$estimate[ii[d$effect[ii] == "NDE"]]
  nie <- d$estimate[ii[d$effect[ii] == "NIE"]]
  if (length(nde) != 1L || length(nie) != 1L) stop("Missing effects")
  # Descriptive ratio; not asserted to match manuscript PM definition.
  if (is.finite(nde + nie) && abs(nde + nie) > .Machine$double.eps)
    d$percent_mediated[ii] <- 100 * nie/(nde + nie)
}
d$se <- NA_real_; d$ci_lower <- NA_real_; d$ci_upper <- NA_real_
d$validation_note <- "Point estimates only; SE/CI not validated; PM exploratory"
cols <- c("estimator","variant","effect","estimate","manuscript_estimate","absolute_difference","percent_mediated","se","ci_lower","ci_upper","validation_note")
write.csv(d[,cols], "results/tables/table2_v10_point_estimates_ONLY.csv", row.names=FALSE, na="")
message("Wrote point-estimate-only output; bootstrap not run. PM descriptive and not validated.")
