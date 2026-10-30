
for (pkg in c("here", "ggplot2", "survival")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package: ", pkg)
}
source(here::here("analysis", "helpers", "frozen_input.R"), local = TRUE)
source(here::here("analysis", "helpers", "day3_window.R"), local = TRUE)
source(here::here("analysis", "helpers", "rna_classification.R"), local = TRUE)

d <- read_application_longitudinal()
a <- select_day3_window(d)
lloq <- 2.0
base_status <- rna_status_at(d, 0)
keep <- base_status == "Observed" & a$status %in% c("Observed", "<LLoQ") &
  !is.na(d$day0_value) & !is.na(d$A)
keep[is.na(keep)] <- FALSE
z <- data.frame(
  id = d$id[keep],
  A = factor(d$A[keep]),
  baseline = as.numeric(d$day0_value[keep]),
  status = a$status[keep],
  day3 = as.numeric(a$selected_value[keep])
)
if (nrow(z) < 20 || length(unique(z$A)) < 2L) stop("Too few records / treatment groups")
z$observed <- z$status == "Observed"
if (sum(z$observed) < 10 || sum(!z$observed) < 2) stop("Insufficient observed/censored data")
if (any(!is.finite(z$day3[z$observed])) ||
    any(z$day3[z$observed] <= lloq)) stop("Unexpected quantifiable Day 3 values")
# For left-censored outcomes, the observed threshold is 2.0, NOT the latent RNA.
z$y <- ifelse(z$observed, z$day3, lloq)
fit <- survival::survreg(
  survival::Surv(z$y, z$observed, type = "left") ~ A + baseline,
  data = z, dist = "gaussian"
)
z$mu <- as.numeric(predict(fit, newdata = z, type = "response"))
sigma <- fit$scale
z$threshold_z <- (lloq - z$mu) / sigma
z$cdf_threshold <- pnorm(z$threshold_z)
z$resid_z <- ifelse(z$observed, (z$y - z$mu) / sigma, NA_real_)

# Diagnostic 1: conditional randomized probability-integral-transform (PIT).
# Observed: U=Phi((Y-mu)/sigma). Censored: U~Uniform(0,Phi((LLoQ-mu)/sigma)).
# Under the fitted model U should be Uniform(0,1), hence qnorm(U) ~ N(0,1).
# Randomization is necessary because censored values are intervals, not points.
set.seed(20261020)
u <- ifelse(z$observed, pnorm(z$resid_z),
            runif(nrow(z)) * z$cdf_threshold)
u <- pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps)
qq <- data.frame(theoretical = qnorm(ppoints(length(u))),
                 residual = sort(qnorm(u)))
p1 <- ggplot2::ggplot(qq, ggplot2::aes(theoretical, residual)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, colour = "steelblue", linewidth = 0.8) +
  ggplot2::geom_point(size = 1.5, alpha = 0.75) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Day 3: censored-normal conditional residual QQ",
    subtitle = sprintf("Randomized PIT; n=%d (observed=%d, censored=%d); baseline quantifiable only",
                       nrow(z), sum(z$observed), sum(!z$observed)),
    x = "Standard normal quantile", y = "Randomized standardized residual"
  ) + ggplot2::theme_bw()

# Diagnostic 2: observed values only, correctly conditional on exceeding LLoQ.
# For each observed Y, the truncated-normal PIT is Uniform(0,1) under the model.
obs <- z[z$observed, , drop = FALSE]
f_y <- pnorm(obs$resid_z)
obs$truncated_pit <- (f_y - obs$cdf_threshold) / (1 - obs$cdf_threshold)
obs$truncated_pit <- pmin(pmax(obs$truncated_pit, .Machine$double.eps),
                          1 - .Machine$double.eps)
qq_obs <- data.frame(theoretical = qnorm(ppoints(nrow(obs))),
                     residual = sort(qnorm(obs$truncated_pit)))
p2 <- ggplot2::ggplot(qq_obs, ggplot2::aes(theoretical, residual)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, colour = "steelblue", linewidth = 0.8) +
  ggplot2::geom_point(size = 1.5, alpha = 0.75) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    title = "Day 3: quantifiable conditional truncated-normal QQ",
    subtitle = sprintf("n=%d observed; accounts for individual-specific LLoQ truncation", nrow(obs)),
    x = "Standard normal quantile", y = "Truncated-normal PIT residual"
  ) + ggplot2::theme_bw()

fig_dir <- here::here("results", "figures")
qc_dir <- here::here("results", "qc")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
for (ext in c("png", "pdf")) {
  ggplot2::ggsave(file.path(fig_dir, paste0("qq_day3_conditional_randomized.", ext)),
                  p1, width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(fig_dir, paste0("qq_day3_conditional_quantifiable.", ext)),
                  p2, width = 8, height = 6, dpi = 300)
}
utils::write.csv(data.frame(id = z$id, treatment = z$A, baseline = z$baseline,
  status = z$status, mu = z$mu, censoring_probability = z$cdf_threshold,
  standardized_residual = z$resid_z),
  file.path(qc_dir, "qq_day3_conditional_diagnostics.csv"), row.names = FALSE)
utils::capture.output(summary(fit), file = file.path(qc_dir, "qq_day3_censored_model.txt"))
message("Saved two conditional QQ plots in ", fig_dir)
message("NOTE: Randomized censored residuals depend on seed and cannot verify the unobserved lower tail.")
message("NOTE: This is NOT the reverse-Kaplan-Meier QQ used in the JID paper; it is a conditional PIT alternative.")
print(p1)
print(p2)
