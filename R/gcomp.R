# AUTO-REFACTORED FROM ORIGINAL APPLICATION CODE.
# Function bodies are preserved unless explicitly noted.

suppressPackageStartupMessages({
 library(data.table)
 library(dplyr)
 library(hal9001)
 library(future.apply)
})

gcomp.nde <- function(data, L_vars = NULL) {
  # Auto-detect L variables if not provided
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Extract variables from the dataset
  y <- data$Y
  m <- data$M
  a <- data$A
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  weights <- data$w_norm
  
  # Build covariates data frame dynamically (handle both data.table and data.frame)
  if (is.data.table(data)) {
    covariates <- data[, ..L_vars]
  } else {
    covariates <- data[, L_vars, drop = FALSE]
  }
  
  # Build outcome model formula dynamically
  formula_y <- paste("y ~", paste(c("m", "a", L_vars, paste("m", "a", sep = "*")), collapse = " + "))
  
  # Fit a weighted model for y including m, a, and covariates
  lm_y <- glm(as.formula(formula_y), data = cbind(data.frame(y = y, m = m, a = a), as.data.frame(covariates)), 
              weights = weights, family = binomial)
  
  # Predict potential outcomes under different treatments
  pred_data_a1 <- cbind(as.data.frame(covariates), a = 1, m = m)
  pred_data_a0 <- cbind(as.data.frame(covariates), a = 0, m = m)
  
  pred_y1 <- predict(lm_y, newdata = pred_data_a1, type = "response")
  pred_y0 <- predict(lm_y, newdata = pred_data_a0, type = "response")
  
  # Build second-stage model formula with interactions (only if L_vars is not empty)
  formula_parts <- c("a", L_vars)
  if (length(L_vars) > 0) {
    interactions <- paste("a", L_vars, sep = "*")
    formula_parts <- c(formula_parts, interactions)
  }
  formula_y1 <- paste("pred_y1 ~", paste(formula_parts, collapse = " + "))
  
  # Fit models using predicted outcomes
  lm_y1 <- glm(as.formula(formula_y1), data = cbind(data.frame(pred_y1 = pred_y1, a = a), as.data.frame(covariates)),
               weights = weights, family = quasibinomial(link = "logit"))
  lm_y0 <- glm(as.formula(gsub("pred_y1", "pred_y0", formula_y1)), 
               data = cbind(data.frame(pred_y0 = pred_y0, a = a), as.data.frame(covariates)),
               weights = weights, family = quasibinomial(link = "logit"))
  
  # Predict the causal effect when a = 0
  pred_data_a0_simple <- cbind(as.data.frame(covariates), a = 0)
  y10 <- predict(lm_y1, newdata = pred_data_a0_simple, type = "response")
  y00 <- predict(lm_y0, newdata = pred_data_a0_simple, type = "response")
  
  # Calculate and return the estimate
  y00[is.na(y00) == TRUE] <- 0
  y10[is.na(y10) == TRUE] <- 0
  nde.rd <- mean(y10 - y00)
  nde.rr <- mean(y10) / mean(y00)
  nde.or <- (mean(y10) / (1 - mean(y10))) / (mean(y00) / (1 - mean(y00)))
  
  return(list(nde.rd = nde.rd, nde.rr = nde.rr, nde.or = nde.or))
}

gcomp.nie <- function(data, L_vars = NULL) {
  # Auto-detect L variables if not provided
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Extract variables from the dataset
  y <- data$Y
  m <- data$M
  a <- data$A
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  weights <- data$w_norm
  
  # Build covariates data frame dynamically (handle both data.table and data.frame)
  if (is.data.table(data)) {
    covariates <- data[, ..L_vars]
  } else {
    covariates <- data[, L_vars, drop = FALSE]
  }
  
  # Build outcome model formula dynamically
  formula_y <- paste("y ~", paste(c("m", "a", L_vars, paste("m", "a", sep = "*")), collapse = " + "))
  
  # Fit a weighted model for y including m, a, and covariates
  lm_y <- glm(as.formula(formula_y), data = cbind(data.frame(y = y, m = m, a = a), as.data.frame(covariates)),
              weights = weights, family = binomial)
  
  # Predict potential outcomes under a = 1
  pred_data_a1 <- cbind(as.data.frame(covariates), a = 1, m = m)
  pred_y1 <- predict(lm_y, newdata = pred_data_a1, type = "response")
  
  # Build second-stage model formula (without interactions for NIE)
  formula_y1 <- paste("pred_y1 ~", paste(c("a", L_vars), collapse = " + "))
  
  # Fit model using predicted outcomes
  lm_y1 <- glm(as.formula(formula_y1), data = cbind(data.frame(pred_y1 = pred_y1, a = a), as.data.frame(covariates)),
               weights = weights, family = quasibinomial(link = "logit"))
  
  # Predict the causal effect for a = 1 vs a = 0
  pred_data_a1_simple <- cbind(as.data.frame(covariates), a = 1)
  pred_data_a0_simple <- cbind(as.data.frame(covariates), a = 0)
  
  y11 <- predict(lm_y1, newdata = pred_data_a1_simple, type = "response")
  y10 <- predict(lm_y1, newdata = pred_data_a0_simple, type = "response")
  
  # Calculate and return the estimate
  y11[is.na(y11) == TRUE] <- 0
  y10[is.na(y10) == TRUE] <- 0
  nie.rd <- mean(y11 - y10)
  nie.rr <- mean(y11) / mean(y10)
  nie.or <- (mean(y11) / (1 - mean(y11))) / (mean(y10) / (1 - mean(y10)))
  
  return(list(nie.rd = nie.rd, nie.rr = nie.rr, nie.or = nie.or))
}

gcomp_analysis <- function(data, k = 5, p_cen = censor_rate, 
                           n_bootstrap = 1000, B1 = 150, B2 = 100, 
                           nominal_coverage = 0.95) {
  
  # NDE Analysis
  nde_gcom_glm <- gcomp.nde(data)
  
  choose_ciNDE <- choose_ci(
    data = data, 
    p_cen = p_cen, 
    k = k, 
    type = "NDE", 
    measure = "RD",
    B1 = B1,
    B2 = B2,
    nominal_coverage = nominal_coverage
  )
  
  # m-out-of-n bootstrap with SE and p-value
  CI.nderes <- gcom_boots_ci_m_out_n_with_se(
    data, 
    m = choose_ciNDE$opt_m, 
    n_bootstrap = n_bootstrap, 
    type = "NDE", 
    measure = "RD"
  )
  
  # Regular bootstrap with SE and p-value
  CI.nderes.reg <- gcom_boots_ci_with_se(
    data, 
    n_bootstrap = n_bootstrap, 
    type = "NDE", 
    measure = "RD"
  )
  
  # NIE Analysis
  nie_gcom_glm <- gcomp.nie(data)
  
  choose_ciNIE <- choose_ci(
    data = data, 
    p_cen = p_cen, 
    k = k, 
    type = "NIE", 
    measure = "RD",
    B1 = B1,
    B2 = B2,
    nominal_coverage = nominal_coverage
  )
  
  # m-out-of-n bootstrap with SE and p-value
  CI.nieres <- gcom_boots_ci_m_out_n_with_se(
    data, 
    m = choose_ciNIE$opt_m, 
    n_bootstrap = n_bootstrap, 
    type = "NIE", 
    measure = "RD"
  )
  
  # Regular bootstrap with SE and p-value
  CI.nieres.reg <- gcom_boots_ci_with_se(
    data, 
    n_bootstrap = n_bootstrap, 
    type = "NIE", 
    measure = "RD"
  )
  
  # build a compact summary table
  summary_table <- data.frame(
    Effect = c("NDE", "NIE"),
    Point_Est = c(nde_gcom_glm$nde.rd, nie_gcom_glm$nie.rd),
    SE_Regular = c(CI.nderes.reg$se, CI.nieres.reg$se),
    SE_m_out_n = c(CI.nderes$se, CI.nieres$se),
    CI_Lower_Regular = c(CI.nderes.reg$ci_quantile[1], CI.nieres.reg$ci_quantile[1]),
    CI_Upper_Regular = c(CI.nderes.reg$ci_quantile[2], CI.nieres.reg$ci_quantile[2]),
    CI_Lower_m_out_n = c(CI.nderes$ci_quantile[1], CI.nieres$ci_quantile[1]),
    CI_Upper_m_out_n = c(CI.nderes$ci_quantile[2], CI.nieres$ci_quantile[2]),
    Pvalue_Regular = c(CI.nderes.reg$p_value, CI.nieres.reg$p_value),
    Pvalue_m_out_n = c(CI.nderes$p_value, CI.nieres$p_value),
    stringsAsFactors = FALSE
  )
  
  # Return comprehensive results, include summary_table
  return(list(
    NDE = list(
      point_estimate = nde_gcom_glm$nde.rd,
      point_estimate_rr = nde_gcom_glm$nde.rr,
      point_estimate_or = nde_gcom_glm$nde.or,
      # m-out-of-n results
      ci_m_out_n = CI.nderes$ci_quantile,
      ci_m_out_n_normal = CI.nderes$ci_normal,
      se_m_out_n = CI.nderes$se,
      variance_m_out_n = CI.nderes$variance,
      p_value_m_out_n = CI.nderes$p_value,
      # Regular bootstrap results
      ci_regular = CI.nderes.reg$ci_quantile,
      ci_regular_normal = CI.nderes.reg$ci_normal,
      se_regular = CI.nderes.reg$se,
      variance_regular = CI.nderes.reg$variance,
      p_value_regular = CI.nderes.reg$p_value,
      # Metadata
      optimal_m = choose_ciNDE$opt_m,
      optimal_c = choose_ciNDE$optimal_c,
      coverage_rates = choose_ciNDE$coverage_rates
    ),
    NIE = list(
      point_estimate = nie_gcom_glm$nie.rd,
      point_estimate_rr = nie_gcom_glm$nie.rr,
      point_estimate_or = nie_gcom_glm$nie.or,
      # m-out-of-n results
      ci_m_out_n = CI.nieres$ci_quantile,
      ci_m_out_n_normal = CI.nieres$ci_normal,
      se_m_out_n = CI.nieres$se,
      variance_m_out_n = CI.nieres$variance,
      p_value_m_out_n = CI.nieres$p_value,
      # Regular bootstrap results
      ci_regular = CI.nieres.reg$ci_quantile,
      ci_regular_normal = CI.nieres.reg$ci_normal,
      se_regular = CI.nieres.reg$se,
      variance_regular = CI.nieres.reg$variance,
      p_value_regular = CI.nieres.reg$p_value,
      # Metadata
      optimal_m = choose_ciNIE$opt_m,
      optimal_c = choose_ciNIE$optimal_c,
      coverage_rates = choose_ciNIE$coverage_rates
    ),
    summary_table = summary_table
  ))
}
