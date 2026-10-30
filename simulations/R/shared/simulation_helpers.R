# Shared simulation helpers extracted verbatim from identical copies in the legacy scenario scripts.
# Statistical function bodies were not rewritten during this refactor.

NaiveImpute <- function(data, LOD = LOD_value) {
  # Identify rows where M == "LOD"
  lod_rows <- data[data$M == "LOD", , drop = FALSE]
  
  # If there are LOD rows, replace M with LOD_value / 2
  if (nrow(lod_rows) > 0) {
    lod_rows$M <- as.numeric(LOD) / 2
  }
  
  # Combine original non-LOD data with imputed LOD rows
  data <- rbind(data[data$M != "LOD", ], lod_rows)
  # Ensure M is numeric
  data$M <- as.numeric(data$M)
  # Ensure outcome and weights are numeric
  data$C <- as.numeric(data$C)
  data$w_norm <- rep(1, nrow(data))
  data <- data %>% arrange(id)
  return(data)
}

compute_eff <- function(data, indices, type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  sample_data <- data[indices, ]  # Resample data
  bres <- if (type == "NDE") gcomp.nde(sample_data) else gcomp.nie(sample_data)
  
  # Extract the corresponding estimate
  if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
  if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
  if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
}

gcom_boots_ci <- function(data, n_bootstrap = 1000, alpha = 0.05, 
                          type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  # Parallel bootstrap estimation using future_lapply()
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = nrow(data), replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    # Extract the corresponding estimate
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    
    return(NA)  # Handle unexpected measure values safely
  })
  
  # Convert to numeric vector and remove NAs
  boots_est <- na.omit(unlist(boots_est))
  
  # Check for empty bootstrap results
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA. Check data or computation methods.")
    return(NULL)
  }
  
  # Compute confidence intervals
  ci_lower <- quantile(boots_est, probs = alpha / 2, na.rm = TRUE)
  ci_upper <- quantile(boots_est, probs = 1 - alpha / 2, na.rm = TRUE)
  
  list(
    point_estimate = mean(boots_est, na.rm = TRUE),
    variance = var(boots_est, na.rm = TRUE),
    ci = c(lower = ci_lower, upper = ci_upper)
  )
}

gcom_boots_ci_m_out_n <- function(data, m, n_bootstrap = 1000, alpha = 0.05, 
                                  type = c("NDE", "NIE"), measure = c("RD", "RR", "OR")) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  
  boots_est <- numeric(n_bootstrap)
  
  # Parallel bootstrap estimation using future_lapply
  boots_est <- future_lapply(seq_len(n_bootstrap), function(i) {
    bootstrap_sample <- data[sample(nrow(data), size = m, replace = TRUE), ]
    
    # Compute bootstrap estimate based on type and measure
    bres <- if (type == "NDE") gcomp.nde(bootstrap_sample) else gcomp.nie(bootstrap_sample)
    
    # Extract the corresponding estimate
    if (measure == "RD") return(if (type == "NDE") bres$nde.rd else bres$nie.rd)
    if (measure == "RR") return(if (type == "NDE") bres$nde.rr else bres$nie.rr)
    if (measure == "OR") return(if (type == "NDE") bres$nde.or else bres$nie.or)
    
    return(NA)  # Handle unexpected measure values safely
  }, future.seed = TRUE)  # Ensure reproducibility
  
  boots_est <- na.omit(unlist(boots_est))  # Flatten and remove NAs
  
  if (length(boots_est) == 0) {
    warning("All bootstrap estimates are NA. Check data or computation methods.")
    return(NULL)
  }
  
  # Scale variance by (n/m) for m-out-of-n bootstrap correction
  point_estimate <- mean(boots_est, na.rm = TRUE)
  scaled_variance <- var(boots_est, na.rm = TRUE) * (m / nrow(data))
  ci <- quantile(boots_est, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  list(
    point_estimate = point_estimate,
    variance = scaled_variance,
    ci = ci
  )
}

gcomp.nde <- function(data) {
  # Extract variables from the dataset
  y <- data$Y        # Assuming Y is the outcome variable
  m <- data$M        # Assuming M is the mediator
  a <- data$A        # Assuming A is the treatment
  l1 <- data$L1      # Assuming L1 is the first covariate
  l2 <- data$L2      # Assuming L2 is the second covariate
  l3 <- data$L3      # Assuming L3 is the third covariate
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  weights <- data$w_norm  # Extract the weights
  
  # Combine L1, L2, L3 into a data frame for modeling convenience
  covariates <- data.frame(l1 = l1, l2 = l2, l3 = l3)
  
  # Fit a weighted linear model for y including m, a, and covariates
  lm_y <- glm(y ~ m + a + l1 + l2 + l3 + m*a, weights = weights, family = binomial)
  
  # Predict potential outcomes under different treatments
  pred_y1 <- predict(lm_y, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 1, m = m), type = "response")
  pred_y0 <- predict(lm_y, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 0, m = m), type = "response")
  
  # Fit a weighted linear model using the predicted outcome
  lm_y1 <- glm(pred_y1 ~ a + l1 + l2 + l3 + a*l1 + a*l2 + a*l3, 
               weights = weights, family = quasibinomial(link = "logit"))
  lm_y0 <- glm(pred_y0 ~ a + l1 + l2 + l3 + a*l1 + a*l2 + a*l3, 
               weights = weights, family = quasibinomial(link = "logit"))
  
  # Predict the causal effect when a = 0
  y10 <- predict(lm_y1, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 0), type = "response")
  y00 <- predict(lm_y0, newdata = data.frame(l1 = l1, l2 = l2, l3 = l3, a = 0), type = "response")
  
  # Calculate and return the estimate
  nde.rd <- mean(y10 - y00); nde.rr <- mean(y10)/mean(y00); nde.or <- (mean(y10)/(1-mean(y10))) / (mean(y00)/(1-mean(y00)))
  return(list(nde.rd = nde.rd, nde.rr = nde.rr, nde.or = nde.or))
}

gcomp.nie <- function(data) {
  # Extract variables from the dataset
  y <- data$Y        # Outcome variable
  m <- data$M        # Mediator
  a <- data$A        # Treatment
  l1 <- data$L1      # Covariate 1
  l2 <- data$L2      # Covariate 2
  l3 <- data$L3      # Covariate 3
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  weights <- data$w_norm  # Weights for weighted regression
  
  # Combine covariates into a data frame for modeling convenience
  covariates <- data.frame(l1 = l1, l2 = l2, l3 = l3)
  
  # Fit a weighted linear model for y including m, a, and covariates
  lm_y <- glm(y ~ m + a + l1 + l2 + l3 + m*a, weights = weights, family = binomial)
  
  # Predict potential outcomes under different treatments
  pred_y1 <- predict(lm_y, newdata = transform(covariates, a = 1, m = m),  type = "response")
  
  # Fit a weighted linear model using the predict potential outcomes
  lm_y1 <- glm(pred_y1 ~ a + l1 + l2 + l3, weights = weights, family = quasibinomial(link = "logit"))
  
  # Predict the causal effect when a = 0
  y11 <- predict(lm_y1, newdata = transform(covariates, a = 1),  type = "response")
  y10 <- predict(lm_y1, newdata = transform(covariates, a = 0),  type = "response")
  
  # Calculate and return the estimate
  nie.rd <- mean(y11 - y10); nie.rr <- mean(y11)/mean(y10); nie.or <- (mean(y11)/(1-mean(y11))) / (mean(y10)/(1-mean(y10)))
  # Return the estimate of the indirect effect
  return(list(nie.rd = nie.rd, nie.rr = nie.rr, nie.or = nie.or))
}

mod_update_hal <- function(data, basis_list = NULL) {
  # Check and assign w_norm if NULL (for the initial value of the weights)
  if (is.null(data$w_norm)) {
    data$w_norm <- 1
  }
  
  x_names <- c("L1", "L2", "L3", "A", "M")
  y_out_names <- "Y"
  weights_names <- "w_norm"
  
  # fit enumerate_basis on full data to get set of basis functions
  if (is.null(basis_list)) {
    basis_list <- hal9001::enumerate_basis(
      as.matrix(data[, ..x_names]),
      max_degree = 2, 
      smoothness_orders = 1, 
      num_knots = c(100, 75, 50)
    )
  }
  # Function to fit HAL model
  out_mod <- hal9001::fit_hal(
    X = as.matrix(data[, ..x_names]),
    Y = as.numeric(data[, get(y_out_names)]),
    weights = as.numeric(data[, get(weights_names)]),
    family = "binomial",
    basis_list = basis_list,
    fit_control = list(
      cv_select = TRUE,
      use_min = TRUE,
      lambda.min.ratio = 1e-5,
      prediction_bounds = "default"
    )
  )
  
  out_probY1 <- predict(out_mod, new_data = as.matrix(data[, ..x_names]), type = "response")
  out_probY0 <- 1 - out_probY1
  
  # Create a new column for predicted probabilities based on observed Y values
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)
  
  # Return predictions, coefficients, and terms
  return(list(
    out_prob = out_prob,
    basis_list = basis_list
  ))
}

os_compute_eff <- function(data, indices, type = c("NDE", "NIE"), fast = FALSE) {
  type <- match.arg(type)
  sample_data <- data[indices, ]  # Resample data
  
  # Choose function based on fast argument
  if (fast) {
    bres <- Onestep_med_effs_fast(dat_aug = sample_data, effect = type)
  } else {
    bres <- Onestep_med_effs(dat_aug = sample_data, effect = type)
  }
  
  # Extract estimate and EIF
  est_name <- if (type == "NDE") "nde_onestep_est" else "nie_onestep_est"
  
  return(list(
    eif = bres$eif,          # Assumes eif is always returned regardless of effect
    est = bres[[est_name]]
  ))
}
