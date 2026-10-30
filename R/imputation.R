# Authoritative application imputation implementation.
# Source R/config.R and R/utils.R first; pass L_vars explicitly for application fits.
# See docs/IMPUTATION.md for comparison with the separately supplied Imputation.R.
# AUTO-REFACTORED FROM ORIGINAL APPLICATION CODE.
# Function bodies are preserved unless explicitly noted.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(hal9001)
  library(nadir)
})

parse_beta_m <- function(beta_m) {
  # Extract SD parameter (assumed to be last element named "sd")
  sd_param <- beta_m[names(beta_m) == "sd"]
  
  # Extract mean coefficients (all except "sd")
  mean_coefs <- beta_m[names(beta_m) != "sd"]
  
  list(
    mean_coefs = mean_coefs,
    sd_param = sd_param
  )
}

compute_mean_M <- function(data, beta_m, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  mean_coefs <- beta_m[names(beta_m) != "sd"]
  means <- numeric(nrow(data))
  
  for (coef_name in names(mean_coefs)) {
    if (coef_name == "b0") {
      # Intercept term
      means <- means + mean_coefs[coef_name]
    } else if (coef_name %in% L_vars) {
      # L variable terms
      means <- means + mean_coefs[coef_name] * data[[coef_name]]
    } else if (coef_name == "A") {
      # A term
      means <- means + mean_coefs[coef_name] * data$A
    } else if (grepl("^[A].*[L0-9]+$", coef_name)) {
      # Interaction terms (A:L interactions)
      interaction_parts <- strsplit(coef_name, "(?=L)", perl = TRUE)[[1]]
      if (length(interaction_parts) == 2) {
        var_name <- trimws(interaction_parts[2])
        if (var_name %in% L_vars) {
          means <- means + mean_coefs[coef_name] * data$A * data[[var_name]]
        }
      }
    } else if (grepl("^M", coef_name)) {
      # M term (if present)
      means <- means + mean_coefs[coef_name] * data$M
    }
  }
  
  return(means)
}

sample_density_proposal <- function(n, data, LOD, beta_m, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Compute means dynamically
  means <- compute_mean_M(data, beta_m, L_vars = L_vars)
  sd_param <- beta_m[names(beta_m) == "sd"]
  
  # Initialize an empty vector to store valid samples
  valid_samples <- numeric(0)
  
  # Log-normal sampling loop until we have enough valid samples
  while (length(valid_samples) < n) {
    remaining <- n - length(valid_samples)
    samples <- rlnorm(remaining, meanlog = means, sdlog = sd_param)
    valid_samples <- c(valid_samples, samples[samples <= LOD])
  }
  
  return(valid_samples[1:n])
}

compute_f_M <- function(data, beta_m, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Compute means dynamically
  means <- compute_mean_M(data, beta_m, L_vars = L_vars)
  sd_param <- beta_m[names(beta_m) == "sd"]
  
  # Compute the log-normal density for each M
  f_M <- dlnorm(data$M, meanlog = means, sdlog = sd_param)
  return(f_M)
}

compute_f_M_HAL <- function(data, L_vars = NULL, K = 3, epsilon = 1e-8) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Filter L_vars to only those that exist in data
  L_vars <- intersect(L_vars, names(data))

  required_vars <- c("M", "A", L_vars)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
  }
  
  if (is.null(data$w_norm)) {
    data$w_norm <- 1
  }
  
  # Create fold assignments for cross-validation
  n <- nrow(data)
  fold_ids <- sample(rep(1:K, length.out = n))
  f_M_all <- numeric(n)
  
  for (k in 1:K) {
    # Split data
    train_data <- data[fold_ids != k, ]
    valid_data <- data[fold_ids == k, ]
    
    # Build formula dynamically
    formula_str <- paste("M ~", paste(c(L_vars, "A"), collapse = " + "))
    
    # Fit heteroscedastic model on training data
    fit_k <- nadir::lnr_heteroskedastic_density(
      data = train_data,
      formula = as.formula(formula_str),
      mean_lnr = lnr_hal,
      var_lnr = lnr_hal,
      mean_lnr_args = list(
        weights = train_data$w_norm,
        smoothness_orders = 1,
        max_degree = 2,
        num_knots = c(100, 75, 50),
        fit_control = list(
          cv_select = TRUE,
          use_min = TRUE,
          lambda.min.ratio = 1e-5,
          prediction_bounds = "default"
        )
      ),
      var_lnr_args = list(
        weights = train_data$w_norm,
        smoothness_orders = 1,
        max_degree = 2,
        num_knots = c(100, 75, 50),
        fit_control = list(
          cv_select = TRUE,
          use_min = TRUE,
          lambda.min.ratio = 1e-5,
          prediction_bounds = "default"
        )
      )
    )
    
    # Predict density on validation data
    f_M_valid <- fit_k(valid_data)
    f_M_valid <- pmax(f_M_valid, epsilon)
    
    # Store results
    f_M_all[fold_ids == k] <- f_M_valid
  }
  
  loss <- mean(-log(f_M_all))
  
  return(list(
    f_M = f_M_all,
    loss = loss,
    fold_ids = fold_ids
  ))
}

Impute_step <- function(data, beta_m, S, LOD, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Duplicate and impute LOD rows
  imputed_rows_list <- lapply(1:nrow(data[data$M == "LOD", ]), function(i) {
    replicated_row <- data[data$M == "LOD", ][rep(i, S), , drop = FALSE]
    replicated_row$M <- sample_density_proposal(
      n = S,
      data = as.data.frame(replicated_row),
      LOD = LOD,
      beta_m = beta_m,
      L_vars = L_vars
    )
    return(replicated_row)
  })
  
  # Combine non-LOD data with imputed LOD rows
  data <- rbind(data[data$M != "LOD", ], do.call(rbind, imputed_rows_list))
  # Convert all character columns to numeric or factor as needed
  data <- data %>%
    mutate(across(where(is.character), as.numeric))
  # Ensure the outcome and weights are numeric
  data$C <- as.numeric(data$C)
  # Drop rows with any missing values, then return ordered by id
  data <- data[complete.cases(data), ]
  return(data[order(data$id), ])
}

mod_update_glm <- function(data, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Filter L_vars to only those that exist in data
  L_vars_present <- intersect(L_vars, names(data))
  
  # Check and assign w_norm if NULL
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  
  # Build formula dynamically with only existing L variables
  formula_str <- paste("Y ~", paste(c(L_vars_present, "A", "M", "A*M"), collapse = " + "))
  out_mod <- glm(as.formula(formula_str), data = data, family = binomial, weights = w_norm)
  out_probY1 <- predict(out_mod, type = "response")
  out_probY0 <- 1 - out_probY1
  
  # Create a new column for predicted probabilities based on observed Y values
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)
  
  # Return predictions and coefficients
  return(list(
    out_prob = out_prob,
    out_coef = out_mod$coefficients
  ))
}

mod_update_hal <- function(data, L_vars = NULL, basis_list = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)

  if (is.null(data$w_norm)) {
    data$w_norm <- 1
  }

  x_names <- c(L_vars, "A", "M")
  y_out_names <- "Y"
  weights_names <- "w_norm"

  # Validate predictors, outcome and weights exist
  missing_x <- setdiff(x_names, names(data))
  if (length(missing_x) > 0) {
    stop("Missing predictor columns: ", paste(missing_x, collapse = ", "))
  }
  if (!y_out_names %in% names(data)) {
    stop("Outcome column '", y_out_names, "' not found in data")
  }
  Y_vec <- as.numeric(data[[y_out_names]])
  if (any(is.na(Y_vec))) {
    stop("NA detected in outcome column '", y_out_names, "': hal9001::fit_hal does not support missingness in Y")
  }
  if (!weights_names %in% names(data) || all(is.na(data[[weights_names]]))) {
    data[[weights_names]] <- 1
  }
  weights_vec <- as.numeric(data[[weights_names]])

  # Generate basis functions if not provided
  if (is.null(basis_list)) {
    basis_list <- hal9001::enumerate_basis(
      as.matrix(data[, x_names]),
      max_degree = 2,
      smoothness_orders = 1,
      num_knots = c(100, 75, 50)
    )
  }

  # Fit HAL model using direct column extraction (no get())
  out_mod <- hal9001::fit_hal(
    X = as.matrix(data[, x_names]),
    Y = Y_vec,
    weights = weights_vec,
    family = "binomial",
    basis_list = basis_list,
    fit_control = list(
      cv_select = TRUE,
      use_min = TRUE,
      lambda.min.ratio = 1e-5,
      prediction_bounds = "default"
    )
  )

  out_probY1 <- predict(out_mod, new_data = as.matrix(data[, x_names]), type = "response")
  out_probY0 <- 1 - out_probY1
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)

  return(list(
    out_prob = out_prob,
    basis_list = basis_list
  ))
}

Wet_step <- function(data_imputed, mod_pred, beta_m_new, beta_m_0, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data_imputed)
  
  # Update probabilities
  data_imputed$out_prob <- mod_pred$out_prob
  
  # Filter censored data
  censored_data <- data_imputed[data_imputed$C == 0, ]
  
  # Compute new and initial weights for censored data
  cenwet_new <- compute_f_M(censored_data, beta_m_new, L_vars = L_vars)
  cenwet_int <- compute_f_M(censored_data, beta_m_0, L_vars = L_vars)
  
  censored_data$w <- censored_data$out_prob * (cenwet_new / cenwet_int)
  censored_data$w_norm <- ave(censored_data$w, censored_data$id, FUN = function(x) x / sum(x))
  
  # For uncensored data
  uncensored_data <- data_imputed[data_imputed$C == 1, ]
  uncensored_data$w_norm <- 1; uncensored_data$w <- 1
  data_imputed <- rbind(censored_data, uncensored_data)
  
  # Remove temporary columns
  data_weted <- data_imputed[, !names(data_imputed) %in% c("out_prob"), with = FALSE]
  return(data_weted[order(data_weted$id),])
}

Wet_step_HAL <- function(data_imputed, mod_pred, beta_m_0, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data_imputed)
  
  # Update probabilities
  data_imputed$out_prob <- mod_pred$out_prob
  
  # Compute HAL-based M density
  Mdens_HAL <- compute_f_M_HAL(data = data_imputed, L_vars = L_vars)
  data_imputed$Mdens_new <- Mdens_HAL$f_M
  Mdens_error <- Mdens_HAL$loss
  
  censored_data <- data_imputed[data_imputed$C == 0, ]
  
  cendens_new <- censored_data$Mdens_new
  cendens_int <- compute_f_M(censored_data, beta_m_0, L_vars = L_vars)
  
  censored_data$w <- censored_data$out_prob * (cendens_new / cendens_int)
  censored_data$w_norm <- ave(censored_data$w, censored_data$id, FUN = function(x) x / sum(x))
  
  uncensored_data <- data_imputed[data_imputed$C == 1, ]
  uncensored_data$w_norm <- 1; uncensored_data$w <- 1
  data_imputed <- rbind(censored_data, uncensored_data)
  
  data_weted <- data_imputed[, !names(data_imputed) %in% c("out_prob", "Mdens_new"), with = FALSE]
  
  return(list(
    wet_dat = data_weted[order(data_weted$id), ],
    dens_errors = Mdens_error
  ))
}

M_stepf_M <- function(data, beta_m, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Define an internal function for optimization
  neg_log_likelihood <- function(beta_m, data, L_vars) {
    # Extract weights
    weights <- data$w_norm
    # Compute means dynamically
    means <- compute_mean_M(data, beta_m, L_vars = L_vars)
    sd_param <- beta_m[names(beta_m) == "sd"]
    
    # Compute log-likelihood for each data point
    log_likelihoods <- log(dlnorm(data$M, meanlog = means, sdlog = sd_param))
    
    # Incorporate weights into the likelihood calculation
    weighted_log_likelihood <- weights * log_likelihoods
    
    # Return the negative of the weighted log-likelihood sum
    return(-sum(weighted_log_likelihood))
  }
  
  optim_results <- optim(par = beta_m, fn = neg_log_likelihood, data = data, L_vars = L_vars)
  
  return(beta_m_new = optim_results$par)
}

EM_algorithm <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 50, tol = 1e-20, verbose = TRUE, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  iter <- 0; diff <- Inf; rel_diff <- Inf; beta_m <- beta_m_0 
  
  data_imputed <- Impute_step(data, beta_m_0, S, LOD, L_vars = L_vars)
  mod_pred <- mod_update_glm(data_imputed, L_vars = L_vars)
  
  start_time <- Sys.time()
  
  while (iter < max_iter && diff > tol && rel_diff > tol) {
    iter_time <- Sys.time()
    
    # W-Step
    data_weted <- Wet_step(data_imputed, mod_pred, beta_m_new = beta_m, beta_m_0 = beta_m_0, L_vars = L_vars)
    
    # M-Step
    beta_m_new <- M_stepf_M(data = data_weted, beta_m, L_vars = L_vars)
    
    # Compute convergence criteria
    diff <- sum((beta_m_new - beta_m)^2)
    out_coef_old <- mod_pred$out_coef
    
    # Parameter update
    beta_m <- beta_m_new
    mod_pred <- mod_update_glm(data_weted, L_vars = L_vars)
    
    # Update diff with coefficient changes
    diff <- diff + sum((mod_pred$out_coef - out_coef_old)^2)
    rel_diff <- diff / (sum(beta_m^2) + sum(out_coef_old^2) + 1e-8)
    iter_time_end <- Sys.time()
    iter_duration <- round(difftime(iter_time_end, iter_time, units = "secs"), 4)
    
    if (verbose) {
      cat("Iteration:", iter, 
          "| Absolute Diff:", diff, 
          "| Relative Diff:", rel_diff,
          "| Time:", iter_duration, "s\n")
    }
    
    iter <- iter + 1
  }
  
  total_runtime <- round(difftime(Sys.time(), start_time, units = "secs"), 4)
  
  if (iter >= max_iter) {
    message("EM Algorithm stopped: Maximum iterations reached (", max_iter, ").")
  } else if (diff <= tol || rel_diff <= tol) {
    message("EM Algorithm converged successfully in ", iter, " iterations.")
  }
  
  cat("Total Runtime:", total_runtime, "seconds\n")
  
  return(list(beta_m = beta_m, data_weted = data_weted))
}

EM_algorithm_HAL <- function(data, beta_m_0, S = 2, LOD = LOD_value, max_iter = 100, tol = 1e-5, verbose = TRUE, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  iter <- 0; diff <- Inf; rel_diff <- Inf
  
  # I-Step: Initial Imputation
  data_imputed <- Impute_step(data, beta_m_0, S, LOD, L_vars = L_vars)
  halmod_pred <- mod_update_hal(data = data_imputed, L_vars = L_vars)
  halbasis <- halmod_pred$basis_list
  
  wet_step <- Wet_step_HAL(data_imputed = data_imputed, mod_pred = halmod_pred, beta_m_0 = beta_m_0, L_vars = L_vars)
  data_weted <- wet_step$wet_dat
  errors <- wet_step$dens_errors
  
  start_time <- Sys.time()
  
  while (iter < max_iter && diff > tol && rel_diff > tol) {
    iter_time <- Sys.time()
    
    # W-Step
    wet_step <- Wet_step_HAL(data_imputed = data_weted, mod_pred = halmod_pred, beta_m_0 = beta_m_0, L_vars = L_vars)
    data_weted <- wet_step$wet_dat
    errors <- wet_step$dens_errors
    
    # M-Step
    errors_new <- compute_f_M_HAL(data = data_weted, L_vars = L_vars)$loss
    
    # Convergence Criteria
    diff <- sum((errors_new - errors)^2)
    
    # Parameter Update
    halmod_pred <- mod_update_hal(data_weted, L_vars = L_vars, basis_list = halbasis)
    
    rel_diff <- diff / (sum(errors^2) + 1e-8)
    iter_time_end <- Sys.time()
    iter_duration <- round(difftime(iter_time_end, iter_time, units = "secs"), 4)
    
    if (verbose) {
      cat("Iteration:", iter, 
          "| Absolute Diff:", diff, 
          "| Relative Diff:", rel_diff,
          "| Time:", iter_duration, "s\n")
    }
    
    iter <- iter + 1
  }
  
  total_runtime <- round(difftime(Sys.time(), start_time, units = "secs"), 4)
  
  if (iter >= max_iter) {
    message("EM Algorithm stopped: Maximum iterations reached (", max_iter, ").")
  } else if (diff <= tol || rel_diff <= tol) {
    message("EM Algorithm converged successfully in ", iter, " iterations.")
  }
  
  cat("Total Runtime:", total_runtime, "seconds\n")
  
  return(list(data_weted = data_weted))
}

fractional_impute <- function(data, beta_m_0, S = 20, LOD, max_iter = 50, tol = 1e-20, verbose = TRUE, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Run EM algorithm
  em_result <- EM_algorithm(data, beta_m_0, S = S, LOD = LOD, max_iter = max_iter, tol = tol, verbose = verbose, L_vars = L_vars)
  
  # Extract components from EM result
  data_aug <- em_result$data_weted
  beta_mres <- em_result$beta_m
  
  # Return comprehensive results
  return(list(
    beta_m = beta_mres,
    data_aug = data_aug,
    L_vars = L_vars,
    em_full_result = em_result
  ))
}

fractional_impute_CDE <- function(data, beta_m_0, S = 20, LOD, max_iter = 25, tol = 1e-5, verbose = TRUE, L_vars = NULL) {
  if (is.null(L_vars)) L_vars <- get_L_vars(data)
  
  # Run EM algorithm
  em_result <- EM_algorithm_HAL(data, beta_m_0, S = S, LOD = LOD, max_iter = max_iter, tol = tol, verbose = verbose, L_vars = L_vars)
  
  # Extract components from EM result
  data_aug <- em_result$data_weted
  
  # Return comprehensive results
  return(list(
    data_aug = data_aug,
    L_vars = L_vars,
    em_full_result = em_result
  ))
}

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
  return(data)
}
