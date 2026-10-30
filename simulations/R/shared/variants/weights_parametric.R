# Shared verbatim implementation of Wet_step() used by 4 simulation scenarios.
# Extracted without rewriting the function body.

Wet_step <- function(data_imputed, mod_pred, beta_m_new, beta_m_0) {
  # Update probabilities
  #data_imputed$cen_prob <- pmin(pmax(mod_pred$cen_prob, 1e-40), 1 - 1e-40)
  data_imputed$out_prob <- mod_pred$out_prob
  
  # Filter censored data
  censored_data <- data_imputed[data_imputed$C == 0, ]
  
  # Compute new and initial weights for censored data
  cenwet_new <- compute_f_M(censored_data, beta_m_new)
  cenwet_int <- compute_f_M(censored_data, beta_m_0)
  
  censored_data$w <- censored_data$out_prob * (cenwet_new / cenwet_int)
  
  # Normalize weights by ID
  censored_data$w_norm <- ave(censored_data$w, censored_data$id, FUN = function(x) x / sum(x))
  
  # For uncensored data
  uncensored_data <- data_imputed[data_imputed$C == 1, ]
  uncensored_data$w_norm <- 1;   uncensored_data$w <- 1
  data_imputed <- rbind(censored_data, uncensored_data)
  
  # Remove temporary columns
  data_weted <- data_imputed[, !names(data_imputed) %in% c("out_prob"), with = FALSE]
  return(data_weted[order(data_weted$id),])
}
