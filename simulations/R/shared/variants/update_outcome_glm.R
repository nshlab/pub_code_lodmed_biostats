# Shared verbatim implementation of mod_update_glm() used by 3 simulation scenarios.
# Extracted without rewriting the function body.

mod_update_glm <- function(data) {
  # Check and assign w_norm if NULL
  if (is.null(data$w_norm)) {data$w_norm <- 1}
  # Censoring Modeling (parametric model, does not include Y)
  #cen_mod <- glm(C ~ L1 + L2 + L3 + A + M + A*M, data = data,family = binomial, weights = w_norm)
  
  #cen_prob <- predict(cen_mod, type = "response")
  # Outcome Modeling
  out_mod <- glm(Y ~ L1 + L2 + L3 + A + M + A*M, data = data, family = binomial, weights = w_norm)
  out_probY1 <- predict(out_mod, type = "response")
  out_probY0 <- 1 - out_probY1
  
  # Create a new column for predicted probabilities based on observed Y values
  out_prob <- ifelse(data$Y == 1, out_probY1, out_probY0)
  
  # Return predictions and coefficients
  return(list(
    #cen_prob = cen_prob, 
    out_prob = out_prob,
    #cen_coef = cen_mod$coefficients,
    out_coef = out_mod$coefficients
  ))
}
