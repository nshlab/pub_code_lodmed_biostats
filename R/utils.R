# Utility helpers. Pass ACTIV2_L_VARS explicitly in application analyses.

get_L_vars <- function(data, L_vars = NULL) {
  if (!is.null(L_vars)) return(intersect(L_vars, names(data)))
  if (exists("L_vars", envir = .GlobalEnv)) {
    lv <- get("L_vars", envir = .GlobalEnv)
    if (is.character(lv)) return(intersect(lv, names(data)))
  }
  # common-name heuristics
  patterns <- c("^age_", "^sex", "^bmi_cat", "^symptom_duration$", "^risk", "^DIAB$", "^CVD$", "^HPTN$", "^CURSMOKR$")
  matched <- unique(unlist(lapply(patterns, function(p) grep(p, names(data), value = TRUE))))
  matched <- intersect(matched, names(data))
  if (length(matched) > 0) return(matched)
  # fallback to old L1,L2,... naming
  names(data)[grepl("^L[0-9]+$", names(data))]
}
get_M_vars <- function(data) {
  names(data)[grepl("^M_?[0-9]*$", names(data))]
}
