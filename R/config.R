# Application configuration -------------------------------------------------
ACTIV2_L_VARS <- c(
  "age_get60", "raceWhite", "enthnicityHispanic",
  "sexM", "bmi_catOverweight", "symptom_duration", "vaccine",
  "DIAB", "HPTN", "CURSMOKR", "M_base"
)
# Historical Table 2 estimator adjustment set (numerically checked in V8).
# Distinct from the 11-variable mediator fractional-imputation specification.
APPLICATION_ESTIMATOR_L_VARS <- c(
  "age_get60", "sexM", "bmi_catOverweight", "symptom_duration",
  "DIAB", "HPTN", "CURSMOKR"
)
ACTIV2_CORE_VARS <- c("A", "Y", "M", "C")
ACTIV2_REQUIRED_VARS <- c(ACTIV2_L_VARS, ACTIV2_CORE_VARS)
LOD_MED_SEED <- 2026L
set.seed(LOD_MED_SEED)

assert_activ2_columns <- function(data, extra = character()) {
  req <- unique(c(ACTIV2_REQUIRED_VARS, extra))
  missing <- setdiff(req, names(data))
  if (length(missing)) stop("Missing required ACTIV-2 variables: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}
