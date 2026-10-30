# Run from repository root: Rscript tests/application/01_audit_inputs.R
files <- c(naive="my_data_naive_AN_imputed.csv", lognormal="my_data_aug_Lognormal_AN_imputed.csv", cde="my_data_aug_CDE_AN_imputed.csv")
base <- file.path("data", "private")
if (!all(file.exists(file.path(base, files)))) stop("Missing restricted application CSV inputs")
vars <- c("age_get60","raceWhite","enthnicityHispanic","sexM","bmi_catOverweight","symptom_duration","vaccine","DIAB","HPTN","CURSMOKR","M_base")
for (variant in names(files)) {
  d <- read.csv(file.path(base, files[[variant]]), check.names=FALSE)
  if (!all(c("A","Y","M","C","w_norm",vars) %in% names(d))) stop("Required columns missing: ", variant)
  id <- if ("PUID" %in% names(d)) "PUID" else if ("id" %in% names(d)) "id" else NA_character_
  cat("\n", variant, ": rows=", nrow(d), ", participants=", if(is.na(id)) "UNKNOWN" else length(unique(d[[id]])), ", missing M_base rows=", sum(is.na(d$M_base)), "\n",sep="")
  cat("missing values by adjustment variable:\n"); print(colSums(is.na(d[vars])))
  if (!is.na(id)) {
    w <- aggregate(d$w_norm, list(d[[id]]), sum)
    cat("weight sums range:", paste(range(w$x), collapse=" to "), "\n")
  }
}
cat("\nAudit only. This script does not execute estimators or validate manuscript Table 2.\n")
