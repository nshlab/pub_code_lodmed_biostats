# V6 fast numerical smoke regression for sim_iteration orchestration.
if (!requireNamespace("here", quietly=TRUE)) stop("Package 'here' is required.")
source(here::here("tests", "simulation", "smoke_helpers.R"))

scenarios <- data.frame(
  id=c("study1_EM1_ture_nuis_Gcom","study1_Hetero_CDE_HAL","study1_LLoQ_Gcom","study1_Wrong_GLM_Gcom","study2_Heterosced_CDE_HAL_1step","study2_LLoQ_1step"),
  study=c(rep("study1",4),rep("study2",2)),
  legacy=c(
    "Study1/EM1_ture_nuis_Gcom.R.txt","Study1/Hetero_CDE_HAL.R.txt","Study1/LLoQ_Gcom.R.txt","Study1/Wrong_GLM_Gcom.R.txt",
    "Study2/Heterosced_CDE_HAL_1step.R.txt","Study2/LLoQ_1step.R.txt"),
  refactored=c(
    "study1/EM1_ture_nuis_Gcom.R","study1/Hetero_CDE_HAL.R","study1/LLoQ_Gcom.R","study1/Wrong_GLM_Gcom.R",
    "study2/Heterosced_CDE_HAL_1step.R","study2/LLoQ_1step.R"),
  stringsAsFactors=FALSE
)

status <- vector("list", nrow(scenarios))
for (j in seq_len(nrow(scenarios))) {
  s <- scenarios[j,]
  cat("\n--- smoke:", s$id, "---\n")
  ans <- tryCatch(
    run_sim_iteration_smoke(
      here::here("tests","fixtures","simulation",s$legacy),
      here::here("simulations","scenarios",s$refactored),
      study=s$study
    ),
    error=function(e) list(passed=FALSE, max_abs_diff=NA_real_, error=conditionMessage(e))
  )
  msg <- if (!is.null(ans$error)) ans$error else if (isTRUE(ans$passed)) "PASS" else "Numerical result vectors differ"
  cat("passed:", isTRUE(ans$passed), "\n")
  if (!is.null(ans$max_abs_diff)) cat("max_abs_diff:", ans$max_abs_diff, "\n")
  status[[j]] <- data.frame(test=s$id, passed=isTRUE(ans$passed), message=msg, stringsAsFactors=FALSE)
}
status <- do.call(rbind,status)
print(status, row.names=FALSE)
out <- here::here("tests","simulation","smoke_test_status.csv")
write.csv(status,out,row.names=FALSE)
if (!all(status$passed)) stop("One or more fast numerical smoke tests failed. See: ", out)
cat("All fast sim_iteration numerical smoke tests passed.\n")
cat("NOTE: expensive estimators/bootstrap calls are stubbed; this validates orchestration equivalence.\n")
