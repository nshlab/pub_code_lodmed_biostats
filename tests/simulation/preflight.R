# Dependency preflight for full/expensive simulation reproduction.
# Rsolnp is known to fail source compilation with some R 4.5.x / current
# RcppArmadillo toolchains. It is not required for the V5 static/DGP tests or
# the V6 fast orchestration smoke tests.
packages <- c("data.table","dplyr","sl3","hal9001","Rcpp","medoutcon","truncnorm","xgboost","zoo","future.apply","ranger","Rsolnp")
status <- data.frame(package=packages, installed=vapply(packages, requireNamespace, logical(1), quietly=TRUE), row.names=NULL)
print(status, row.names=FALSE)
if (!all(status$installed)) {
  missing <- status$package[!status$installed]
  message("Missing packages affect full/expensive numerical runs only: ", paste(missing, collapse=", "))
  if ("Rsolnp" %in% missing) message("Known local compatibility note: CRAN Rsolnp 2.0.1 may fail to compile under R 4.5.x with current Rcpp/RcppArmadillo. Do not downgrade working core dependencies solely for this smoke test.")
}
