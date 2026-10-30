# Compute point estimates for Table 2.
options(warn = 1)
source(file.path("R", "config.R"))
source(file.path("R", "utils.R"))
source(file.path("R", "gcomp.R"))
source(file.path("R", "medoutcon_estimators.R"))
inputs <- c(naive="my_data_naive_AN_imputed.csv",
            lognormal="my_data_aug_Lognormal_AN_imputed.csv",
            cde="my_data_aug_CDE_AN_imputed.csv")
base <- file.path("data", "private")
if (!all(file.exists(file.path(base, inputs)))) stop("Missing restricted inputs in data/private")
L <- APPLICATION_ESTIMATOR_L_VARS
results <- list()
for (variant in names(inputs)) {
  dat <- read.csv(file.path(base, inputs[[variant]]), check.names=FALSE)
  req <- c("A","Y","M","w_norm",L)
  if (!all(req %in% names(dat))) stop("Missing columns for ", variant, ": ", paste(setdiff(req,names(dat)),collapse=","))
  dat <- dat[complete.cases(dat[,req,drop=FALSE]),,drop=FALSE]
  # Required by medoutcon estimators; matches the validated Stage 4 runner.
  dat$Y_scale <- as.numeric(dat$Y)
  id <- if ("PUID" %in% names(dat)) "PUID" else if ("id" %in% names(dat)) "id" else NA_character_
  n_participants <- if (is.na(id)) NA_integer_ else length(unique(dat[[id]]))
  for (estimator in c("gcomp","onestep","tmle")) for (effect in c("NDE","NIE")) {
    message(variant, " / ", estimator, " / ", effect)
    fit <- tryCatch({
      if (estimator == "gcomp") {
        if (effect == "NDE") gcomp.nde(dat,L_vars=L)$nde.rd else gcomp.nie(dat,L_vars=L)$nie.rd
      } else {
        x <- Onestep_med_effs(dat_aug=dat,effect=effect,L_vars=L,estimator=estimator)
        x[[if (effect=="NDE") "nde_onestep_est" else "nie_onestep_est"]]
      }
    },error=function(e)e)
    error <- if (inherits(fit,"error")) conditionMessage(fit) else ""
    value <- if (nzchar(error)) NA_real_ else as.numeric(fit)
    if (!nzchar(error) && (length(value)!=1L || !is.finite(value))) {value<-NA_real_;error<-"Non-finite or non-scalar estimate"}
    results[[length(results)+1L]] <- data.frame(estimator=estimator,variant=variant,effect=effect,
      estimate=value,error=error,n_rows=nrow(dat),n_participants=n_participants,L_vars=paste(L,collapse=";"))
    dir.create("results/tables",recursive=TRUE,showWarnings=FALSE)
    write.csv(do.call(rbind,results),"results/tables/table2_v9_point_estimates.csv",row.names=FALSE)
  }
}
source("tests/application/06_compare_v9.R")
