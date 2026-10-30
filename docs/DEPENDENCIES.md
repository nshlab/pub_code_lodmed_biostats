# R dependencies

Packages referenced by active code are listed below. No automatic installation occurs during analysis. 

`MetBrewer`, `R6`, `Rcpp`, `Rsolnp`, `caret`, `data.table`, `dplyr`, `earth`, `future`, `future.apply`, `ggplot2`, `gt`, `gtsummary`, `hal9001`, `haven`, `here`, `latex2exp`, `medoutcon`, `nadir`, `parallelly`, `purrr`, `ranger`, `readr`, `sl3`, `stringr`, `survival`, `tibble`, `tidyr`, `tidyverse`, `tlverse`, `truncnorm`, `xgboost`, `zoo`


## Workflow-specific requirements

- Minimum syntax: R >=4.1 for the native pipe; no tested minimum package versions are established.
- Entry points/unit controls: `here`.
- Table S1: `here` and base/recommended R; Figure 2 additionally `ggplot2`; Q-Q additionally `survival`.
- Simulation runtime: `here`, `future`, `future.apply`, `parallelly`, plus `packages` in `simulations/config/scenarios.R` for the selected scenario. Package preflight is necessary but cannot prove learner API compatibility.
- Table 2: imports in `R/gcomp.R` and `R/medoutcon_estimators.R`, including `nadir`, `sl3`, `R6`, `data.table`, `dplyr`, `hal9001`, `medoutcon`, `truncnorm`, `xgboost`, `future.apply`, `ranger`, `earth`; `caret` is optional in the existing implementation. Learner backends may require additional packages such as `glmnet`, `gam`, and `nnls`.
- Historical plotting: `here`, `dplyr`, `purrr`, `readr`, `stringr`, `tibble`, `tidyr`, `ggplot2`, `latex2exp`, `MetBrewer` and separately supplied summary CSVs.
