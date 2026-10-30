# Imputation implementation

`R/imputation.R` is the sole application module. 

Load in order:

```r
source("R/config.R")
source("R/utils.R")
source("R/imputation.R")
```

Pass the intended `L_vars` explicitly. `get_L_vars()` is defined in `R/utils.R`; importing the comparison file's narrower helper would replace its application-aware behavior. `R/config.R` also sets the existing seed to 2026. Functions whose default LOD references `LOD_value` should be called with an explicit authorized LOD.

## Differences from the separately supplied file

| Component | Existing authoritative module | Separately supplied file | Release decision |
|---|---|---|---|
| Covariate helper | Utility supports explicit/global lists, application-name heuristics, then L1/L2/etc. | Defines a helper selecting only L-number names | Keep existing helper and loading order |
| HAL density / GLM covariates | Intersects specified covariates with available columns | Does not intersect in those functions | Preserve existing behavior; users must check explicit specifications |
| `Impute_step` | Removes incomplete rows before ordering by ID | Retains incomplete rows | Meaningful difference in analysis rows; do not merge |
| `mod_update_hal` | Ordinary column selection, explicit outcome/weight extraction and missing-outcome checks | data.table-specific `..x_names` / `get()` selection; fewer checks | Preserve existing implementation; container-class differences require dedicated runtime study |
| `fractional_impute` | Default S=20 | Default S=2 | Preserve S=20 |
| `fractional_impute_CDE` | Default S=20, max_iter=25 | Default S=2, max_iter=100 | Preserve existing defaults |
| `NaiveImpute` | Returns existing assembled order | Sorts by ID before returning | Preserve row order |
| Top-level demonstrations | None | Uses external `data`, runs fits and reads a CSV | Exclude comparison file and demonstrations |