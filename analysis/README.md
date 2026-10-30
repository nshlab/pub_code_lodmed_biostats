# ACTIV-2 application

Run from the repository root after obtaining the authorized inputs in [data instructions](../data/README.md).

| Task | Command | Outputs |
|---|---|---|
| Table S1 draft | `Rscript analysis/run.R descriptive` | `results/tables/tableS1.csv`, baseline/category/missingness QC and notes |
| Figure 2 | `Rscript analysis/run.R figure2` | `results/figures/figure2a.pdf/png`, `figure2b.pdf/png`, aggregate counts |
| Supplemental diagnostics | `Rscript analysis/run.R qqplot` | Conditional Gaussian Q-Q plots and diagnostic summaries in `results/` |
| Table 2 point estimates | `Rscript analysis/run.R table2` | `results/tables/table2_v9_point_estimates.csv`, comparison and `table2_v10_point_estimates_ONLY.csv` |

`Rscript main.R TASK` is a compatibility alias. Output basenames containing v9/v10 are retained to preserve existing consumers; these are not alternate implementations. `05_record_session.R` runs before Table 2, and `04_export_table2.R` consumes its benchmark comparison.

Day 3 uses assigned visits 2/3/4, preferring Day 3 and otherwise the smaller backup value. Selected-record status is paired to the measurement. Other visits receive no window. Seven RNA visits are available; the full historical panel (b) requires additional unavailable visits. Q-Q diagnostics use quantifiable baseline measurements and fit a left-censored Gaussian model adjusted for treatment and baseline RNA.

Table S1 preserves observed numeric baseline codes, including below-LLoQ codes, without historical replacement-to-1/BASE fallback. 

Table 2 reads the three prepared inputs named in [data instructions](../data/README.md). It sources `R/config.R`, `R/utils.R`, `R/gcomp.R`, and `R/medoutcon_estimators.R` in that order. It does not call imputation. The seven-variable estimator set remains distinct from the 11-variable imputation specification. Complete-case filtering is preserved. 

Reusable imputation: source `R/config.R`, then `R/utils.R`, then `R/imputation.R`; pass `L_vars` explicitly. Bootstrap routines require their corresponding estimator functions and attached packages to be loaded first. No raw-to-Table-2 generation script is asserted. See [imputation comparison](../docs/IMPUTATION.md).
