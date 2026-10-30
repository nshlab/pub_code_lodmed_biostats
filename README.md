# Evaluating causal indirect effects when mediators are left-censored by assay limit of quantification

Research code for mediation analysis with a left-censored mediator, including parametric and HAL-based fractional imputation, G-computation, one-step and targeted estimation, bootstrap routines, six simulation scenarios, and an ACTIV-2 application.

The repository includes an ACTIV-2 participant-level analysis dataset and its data dictionary under `data/release_data/`. 

## Structure

| Directory | Contents |
|---|---|
| `R/` | Application estimation, imputation, bootstrap, configuration and utilities |
| `simulations/` | Scenario registry, data generating process, shared implementations, runner, summaries and plotting |
| `analysis/` | Table S1, Figure 2, supplemental Q-Q diagnostics and Table 2 point-estimate workflow |
| `tests/` | Application controls, simulation regression checks and historical reference fixtures |
| `data/` | ACTIV-2 analysis dataset and data dictionary |
| `docs/` | Dependencies and imputation implementation notes |

## Software

Use R >= 4.1 and the author's tested package environment where available.

`here` is required by the entry points. Table S1 uses base R plus `here`; Figure 2 additionally requires `ggplot2`; Q-Q diagnostics also require `survival`. Estimation and simulation dependencies include `data.table`, `dplyr`, `future`, `future.apply`, `parallelly`, `hal9001`, `nadir`, `sl3`, `medoutcon`, `truncnorm` and learner backends.

See [dependency guidance](docs/DEPENDENCIES.md) and scenario-specific preflight checks. No validated lockfile or exact historical package versions are included. Install the required packages before running the scripts; scripts do not install packages automatically.

## Simulations

Run all commands from the repository root.

### List available scenarios

```bash
Rscript simulations/run.R --list
```

### Check dependencies

Check the dependencies required by the selected scenario before running it:

```bash
Rscript simulations/run.R study1_EM1_ture_nuis_Gcom --preflight=true
Rscript simulations/run.R study1_Hetero_CDE_HAL --preflight=true
Rscript simulations/run.R study2_Heterosced_CDE_HAL_1step --preflight=true
```

### Run a small pilot

Use two Monte Carlo replicates to check execution and inspect outputs:

```bash
Rscript simulations/run.R study1_EM1_ture_nuis_Gcom --replicates=10 --workers=1 --seed=2026 --run-id=pilot
```

The pilot reduces the Monte Carlo replicate count only. Nested estimation and bootstrap computations remain unchanged and may still be expensive.

### Run 300 Monte Carlo replicates

Examples for all six configured scenarios:

```bash
# Study 1: G-computation scenarios
Rscript simulations/run.R study1_EM1_ture_nuis_Gcom --replicates=300 --workers=1 --seed=2026 --run-id=mc300
Rscript simulations/run.R study1_Hetero_CDE_HAL --replicates=300 --workers=1 --seed=2026 --run-id=mc300
Rscript simulations/run.R study1_LLoQ_Gcom --replicates=300 --workers=1 --seed=2026 --run-id=mc300
Rscript simulations/run.R study1_Wrong_GLM_Gcom --replicates=300 --workers=1 --seed=2026 --run-id=mc300

# Study 2: one-step scenarios
Rscript simulations/run.R study2_Heterosced_CDE_HAL_1step --replicates=300 --workers=1 --seed=2026 --run-id=mc300
Rscript simulations/run.R study2_LLoQ_1step --replicates=300 --workers=1 --seed=2026 --run-id=mc300
```

These commands override the Monte Carlo replicate count only. Sample size, censoring rate, imputation count and estimator settings remain those configured for each scenario. 

### Use parallel workers

Where sufficient CPU cores and memory are available, parallelize the outer Monte Carlo replicates:

```bash
Rscript simulations/run.R study1_EM1_ture_nuis_Gcom --replicates=300 --workers=4 --seed=2026 --run-id=mc300_parallel
```

The worker count must not exceed the available cores. Inner bootstrap computations remain sequential.


### Summarize results

After a completed run, supply independently verified NDE and NIE truths:

```bash
Rscript simulations/summarize.R simulations/results/raw/study1_EM1_ture_nuis_Gcom/mc300 TRUE_NDE TRUE_NIE
```

Replace `TRUE_NDE` and `TRUE_NIE` with numerical values for the intended DGP. Inspect `status.csv` and replicate warnings before interpreting the summaries.

See [simulation instructions](simulations/README.md) for scenario defaults, RNG policy, checkpointing and output formats. These command examples have not been runtime-validated during the repository audit.


## ACTIV-2 application

The dataset and dictionary are located at:

- `data/release_data/ACTIV2_PUBLIC.csv`
- `data/release_data/ACTIV2_PUBLIC_DICTIONARY.csv`

Consult the dictionary for variable definitions, coding and measurement units.

```bash
Rscript analysis/run.R --list
Rscript analysis/run.R descriptive
Rscript analysis/run.R figure2
Rscript analysis/run.R qqplot
```

These tasks write outputs under `results/tables/`, `results/figures/` and `results/qc/`. 

Table S1 and Figure 2 use observed data without model-based imputation. The Day 3 assessment prefers the Day 3 measurement; otherwise, it selects the lower available Day 2/4 value, with ties assigned to Day 2. This rule uses assigned visit labels; actual collection dates are unavailable.

RNA values >2.0 log10 copies/mL are classified as quantifiable; values <=2.0 are classified as at/below LLoQ. Figure 2 panel (b) displays the seven available visits: Days 0, 2, 3, 4, 7, 14 and 28. Table S1 summarizes recorded numeric Day 0 values, including below-LLoQ codes.

### Table 2

Before estimating the effects in Table 2, first generate the imputed datasets using the imputation functions in `R/imputation.R`. Apply each of the three strategies—LLoQ/2 substitution, parametric log-normal fractional imputation, and semi-parametric conditional-density fractional imputation—and save the resulting datasets as:

- `data/private/my_data_naive_AN_imputed.csv`
- `data/private/my_data_aug_Lognormal_AN_imputed.csv`
- `data/private/my_data_aug_CDE_AN_imputed.csv`

See [imputation instructions](docs/IMPUTATION.md) for function-loading order and implementation details. The Table 2 estimation script reads these prepared datasets; it does not perform imputation automatically.


```bash
Rscript tests/application/audit_inputs.R
Rscript analysis/run.R table2
```

## Validation

```bash
Rscript tests/application/test_rna_classification.R
Rscript tests/application/test_day3_window.R
Rscript tests/simulation/run_all.R
Rscript tests/simulation/run_smoke_all.R
Rscript tests/simulation/check_workflow.R
```

The first two tests use synthetic controls. Simulation regression tests compare preserved function definitions and DGP behavior against historical fixtures. Smoke tests stub expensive fitting and bootstrap computations.

## Data availability, citation and license

The ACTIV-2 analysis dataset and dictionary are included under `data/release_data/`. Their inclusion and filenames do not independently establish authorization for public redistribution. 

# Citation

If you use this code and the ACTIV-2 data, please cite the accompanying manuscript:

Jiang, C., Hughes, M. D., and Hejazi, N. S. *Evaluating Causal Indirect Effects When Mediators Are Left-Censored by Assay Limit of Quantification*.  [arXiv:2605.20615](https://arxiv.org/abs/2605.20615).

Contributed code and documentation use the existing MIT [license](LICENSE). ACTIV-2 data and external packages are subject to their own applicable terms.