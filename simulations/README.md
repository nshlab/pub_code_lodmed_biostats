# LLoQ mediation simulations

Six scenarios and their supplied driver defaults are retained.

## Folder responsibilities

| Location | Purpose |
|---|---|
| `config/scenarios.R` | Authoritative settings and dependency registry |
| `R/shared/` | Existing DGP and reusable statistical functions |
| `R/shared/variants/` | Distinct historical implementations, with a usage index |
| `scenarios/study1/` | Four G-computation definitions |
| `scenarios/study2/` | Two one-step definitions |
| `R/workflow.R` | Definition loading, preflight, RNG, checkpointing and resume |
| `run.R` | Unified command-line entry point |
| `summarize.R` | Replicate exports and Monte Carlo metrics |
| `plots/` | Historical plotting scripts requiring separately supplied CSVs |
| `../tests/simulation/` | Regression, smoke and workflow checks |

## Commands from the repository root

```bash
Rscript simulations/run.R --list
Rscript simulations/run.R study1_LLoQ_Gcom --preflight=true
Rscript simulations/run.R study1_LLoQ_Gcom --replicates=2 --workers=1 --seed=2026 --run-id=pilot
Rscript simulations/run.R study1_LLoQ_Gcom --workers=2 --run-id=full
```

A pilot reduces Monte Carlo replicates only. Estimator/bootstrap settings remain unchanged, so even two replicates may be expensive. Dependencies must already be installed; scripts do not install packages.

## RNG, outputs and recovery

Each replicate uses L'Ecuyer-CMRG with a seed derived from the base seed and replicate ID. The manifest records settings, R/package versions and source hashes. Successful checkpoints are reused and missing/failed replicates are retried only when the manifest matches. Choose a new run ID after changing source, settings, replicate count or dependency versions.

Workers intercept the historical append-write and atomically save separate checkpoints. Only the outer replicate loop is parallel; inner bootstrap tasks are sequential. This is a new deterministic execution policy, not recovery of unspecified historical Monte Carlo draws.

Outputs include `manifest.rds`, `sessionInfo.txt`, `replicate_000001.rds` etc., `results.rds`, and `status.csv`. Errors and warnings are retained in replicate records. Inspect them before accepting a run. Do not launch the same scenario/run ID concurrently. After an interruption, remove a stale `.running` lock only after confirming no process still uses the run.

## Summaries

Supply approved numerical NDE/NIE truths in place of the placeholders:

```bash
Rscript simulations/summarize.R simulations/results/raw/study1_LLoQ_Gcom/full TRUE_NDE TRUE_NIE
```

This writes `estimates.csv` and `summary.csv`: bias, empirical variance, MSE, adjusted/standard coverage, Monte Carlo SEs and mean standard variance estimate. Incomplete/failed or nonfinite runs are rejected. Full named vectors, including nuisance parameters, remain in `results.rds`. Study 2's `selected_c` is the historical `optimal_c`, not `opt_m`.


The new summary format does not feed the historical plotting scripts automatically. Those scripts require their original two-row CSV schema and filenames. Precomputed CSVs are absent; the directories are placeholders for separately supplied validated files.

