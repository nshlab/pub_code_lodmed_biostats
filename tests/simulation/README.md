# Simulation checks

Run from the repository root with `here` and the required test dependencies installed.

```bash
Rscript tests/simulation/run_all.R
Rscript tests/simulation/run_smoke_all.R
Rscript tests/simulation/check_workflow.R
Rscript tests/simulation/preflight.R
```

`run_all.R` compares common function bodies and fixed-seed DGP outputs with historical function-definition fixtures. `sim_iteration` is intentionally excluded from body equality because its output path/orchestration was refactored previously. `run_smoke_all.R` compares result assembly using real DGP draws and deterministic stubs for expensive fitting/bootstrap. It exercises the EM-MLE branch, not every branch or full scientific computation. `check_workflow.R` uses a stochastic stand-in to test checkpoint/resume/seed behavior and worker invariance where cores permit. `preflight.R` checks full-run dependencies.

Historical sources required by those comparisons are reduced to non-executable `.R.txt` fixtures under `tests/fixtures/simulation/`. Package loads, installation instructions, exploratory fits, parallel plans and driver tails were removed. The intercepted legacy `sim_iteration` output path is neutralized. Statistical function text is preserved; tests read definitions only. These are frozen test references, not duplicate supported analysis modules. No `legacy/` directory is needed.
