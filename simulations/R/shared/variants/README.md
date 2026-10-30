# Shared implementation variants

Filenames describe each implementation. Same-named functions are not interchangeable. Each scenario sources only its original variant; all 12 files are preserved byte-for-byte. This index makes the otherwise opaque filenames traceable.

| Implementation | Used by |
|---|---|
| `em_parametric.R` | `study1/EM1_ture_nuis_Gcom`, `study1/LLoQ_Gcom` |
| `em_hal.R` | `study1/Hetero_CDE_HAL`, `study2/Heterosced_CDE_HAL_1step` |
| `impute_integer_indicator.R` | `study1/Hetero_CDE_HAL`, `study2/Heterosced_CDE_HAL_1step` |
| `impute_numeric_indicator.R` | `study1/EM1_ture_nuis_Gcom`, `study1/LLoQ_Gcom`, `study1/Wrong_GLM_Gcom`, `study2/LLoQ_1step` |
| `mediator_m_step.R` | `study1/EM1_ture_nuis_Gcom`, `study2/LLoQ_1step` |
| `onestep_mediation.R` | `study1/Hetero_CDE_HAL`, `study1/EM1_ture_nuis_Gcom`, `study1/LLoQ_Gcom`, `study1/Wrong_GLM_Gcom` |
| `weights_hal.R` | `study1/Hetero_CDE_HAL`, `study2/Heterosced_CDE_HAL_1step` |
| `weights_parametric.R` | `study1/EM1_ture_nuis_Gcom`, `study1/LLoQ_Gcom`, `study1/Wrong_GLM_Gcom`, `study2/LLoQ_1step` |
| `select_gcomp_ci.R` | `study1/Hetero_CDE_HAL`, `study1/EM1_ture_nuis_Gcom` |
| `mediator_lognormal_density.R` | `study1/Hetero_CDE_HAL`, `study1/EM1_ture_nuis_Gcom`, `study2/Heterosced_CDE_HAL_1step`, `study2/LLoQ_1step` |
| `update_outcome_glm.R` | `study1/EM1_ture_nuis_Gcom`, `study1/LLoQ_Gcom`, `study2/LLoQ_1step` |
| `sample_truncated_lognormal.R` | `study1/Hetero_CDE_HAL`, `study1/EM1_ture_nuis_Gcom`, `study2/Heterosced_CDE_HAL_1step`, `study2/LLoQ_1step` |

The two imputation variants differ in whether `C` is converted to integer or numeric. They remain separate to preserve original implementations. Filenames changed; function names and bodies did not.
