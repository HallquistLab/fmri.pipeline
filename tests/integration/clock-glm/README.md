# End-to-end clock-task GLM tests

Run this suite on a Slurm system with access to the clock fixture. It builds a fresh pipeline for each engine, runs the native estimators through `run_glm_pipeline()`, waits for cleanup, and reads the resulting images. Each run uses its own SQLite databases, caches and model directories.

| Profile | Participants and runs | Images |
| --- | --- | --- |
| `full` (default) | All 10 participants and all eight runs each | Original whole-brain grid and full time series |
| `smoke` | First four participants by ID, runs 1 and 2 | Central 17 × 17 × 17 voxels of each real acquisition, full time series |

Both profiles use the same clock and feedback events, reaction-time durations, clock-aligned entropy modulator, four CSF/white-matter nuisance regressors, and two leading-volume drops. The suite builds these from `gpa_base.rds` and `sample_2_L1_spec.yaml`. It retains the selected cohort so an unexpectedly dropped subject or run fails validation. Each group analysis uses the fixture’s mean-centered age model.

| Engine | Estimation and assertions |
| --- | --- |
| FSL | FEAT for each run, fixed-effects combination within participant, then a group GLM. Checks participant/run coverage, estimator provenance, COPE and T images, and the subject COPEs used in the group design. |
| SPM | A multi-run first-level model per participant followed by SPM group regression. Checks named contrast and statistic images and the specific subject contrast files supplied to the group fit. |
| AFNI | Multi-run REML per participant followed by group OLS. Checks run boundaries, the estimator, named coefficient and T images, and group input identities. |

## Run on Longleaf

Install the checkout under test in an R library visible to compute nodes. Load R, FSL, AFNI and licensed MATLAB in the submission shell. Set `SPM_PATH` to an SPM installation containing `spm.m` (the default on Longleaf is `/proj/mnhallqlab/lab_resources/spm12`). The engine workers inherit the submission R library paths. The default fixture is `/proj/mnhallqlab/projects/fmri.pipeline_test_data`. Set `FMRI_TEST_DATA` to use a copy with the same layout.

Run the driver in a compute allocation when processing whole-brain images. Finalization trims time series synchronously before submitting estimation jobs. From the repository root, use a **new directory on shared storage** for every invocation. The parent directory must exist. Allow space for trimmed time series and native fit outputs for each engine. The full profile processes 240 run-equivalents across the three engines.

```bash
test_lib="$SCRATCH/clock-e2e-r-lib"
mkdir -p "$test_lib"
export R_LIBS="$test_lib${R_LIBS:+:$R_LIBS}"
R CMD INSTALL --library="$test_lib" .
Rscript --vanilla tests/integration/clock-glm/run.R "$SCRATCH/clock-smoke-01" smoke
Rscript --vanilla tests/integration/clock-glm/run.R "$SCRATCH/clock-full-01" full
```

Select one or more engines with the third argument. Use `prepare` as a fourth argument to check the complete fixture, generate confounds and save finalized GPAs without submitting estimation jobs.

```bash
Rscript --vanilla tests/integration/clock-glm/run.R "$SCRATCH/clock-afni-01" smoke afni
Rscript --vanilla tests/integration/clock-glm/run.R "$SCRATCH/clock-prep-01" full fsl,spm,afni prepare
```

The driver writes the selected inputs to `selected_runs.csv`, the effective model to `clock.yaml`, and one result per backend to `results.csv`. A failed assertion or failed job makes the command exit nonzero. The default wait limit is 24 hours per backend, configurable with `GLM_TEST_TIMEOUT_HOURS`. Slurm jobs continue if the driver exits or times out. Each backend’s project configuration identifies its batch directory, scripts and logs.

To recheck completed outputs without resubmission, use the same installed package and pass the backend directories to the validator.

```bash
Rscript --vanilla tests/integration/clock-glm/validate.R "$SCRATCH/clock-full-01"/{fsl,spm,afni}
```

This external-data suite runs on explicit invocation. Ordinary package checks need neither the clock data nor the three estimators.
