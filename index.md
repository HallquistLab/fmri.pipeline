# fmri.pipeline

`fmri.pipeline` is an R package for specifying, running, and tracking
task fMRI general linear models with FSL, AFNI, and SPM. Use the same
trial, run, and participant tables to generate suite-specific designs
and contrasts, then run them locally or on a cluster.

The package is built around a single analysis object, usually called
`gpa`, that stores the input tables, model specifications, output
locations, backend choices, scheduler settings, and job-tracking state
for one analysis.

Current package version: `0.5`.

## Why this package exists

Large fMRI GLM projects have a familiar failure mode: the actual
analysis is a constellation of design files, shell scripts, scheduler
jobs, and partially remembered choices. `fmri.pipeline` tries to make
that constellation explicit. The package helps you:

- define run-level, subject/session-level, and group-level models in one
  place;
- generate backend-specific analysis files for FSL, AFNI, and SPM;
- submit staged analyses through Slurm, TORQUE, or local execution;
- keep track of job dependencies, logs, and analysis status;
- reuse the same input tables and YAML model specification across model
  revisions;
- support longitudinal analyses where runs are nested within sessions
  and sessions are nested within subjects.

Choose an engine with `glm_software = "fsl"`, `"afni"`, or `"spm"` in
[`setup_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md).
FSL fits runs separately and can carry their variances into higher-level
models. AFNI estimates shared task effects across runs. SPM models runs
as separate SPM sessions by default and can contrast their coefficients
within one first-level model. The [GLM engine
guide](https://hallquistlab.github.io/fmri.pipeline/articles/glm-engines.html)
follows a two-run task through each suite and shows where to inspect its
design and contrast maps.

`fmri.pipeline` is not a preprocessing pipeline. It expects preprocessed
4D BOLD NIfTI files, event timing, and optional confound files from
tools such as fMRIPrep, custom preprocessing workflows, or lab-specific
derivatives.

## Analysis structure

The package uses a stable three-level vocabulary:

- **Level 1** fits the BOLD time series using trial events, signals, and
  nuisance regressors.
- **Level 2** combines run-level contrasts in multi-run FSL analyses.
  SPM and AFNI handle the runs within L1.
- **Level 3** models subject or subject-session contrast maps across
  participants.

Models are organized as sets at each level. For a multi-run FSL
analysis, two compatible models at each of the three levels yield eight
combinations. For native AFNI or SPM, select the L1 and L3 models to
run.

## Input data

The preferred setup uses three explicit long-form data frames:

| Table | One row per | Purpose |
|----|----|----|
| `trial_data` | trial/event | Event onsets, durations, conditions, and trial-level modulators for Level 1 regressors. |
| `run_data` | subject/session/run | Paths to preprocessed run NIfTIs, run-level covariates, TR, confounds, and run exclusion flags. |
| `subject_data` | subject or subject/session | Level 3 covariates. For longitudinal analyses, use one row per `id` + `session`. |

Default identifier columns are `id`, `session`, and `run_number`. If
your project uses different names, pass a named `vm` vector to
[`setup_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md).

For BIDS-like derivatives, helper functions can create starter tables:

``` r

run_df <- generate_run_data_from_bids(
  bids_dir = "/proj/lab/derivatives/fmriprep",
  task_name = "taskname",
  desc = "preproc",
  suffix = "bold"
)

subject_df <- generate_subject_data_from_bids(
  bids_dir = "/proj/lab/derivatives/fmriprep"
)

trial_df <- generate_trial_data_from_bids(
  bids_dir = "/proj/lab/derivatives/fmriprep",
  task_name = "taskname"
)
```

Treat these helpers as a starting point. In most real analyses, you will
still join in task variables, exclusions, and study-specific covariates
before submitting the pipeline.

## Basic usage

Install from GitHub:

``` r

install.packages("remotes")
remotes::install_github("HallquistLab/fmri.pipeline")
```

Create the analysis object from explicit input tables:

`session_data` is not a user-provided input;
[`setup_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md)
derives `gpa$session_data` internally from `run_data` and
`subject_data`.

``` r

library(fmri.pipeline)

gpa <- setup_glm_pipeline(
  analysis_name = "facehouse_longitudinal_demo",
  scheduler = "slurm",
  output_directory = "/proj/lab/glm_outputs",
  subject_data = subject_df,
  run_data = run_df,
  trial_data = trial_df,
  glm_software = c("fsl", "afni"),
  n_expected_runs = 4L,
  l1_models = NULL, l2_models = NULL, l3_models = NULL,
  confound_settings = list(
    motion_params_file = NULL,
    confound_input_file = NULL,
    spike_volumes = NULL
  )
)
gpa$level_backends <- list(l1 = "fsl", l2 = "fsl", l3 = c("fsl", "afni"))
```

Define models interactively:

``` r

gpa <- build_l1_models(gpa)
gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)

export_glm_config(gpa, file = "task_glm_spec.yaml")
```

Or build the model sets from a YAML specification:

``` r

gpa <- build_l1_models(gpa, from_spec_file = "task_glm_spec.yaml")
gpa <- build_l2_models(gpa, from_spec_file = "task_glm_spec.yaml")
gpa <- build_l3_models(gpa, from_spec_file = "task_glm_spec.yaml")
```

Run selected model combinations:

``` r

run_glm_pipeline(
  gpa,
  l1_model_names = "facehouse",
  l2_model_names = "l2_id_session",
  l3_model_names = c("l3_pooled_subject_ev", "l3_3dlmer_session")
)
```

Check status and inspect logs:

``` r

diagnose_pipeline(gpa)
view_log(gpa, level = 1, id = "01")
```

## Estimation by engine

### FSL

`glm_software = "fsl"` is the default. FEAT fits one time-series model
per run. For a multi-run analysis, build an L2 model to combine the run
COPEs and their VARCOPEs. Set `l2_scope = "id_session"` when the group
model needs one estimate per visit. With a single run, L1 feeds L3
directly. Inspect the generated `.fsf` files, the run’s `.feat/stats/`
directory, and the higher-level `.gfeat` outputs.

### AFNI

Select `glm_software = "afni"` in
[`setup_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md)
to use AFNI throughout the analysis. Build L1 and L3 models with the
usual YAML specification, then submit them with:

``` r

run_glm_pipeline(
  gpa,
  glm_software = "afni",
  l1_model_names = "all",
  l2_model_names = "none",
  l3_model_names = "all"
)
```

AFNI estimates one model per subject/session. All retained runs share
task coefficients; each run has its own polynomial baseline and nuisance
coefficients. It always uses the pipeline’s convolved regressors,
preserving event durations, parametric modulators, and signal
normalizations. Runs with `exclude_run = TRUE` and excluded subjects do
not enter the fit. This is a joint time-series fit: AFNI does not run a
separate L2 job or apply L2 model specifications. Use FSL’s explicit L2
path when you need run-level COPE/VARCOPE combination or L2 contrasts.

Configure AFNI through `glm_settings` when creating the pipeline:

``` r

glm_settings = list(afni = list(
  estimation = "REML",        # default; "OLS" uses 3dDeconvolve alone
  polort = "A",               # AFNI automatic run drift order, or an integer >= 0
  mask = NULL,                # optional mask on the input image grid
  censor_column = "afni_censor_file",
  write_residuals = TRUE,
  force_l1_creation = FALSE,
  force_l3_creation = FALSE
))
```

When `run_data$afni_censor_file` is supplied, each path points to a
single-column binary vector (`1 = keep`, `0 = censor`). Censor files and
nuisance files may contain either retained volumes or original volumes
with the specified leading drops. Length mismatches fail setup. Spikes
on retained volumes remain nuisance regressors. The AFNI design drops
spike columns that become empty after censoring. Without an explicit
mask, L1 intersects masks derived from the retained runs. Image geometry
must match; AFNI estimation does not perform registration or intensity
normalization. Supply preprocessed images in a common space and on a
comparable intensity scale.

L3 uses ordinary least squares for the pipeline’s group design and
contrasts, through `3dDeconvolve` with no temporal baseline or
prewhitening. `per_session` analyses stay separate by session;
`pooled_sessions_subject_ev` retains the pipeline’s subject effects. For
longitudinal mixed models, specify `l3_input_mode: 3dlmer`. AFNI-only
analyses use AFNI contrasts as group inputs. To retain a hybrid
workflow, set `execution_backend: afni` and `producer_backend: fsl` on
the L3 model and select FSL at L1/L2. Existing mixed FSL/AFNI `3dlmer`
configurations continue to default to FSL inputs.

Outputs live under `afni_l1/` and `afni_l3/` (mixed models use
`afni_3dlmer/`). Each ordinary GLM directory contains `run_afni_glm.sh`,
`afni_glm.log`, `design.xmat.1D`, `stats.nii.gz`, a mask, and optional
residuals. `contrasts.tsv` maps original contrast names to
`con_####.nii.gz` and `tstat_####.nii.gz`; `regressors.tsv` maps names
to AFNI labels. Completion requires both the success marker and every
expected output. Failed jobs return a nonzero status to the scheduler.
After changing event/model specifications, set
`force_l1_creation = TRUE` to rebuild the cached design, then rerun
setup and estimation.

Resource settings are under `gpa$parallel$afni`: `l1_glm_time`,
`l1_glm_memgb`, `l1_glm_cpus`, `l1_glm_alljobs_time`, and the
corresponding `l3_glm_time`, `l3_glm_memgb`, `l3_glm_cpus`. Existing
`l3_lmer_*` settings continue to control mixed models. Configure AFNI
modules or PATH through `parallel$compute_environment$afni`.

AFNI’s [3dDeconvolve
documentation](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dDeconvolve.html)
and [3dREMLfit
documentation](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dREMLfit.html)
describe the OLS and REML estimators.

### SPM

Set `glm_software = "spm"` and give `glm_settings$spm$spm_path` the path
to your SPM installation. The default first-level design treats each run
as an SPM session and creates one `SPM.mat` per participant, study
visit, and task model. Build the L1 and L3 models and submit with
`l2_model_names = "none"`. Pipeline L3 is the SPM second-level analysis.

SPM L1 contrasts are matched to the fitted design by regressor name. A
missing required regressor or a non-estimable contrast stops the job.
Group inputs are selected using the actual `SPM.xCon` names recorded in
`spm_contrasts.tsv`, so a change in contrast order cannot silently
select a different effect.

SPM setup checks retained scan counts and nuisance rows against each
run’s `run_volumes` and `drop_volumes`. Original images use the
appropriate leading scan offset; already-trimmed images start at scan 1.
Nuisance files may contain original or retained rows. Missing files and
inconsistent lengths are errors. Compressed images use separate
temporary paths for each analysis and run.

A complete SPM fit needs a successful full run and validated, nonempty
model, beta, contrast, statistic, and requested residual files. Setup
records the model specification and input signatures in `spm_setup.rds`;
images use path, size, and modification time, while nuisance files are
content-hashed. Rerun
[`setup_l1_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_l1_models.md)
or
[`setup_l3_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_l3_models.md)
after changing inputs or model settings. Existing fits without the new
validation metadata require setup and estimation once before they can be
reused. Runs with disabled setup, estimation, or contrast stages receive
a stage-completion marker and remain incomplete for group input
selection.

## YAML model specifications

The model builders read the model sections of the YAML file. Set the
analysis-wide `glm_software` and `level_backends` in R as shown above.
The `analysis` block below records those choices alongside the
face/house model definitions:

``` yaml
analysis:
  name: facehouse_longitudinal_demo
  scheduler: slurm
  glm_software:
    - fsl
    - afni
  level_backends:
    l1: fsl
    l2: fsl
    l3:
      - fsl
      - afni
  n_expected_runs: 4

onsets:
  - onset

durations:
  - duration

events:
  stim:
    onset: onset
    duration: duration

signals:
  face:
    event: stim
    trial_subset_expression: object == 'face'
    normalization: none
    value_fixed: 1

  house:
    event: stim
    trial_subset_expression: object == 'house'
    normalization: none
    value_fixed: 1

l1_models:
  facehouse:
    signals:
      - face
      - house
    contrasts:
      diagonal: true

l2_models:
  l2_id_session:
    level: 2
    l2_scope: id_session
    model_formula: ~ stress_label
    reference_level:
      stress_label: stress
    contrasts:
      diagonal: true
      cond_means:
        - stress_label
      pairwise_diffs:
        - stress_label
      overall_response: true

l3_models:
  l3_pooled_subject_ev:
    level: 3
    l3_input_mode: pooled_sessions_subject_ev
    execution_backend: fsl
    producer_backend: fsl
    model_formula: ~ session_label
    reference_level:
      session_label: analgesic
    diagonal: true
    cond_means:
      - session_label
    pairwise_diffs:
      - session_label
    overall_response: true

  l3_3dlmer_session:
    level: 3
    l3_input_mode: 3dlmer
    execution_backend: afni
    producer_backend: fsl
    model_formula: ~ session_label
    random_effects: (1 | Subj)
    reference_level:
      session_label: analgesic
    diagonal: true
    cond_means:
      - session_label
    pairwise_diffs:
      - session_label
    overall_response: true
```

This example demonstrates the current longitudinal pattern: keep runs
nested within subject-session at Level 2 (`l2_scope: id_session`), then
fit either an FSL-native session model with subject explanatory
variables or an AFNI `3dLMEr` model using FSL-produced Level 2 inputs.

## Output layout

Each analysis is written under `output_directory/analysis_name`. The
exact layout depends on selected backends, but typical outputs include:

``` text
analysis_name/
  configuration_files/     # exported specs, caches, and finalized pipeline objects
  scheduler_scripts/       # one batch_<uuid>/ folder per run_glm_pipeline() call
  feat_l1/                 # FSL Level 1 setup and outputs
  feat_l2/                 # FSL Level 2 setup and outputs
  feat_l3/                 # FSL Level 3 setup and outputs
  spm_l1/                  # SPM first-level analyses
  spm_l3/                  # SPM group models
  afni_l1/                 # AFNI subject/session fits
  afni_l3/                 # AFNI ordinary group GLMs
  afni_3dlmer/             # AFNI 3dLMEr data tables, scripts, and outputs
  logs/                    # text and JSON logs
  glm_pipeline.sqlite      # job tracking database
```

The pipeline writes backend preflight information before submission,
then tracks submitted jobs through the SQLite database.
[`diagnose_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/diagnose_pipeline.md)
summarizes the job tree and helps locate failed or incomplete stages.

FSL group models can include an `inference` configuration to run pTFCE
automatically after L3 estimation:

``` r

gpa$l3_models$models$group$inference <- fwe_spec(
  "group_ptfce", targets = list(l3_model = "group"),
  fwe_alpha = 0.05, two_sided = TRUE
)
run_glm_pipeline(gpa, l3_model_names = "group")
```

The L3 model builder offers the same settings, and
[`export_glm_config()`](https://hallquistlab.github.io/fmri.pipeline/reference/export_glm_config.md)
preserves them in YAML. FWE is controlled separately within each
selected map. Use `run_fwe_correction(gpa)` to retry correction without
fitting models again, and `refresh_fwe_correction(gpa)` to collect
scheduled results. See the FWE vignette for custom masks, resource
settings, and the standalone planning workflow.

## Longitudinal workflows

Longitudinal fMRI analyses require special care because runs, sessions,
and subjects represent distinct sources of variation. The current
recommended patterns are:

- For multi-run FSL, use `l2_scope = "id_session"` to keep one L2
  estimate per visit. Native AFNI L1 already produces one fit per visit.
  SPM preserves that distinction with its default
  `l1_session_mode = "separate"`.
- Use `l3_input_mode = "pooled_sessions_subject_ev"` for balanced
  categorical visit comparisons that should remain inside FSL.
- Use `l3_input_mode = "3dlmer"` for irregular visits, incomplete
  follow-up, continuous time effects, treatment-by-time interactions, or
  random-effects structures that are better handled by AFNI `3dLMEr`.
- In FSL, reserve `l2_scope = "id"` for cases where sessions are
  essentially pseudo-runs and longitudinal change is not the target.

The [longitudinal
vignette](https://hallquistlab.github.io/fmri.pipeline/articles/longitudinal.html)
walks through a face/house study with repeated visits.

## Documentation

The pkgdown site is the primary documentation entry point:

- Package website: <https://hallquistlab.github.io/fmri.pipeline/>
- Quick start:
  <https://hallquistlab.github.io/fmri.pipeline/articles/quickstart.html>
- Choosing a GLM engine:
  <https://hallquistlab.github.io/fmri.pipeline/articles/glm-engines.html>
- Design concepts:
  <https://hallquistlab.github.io/fmri.pipeline/articles/design.html>
- Longitudinal analyses:
  <https://hallquistlab.github.io/fmri.pipeline/articles/longitudinal.html>
- Pipeline diagnosis:
  <https://hallquistlab.github.io/fmri.pipeline/articles/diagnosis.html>
- Reference index:
  <https://hallquistlab.github.io/fmri.pipeline/reference/index.html>

MEDuSA deconvolution, event-locked alignment, and multilevel trajectory
modeling (`mixed_by`, `jlmer`, `meta_mixed_by`) now live in the
standalone [`medusa`](https://github.com/HallquistLab/medusa) package.
Legacy pre-0.4-1 MEDuSA entry points in `fmri.pipeline` emit deprecation
warnings and forward to `medusa` when installed.

## End-to-end tests

The clock-task suite runs FSL, SPM and AFNI from trial data through
group contrasts on Slurm. It includes a four-participant smoke profile
and the complete 80-run fixture. [Run the clock GLM
suite](https://hallquistlab.github.io/fmri.pipeline/tests/integration/clock-glm/README.md).

## Development status

`fmri.pipeline` is under active development. For new projects, prefer
YAML model specifications, explicit
`subject_data`/`run_data`/`trial_data` tables, and small test runs
before launching the full model grid.

Issues and feature requests are tracked at
<https://github.com/HallquistLab/fmri.pipeline/issues>.
