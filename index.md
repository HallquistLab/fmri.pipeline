# fmri.pipeline

`fmri.pipeline` is an R package for setting up, submitting, and tracking
hierarchical fMRI GLM analyses on high-performance computing systems. It
is designed for studies where the scientifically interesting model is
larger than a single FEAT or SPM job: many subjects, many runs, multiple
model variants, longitudinal visits, and downstream group models that
need to be reproducible months or years later.

The package is built around a single analysis object, usually called
`gpa`, that stores the input tables, model specifications, output
locations, backend choices, scheduler settings, and job-tracking state
for one analysis.

Current package version: `0.4` (May 2026).

## Why this package exists

Large fMRI GLM projects have a familiar failure mode: the actual
analysis is a constellation of design files, shell scripts, scheduler
jobs, and partially remembered choices. `fmri.pipeline` tries to make
that constellation explicit. The package helps you:

- define run-level, subject/session-level, and group-level models in one
  place;
- generate backend-specific analysis files for FSL, AFNI, and selected
  SPM workflows;
- submit staged analyses through Slurm, TORQUE, or local execution;
- keep track of job dependencies, logs, and analysis status;
- reuse the same input tables and YAML model specification across model
  revisions;
- support longitudinal analyses where runs are nested within sessions
  and sessions are nested within subjects.

The most mature end-to-end backend is FSL FEAT/FLAME. Current
development also supports backend selection and hybrid workflows,
especially FSL at Levels 1-2 followed by AFNI `3dLMEr` at Level 3 for
longitudinal mixed-effects models.

`fmri.pipeline` is not a preprocessing pipeline. It expects preprocessed
4D BOLD NIfTI files, event timing, and optional confound files from
tools such as fMRIPrep, custom preprocessing workflows, or lab-specific
derivatives.

## Analysis structure

The package uses a stable three-level vocabulary:

- **Level 1**: run-level time-series GLM. Events and parametric
  modulators are converted into design matrices or timing files and
  estimated separately for each run.
- **Level 2**: within-subject or within-subject-session run combination.
  Multiple Level 1 COPEs are combined into subject-level or
  session-level estimates.
- **Level 3**: sample-level, group-level, or longitudinal model. Subject
  or subject-session maps are modeled across participants.

Model definitions are organized as sets at each level. If you define two
L1 models, two L2 models, and two L3 models, the pipeline can run the
full grid of eight model combinations. This makes systematic model
comparison possible, but it also means you should start small, verify
the analysis path, and then expand.

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
  confound_settings = list(
    motion_params_file = NULL,
    confound_input_file = NULL,
    spike_volumes = NULL
  )
)
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

## YAML model specifications

YAML specifications are the recommended way to make analysis decisions
reviewable and repeatable. A compact face/house example looks like this:

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
  afni_3dlmer/             # AFNI 3dLMEr data tables, scripts, and outputs
  logs/                    # text and JSON logs
  glm_pipeline.sqlite      # job tracking database
```

The pipeline writes backend preflight information before submission,
then tracks submitted jobs through the SQLite database.
[`diagnose_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/diagnose_pipeline.md)
summarizes the job tree and helps locate failed or incomplete stages.

## Longitudinal workflows

Longitudinal fMRI analyses require special care because runs, sessions,
and subjects represent distinct sources of variation. The current
recommended patterns are:

- Use `l2_scope = "id_session"` for most longitudinal studies. This
  creates one Level 2 output per subject-session.
- Use `l3_input_mode = "pooled_sessions_subject_ev"` for balanced
  categorical visit comparisons that should remain inside FSL.
- Use `l3_input_mode = "3dlmer"` for irregular visits, incomplete
  follow-up, continuous time effects, treatment-by-time interactions, or
  random-effects structures that are better handled by AFNI `3dLMEr`.
- Reserve `l2_scope = "id"` for cases where sessions are essentially
  pseudo-runs and longitudinal change is not the target.

The repository includes a worked face/house demonstration in
`longitudinal_fmri_demo_jun2026/`, including a walkthrough script and
YAML model file.

## Documentation

The pkgdown site is the primary documentation entry point:

- Package website: <https://hallquistlab.github.io/fmri.pipeline/>
- Quick start:
  <https://hallquistlab.github.io/fmri.pipeline/articles/quickstart.html>
- Design concepts:
  <https://hallquistlab.github.io/fmri.pipeline/articles/design.html>
- Longitudinal analyses:
  <https://hallquistlab.github.io/fmri.pipeline/articles/longitudinal.html>
- Pipeline diagnosis:
  <https://hallquistlab.github.io/fmri.pipeline/articles/diagnosis.html>
- Reference index:
  <https://hallquistlab.github.io/fmri.pipeline/reference/index.html>

MEDuSA documentation now follows a three-part learning sequence:

1.  [Conceptual
    foundations](https://hallquistlab.github.io/fmri.pipeline/vignettes/medusa-concepts.Rmd):
    convolution, regularization, HRF learning, event mixing, units,
    duration, and precision-weight propagation, with an annotated
    reading guide.
2.  [Worked
    workflows](https://hallquistlab.github.io/fmri.pipeline/vignettes/medusa-workflows.Rmd):
    executable simulation and recipes for NIfTI inputs, file alignment,
    and downstream precision-weighted / ZCA-whitened multilevel models
    (`mixed_by`).
3.  [Technical companion (PRs
    \#42–#55+)](https://hallquistlab.github.io/fmri.pipeline/medusa_deconvolution_review_prs_42_53.qmd):
    equations, simulation demonstrations, implementation details, and
    validation limits.

The first two are standard package vignettes
([`vignette("medusa-concepts")`](https://hallquistlab.github.io/fmri.pipeline/articles/medusa-concepts.md)
and
[`vignette("medusa-workflows")`](https://hallquistlab.github.io/fmri.pipeline/articles/medusa-workflows.md)
after building vignettes). The technical companion is a standalone
Quarto report at the repository root, excluded from package and pkgdown
builds; its historical basename is retained. See [rendering
instructions](https://hallquistlab.github.io/fmri.pipeline/vignettes/medusa-docs.md)
for producing all three HTML documents from the same checkout.

## Development status

`fmri.pipeline` is under active development. FSL-centered analyses are
the most established path. Backend resolution, SPM support, AFNI Level 3
workflows, longitudinal model specification, and diagnostic tooling are
current areas of active work. For new projects, prefer YAML model
specifications, explicit `subject_data`/`run_data`/`trial_data` tables,
and small test runs before launching the full model grid.

Issues and feature requests are tracked at
<https://github.com/HallquistLab/fmri.pipeline/issues>.
