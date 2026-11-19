# Choosing a GLM engine: FSL, AFNI, or SPM

Suppose each participant completes two runs at each visit and you want
the group mean response to a cue. All three engines can estimate it.
Your choice determines how the two runs contribute to each person’s cue
contrast and what information reaches the group model. The pipeline
calls the time-series fit **L1** and the group fit **L3** in all three
cases.

| Choose this route when you need… | How the runs enter the fit | What reaches the group model |
|----|----|----|
| **FSL** for run-wise estimates, run-level predictors or precision weighting | FEAT fits each run. A separate L2 model combines its COPEs (contrast estimates) using their VARCOPEs (estimated variances) in the default fixed-effects fit. | FLAME can use the combined effect and its variance. A single-run analysis takes its inputs directly from L1. |
| **AFNI** for a shared task effect across runs, run-wise nuisance terms or explicit volume censoring | One L1 model per participant and visit shares the cue coefficient across retained runs. Each run gets its own baseline and nuisance coefficients. | An ordinary group GLM fits the contrast estimates by OLS. For repeated visits, `3dLMEr` fits fixed and random effects to visit-specific estimates. |
| **SPM** for SPM condition and modulator designs or contrasts across run-specific coefficients | By default each run is an SPM session in one first-level model per participant and visit. The default cue contrast averages the matching cue coefficients across runs. | An SPM second-level analysis (pipeline L3) uses first-level contrast images. The wrapper generates a one-sample test or a multiple-regression design. |

## Run the same task with one engine

The [quick
start](https://hallquistlab.github.io/fmri.pipeline/articles/quickstart.md)
defines the three input tables and a `clock_glm_spec.yaml` with models
named `task`, `run_average`, and `group_mean`. The example below uses
that two-run design with a TR of 0.9 seconds. Set `engine` to one value
and supply paths and nuisance settings for your study. If your cluster
uses modules, configure the suite in `parallel$compute_environment` so
the scheduled job can find the same installation.

``` r

library(fmri.pipeline)

engine <- "afni"  # "fsl", "afni", or "spm"
settings <- if (engine == "spm") {
  list(spm = list(spm_path = "/path/to/spm12", matlab_cmd = "matlab"))
} else {
  "default"
}

gpa <- setup_glm_pipeline(
  analysis_name = paste0("clock_", engine),
  scheduler = "slurm",
  trial_data = trial_df, run_data = run_df, subject_data = subject_df,
  tr = 0.9, n_expected_runs = 2L,
  glm_software = engine, glm_settings = settings,
  l1_models = NULL, l2_models = NULL, l3_models = NULL
)

spec <- "clock_glm_spec.yaml"
gpa <- build_l1_models(gpa, from_spec_file = spec)
if (engine == "fsl") {
  gpa <- build_l2_models(gpa, from_spec_file = spec)
}
gpa <- build_l3_models(gpa, from_spec_file = spec)

run_glm_pipeline(
  gpa, l1_model_names = "task",
  l2_model_names = if (engine == "fsl") "run_average" else "none",
  l3_model_names = "group_mean"
)
```

Here `group_mean` uses `l3_input_mode: per_session`, so each visit gets
its own group analysis. The SPM and AFNI routes take their group inputs
directly from L1. When building the SPM group model without an L2
specification, the builder can warn that L2/L3 compatibility has not
been checked. This example needs no projected run-level predictors.

## Coming from FSL

A pipeline *signal* supplies a task EV, and an L1 model selects signals
and contrasts. The default `use_preconvolve = TRUE` passes sampled
regressors to FEAT. Set it to `FALSE` at setup to have FEAT convolve
onset/duration/value files. Check the EVs and contrasts in the generated
`.fsf` and the fitted `design.mat` and `design.con` under `feat_l1/`. A
completed run has a familiar `.feat/stats/` directory with
`cope1.nii.gz` and `varcope1.nii.gz` for the first contrast. [FSL FEAT
guide](https://fsl.fmrib.ox.ac.uk/fsl/docs/task_fmri/feat/index.html)

If the second run was performed under a different condition, put that
run-level predictor in `run_data` and include it in the L2 formula. Use
`l2_scope: id_session` to keep one estimate per visit. Look in
`feat_l2/` and `feat_l3/` for the higher-level `.gfeat` results.
Choosing AFNI for L3 sends the FSL contrast estimates forward without
their VARCOPE weights. For automatic pTFCE correction after L3
estimation, use the FSL group route described in the [FWE
guide](https://hallquistlab.github.io/fmri.pipeline/articles/fwe.md).

## Coming from AFNI

Start with the preprocessed runs you would normally pass to the
regression step. Supply images on a common grid with comparable
intensity scaling. Define timing and amplitudes in `trial_data` and the
L1 model. The wrapper writes convolved task columns as `-stim_file`
inputs to `3dDeconvolve`. It records run boundaries and builds
run-specific polynomial baselines (`polort = "A"` by default). For
volume censoring, point `run_data$afni_censor_file` to a binary file per
run, with 1 to keep and 0 to censor. Spike regressors in the confounds
model are a separate choice.

By default `3dDeconvolve` builds the design and `3dREMLfit` estimates
the time series with temporal prewhitening. Set
`glm_settings = list(afni = list(estimation = "OLS"))` to fit by
ordinary least squares in `3dDeconvolve`. Inspect `run_afni_glm.sh` and
`design.xmat.1D` under `afni_l1/`. Use `contrasts.tsv` to find a named
GLT’s coefficient and t-statistic images. The coefficient image feeds
L3. Ordinary group fits live in `afni_l3/` and mixed models in
`afni_3dlmer/`. [AFNI 3dREMLfit
help](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dREMLfit.html)

## Coming from SPM

A *session* column in `run_data` identifies a study visit. With the
defaults `l1_session_mode = "separate"` and `concatenate_runs = FALSE`,
two visits with two runs each produce two `SPM.mat` files per
participant, each with two SPM sessions. Conditions and parametric
modulators are generated from pipeline events and signals. The generated
condition batch enables SPM orthogonalization, so modulator order can
change what a coefficient measures. Inspect the fitted columns in
`SPM.xX.name` and contrast weights in `SPM.xCon`. Pipeline defaults
include a 100-second high-pass cutoff and `cvi = "AR(1)"`. [SPM model
specification
manual](https://www.fil.ion.ucl.ac.uk/spm/docs/manual/fmri_spec/fmri_spec/)

Set options inside `glm_settings$spm`. Two choices change different
aspects of the first-level design:

- `concatenate_runs = TRUE` puts runs in one SPM session, sharing task
  columns, and invokes `spm_fmri_concatenate` with the run lengths.
- `l1_session_mode = "pooled"` places runs from multiple study visits in
  one participant model. Keep the default for a separate model at each
  visit.

To model a run-varying predictor, build its L2 specification **before L1
setup** and select it with `l2_projection_model`. SPM can project its
varying terms and requested task interactions into the first-level
analysis. A predictor constant within an SPM session is aliased with the
session constant as a main-effect column. Check the generated design and
contrasts for the run comparison you intend.

Find `SPM.mat`, beta images and `spm_contrasts.tsv` under `spm_l1/`. The
group setup matches contrast names to fitted `con_*.nii` files. Leave
`producer_backend` unset for native SPM L3. The submission interface
accepts explicit producer overrides only for FSL and AFNI. Inspect the
group design under `spm_l3/`. The wrapper implements the
`flexible_factorial` setting through multiple regression for
non-intercept designs, preserving the pipeline design matrix and
contrasts.

## Use FSL estimates in an AFNI group model

A vector such as `glm_software = c("fsl", "afni")` enables both engines.
For a longitudinal analysis, set that vector at setup. Specify which
engine handles each stage on the resulting object:

``` r

gpa$level_backends <- list(l1 = "fsl", l2 = "fsl", l3 = "afni")
```

Build the FSL L2 model with `l2_scope: id_session`. Add the following L3
model to the YAML specification, with a visit covariate named
`session_label` in your input data:

``` yaml
l3_models:
  visit_effect:
    level: 3
    model_formula: ~ session_label
    l3_input_mode: 3dlmer
    random_effects: "(1 | Subj)"
    execution_backend: afni
    producer_backend: fsl
    contrasts:
      diagonal: true
```

Here `producer_backend` identifies the suite that supplies the contrast
maps. `execution_backend` identifies the suite fitting the group model.
For an AFNI-only mixed model, choose AFNI at setup, use
`producer_backend: afni` and skip L2. AFNI group models support one FSL
or AFNI producer per model. The native SPM route uses its own
first-level contrasts. See the [longitudinal
examples](https://hallquistlab.github.io/fmri.pipeline/articles/longitudinal.md)
for visit comparisons and the [design
guide](https://hallquistlab.github.io/fmri.pipeline/articles/design.html#choosing-model-backends)
for run-time overrides.

At submission, check the backend preflight report in the log for the
execution engine and the level supplying each group contrast. When
comparing fits, check the retained volumes, mask, task columns, contrast
weights and scaling for one participant first. For the cue example, a
shared AFNI coefficient, an SPM contrast across run-specific
coefficients and a precision-weighted FSL estimate can differ even with
identical event onsets.
