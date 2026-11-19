# Main worker function for setting up an analysis pipeline

Main worker function for setting up an analysis pipeline

## Usage

``` r
setup_glm_pipeline(
  analysis_name = "glm_analysis",
  scheduler = "slurm",
  output_directory = file.path(getwd(), analysis_name),
  subject_data = NULL,
  run_data = NULL,
  trial_data = NULL,
  ppi_data = NULL,
  group_output_directory = "default",
  output_locations = "default",
  vm = c(id = "id", session = "session", run_number = "run_number", block_number =
    "block_number", trial = "trial", run_trial = "run_trial", mr_dir = "mr_dir",
    run_nifti = "run_nifti", exclude_subject = "exclude_subject"),
  bad_ids = NULL,
  tr = NULL,
  fmri_file_regex = ".*\\.nii(\\.gz)?",
  fmri_path_regex = NULL,
  run_number_regex = ".*run-*([0-9]+).*",
  drop_volumes = 0L,
  l1_models = "prompt",
  l2_models = "prompt",
  l3_models = "prompt",
  glm_software = "fsl",
  n_expected_runs = 1L,
  use_preconvolve = TRUE,
  glm_settings = "default",
  confound_settings = list(motion_params_file = NULL, motion_params_colnames = NULL,
    confound_input_file = NULL, confound_input_colnames = NULL, l1_confound_regressors =
    NULL, exclude_run = NULL, truncate_run = NULL, exclude_subject = NULL, spike_volumes
    = NULL),
  parallel = list(l1_setup_cores = 4L, pipeline_cores = "default"),
  additional = list(feat_l1_args = list(z_thresh = 1.96, prob_thresh = 0.05)),
  lgr_threshold = "info"
)
```

## Arguments

- analysis_name:

  A character string providing a useful name for identifying this
  analysis. Practically, this influences the top-level folder name of
  the group analysis outputs, as well as the name of .RData objects
  saved by this function to the output directory for the analysis.

- scheduler:

  Which HPC scheduler system should be used for queueing jobs. Options
  are 'slurm', 'torque', or 'local'.

- output_directory:

  The output directory for all configuration and job submission files
  for this analysis. Default is a subfolder called "glm_out" in the
  current working directory.

- subject_data:

  A data.frame containing Level 3 covariates and subject/session
  metadata. For cross-sectional analyses this is typically one row per
  subject. For longitudinal analyses, this must contain one row per
  `id`/`session`; person-level covariates such as sex, baseline age, or
  treatment arm should be repeated across that subject's sessions.
  Columns from `subject_data` can be used as covariates in group (aka
  'level 3') analyses. If `NULL`, then this will be distilled from
  `trial_data` by looking for variables that vary at the same level as
  `vm["id"]` and `vm["session"]`. Do not provide `session_data`;
  `setup_glm_pipeline()` derives `gpa$session_data` internally from
  `run_data` and `subject_data`.

- run_data:

  A data.frame containing all run-level data such as run condition or
  run number. Columns from `run_data` can be used as covariates in
  subject (aka 'level 2') analyses. If `NULL`, this will be distilled
  from `trial_data` by looking for variables that vary at the same level
  as `vm["run_number"]`. Columns that are constant within each
  `id`/`session` are copied into the internally derived
  `gpa$session_data`.

- trial_data:

  A data.frame containing trial-level statistics for all subjects. Data
  should be stacked in long form such that each row represents a single
  trial for a given subject and the number of total rows is subjects x
  trials. If you wish, you can pass a single trial-level data frame that
  also contains all run-level and subject-level covariates (i.e., a
  combined long format, where variables at different levels are all
  included as columns). In this case, `setup_glm_pipeline` will detect
  which variables occur at each level and parse these accordingly into
  `subject_data` and `run_data`; `session_data` is still constructed
  internally.

- ppi_data:

  An optional \`data.frame\` containing physiological signals (e.g.,
  deconvolved time series in seed ROIs) that can be used in setting up
  the level 1 models for a PPI analysis. Specifically, if this
  \`data.frame\` is provided, you will be able to create interaction
  regressors for any convolved event-related signal and the columns of
  \`ppi_data\`. This data.frame must include subject id, session id (if
  relevant), and run_number as columns. All additional columns should be
  eligible physiological signals for PPI analysis.

- group_output_directory:

  Group-level output directory policy.

- output_locations:

  Output location template settings.

- vm:

  A named character vector containing key identifying columns in
  `subject_data` and `trial_data`. Minimally, this vector should contain
  the elements 'id'

- bad_ids:

  An optional vector of ids in `subject_data` and `trial_data` that
  should be excluded from analysis.

- tr:

  The repetition time of the scanning sequence in seconds. Used for
  setting up design matrices. If this is NULL, the function will look
  for a `tr` field in the `run_data` object, which specifies TR at run
  level (e.g., if varies).

- fmri_file_regex:

  A character string containing a Perl-compatible regular expression for
  the subfolder and filename within the `mr_dir` field in
  `subject_data`.

- fmri_path_regex:

  Optional regular expression for extracting fMRI paths.

- run_number_regex:

  Optional regular expression for extracting run numbers.

- drop_volumes:

  The number of volumes to drop from the fMRI data and convolved
  regressors prior to analysis. Default is 0.

- l1_models:

  An `l1_model_set` object containing all level 1 (run) models to be
  included in GLM pipeline. If "prompt" is passed, `setup_glm_pipeline`
  will call `build_l1_models` to setup l1 models interactively.
  Optionally, this argument can be `NULL` if you want to setup l1 models
  later, though the resulting object will not be functional within the
  pipeline until l1 models are provided.

- l2_models:

  An `hi_model_set` object containing all level 2 (subject) models to be
  included in GLM pipeline. If "prompt" is passed, `setup_glm_pipeline`
  will call `build_l2_models` to setup l2 models interactively. Use
  `NULL` to defer setup or omit L2 for a native AFNI or SPM fit. SPM can
  optionally use an L2 specification to project run-level predictors
  into L1.

- l3_models:

  An `hi_model_set` object containing all level 3 (sample) models to be
  included in GLM pipeline. If "prompt" is passed, `setup_glm_pipeline`
  will call `build_l3_models` to setup l3 models interactively.
  Optionally, this argument can be `NULL` if you want to setup l3 models
  later, though the resulting object will not be functional within the
  pipeline until l3 models are provided.

- glm_software:

  Default GLM engine(s), case insensitive. Use `"fsl"` (default) for
  run-wise FEAT and an explicit multi-run L2, `"afni"` for shared task
  effects across runs within a visit, or `"spm"` for an SPM first-level
  model and group analysis. A vector enables multiple engines. Set
  `gpa$level_backends` to assign them to levels.

- n_expected_runs:

  Number of expected runs per subject. Used to determine 2- versus
  3-level analysis (for FSL), and for providing feedback about subjects
  who have unexpected numbers of runs.

- use_preconvolve:

  For FSL, pass pipeline-convolved EVs (`TRUE`, the default) or let FEAT
  convolve onset/duration/value files (`FALSE`). AFNI uses
  pipeline-convolved task regressors. SPM builds conditions and
  parametric modulators from event timing.

- glm_settings:

  A named list of engine settings, e.g.,
  `list(afni = list(estimation = "OLS"))` or
  `list(spm = list(spm_path = "/path/to/spm12"))`. The default AFNI
  estimator is REML. Configure SPM with your local installation path.

- confound_settings:

  Confound-file and censoring settings.

- parallel:

  Parallelization and scheduler settings.

- additional:

  A list of additional metadata that will be added to the `glm.pipeline`
  object returned by the function. This can be useful if there are other
  identifiers that you want for long-term storage or off-shoot
  functions.

- lgr_threshold:

  The logging threshold used to determine whether to output messages of
  different severity to the screen and to log files. Default is "info",
  which produces all messages, warnings, and errors, but not debug or
  trace statements. To output only concerning errors, change to "error".
  See: <https://s-fleck.github.io/lgr/articles/lgr.html> for details.

## See also

[`run_glm_pipeline`](https://hallquistlab.github.io/fmri.pipeline/reference/run_glm_pipeline.md),
[`vignette("glm-engines", package = "fmri.pipeline")`](https://hallquistlab.github.io/fmri.pipeline/articles/glm-engines.md)
