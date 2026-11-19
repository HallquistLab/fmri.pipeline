# Set up an AFNI first-level GLM

Runs within a subject/session share task coefficients and have
independent polynomial baselines and nuisance coefficients. The
pipeline's convolved design is passed to 3dDeconvolve; 3dREMLfit is the
default estimator.

## Usage

``` r
afni_l1_model(
  id,
  session,
  l1_confound_files = NULL,
  d_obj,
  gpa,
  model_name,
  run_nifti = NULL,
  censor_files = NULL
)
```

## Arguments

- id, session:

  Subject and session identifiers.

- l1_confound_files:

  Optional nuisance files, one per BDM run.

- d_obj:

  A build_design_matrix (bdm) object.

- gpa:

  A glm_pipeline_arguments object.

- model_name:

  Name of the L1 model.

- run_nifti:

  Optional input images when absent from d_obj.

- censor_files:

  Optional binary AFNI censor files, one per BDM run.

## Value

A one-row AFNI setup/status data.frame.
