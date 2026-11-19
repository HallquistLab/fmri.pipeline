# This function generates the inputs for an FSL level 2 analysis, where multiple runs for a subject are combined using fixed effects estimation.

This function generates the inputs for an FSL level 2 analysis, where
multiple runs for a subject are combined using fixed effects estimation.

## Usage

``` r
setup_l2_models(
  gpa,
  l1_model_names = NULL,
  l2_model_names = NULL,
  backend = NULL
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing the analysis
  specification

- l1_model_names:

  a subset of L1 models to be passed to L2 by this function. If not
  specified, all models in gpa\$l1_models will be included

- l2_model_names:

  a subset of L2 models to be setup by this function. If not specified,
  all models in gpa\$l2_models will be included

- backend:

  optional backend filter (e.g., "fsl"). If supplied, only those
  backends will be processed.

## Details

This function will setup FSL level 2 (subject) .fsf files for all
combinations of `l2_model_names` and `l1_model_names`. Native AFNI and
SPM analyses combine runs within L1 and proceed directly to group setup.

## Author

Michael Hallquist
