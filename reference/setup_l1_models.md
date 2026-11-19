# Prepare first-level designs and files for the selected GLM engines

Prepare first-level designs and files for the selected GLM engines

## Usage

``` r
setup_l1_models(gpa, l1_model_names = NULL)
```

## Arguments

- gpa:

  A `glm_pipeline_arguments` object from
  [`setup_glm_pipeline`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md).

- l1_model_names:

  Names of the L1 models to prepare. The default selects all L1 models.

## Value

The updated `gpa`. Its `l1_model_setup` list holds run metadata and
setup tables named `fsl`, `afni`, and `spm`. Inspect the table for the
engine you selected.

## Details

FSL generates a FEAT design per run. AFNI prepares one joint time-series
fit per participant and visit. SPM prepares a first-level model per
participant and visit by default, with separate SPM sessions for the
retained runs. Set `gpa$parallel$l1_setup_cores` to control parallel
setup across participants.

## Author

Michael Hallquist
