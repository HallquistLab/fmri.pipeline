# Prepare group models from completed lower-level contrast maps

Prepare group models from completed lower-level contrast maps

## Usage

``` r
setup_l3_models(
  gpa,
  l3_model_names = NULL,
  l2_model_names = NULL,
  l1_model_names = NULL,
  backend = NULL
)
```

## Arguments

- gpa:

  A `glm_pipeline_arguments` object containing the model specifications.

- l3_model_names:

  Names of the group models to prepare. Defaults to all defined L3
  models.

- l2_model_names:

  Names of upstream FSL run-combination models when inputs come from L2.
  Defaults to all defined L2 models in that case.

- l1_model_names:

  Names of the upstream task models. Defaults to all defined L1 models.

- backend:

  Optional filter on the resolved group engines: `"fsl"`, `"afni"`, or
  `"spm"`.

## Value

The updated `gpa`. Group setup records for the selected engines are
stored under `l3_model_setup`.

## Details

FSL uses FEAT contrasts from L1 for single-run analyses and from L2 for
multi-run analyses. SPM uses named contrast images from its first-level
models. AFNI ordinary group GLMs and longitudinal `3dLMEr` models accept
AFNI subject-session contrast estimates or FSL COPEs, according to the
L3 model's producer setting.

## See also

[`setup_l1_models`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_l1_models.md),
[`refresh_glm_status`](https://hallquistlab.github.io/fmri.pipeline/reference/refresh_glm_status.md)

## Author

Michael Hallquist
