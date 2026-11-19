# primary function for running a GLM analysis pipeline

primary function for running a GLM analysis pipeline

## Usage

``` r
run_glm_pipeline(
  gpa,
  l1_model_names = "prompt",
  l2_model_names = "prompt",
  l3_model_names = "prompt",
  glm_software = NULL,
  level_backends = NULL,
  backend_overrides = NULL,
  run_fwe = TRUE
)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object containing a model specification
  (created by setup_glm_pipeline)

- l1_model_names:

  a character vector of level 1 model names (specified during
  build_l1_models) that should be executed

- l2_model_names:

  Names of standalone L2 models to run for multi-run FSL. Use `"none"`
  for a native AFNI or SPM L1-to-L3 analysis.

- l3_model_names:

  a character vector of level 3 model names (specified during
  build_l3_models) that should be executed

- glm_software:

  Optional global engine override for this submission. A non-NULL value
  resets the level defaults to the supplied engine(s): `"fsl"`,
  `"afni"`, or `"spm"`.

- level_backends:

  Optional list overriding specific level defaults for this run, e.g.,
  `list(l1 = "fsl", l2 = "fsl", l3 = "afni")`. Model-specific settings
  take precedence.

- backend_overrides:

  Optional overrides keyed by level and model name. For example,
  `list(l3 = list(execution = c(group = "afni"), producer = c(group = "fsl")))`
  fits the named group model with AFNI using FSL contrast maps. These
  overrides take precedence over the model's `execution_backend` and
  `producer_backend`.

- run_fwe:

  run configured FWE correction after the selected L3 models finish.

## See also

[`setup_glm_pipeline`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md),
[`vignette("glm-engines", package = "fmri.pipeline")`](https://hallquistlab.github.io/fmri.pipeline/articles/glm-engines.md)
