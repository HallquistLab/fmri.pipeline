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

  a character vector of level 2 model names (specified during
  build_l2_models) that should be executed

- l3_model_names:

  a character vector of level 3 model names (specified during
  build_l3_models) that should be executed

- glm_software:

  which GLM software should be used for model estimation (fsl, spm, or
  afni)

- level_backends:

  optional per-level backend override list keyed by \`l1\`, \`l2\`,
  and/or \`l3\`

- backend_overrides:

  optional model-specific backend override list

- run_fwe:

  run configured FWE correction after the selected L3 models finish.
