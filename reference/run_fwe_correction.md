# Run configured group inference without refitting GLMs

Runs the declarative corrections stored in L3 models' `inference` fields
and in `gpa$fwe_correction`. Scheduler runs return after submission; use
[`refresh_fwe_correction()`](https://hallquistlab.github.io/fmri.pipeline/reference/refresh_fwe_correction.md)
to collect them.

## Usage

``` r
run_fwe_correction(
  gpa,
  l3_model_names = names(gpa$l3_models$models),
  force = FALSE,
  dry_run = FALSE,
  ...
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object.

- l3_model_names:

  models whose configured corrections should run; defaults to all L3
  models.

- force:

  rerun corrections even when verified outputs exist.

- dry_run:

  preview plans and commands without writing or submitting.

- ...:

  additional arguments to
  [`run_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_fwe_plan.md).

## Value

Updated `gpa`, with `fwe_runs` and `fwe_status`.
