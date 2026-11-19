# Submit AFNI GLMs to the configured scheduler

Submit AFNI GLMs to the configured scheduler

## Usage

``` r
run_afni_sepjobs(
  gpa,
  level = 1L,
  model_names = NULL,
  rerun = FALSE,
  wait_for = ""
)
```

## Arguments

- gpa:

  A glm_pipeline_arguments object.

- level:

  Analysis level (1 or 3). Runs are combined within L1.

- model_names:

  Optional model subset.

- rerun:

  Whether to rerun complete models.

- wait_for:

  Optional parent job IDs.

## Value

Submitted job IDs.
