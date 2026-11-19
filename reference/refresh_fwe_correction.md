# Refresh configured group inference results

Refresh configured group inference results

## Usage

``` r
refresh_fwe_correction(
  gpa,
  correction_names = names(gpa$fwe_runs),
  require_complete = FALSE
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` with `fwe_runs`.

- correction_names:

  corrections to refresh; defaults to all recorded runs.

- require_complete:

  error after persisting state if any correction is incomplete.

## Value

Updated `gpa`, including result manifests and separate FWE status.
