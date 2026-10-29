# Identify the last valid volume acquired in a given run. Subjects often exhibit head movement after run ends (MATLAB closes), but scan hasn't stopped This occurs because the MB raw transfer of the prior run is occurring, but does not finish before the current run Thus, truncate mr files to be 12 seconds after final feedback presentation, which is how the paradigm timing files are setup

Identify the last valid volume acquired in a given run. Subjects often
exhibit head movement after run ends (MATLAB closes), but scan hasn't
stopped This occurs because the MB raw transfer of the prior run is
occurring, but does not finish before the current run Thus, truncate mr
files to be 12 seconds after final feedback presentation, which is how
the paradigm timing files are setup

## Usage

``` r
truncate_runs(
  mr_df,
  gpa = NULL,
  subj_outdir = NULL,
  truncation_data = NULL,
  lg = NULL,
  force = FALSE
)
```

## Arguments

- force:

  whether to replace an existing truncated image when its source has
  changed.
