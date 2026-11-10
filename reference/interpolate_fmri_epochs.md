# Interpolate an aligned fmri_ts object onto a consistent time grid with respect to a target event

Interpolate an aligned fmri_ts object onto a consistent time grid with
respect to a target event

## Usage

``` r
interpolate_fmri_epochs(
  a_obj,
  evt_time = "evt_time",
  time_before = -3,
  time_after = 3,
  output_resolution = 1,
  group_by = NULL,
  time_audit = FALSE,
  max_gap = NULL,
  readout_kernel = "linear",
  readout_width = NULL
)
```

## Arguments

- a_obj:

  an fmri_ts object that has data aligned to an event of interest

- evt_time:

  column name in `a_obj` event data to which time series should be
  aligned and interpolated

- time_before:

  The earliest time point (in seconds) relative to the event to be
  output in interpolation

- time_after:

  The latest time point (in seconds) relative to the event to be output
  in interpolation

- output_resolution:

  The timestep (in seconds) used for interpolation

- group_by:

  a character vector of keying variables used for aggregation of data
  prior to interpolation

- time_audit:

  If TRUE, additional columns will be added to the output showing how
  alignment is calculated vis-a-vis event timing

- max_gap:

  Maximum source bracket duration (default acquisition TR).

## Author

Michael Hallquist
