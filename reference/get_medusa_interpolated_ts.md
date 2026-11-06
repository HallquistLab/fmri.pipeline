# front-end function for taking a list of windowed time series by mask value, interpolating them onto a time grid, and (optionally) averaging across voxels/units within a value to derive the mean interpolated time series

front-end function for taking a list of windowed time series by mask
value, interpolating them onto a time grid, and (optionally) averaging
across voxels/units within a value to derive the mean interpolated time
series

## Usage

``` r
get_medusa_interpolated_ts(
  fmri_obj,
  event = NULL,
  time_before = -3,
  time_after = 3,
  collide_before = NULL,
  collide_after = NULL,
  pad_before = -1.5,
  pad_after = 1.5,
  output_resolution = NULL,
  group_by = "trial",
  time_audit = FALSE
)
```

## Arguments

- fmri_obj:

  an fmri_ts object containing a single run of data with corresponding
  events

- event:

  the event to which the fmri time series should be aligned (column in
  `fmri_obj$event_data`)

- time_before:

  How many seconds before the `event` do we want data. Default: -3

- time_after:

  How many seconds after the `event` do we want data. Default: 3

- collide_before:

  An optional vector of column names in `trial_df` that set a boundary
  on the earliest time point used in interpolation. This effectively
  truncates the time series for interpolation to a smaller window than
  specified by `time_before`.

- collide_after:

  An optional vector of column names in `trial_df` that set a boundary
  on the latest time point used in interpolation. This effectively
  truncates the time series for interpolation to a smaller window than
  specified by `time_after`.

- pad_before:

  Number of seconds to include in the epoch time window before the event
  of interest. Interpolation spans the window from `time_before` to
  `time_after`, but padding includes data points at the boundary that
  can help provide sufficient data to interpolate early and late times
  within the epoch.

- pad_after:

  Number of seconds to include in the epoch time window after the event
  of interest.

- output_resolution:

  Time step, in seconds, for the interpolated data. Defaults to the
  repetition time stored in `fmri_obj`.

- group_by:

  return interpolated time series for each combination of group_by
  variables. Default is to provide one interpolated time series per
  trial.

- time_audit:

  If TRUE, additional columns will be added to the output showing how
  alignment is calculated vis-a-vis event timing

## Author

Michael Hallquist

## Details

The interpolation grid includes zero and remains within the requested
interval, including zero-width and one-sided windows. When trial is a
grouping variable, missing and unsupported epochs retain missing
estimates and an `epoch_status` column: `"ok"`, `"partial"`,
`"insufficient_data"`, `"missing_onset"`, `"outside_run"`, `"censored"`,
or `"no_data"`. These describe timing and available output support, not
physiological isolation from other trials.
