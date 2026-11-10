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
  time_audit = FALSE,
  overlap_policy = NULL,
  valid_col = NULL,
  neighbor_events = NULL,
  max_gap = NULL,
  readout_kernel = "linear",
  readout_width = NULL
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

- overlap_policy:

  One of `"allow"`, `"truncate"`, or `"exclude"`. NULL uses truncation
  when collision columns are supplied and allows overlap otherwise.
  Explicit allow retains full windows. Exclude removes whole epochs
  whose closed requested windows overlap another target epoch, including
  ties.

- valid_col:

  Optional logical source-data column. FALSE or NA censors that row for
  all value columns without deleting its acquisition time.

- neighbor_events:

  Event columns used for nearest-event annotations; defaults to the
  target event column. Does not itself censor observations.

- max_gap:

  Maximum interpolation bracket duration in seconds, defaulting to the
  acquisition TR. Explicit invalid rows always block interpolation; this
  tolerance only permits wider brackets across deleted rows.

- readout_kernel:

  Readout kernel applied at each output grid time: `"linear"` (default
  two-point bracket interpolation), `"boxcar"`, `"triangle"`, or
  `"gaussian"` local window integration within each contiguous
  uncensored segment.

- readout_width:

  Full window width (in seconds) for non-linear `readout_kernel` modes.
  Defaults to `1.5 * max_gap`.

## Details

The interpolation grid includes zero and remains within the requested
interval, including zero-width and one-sided windows. When trial is a
grouping variable, missing and unsupported epochs retain missing
estimates and an `epoch_status` column: `"ok"`, `"partial"`,
`"insufficient_data"`, `"missing_onset"`, `"outside_run"`, `"censored"`,
or `"no_data"`. These describe timing and available output support, not
physiological isolation from other trials. Interpolation occurs within
each unit before aggregation, preserving explicit missing rows and
acquisition gaps. Off-grid median/SD now summarize the interpolated
units. Each value has contributor/reason counts, source IDs, and the
maximum distinct-epoch reuse count. Source IDs index observed
acquisition times within this run; reuse is computed within the returned
grouping. The `source_clock` attribute maps IDs to absolute acquisition
times. Epoch annotations beginning with `align_` are always included
when trial is a grouping variable. They describe timing, not
physiological attribution. Whole-window exclusion adds
`"overlap_excluded"` to epoch statuses.

## Author

Michael Hallquist
