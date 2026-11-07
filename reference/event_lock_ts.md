# function to get interpolated event locked data

function to get interpolated event locked data

## Usage

``` r
event_lock_ts(
  fmri_obj,
  event = NULL,
  time_before = -3,
  time_after = 3,
  collide_before = NULL,
  collide_after = NULL,
  pad_before = -1,
  pad_after = 1,
  time_audit = FALSE,
  overlap_policy = NULL,
  valid_col = NULL,
  neighbor_events = NULL
)
```

## Arguments

- fmri_obj:

  An fmri_ts object containing time series data that are optionally
  keyed on one or more grouping variables (e.g., ROI). The fmri_obj must
  also contain trial-level data (in the \$event_data field) for
  event-locking to proceed.

- event:

  Name of column in `fmri_obj$event_data` that identifies onset for
  event of interest

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
  can help to have sufficient data to interpolate early and late times
  within the epoch.

- pad_after:

  Number of seconds to include in the epoch time window after the event
  of interest.

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

## Details

Nonfinite target onsets and unsupported windows retain placeholder rows
with missing time and values and an `epoch_status` column. Nonfinite
collision candidates are ignored. Nearest-event boundaries remain strict
and do not guarantee disjoint trial windows.
