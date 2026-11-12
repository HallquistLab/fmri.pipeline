# Diagnose temporal spread of a fixed linear deconvolution operator

Generate isolated unit impulses at actual event times, convolve them
with a known HRF, and measure their recovery and entry into other event
windows.

## Usage

``` r
deconvolution_impulse_diagnostics(
  event_data,
  kernel,
  n_time,
  tr,
  lambda,
  onset = "onset",
  trial = "trial",
  time_offset = 0,
  time_before = -3,
  time_after = 3,
  output_resolution = tr,
  backend = c("cpp", "r", "wiener"),
  penalty = c("diff2", "diff1", "ridge", "matern1", "sobolev"),
  forward_kernel = NULL,
  prewhiten_rho = 0,
  ridge_floor = 1e-08,
  input_scale = 1,
  input_demean = FALSE,
  tail_padding = 0L,
  source_trials = NULL,
  batch_size = 32L,
  return_timeseries = FALSE,
  kernel_normalize = c("auto", "signed_sum", "positive_lobe", "peak")
)
```

## Arguments

- event_data:

  Data.frame with unique trial labels and numeric onsets.

- kernel:

  Frozen inverse-model HRF, normalized according to `kernel_normalize`.

- n_time:

  Number of measured acquisition samples (at least three).

- tr:

  Acquisition interval in seconds.

- lambda:

  Fixed regularization strength; automatic tuning is not performed.

- onset, trial:

  Names of onset and trial columns in event_data.

- time_offset:

  Absolute time assigned to the first acquisition.

- time_before, time_after:

  Requested window limits relative to each event.

- output_resolution:

  Spacing of the diagnostic window grid, in seconds.

- backend:

  Linear backend: cpp, r, or scalar wiener. Adaptive PSD and nonlinear
  estimators are not supported.

- penalty:

  Fixed quadratic penalty for cpp/r.

- forward_kernel:

  Generating HRF, defaulting to kernel. Use a different shape for a
  specified HRF-mismatch experiment, normalized with the same
  `kernel_normalize` convention as `kernel`.

- prewhiten_rho:

  Frozen AR(1) coefficient; zero disables whitening.

- ridge_floor:

  Fixed numerical ridge floor for cpp/r.

- input_scale:

  Fixed positive BOLD divisor from the fitted run. Recovery is
  multiplied back by this divisor to report activity in injected-impulse
  units.

- input_demean:

  Whether to subtract each simulated measured BOLD mean before padding.
  This applies a fixed linear centering operation, not z-scoring.

- tail_padding:

  Number of zero BOLD samples appended as observations. Use
  length(kernel) to reproduce the current voxelwise wrapper's padding.

- source_trials:

  Optional trial labels to inject; all are used by default. All
  event_data rows remain target windows. Subset sources for large
  schedules.

- batch_size:

  Maximum number of impulse columns fitted together.

- return_timeseries:

  Include full acquisition-grid truth, BOLD, and recovered activity for
  each selected source (default FALSE).

- kernel_normalize:

  Normalization applied to both `kernel` and `forward_kernel`: `"auto"`
  (default; honors `attr(..., "normalize")` when present, switches to
  `"positive_lobe"` when `kernel` and `forward_kernel` differ in
  undershoot ratio, and uses `"signed_sum"` otherwise), `"signed_sum"`,
  `"positive_lobe"`, or `"peak"`.

## Value

A list with impulses (one status/summary per selected source), windows
(directed source/target metrics), optional time_series, and frozen
settings.

## Details

Off-grid events deposit mass on the two surrounding acquisition bins;
this is a grid approximation, not a sub-TR neural-resolution claim. The
injected samples sum to one (area tr). Missing/outside onsets and
impulses producing no observed BOLD variation retain explicit statuses
and missing recovery metrics. Windows use overlapping support without
collision censoring. window_l1_ratio compares absolute recovered area in
a target window with that in the source's own window. Windows may
overlap, so ratios are not fractions of neural attribution and may
exceed one. Metrics condition on fixed model settings and the specified
generating HRF. No noise, retuning, observation- mask refitting, or
physiological uncertainty estimation is performed.
