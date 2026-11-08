# Operator-based epoch de-mixing and gain-corrected trial regression across lags

Estimate trial-by-trial correlations and regression slopes across epoch
lags relative to a target event while removing preceding/neighboring
event leakage and correcting regularization-induced amplitude
attenuation.

## Usage

``` r
decon_epoch_regression(
  deconvolution,
  covariates,
  events = NULL,
  target_condition = NULL,
  lags = seq(-2, 4, by = 1),
  readout_kernel = c("triangle", "linear", "boxcar", "gaussian"),
  readout_width = NULL,
  demix_other_events = TRUE,
  demix_target_trials = TRUE,
  demix_ridge = 0.05,
  gain_floor = 0.05,
  restore_scale = FALSE,
  source_model = c("impulse", "window"),
  return_epochs = FALSE,
  return_panel = FALSE
)
```

## Arguments

- deconvolution:

  Result of
  [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_crossrun_hrf.md).

- covariates:

  Numeric vector, matrix, or data.frame of trial-level predictors for
  target trials (concatenated across runs or supplied as a run-named
  list matching `deconvolution$hrf_fit$run_ids`).

- events:

  Optional run-named list of event data.frames with `onset` and optional
  `condition` and `duration`. Defaults to the HRF fit's event tables.

- target_condition:

  Optional character scalar naming the target event condition in
  `events[[r]]$condition` (e.g., `"feedback"`).

- lags:

  Numeric vector of epoch readout lags (in seconds) relative to each
  target event onset.

- readout_kernel:

  Epoch readout kernel: `"triangle"`, `"linear"`, `"boxcar"`, or
  `"gaussian"`.

- readout_width:

  Optional kernel full width in seconds for non-linear readouts
  (defaults to `1.5 * tr`).

- demix_other_events:

  Logical; whether to estimate non-target event amplitudes at their own
  onsets and subtract their operator-predicted leakage into each target
  epoch lag. Default: `TRUE`.

- demix_target_trials:

  Logical; whether to use the stratum resolution matrix at each lag to
  de-mix adjacent target-trial overlap and correct amplitude
  attenuation. Default: `TRUE`.

- demix_ridge:

  Positive ridge stabilization parameter for unit-diagonal cross-talk
  inversion. Default: `0.05`.

- gain_floor:

  Minimum absolute operator gain used when correcting amplitude slopes.
  Default: `0.05`.

- restore_scale:

  Logical; if `TRUE`, multiplies each run and voxel by its recorded
  `deconvolution$scaling[[r]]$scale` divisor. Default: `FALSE`.

- source_model:

  Forward neural source model at each target lag: `"impulse"` (default)
  or `"window"`.

- return_epochs:

  Logical; whether to include raw and de-mixed trial-by-voxel epoch
  matrices per lag. Default: `FALSE`.

- return_panel:

  Logical; whether to return a tidy MEDUSA-compatible `data.table`
  (`panel`) for downstream
  [`mixed_by()`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md)
  modeling. Default: `FALSE`.

## Value

A list with `coefficients` (per voxel, lag, and predictor term),
`summary` (averaged across usable voxels per lag and term), optional
`epochs`, optional `panel`, and `settings`.
