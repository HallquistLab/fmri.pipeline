# Diagnose amplitude attenuation and event leakage for a frozen inverse

Diagnose amplitude attenuation and event leakage for a frozen inverse

## Usage

``` r
diagnose_crossrun_recovery(
  deconvolution,
  events = NULL,
  lags = 0,
  weights = NULL,
  source_kernel_grid = NULL,
  return_matrices = FALSE,
  max_events = 2000L
)
```

## Arguments

- deconvolution:

  Result of deconvolve_crossrun_hrf().

- events:

  Optional run-named event tables with onset and optional condition and
  duration. Defaults to the HRF training event records. Every record is
  a separate unit-amplitude source, regardless of its condition label.

- lags:

  Readout times relative to each event onset, in seconds. Default 0.

- weights:

  Linear readout weights, default one per lag. For example, lags = -1:3
  produces a five-sample window sum at TR = 1.

- source_kernel_grid:

  Optional normalized HRF candidates used to generate synthetic BOLD.
  Default uses each stratum's fitted HRF. Honors
  `attr(source_kernel_grid, "normalize")` when present, otherwise
  matching the fitted normalization mode. Supplying alternatives
  evaluates sensitivity to HRF mismatch; it does not estimate the true
  HRF.

- return_matrices:

  Return event-by-event recovered, ideal, and excess resolution
  matrices, default FALSE.

- max_events:

  Explicit limit per run on the quadratic diagnostic matrix,
  default 2000. Use a smaller prespecified subset or explicitly raise
  the limit.

## Value

A list with per-event diagnostics, stratum summaries, voxel-to-stratum
mappings, optional resolution matrices, and the conditional assumptions.

## Details

For source events B, convolution H, frozen inverse D, and readout S,
recovered = S D H B and ideal = S B. Gain is their diagonal ratio.
Excess leakage is recovered minus ideal, normalized by each receiving
event's ideal diagonal. Subtracting ideal retains natural overlap of
neural readout windows. Previous/next refer to chronological source
records, including other event types. Zero ideal diagonals and readouts
outside the acquisition grid are NA.

The calculation reuses one operator per fitted HRF/noise stratum and
uses BLAS across source events, not a simulation or fit per voxel. It
includes the actual masks, nuisance projection, padding, centering, and
inverse penalty. Scaling divisors cancel in the gain ratios and are
omitted. Their estimated values, noise parameters, selected HRFs, and
lambda are held fixed: these are conditional linear diagnostics, not
uncertainty intervals or a claim about physiological neural units.
Source events of positive duration have unit height (fractional
acquisition-bin occupancy), not unit integrated area. No automatic
correction is applied; a scalar cannot undo a full matrix of attenuation
and cross-talk. Wiener diagnostics are not supported here.
