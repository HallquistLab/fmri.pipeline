# Regularized linear deconvolution of BOLD time series

Estimate a continuous latent activity time series by solving a
regularized inverse convolution problem for one or more BOLD signals.

## Usage

``` r
deconvolve_reglin(
  BOLDobs,
  kernel,
  lambda = NULL,
  lambda_grid = 10^seq(-4, 2, length.out = 20),
  penalty = c("diff2", "diff1", "ridge", "matern1", "sobolev"),
  kernel_grid = NULL,
  hrf_lags = 0L,
  tune_by = c("global", "signal"),
  normalize = TRUE,
  demean = TRUE,
  trim_kernel = TRUE,
  ridge_floor = 1e-08,
  return_diagnostics = FALSE,
  prewhiten_gcv = FALSE,
  prewhiten_rho = NULL,
  backend = c("cpp", "r", "wiener"),
  wiener_shrinkage = c("scalar", "psd"),
  wiener_psd_smooth = 7L,
  wiener_psd_floor = 1e-08,
  gcv_rule = c("min", "1se"),
  scaling = NULL,
  kernel_normalize = c("auto", "signed_sum", "positive_lobe", "peak")
)
```

## Arguments

- BOLDobs:

  matrix of observed BOLD time series (n_timepoints x n_signals) or a
  numeric vector.

- kernel:

  assumed HRF kernel. A vector or single-column matrix.

- lambda:

  regularization strength. If `NULL`, chosen from `lambda_grid` by
  generalized cross-validation.

- lambda_grid:

  candidate regularization strengths used when `lambda` is `NULL`.

- penalty:

  regularization penalty: `"diff2"` for smooth second differences,
  `"diff1"` for smooth first differences, `"ridge"` for amplitude
  shrinkage, `"matern1"` for a damped first-order
  Matern/Ornstein-Uhlenbeck penalty (`0.10 * I + 0.90 * crossprod(D1)`),
  or `"sobolev"` for a critically damped second-order Sobolev penalty
  (`0.05 * I + 0.45 * crossprod(D1) + 0.50 * crossprod(D2)`) that
  suppresses negative sidelobe ringing and stabilizes low frequencies.

- kernel_grid:

  optional list or matrix of candidate HRF kernels. If supplied,
  `kernel` and `hrf_lags` are ignored for tuning.

- hrf_lags:

  integer shifts, in TR bins, applied to `kernel`. Negative values make
  the candidate HRF peak earlier; positive values make it later.

- tune_by:

  whether to tune `lambda`/HRF globally across all signals or
  independently for each signal.

- normalize:

  whether to z-score `BOLDobs` before deconvolution. Default: `TRUE`.

- demean:

  whether to subtract column means before deconvolution. Ignored when
  `normalize = TRUE`, because z-scoring already demeans.

- trim_kernel:

  whether to remove the initial `length(kernel) - 1` latent samples,
  corresponding to pre-observation HRF support.

- ridge_floor:

  small diagonal ridge added for numerical stability.

- return_diagnostics:

  whether to return a list containing diagnostics rather than only the
  deconvolved matrix.

- prewhiten_gcv:

  whether to estimate AR(1) autocorrelation and use a prewhitened
  objective for lambda/HRF selection and deconvolution.

- prewhiten_rho:

  optional AR(1) coefficient used when `prewhiten_gcv = TRUE`. If
  `NULL`, a pooled median AR(1) coefficient is estimated from first-pass
  reconstruction residuals.

- backend:

  implementation backend. `"cpp"` uses the RcppArmadillo regularized
  inverse implementation; `"r"` uses the reference R implementation;
  `"wiener"` uses a frequency-domain Wiener-style deconvolution for
  comparison.

- wiener_shrinkage:

  shrinkage strategy for `backend = "wiener"`. `"scalar"` uses
  `Conj(H) / (|H|^2 + lambda)`. `"psd"` estimates smoothed signal and
  residual noise spectra from a scalar Wiener first pass, then uses
  frequency-dependent shrinkage.

- wiener_psd_smooth:

  odd integer smoothing width for the estimated signal/noise spectra in
  `wiener_shrinkage = "psd"`.

- wiener_psd_floor:

  lower bound for estimated spectra.

- gcv_rule:

  lambda/HRF selection rule. `"min"` selects the minimum mean GCV
  candidate. `"1se"` selects the largest lambda whose mean GCV is within
  one standard error of the minimum. The standard error is computed
  across signal columns for `tune_by = "global"`.

- scaling:

  Explicit input scaling: `"zscore"`, `"noise_sd"`, `"demean"`,
  `"none"`, or `"psc"` (percent signal change). `"noise_sd"` centers
  each signal and divides by a high-frequency second-difference estimate
  of observation noise SD (adjusted for AR(1) autocorrelation when
  available), standardizing SNR without conflating neural task variance
  with noise variance. When supplied, overrides `normalize` and
  `demean`. Otherwise those flags retain their existing behavior; the
  default is z-scoring.

- kernel_normalize:

  Normalization applied to `kernel` or `kernel_grid`: `"auto"` (default;
  honors `attr(..., "normalize")` when present, switches to
  `"positive_lobe"` when candidate kernels span heterogeneous undershoot
  ratios, and uses `"signed_sum"` otherwise), `"signed_sum"`
  (`sum(h) == 1`), `"positive_lobe"` (`sum(pmax(h, 0)) == 1`), or
  `"peak"` (`max(h) == 1`).

## Value

A matrix of deconvolved time series, or a list with `activity` and
tuning diagnostics when `return_diagnostics = TRUE`.

## Details

Constant signals and signals containing nonfinite observations are
excluded from fitting and tuning. Their activity columns are returned as
`NA`; diagnostics include a per-column `signal_status` vector (`"ok"`,
`"constant"`, or `"nonfinite"`). Partial missingness excludes the entire
signal; observations are not interpolated or deleted. Z-scoring remains
the default and is recommended for standardized within-run trial
responses. For fixed HRF/regularization, the linear inverse is
homogeneous in its input: rescaling the observations rescales activity
without changing the conditioning of the inverse matrix. Without
z-scoring, global GCV gives signals with larger variance more weight.
Demeaning retains input units; `"none"` retains offsets too and is
mainly useful for already-centered data or simulations. Percent signal
change uses `100 * (y - mean(y)) / mean(y)` and requires a positive
baseline well away from zero. It is unsuitable for zero-mean residuals.
No scaling choice removes HRF uncertainty or regularization-induced
attenuation.

Scaling is applied once. The matrix attributes `scaled:center`,
`scaled:scale`, and `scaling` record the actual BOLD transformation
`(y - center) / scale`; diagnostics expose the same values. With
`scaling = "none"`, the recorded center is zero. The activity scale
depends on this transformation and the HRF normalization convention
(`kernel_normalize`); it is not a calibrated physiological
neural-amplitude scale. The removed BOLD mean is not an identified
neural baseline.

Scalar Wiener GCV uses the exact trace of the finite observed smoother,
including zero extension, latent truncation, and the implemented
filtered HRF when whitening. This trace is conditional on preprocessing
and the AR(1) coefficient; it does not include derivatives of estimated
preprocessing/noise parameters. Adaptive PSD shrinkage requires a fixed
lambda and a single HRF. Its diagnostics label the filter trace
`conditional_fixed_psd` and return missing GCV values, because the trace
holding spectra fixed is not the degrees of freedom of the full
nonlinear, data-adaptive estimator. Scalar diagnostics use
`df_method = "finite_observation_trace"`.
