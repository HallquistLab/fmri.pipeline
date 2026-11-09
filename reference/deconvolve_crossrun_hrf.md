# Deconvolve runs using frozen HRFs learned from other runs

Deconvolve runs using frozen HRFs learned from other runs

## Usage

``` r
deconvolve_crossrun_hrf(
  bold,
  hrf_fit,
  lambda = NULL,
  scaling = "zscore",
  penalty = c("diff2", "diff1", "ridge", "matern1", "sobolev"),
  backend = c("cpp", "r", "wiener"),
  batch_size = 2048L,
  tail_padding = 0L,
  ridge_floor = 1e-08,
  operator_cache_mb = 128,
  nuisance_mode = NULL,
  observation_mask = NULL,
  lambda_grid = 10^seq(-2, 1, length.out = 16),
  gcv_rule = c("min", "1se")
)
```

## Arguments

- bold:

  Named time-by-voxel run matrices with IDs/order/dimensions matching
  the HRF fit. These should be the same preprocessed BOLD used for
  estimation.

- hrf_fit:

  Object returned by estimate_crossrun_hrf().

- lambda:

  Nonnegative regularization strength. If `NULL` (the default for modern
  `ar1`/`calibrated` fits), selected automatically from `lambda_grid` by
  gain-calibrated, task-weighted adaptive prewhitened generalized
  cross-validation using `gcv_rule`.

- scaling:

  Input scaling: zscore (default), noise_sd, demean, none, or psc.

- penalty:

  Quadratic penalty: diff2, diff1, ridge, matern1, or sobolev.

- backend:

  Fixed linear inverse backend: cpp, r, or scalar wiener.

- batch_size:

  Maximum voxel columns per scaling/multiplication block.

- tail_padding:

  Number of appended zero BOLD observations, default zero. Use
  nrow(hrf_fit\$kernels) to reproduce the current voxelwise wrapper
  assumption.

- ridge_floor:

  Numerical diagonal ridge for cpp/r.

- operator_cache_mb:

  Maximum inverse-operator cache size in MiB, default 128. Zero disables
  reuse across runs. A single operator may exceed the cache budget while
  it is being applied, but is not retained in the cache.

- nuisance_mode:

  Jointly fit unpenalized intercept, trend, and supplied nuisance terms,
  or ignore them (default, preserving prepared input signals). Joint
  projection is applied to both data and design, but can amplify weakly
  identified low-frequency activity in a flexible continuous inverse.
  Validate it for the intended design rather than assuming it improves
  recovery.

- observation_mask:

  Use fit's valid_rows in the inverse (fit), or all observations (all).
  Defaults to fit for AR/calibrated, all for OLS/legacy.

- lambda_grid:

  Candidate regularization values evaluated when `lambda = NULL`.

- gcv_rule:

  Selection rule (`"min"` or `"1se"`) used when `lambda = NULL`.

## Value

A list with activity matrices, per-run scaling/status tables, the HRF
fit and frozen inverse settings. Activity matrices retain voxel IDs and
carry tr and time_offset attributes for subsequent event alignment.

## Details

One inverse operator is constructed per selected HRF/run length, cached
for consecutive runs with identical masks/nuisance designs, and applied
using block matrix products for cpp/r. Scalar Wiener uses direct block
FFT fits. HRFs are frozen; no voxelwise factorization or R fitting loop
is required. AR/calibrated fits use a masked GLS inverse with the frozen
run-specific noise covariance. A gap of d acquisitions uses rho^d and
the corresponding accumulated innovation variance. Censored rows
contribute no observations; latent activity is still estimated on the
complete original grid. Its observation_mask attribute identifies
modeled versus observed samples. Nonfinite/constant signals on usable
observations remain NA. Nonfinite censored observations do not
invalidate a run. Joint nuisance fitting avoids projecting only the BOLD
and inadvertently changing the task model. Wiener requires OLS,
nuisance_mode = "ignore", observation_mask = "all". Scaling is applied
once before optional padding, and its centers/divisors are returned.
Different run z-scores do not preserve absolute between-run units.
