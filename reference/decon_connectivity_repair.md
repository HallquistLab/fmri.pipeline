# Measurement-error-corrected inter-region connectivity and PPI via REPAIR

Estimates trial-by-trial seed-to-target neural connectivity, directed
lagged connectivity (controlling for baseline target activity at the
seed lag), and task-modulated psychophysiological interactions (PPI)
from a
[`decon_epoch_regression`](https://hallquistlab.github.io/fmri.pipeline/reference/decon_epoch_regression.md)
trial panel using the REPAIR modified method-of-moments (`M-MOMS` /
`M-MOM`) framework of Vig & Hallquist (2026). Propagates per-trial
operator standard errors (`se_demixed`) into exact heteroscedastic
predictor and outcome error covariances,including the rank-1 interaction
error block for PPI (`seed * x`) and dual predictor error variances for
Granger-style auto-lag controls (`targ_t0`).

## Usage

``` r
decon_connectivity_repair(
  panel_dt,
  seed_region,
  target_region,
  seed_lag = 0,
  target_lag = 0,
  covariates = NULL,
  test_ppi = FALSE,
  control_target_seed_lag = FALSE,
  orthogonalize_seed_x = TRUE,
  method = c("mmoms", "mmom", "eb_predictor", "naive_gls", "naive_unweighted"),
  use_whitened = NULL,
  region_col = NULL,
  id_col = "id",
  return_augmented_panel = FALSE
)
```

## Arguments

- panel_dt:

  A `data.frame` or `data.table` returned in `$panel` by
  [`decon_epoch_regression`](https://hallquistlab.github.io/fmri.pipeline/reference/decon_epoch_regression.md).

- seed_region:

  Character scalar naming the seed region/parcel/voxel.

- target_region:

  Character scalar naming the target region/parcel/voxel.

- seed_lag:

  Numeric lag (in seconds) for the seed readout (default `0`).

- target_lag:

  Numeric lag (in seconds) for the target readout (default `0`).

- covariates:

  Optional character vector of error-free trial covariate columns in
  `panel_dt`. Defaults to `"x"` when present in `panel_dt`, or
  `character(0)` otherwise.

- test_ppi:

  Logical; if `TRUE`, includes the interaction between the seed readout
  and the first covariate in `covariates` (`seed_c:x_c`) with exact
  per-trial interaction error covariance \\\tilde{\psi}\_{S,j} \[1,
  x\_{c,j}; x\_{c,j}, x\_{c,j}^2\]\\. Default: `FALSE`.

- control_target_seed_lag:

  Logical; if `TRUE` and `abs(target_lag - seed_lag) > 1e-6`, includes
  the target region's readout at `seed_lag` (`targ_t0_val`) as an
  error-prone control covariate with measurement-error variance
  `se_demixed^2`. Default: `FALSE`.

- orthogonalize_seed_x:

  Logical; if `TRUE` (default), residualizes the run-centered seed
  readout on the error-free covariates within each subject-run and
  scales the seed error variance by \\1 - h\_{x,r,j}\\ prior to forming
  PPI interactions, preventing quadratic task drive (\\x_c^2\\) from
  confounding physiological-by-psychological interactions.

- method:

  Estimation method: `"mmoms"` (default; adaptive step-up REPAIR
  `M-MOMS`), `"mmom"` (fixed minimum-regularization REPAIR `M-MOM`),
  `"eb_predictor"` (Skrondal & Laake 2001 role-specific Empirical-Bayes
  shrunk predictor on likelihood-only outcome), `"naive_gls"`
  (precision-weighted WLS without predictor error subtraction), or
  `"naive_unweighted"` (unweighted OLS).

- use_whitened:

  Logical; if `TRUE`, uses `decon_demixed_whitened` and `se_whitened`
  when present. Defaults to `FALSE` for `"mmoms"` and `"mmom"` (so
  elementwise PPI products are formed on unmixed trial states before
  error correction) and `TRUE` for `"naive_gls"`.

- region_col:

  Optional column name identifying regions/parcels/voxels in `panel_dt`
  (auto-detected from `"region"`, `"parcel"`, `"voxel"`).

- id_col:

  Column name identifying subjects in `panel_dt` (default `"id"`).

- return_augmented_panel:

  Logical; if `TRUE`, returns the merged trial panel augmented with
  moment-preserving REPAIR predictor columns (`seed_repair`,
  `ppi_repair`, `repair_weight`) that reproduce subject-level `M-MOMS`
  slopes in weighted least squares / `lmer()`.

## Value

A list of class `"medusa_connectivity_repair"` containing `group`
(group-level random-effects meta-analytic inference table), `subjects`
(per-subject estimates, SEs, and REPAIR diagnostics), optional
`panel_repair`, and `settings`.
