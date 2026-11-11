# Joint Multi-Lag Trajectory Modeling and Operator-Informed Inference for MEDuSA

Replaces independent lag-by-lag post-hoc testing in MEDuSA Stage 4 with
a unified multi-lag estimation and inference framework that exploits the
exact Stage-2/3 cross-lag point-spread resolution matrix
\\C\_{\mathrm{sym}}\\ and cross-lag operator error correlation
\\R\_{\mathrm{op}}\\:

1.  **Lobe-Normalized Spectral Lavrentiev–Tikhonov Cross-Lag
    De-Blurring + Correlated-SURE 2nd-Order Boosted P-Spline
    Smoothing**: sharpens multiphasic and ramped within-trial neural
    waveforms without boundary distortion or interior gain inflation.

2.  **Multivariate Empirical Bayes (EB) Subject Trajectory Shrinkage**:
    borrows strength across all \\L\\ lags simultaneously via the
    between-subject trajectory covariance \\\hat{\Sigma}\_b\\ and
    subject-specific cross-lag error covariance \\\hat{\Sigma}\_{e,s}\\.

3.  **Operator-Informed Simultaneous Cross-Lag Inference**: calibrates
    epoch-wide family-wise error rate (FWER) via exact correlated
    multivariate-\\t\\ Max-\\\|T\|\\ simultaneous confidence bands
    (`p_maxT`, `conf.low.sim`, `conf.high.sim`) and 1D Threshold-Free
    Cluster Enhancement (`tfce_stat`, `p_tfce`), avoiding the severe
    conservatism of Bonferroni/BY-FDR under strong positive cross-lag
    correlation.

4.  **Continuous Sub-Grid Spline Latency & Shape Estimation**: estimates
    continuous peak latency (`peak_lag`), center-of-mass latency
    (`com_lag`), and full-width at half-maximum (`fwhm`) with
    parametric-bootstrap SEs and confidence intervals.

## Usage

``` r
medusa_multilag_infer(
  data,
  covariates = NULL,
  lag_resolution = NULL,
  lag_error_cor = NULL,
  trajectory_mode = c("auto", "deblur_pspline", "pspline", "deblur", "none"),
  use_whitened = TRUE,
  tr = 1,
  readout_kernel = "triangle",
  readout_width = 1.5,
  decon_lambda = 0.8,
  rho_ar1 = 0.3,
  demix_within_trial = FALSE,
  recalibrate_wt_lobe = FALSE,
  contrasts = NULL,
  n_perm = 4000L,
  n_boot_latency = 400L,
  fine_dt = 0.01,
  conf_level = 0.95,
  polarity = c("positive", "negative", "abs"),
  seed = NULL
)
```

## Arguments

- data:

  Either: (a) a
  [`decon_epoch_regression`](https://hallquistlab.github.io/fmri.pipeline/reference/decon_epoch_regression.md)
  result list (or a list of per-subject `decon_epoch_regression`
  results), (b) a trial-level MEDuSA panel `data.frame`/`data.table`
  (single-subject or multi-subject), or (c) a lag-level coefficient
  table with columns `lag` (or `evt_time`), `estimate` (or
  `beta_corrected`), and `std.error` (or `se_corrected`).

- covariates:

  Optional character vector of covariate names when `data` is a
  trial-level panel. If `NULL`, inferred automatically from panel
  columns.

- lag_resolution:

  Optional named list (by region/parcel/voxel) of \\L \times L\\
  cross-lag resolution matrices \\C\_{\mathrm{sym}}\\, or a single \\L
  \times L\\ matrix.

- lag_error_cor:

  Optional named list (by region/parcel/voxel) of \\L \times L\\
  operator cross-lag error correlation matrices \\R\_{\mathrm{op}}\\, or
  a single matrix.

- trajectory_mode:

  Trajectory estimation mode: `"auto"` (default; applies spectral
  cross-lag de-blurring when `!demix_within_trial` followed by
  correlated-SURE boosted P-spline smoothing), `"deblur_pspline"`,
  `"pspline"`, `"deblur"`, or `"none"`.

- use_whitened:

  Logical; if `TRUE` (default) and whitened columns
  (`decon_demixed_whitened`, `*_whitened`) are present in a trial panel,
  uses them for within-subject lag regressions.

- tr:

  Acquisition TR in seconds (default `1.0`), used when constructing
  reference cross-lag operators if `lag_resolution` is not attached.

- readout_kernel:

  Readout kernel (`"triangle"`, `"boxcar"`, `"gaussian"`, `"linear"`)
  for reference operator construction.

- readout_width:

  Readout kernel width in seconds (default `1.5`).

- decon_lambda:

  Reference deconvolution penalty \\\lambda\\ (default `0.8`).

- rho_ar1:

  Reference AR(1) parameter (default `0.30`).

- demix_within_trial:

  Logical indicating whether Stage-3 within-trial cross-lag de-mixing
  was already applied to the input panel/coefficients (default `FALSE`).

- recalibrate_wt_lobe:

  Logical; if `TRUE` and `demix_within_trial = TRUE`, recalibrates
  interior vs. boundary lobe gain on `W_lag`.

- contrasts:

  Optional named list of lag contrast specifications. Each entry can be
  either a numeric weight vector of length \\L\\ or a list
  `list(early = c(t1, t2), late = c(t3, t4))` specifying two lag windows
  to contrast (`late - early`).

- n_perm:

  Number of correlated multivariate-\\t\\ Monte Carlo draws for
  simultaneous Max-\\\|T\|\\ and 1D TFCE calibration (default `4000L`).

- n_boot_latency:

  Number of parametric bootstrap draws for continuous spline latency SEs
  and confidence intervals (default `400L`).

- fine_dt:

  Sub-grid spline evaluation step in seconds (default `0.01`).

- conf_level:

  Simultaneous and pointwise confidence level (default `0.95`).

- polarity:

  Peak polarity for latency extraction (`"positive"`, `"negative"`, or
  `"abs"`).

- seed:

  Optional random seed for reproducible Monte Carlo / bootstrap draws.

## Value

An object of class `"medusa_multilag"` containing:

- `coefficients`:

  `data.table` of group/unit lag estimates (`estimate`, `std.error`,
  `estimate_multilag`, `std.error_multilag`), pointwise
  \\t\\/\\p\\-values, FDR-adjusted \\p\\-values (`p_fdr_bh`,
  `p_fdr_by`), simultaneous Max-\\\|T\|\\ \\p\\-values and confidence
  bands (`p_maxT`, `conf.low.sim`, `conf.high.sim`), 1D TFCE statistics
  and family-wise \\p\\-values (`tfce_stat`, `p_tfce`), and epoch
  omnibus tests (`omnibus_stat`, `omnibus_p`).

- `latency`:

  `data.table` of continuous sub-grid spline peak latency (`peak_lag`,
  `peak_lag_se`, `peak_lag_conf.low`, `peak_lag_conf.high`), discrete
  grid peak (`peak_lag_grid`), center-of-mass latency (`com_lag`,
  `com_lag_se`), and FWHM (`fwhm`, `fwhm_se`) per unit and term.

- `contrasts`:

  `data.table` of covariance-propagated multi-lag window contrasts
  (`estimate`, `std.error`, `statistic`, `df`, `p.value`, `conf.low`,
  `conf.high`).

- `subjects`:

  `data.table` of subject-specific lag trajectories (`estimate_raw`,
  `estimate_smooth`, and Multivariate Empirical Bayes `estimate_eb`)
  when multi-subject data are supplied.

- `subject_latency`:

  `data.table` of subject-specific continuous spline peak and
  center-of-mass latencies extracted from `estimate_eb`.

- `covariance`:

  Nested list of \\L \times L\\ raw and multi-lag trajectory covariance
  and correlation matrices per unit and term.
