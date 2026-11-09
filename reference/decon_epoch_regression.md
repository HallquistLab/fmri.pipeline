# Operator-based epoch de-mixing and gain-corrected trial regression across lags

Estimate trial-by-trial correlations and regression slopes across epoch
lags relative to a target event while removing preceding/neighboring
event leakage, correcting regularization-induced amplitude attenuation,
and propagating per-trial operator precision weights and within-run
error covariance.

## Usage

``` r
decon_epoch_regression(
  deconvolution,
  covariates = NULL,
  events = NULL,
  target_condition = NULL,
  lags = seq(-2, 4, by = 1),
  readout_kernel = c("triangle", "linear", "boxcar", "gaussian"),
  readout_width = NULL,
  demix_other_events = TRUE,
  demix_target_trials = TRUE,
  demix_within_trial = FALSE,
  demix_ridge = 0.05,
  within_trial_ridge = 0.06,
  gain_floor = 0.05,
  restore_scale = FALSE,
  source_model = c("impulse", "window", "duration"),
  trial_durations = NULL,
  duration_scale = c("rate", "integral"),
  weighting = c("precision", "gls", "none"),
  tau2_mode = c("auto", "residual", "marginal"),
  aggregate_by_parcel = FALSE,
  parcel_groups = NULL,
  whiten_panel = FALSE,
  return_epochs = FALSE,
  return_panel = FALSE
)
```

## Arguments

- deconvolution:

  Result of
  [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_crossrun_hrf.md).

- covariates:

  Optional numeric vector, matrix, or data.frame of trial-level
  predictors for target trials (concatenated across runs or supplied as
  a run-named list matching `deconvolution$hrf_fit$run_ids`). When
  `NULL` (default), performs covariate-agnostic epoch readout,
  de-mixing, biological variance estimation (`tau2_mode = "marginal"`),
  Empirical-Bayes precision weighting, and optional within-run ZCA
  whitening (`whiten_panel = TRUE`), returning the trial panel
  (`return_panel` defaults to `TRUE` when `covariates = NULL`) for
  downstream two-stage
  [`mixed_by`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md)
  behavioral modeling or inter-region connectivity (`R2 ~ R1_seed * x`).

- events:

  Optional run-named list of event data.frames with `onset` and optional
  `condition` and `duration`. Defaults to the HRF fit's event tables.

- target_condition:

  Optional character scalar naming the target event condition in
  `events[[r]]$condition` (e.g., `"feedback"`). When supplied,
  non-target events in `events[[r]]` are treated as interfering events
  whose leakage into the target epoch is removed when
  `demix_other_events = TRUE`.

- lags:

  Numeric vector of epoch readout lags (in seconds) relative to each
  target event onset.

- readout_kernel:

  Epoch readout kernel: `"linear"`, `"boxcar"`, `"triangle"`, or
  `"gaussian"`. Default: `"triangle"`.

- readout_width:

  Optional kernel full width in seconds for non-linear readouts
  (defaults to `1.5 * tr`).

- demix_other_events:

  Logical; whether to estimate non-target event contributions at their
  own onsets and subtract their operator-predicted leakage into each
  target epoch lag. When `demix_ridge >= 0.03` (default `0.05`), uses a
  condition-mean non-target de-mixing operator that removes preceding
  cue/event bleed without inflating target-trial variance; when
  `demix_ridge < 0.03` (e.g. `0.02`), activates a 2nd-order iterated
  Tikhonov joint \\(T, O)\\ operator to de-mix strongly trial-varying
  non-target confounds. Default: `TRUE`.

- demix_target_trials:

  Logical; whether to use the stratum resolution matrix at each lag to
  de-mix adjacent target-trial overlap and restore
  regularization-attenuated amplitude gain via 2nd-order iterated
  Tikhonov inversion (`solve_duration_demix_operator()`). Default:
  `TRUE`.

- demix_within_trial:

  Logical; whether to apply an additional multi-lag point-spread
  de-mixing step across the requested `lags` (when at least 3 lags are
  evaluated) to sharpen temporal localization across lags within each
  target trial. Recommended (`TRUE`) when estimating within-trial
  temporal waveforms \\\beta(\tau)\\ across a multi-lag grid. Default:
  `FALSE`.

- demix_ridge:

  Positive ridge stabilization parameter for target and non-target
  cross-talk inversion. Values `>= 0.03` (default `0.05`) use
  condition-mean non-target de-mixing; values `< 0.03` activate joint
  trial-varying non-target de-mixing. Default: `0.05`.

- within_trial_ridge:

  Positive ridge stabilization parameter for multi-lag within-trial
  point-spread inversion when `demix_within_trial = TRUE`. Default:
  `0.06`.

- gain_floor:

  Minimum absolute operator gain used when correcting amplitude slopes.
  Default: `0.05`.

- restore_scale:

  Logical; if `TRUE`, multiplies each run and voxel by its recorded
  `deconvolution$scaling[[r]]$scale` divisor so that `beta_raw` and
  `beta_corrected` (and panel readouts) are reported in unscaled
  physical input units rather than standardized (`zscore`/`noise_sd`)
  units. Strongly recommended (`TRUE`) when pooling across multiple runs
  in
  [`mixed_by`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md).
  Default: `FALSE`.

- source_model:

  Forward neural source model at each target lag: `"impulse"` (default
  unit impulse at `onset + lag`), `"window"` (local sustained boxcar of
  width `readout_width` centered at `onset + lag`), or `"duration"`
  (trial-varying sustained boxcar over
  `[onset + lag, onset + lag + d_i]` using each target trial's duration
  `d_i` from `trial_durations` or `events[[r]]$duration`, with
  duration-matched readout and unpenalized baseline-rate
  de-confounding).

- trial_durations:

  Optional numeric vector (concatenated across runs) or run-named list
  of non-negative target trial durations (in seconds), such as per-trial
  reaction times. Overrides `events[[r]]$duration` for target trials and
  defaults `source_model` to `"duration"` when `source_model` is not
  explicitly specified.

- duration_scale:

  Scaling convention for positive-duration target trials when
  `source_model = "duration"`: `"rate"` (default; trial amplitude
  represents neural firing rate / boxcar height per TR bin, so
  integrated neural activity scales with duration and constant baseline
  firing rates are de-confounded from duration variability) or
  `"integral"` (each trial boxcar is normalized to unit time-integrated
  neural mass across its duration).

- weighting:

  Trial regression weighting scheme across target trials and runs:
  `"precision"` (default; Empirical-Bayes / floor-stabilized inverse
  error-variance weights propagating per-trial operator standard errors
  and run noise scale), `"gls"` (null-mode-damped, Bartlett-adaptive
  within-run ZCA correlation whitening plus precision weighting
  accounting for off-diagonal trial-to-trial error correlations), or
  `"none"` (unweighted OLS). All three slope estimates
  (`beta_unweighted`, `beta_weighted`, `beta_gls`) and their standard
  errors are always included in `coefficients`.

- tau2_mode:

  Estimator for biological trial-to-trial variance
  \\\hat{\tau}^2\_{\mathrm{bio}}\\ used in Empirical-Bayes trial
  precision weights and within-run ZCA whitening: `"auto"` (default;
  uses `"marginal"` when `covariates = NULL` and an Empirical-Bayes
  blend of residual and marginal estimators when covariates are
  supplied), `"residual"` (method-of-moments estimator from OLS
  covariate regression residuals minus projected operator variance
  \\\mathrm{tr}((\mathbf{I}-\mathbf{P}\_X)\boldsymbol{\Omega})\\), or
  `"marginal"` (covariate-agnostic estimator combining run-centered
  marginal variance, lag-1 trial difference variance, and robust MAD
  scale minus analytical operator variance expectations).

- aggregate_by_parcel:

  Logical; if `TRUE`, aggregates deconvolved signals across voxels
  within each spatial parcel (using inverse-variance voxel precision
  weights and empirical spatial noise correlation propagation) prior to
  Stage-2 epoch de-mixing and Stage-3 regression/ZCA whitening. Replaces
  voxel-by-voxel Stage-3/4 fitting with parcel-level inference while
  calibrating `se_demixed`, `tau2_bio`, `precision_weight`, and ZCA
  whitening directly on the high-SNR parcel series. Default: `FALSE`.

- parcel_groups:

  Optional character or factor vector of parcel/region labels for each
  voxel (either unnamed of length `n_voxels` or named by `voxel_ids`).
  Defaults to `deconvolution$hrf_fit$specification$spatial_groups` when
  `aggregate_by_parcel = TRUE`. Supplying a non-`NULL` `parcel_groups`
  automatically enables `aggregate_by_parcel = TRUE`.

- whiten_panel:

  Logical; if `TRUE` and `return_panel = TRUE`, adds within-run
  null-mode-damped, Bartlett-adaptive ZCA-whitened outcome
  (`decon_demixed_whitened`) and covariate (`<var>_whitened`) columns to
  `panel` so downstream `mixed_by(..., weights = "precision_weight")`
  models account for both heteroscedastic trial precision and
  off-diagonal within-run trial error correlations. Recommended (`TRUE`)
  for group mixed-effects inference. Default: `FALSE`.

- return_epochs:

  Logical; whether to include raw and de-mixed trial-by-voxel (or
  trial-by-parcel) epoch matrices per lag. Default: `FALSE`.

- return_panel:

  Logical; whether to return a tidy MEDUSA-compatible `data.table`
  (`panel`) with columns `run`, `trial`, `event_index`, `onset`,
  `evt_time`, `voxel` (and `parcel` when `aggregate_by_parcel = TRUE`),
  `decon_raw`, `decon_demixed`, `se_raw`, `se_demixed`,
  `precision_weight`, `precision_weight_op`, `precision_weight_raw`,
  `error_ac1`, `tau2_bio`, and the covariate columns for downstream
  [`mixed_by()`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md)
  modeling. Defaults to `TRUE` when `covariates = NULL` and
  `return_epochs = FALSE`.

## Value

A list with `coefficients` (per voxel/parcel, lag, and predictor term),
`summary` (averaged across usable units per lag and term), optional
`epochs`, optional `panel`, and `settings`.

## Details

`decon_epoch_regression()` bridges continuous regularized deconvolution
and trial-level cognitive neuroscience models in five stages:

1.  **Kernel Quadrature Readout:** At each requested lag \\\tau\\,
    continuous deconvolved activity is integrated around each target
    event using `readout_kernel`. For `source_model = "duration"`, the
    readout window centers on the trial's midpoint \\t_i + \tau + 0.5
    d_i\\ and expands to \\\max(w\_{\mathrm{src}}, d_i)\\.

2.  **Non-Target Event De-Mixing (`demix_other_events`):** Using the
    frozen inverse operator \\\mathbf{D}\_r\\ and forward HRF design
    \\\mathbf{G}\_{r,0}\\, the function subtracts resolution-predicted
    leakage from non-target events (e.g., cues or motor responses) using
    a condition-mean operator (`demix_ridge >= 0.03`) or a 2nd-order
    Tikhonov joint \\(T, O)\\ operator (`demix_ridge < 0.03`).

3.  **Target Cross-Talk & Duration De-Mixing (`demix_target_trials`):**
    Inverts the lag-specific target resolution matrix via a 2nd-order
    iterated Tikhonov operator (`solve_duration_demix_operator()`) with
    per-trial effective gain restoration to remove trial-to-trial
    overlap and divide out regularization-induced amplitude shrinkage.
    In `source_model = "duration"`, a centered-subspace penalty and
    rank-1 baseline adjustment keep constant baseline firing rates
    invariant to trial duration (preventing reaction-time confounding).

4.  **Optional Multi-Lag Point-Spread De-Mixing
    (`demix_within_trial`):** When enabled across \\\ge 3\\ lags,
    sharpens temporal specificity across lags using a 2nd-order Tikhonov
    lag-by-lag resolution inverse while preserving cross-trial operator
    covariance structure.

5.  **Analytical Uncertainty, Precision Weights & Null-Mode-Damped ZCA
    Whitening:** Propagates observation-level AR(1) noise covariance
    through the composite linear readout/de-mixing operator
    \\\mathbf{B}\_{r,\tau}\\ to obtain exact per-trial standard errors
    (`se_demixed`), adjacent-trial error autocorrelation (`error_ac1`),
    Empirical-Bayes precision weights (`precision_weight`),
    covariance-propagated sandwich SEs (`se_weighted`), and
    null-mode-damped, Bartlett-adaptive within-run ZCA-whitened
    variables for
    [`mixed_by`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md).

## Recommended Best-Practice Configuration

Across the 5-pillar MEDuSA Monte Carlo validation suite (amplitude
recovery, temporal alignment, inter-region connectivity, factorial power
surfaces, and biophysical Windkessel-Balloon + `fmrisim` scanner
realism), the recommended end-to-end configuration is:

1.  **Stage 0–1 (HRF Estimation):** Use
    `hrf_library(tr, preset = "medusa")` (0.5s peak-delay grid over 4–8s
    with dispersions 0.9, 1.0, 1.1; 27 candidates) and
    [`estimate_crossrun_hrf`](https://hallquistlab.github.io/fmri.pipeline/reference/estimate_crossrun_hrf.md)
    with `noise_model = "ar1"`, `evidence = "calibrated"`,
    `selection_mode = "soft_average"` (evidence-weighted HRF kernel
    averaging across locally supported shapes), and parcel labels in
    `spatial_groups`.

2.  **Stage 2 (Deconvolution):** Use
    [`deconvolve_crossrun_hrf`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_crossrun_hrf.md)
    with `scaling = "noise_sd"` (AR(1) innovation noise SD rather than
    total `"zscore"` variance) and `penalty = "sobolev"`.

3.  **Stage 3 (Epoch De-Mixing & Whitening):** Use
    `decon_epoch_regression(..., demix_other_events = TRUE, demix_target_trials = TRUE, restore_scale = TRUE, weighting = "gls", whiten_panel = TRUE, return_panel = TRUE)`.
    Enable `demix_within_trial = TRUE` with `lags` spaced at
    \\1\times\mathrm{TR}\\ when estimating multi-lag waveforms
    \\\beta(\tau)\\ (if `lags` is sampled finer than
    \\1\times\mathrm{TR}\\, `within_trial_ridge` is automatically scaled
    by \\(\mathrm{TR}/\Delta\tau)^2\\ to prevent sub-TR collinearity
    amplification), and set `source_model = "duration"`
    (`duration_scale = "rate"`) when modeling variable-duration (RT)
    epochs.

4.  **Stage 4 (Group Mixed-Effects Inference):** Pass both the
    ZCA-whitened panel columns and precision weights to
    [`mixed_by`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md):
    `mixed_by(panel, decon_demixed_whitened ~ x_whitened + (1 + x_whitened | id), weights = "precision_weight")`.
