# Estimate HRF shapes from other runs, with optional spatial regularization

Estimate HRF shapes from other runs, with optional spatial
regularization

## Usage

``` r
estimate_crossrun_hrf(
  bold,
  events,
  tr,
  kernel_grid,
  nuisance = NULL,
  valid_rows = NULL,
  time_offset = 0,
  spatial_groups = NULL,
  pooling = c("voxel", "parcel"),
  spatial_strength = 0,
  fallback_kernel = 1L,
  min_r2 = 0.01,
  min_training_runs = 1L,
  min_parcel_voxels = 5L,
  batch_size = 2048L,
  return_scores = FALSE,
  amplitude_model = c("condition", "trial"),
  noise_model = c("ar1", "ols"),
  noise_rho = NULL,
  noise_grid = seq(-0.8, 0.9, by = 0.1),
  evidence = c("calibrated", "legacy"),
  evidence_alpha = 0.05,
  canonical_strength = 0.01
)
```

## Arguments

- bold:

  Named list of time-by-voxel numeric matrices for one participant.
  Voxel column IDs must be unique and identical in all runs, in the same
  order.

- events:

  Named list of data.frames with onset (seconds), condition, and
  optional duration (seconds, default zero). Names/order must match
  bold.

- tr:

  Common acquisition interval in seconds. All runs use the same grid.

- kernel_grid:

  Matrix or list of equal-length candidate causal HRFs, sampled at tr
  starting at lag zero. Use hrf_library() or supply explicit candidates.

- nuisance:

  Optional matching list of time-by-regressor matrices. An intercept and
  linear trend are always included separately for each run.

- valid_rows:

  Optional matching list of logical observation masks; FALSE and NA omit
  rows from training fits without deleting their acquisition times.

- time_offset:

  First acquisition time, scalar or named numeric run vector.

- spatial_groups:

  Optional voxel parcel labels, named in voxel column order. Labels must
  be nonmissing and nonempty. Parcels must be defined independently of
  target-run responses; no implicit whole-brain parcel is constructed.

- pooling:

  Select HRFs per voxel or share one HRF per supplied parcel.

- spatial_strength:

  Nonnegative fixed penalty on voxel candidate distance from the
  training-derived parcel HRF; requires spatial_groups. Zero disables
  the soft penalty. Only used with voxel pooling; not automatically
  tuned.

- fallback_kernel:

  Index of the prespecified fallback candidate, default 1.

- min_r2:

  Minimum mean adjusted partial R-squared for a supported selection,
  default 0.01. This is a configurable heuristic, not a significance
  threshold.

- min_training_runs:

  Minimum usable other runs per voxel, default 1.

- min_parcel_voxels:

  Minimum eligible training voxels per parcel, default 5.

- batch_size:

  Maximum voxel columns per native scoring block.

- return_scores:

  Retain candidate-by-voxel training score matrices per fold.

- amplitude_model:

  Training amplitudes shared within condition (fast default), or free
  for each event record (trial). Trial amplitudes can support HRF
  estimation when signed trial responses average to zero, but cost more
  and require a design with sufficient residual degrees of freedom.

- noise_model:

  Default ar1 uses pooled, grid-REML covariance estimates; ols disables
  whitening. Pooling uses supplied spatial_groups or the run.

- noise_rho:

  Optional fixed AR(1) coefficient, scalar or named in run order.
  Bypasses noise estimation. Requires noise_model = "ar1".

- noise_grid:

  Ordered AR(1) grid within \[-0.95, 0.95\], default -0.8 to 0.9 in
  steps of 0.1. Boundary estimates are explicitly labeled.

- evidence:

  Default calibrated uses conditional GLS F tests with library
  corrections, prespecified-shape comparison, and training-run
  compatibility. Legacy retains the adjusted-R-squared heuristic.
  Conditional tests are exact for known Gaussian AR covariance and a
  correct mean model; estimated noise, non-Gaussian errors, and model
  misspecification require empirical validation.

- evidence_alpha:

  Per-voxel evidence threshold, default 0.05. Corrects the candidate
  search and compatibility comparisons across training runs, not
  multiple testing across the brain. Not a posterior probability.

- canonical_strength:

  Fixed shape-distance penalty toward fallback_kernel in calibrated
  mode, default 0.01. Zero removes this additional shrinkage.

## Value

An object of class crossrun_hrf with normalized kernels, per-target-run
selections and training manifests, parcel anchors, run status/design
summaries, optional candidate scores, voxel/run IDs, and frozen
settings.

## Details

Each target run is excluded from its own HRF selection, including parcel
anchors and fallback decisions. Training fits have run-specific
condition or event amplitudes and nuisance coefficients, with one shared
candidate HRF across conditions. They do not fit unrestricted trial
timecourses or use RPE. Native SVD projections are shared across voxel
blocks. Scores are adjusted partial R-squared after removing nuisance,
averaged equally across usable training runs. Parcel scores average
these voxel scores equally, retaining voxel-specific amplitudes rather
than averaging BOLD first. In voxel mode, selection maximizes score
minus spatial_strength times normalized squared distance to the reliable
parcel anchor. Weak/insufficient voxel fits fall back to a reliable
parcel anchor, then the prespecified fallback kernel. This is a discrete
shape penalty, not Bayesian posterior inference, a local adjacency
prior, or spatial smoothing of neural activity. Library, thresholds,
parcels, and penalty must be fixed externally or tuned within training
data. Calibrated mode additionally penalizes distance from the
prespecified shape. The legacy gate accepted null free-event fits in
colored-noise simulations. Calibrated mode tests task support,
departures from the prespecified shape using pairwise encompassing
models, and compatibility with every usable training run. Fisher
combination uses independent runs; Bonferroni corrections account for
candidate comparisons. Unsupported or inconsistent shapes fall back.
Parcel anchors require min_parcel_voxels individually supported voxels;
voxel p-values are not combined as independent spatial observations. Run
compatibility cannot be assessed with a single training run. There is no
access to target-run vascular shape or guarantee of transfer to that
run. AR covariance is estimated per run using a fixed event-amplitude
mean model and, when residual degrees of freedom allow, its first two
temporal derivatives. Noise likelihood profiles are pooled, not the BOLD
signals. All whitening uses the original acquisition clock, including
gaps. No automatic nested tuning or physiological amplitude calibration
is included. Nonfinite/constant/nuisance-only voxel runs are omitted;
observation masks are common across voxels. Candidate task ranks must
match within a run and leave residual degrees of freedom, otherwise that
run cannot inform HRF selection.
