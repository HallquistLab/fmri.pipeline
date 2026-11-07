# HRF estimation across runs and spatial regularization

This PR adds an opt-in matrix interface for estimating a participant's HRF from
other runs and applying the frozen shapes to continuous linear deconvolution.
It is independent of the alignment changes in PR #44. Existing deconvolution
voxelwise deconvolution defaults and the NIfTI/file wrapper are unchanged. The
matrix interface's noise/evidence defaults are described below. Run this stage before alignment;
it neither truncates epochs nor smooths recovered neural activity across voxels.

## Recommended starting point

Retain a fixed canonical HRF as the baseline and treat learned HRFs as an opt-in
comparison. The repeated validation in [HRF_ESTIMATION_VALIDATION.md](HRF_ESTIMATION_VALIDATION.md)
found benefits with stable HRF mismatch, and harms with low SNR, run-to-run HRF
changes, and neural-duration mismatch. Free-event fitting did not improve noisy
signed-RPE recovery in that experiment, despite succeeding in a noiseless test.

Use a small prespecified library and independently defined parcels when testing
spatial regularization. Hard parcel pooling harmed recovery in heterogeneous
parcels. Neither a universal spatial strength nor parcel size is validated. A
functional parcel need not be a homogeneous vascular territory. Inspect parcel
boundaries, fold stability, and disagreements rather than assuming homogeneity.
The original adjusted-R-squared gate accepted HRFs in every pure colored-noise
voxel/run in that experiment. The defaults now use pooled AR(1) covariance,
conditional GLS evidence tests, a penalty toward the prespecified shape, and
training-run compatibility. See [NOISE_AWARE_HRF_RECOVERY.md](NOISE_AWARE_HRF_RECOVERY.md)
for the change, validation, and limitations. These are conditional evidence tests,
not posterior HRF confidence or whole-brain multiple-testing correction.

```r
library(fmri.pipeline)

# One participant. All runs share the acquisition TR and voxel grid.
# BOLD matrices: time x voxel, retaining all acquisition rows.
# Each matrix needs identical unique voxel column names, in identical order.
bold <- list(runA = bold_A, runB = bold_B, runC = bold_C)
events <- list(runA = events_A, runB = events_B, runC = events_C)
# Each event table: onset (seconds), condition (e.g., cue/response/feedback),
# optional duration (seconds, zero for impulses).
# First acquisition is time zero unless time_offset is supplied.
parcel_ids <- setNames(parcel_for_each_voxel, colnames(bold$runA))
candidates <- hrf_library(tr = TR)

fit <- estimate_crossrun_hrf(
  bold, events, tr = TR, kernel_grid = candidates,
  spatial_groups = parcel_ids, pooling = "parcel",
  nuisance = nuisance_by_run, valid_rows = usable_by_run,
  time_offset = first_acquisition_time_by_run,
  return_scores = TRUE
)

# This is an alternative fit, not another smoothing pass on the parcel fit.
voxel_fit <- estimate_crossrun_hrf(
  bold, events, tr = TR, kernel_grid = candidates,
  spatial_groups = parcel_ids, pooling = "voxel", spatial_strength = 0.02,
  amplitude_model = "trial",
  nuisance = nuisance_by_run, valid_rows = usable_by_run,
  time_offset = first_acquisition_time_by_run
)

# Example lambda and strength are not recommendations for every dataset.
# Choose externally or within training data, without target-run RPE optimization.
result <- deconvolve_crossrun_hrf(
  bold, fit, lambda = 0.01, scaling = "zscore", backend = "cpp",
  batch_size = 2048, operator_cache_mb = 128
)
saveRDS(result, "crossrun_hrf_deconvolution.rds")
activity_A <- result$activity$runA
fit$training_runs$runA      # runB, runC
fit$selections$runA         # which HRF was applied, evidence, and fallback status
fit$parcels$runA            # training-derived parcel anchors and support
```

`nuisance`, `valid_rows`, and `time_offset` are optional. Omit undefined example
objects when not needed. Supplied lists must match BOLD run names and order;
run-specific offsets must be a named numeric vector in that order. A scalar
offset applies to all runs. Native extraction via `extract_nifti_voxel_ts()` returns
voxel-by-time matrices; transpose these for this interface and assign stable
coordinate-derived voxel IDs. Matching dimensions alone cannot establish matching
brain locations. The caller must ensure the runs are registered on the same grid.

Each target run is excluded from its own HRF selection, including parcel anchors,
support gates, and fallbacks. With two runs, A uses B and B uses A; with more runs,
each uses all remaining runs. Runs may have different lengths and event counts.
They must have the same TR and voxel IDs. Run boundaries are never concatenated
for fitting or convolution. Cross-fitted outputs are still statistically dependent
because they share training data; folds are not independent experimental units.

## What the training model estimates

For each candidate HRF, events generate a causal convolution design with one
column per condition by default. Each training run has its own voxel-specific
condition amplitudes and nuisance coefficients. The HRF shape is shared across conditions;
the amplitudes can differ in magnitude and sign. This is not a fit of freely
varying trial timecourses. RPE values are not used by this interface. Include the
relevant event types and durations so an omitted response is less likely to be
absorbed into the estimated HRF. The model assumes amplitudes within a condition
are adequately summarized for estimating its shared response shape; graded trial
variation remains residual variation during this training step.

An important limitation for RPE analyses: a voxel whose response is approximately
proportional to a zero-mean signed RPE may have little average feedback response.
The default condition-mean training model may then provide weak HRF evidence even
when the voxel carries substantial trial variation. `amplitude_model = "trial"`
instead uses one free amplitude per event record and still shares the HRF across
events. This does not force amplitudes to follow RPE or fit a free neural
timecourse within each trial. It can recover HRF evidence when signed responses
average away; a dedicated test covers this case. It also costs more, can overfit
training noise, and can leave too few residual degrees of freedom in dense
designs. Inspect ranks, support and fold stability before preferring it. With
multiple event types per trial, each record has its own amplitude. There is no
automatic switch between the two training models. An independent localizer is
another option; weak condition-mean evidence alone does not imply absent neural
information.

All runs include an intercept and linear trend, plus supplied nuisance columns.
Use consistent preprocessing and nuisance definitions. No task regressor is
silently orthogonalized against other task columns; nuisance is projected out of
both BOLD and task designs. Supplied nuisance can be collinear; SVD determines its
rank after column normalization. Candidate residualized task ranks must be equal
within a run, nonzero, and leave at least two residual degrees of freedom. A run
failing these checks cannot inform selection and is recorded explicitly.

Off-grid impulses deposit unit mass between their two surrounding acquisition
bins. Positive durations use fractional occupancy of bins `[time, time + TR)`.
They represent unit-height inputs; impulse mass and duration occupancy have
different conventions, while condition amplitudes are freely estimated. This is
an acquisition-grid approximation, not evidence of sub-TR neural resolution.
Events before the run are included when their modeled causal response reaches
observations; zero neural history is assumed before the kernel's available
pre-run support. Candidate kernels all have equal lengths, begin at lag zero,
and are normalized to unit signed sum. The built-in nine-candidate library varies
the response-delay and dispersion parameters of the SPM difference of gammas;
the canonical candidate comes first. It does not estimate vascular gain.

For a run, let `SSE0` be residual sum of squares after nuisance, `SSEk` after
nuisance plus the task design for candidate k, `d0` the nuisance-residual degrees
of freedom, and `pk` the added task rank. The training score is:

```
adjusted_partial_R2(k) = 1 - (SSEk / (d0 - pk)) / (SSE0 / d0)
```

Scores are averaged equally over each voxel's usable training runs, rather than
letting a run with larger BOLD variance dominate. Negative scores are retained.
With `noise_model = "ar1"` (default), these sums of squares use whitened data,
task designs, and nuisance designs. Covariance is estimated by pooled grid REML
on each run, using free event amplitudes and temporal derivatives when residual
degrees of freedom permit. Profiles are pooled within supplied parcels, or over
the run if no parcels are supplied. `noise_rho` supplies a fixed coefficient for
validation or externally estimated noise. Inspect `fit$noise`, including boundary
and unidentifiable estimates. No voxel-specific optimizer is used.

The score itself remains a training-fit criterion, not a p-value or held-out
prediction score. Calibrated mode adds separate tests: task versus nuisance,
each alternative's encompassing-model comparison against the prespecified HRF,
and compatibility with each training run. Independent-run p-values are Fisher
combined; Bonferroni corrections account for library comparisons. Compatibility
corrects over competitors and usable training runs. Tests assume a correct mean
model and Gaussian AR covariance; estimated covariance makes their calibration
approximate. `evidence_alpha` is per voxel, not across the brain. One training run
cannot establish stability. Neither target-run HRF stability nor absence of task
signal is established by a fallback. No RPE is used to choose shapes or settings.

For the original scoring/application behavior, explicitly request
`noise_model = "ols", evidence = "legacy"`. `canonical_strength` applies only in
calibrated mode. There is no automatic nested tuning.

`valid_rows` is a logical mask per training run. FALSE and NA omit acquisitions
from training fits while retaining their original times. Nonfinite nuisance is
allowed only on omitted rows. A voxel with nonfinite BOLD on a retained row,
constant BOLD, or no remaining variance after nuisance cannot inform that run's
HRF fit. It may still borrow a supported parcel estimate from other voxels. See
`signal_status`, `designs`, and per-selection training counts.
Runs with fewer than three usable rows are recorded as
`insufficient_observations` and omitted from training evidence; they do not
prevent other runs from being fitted or their own HRFs from using other runs.

## Spatial regularization and fallback

Parcel pooling aggregates candidate scores, not BOLD timeseries. Consequently,
oppositely signed or differently scaled voxel amplitudes do not cancel before
fitting. Each eligible voxel contributes its mean training score equally. Only
voxels with at least `min_training_runs` usable other runs enter the parcel score.
In calibrated mode, an anchor additionally needs at least `min_parcel_voxels`
individually supported voxels for that shape. Their p-values are not combined as
independent spatial observations. Strong evidence against an anchor in a voxel's
training runs prevents borrowing it.

- `pooling = "parcel"` applies that anchor to the whole parcel. This is a hard
  shared-shape constraint, with freely varying voxel amplitudes during training.
- `pooling = "voxel", spatial_strength = 0` selects each voxel's supported HRF
  independently. If parcel labels are supplied, they also enable parcel fallback.
- `pooling = "voxel", spatial_strength > 0` selects the candidate maximizing:

```
mean_training_score(k, voxel)
  - canonical_strength * distance(k, fallback_kernel)
  - spatial_strength * distance(k, parcel_anchor)
distance(k, a) = sum((h_k - h_a)^2) / (sum(h_k^2) + sum(h_a^2))
```

Distance compares actual unit-sum shapes, not arbitrary candidate indices. This
is a discrete penalized estimator centered on a training-derived parcel anchor.
It is not a posterior distribution or an implementation of PyHRF/JDE/HemoLearn.
It has no adjacency graph, cross-parcel coupling, or HRF averaging after selection.
Increasing the strength favors the parcel shape; it never averages the target
voxel's BOLD or neural estimate with neighboring voxels. With no reliable anchor,
the penalty is disabled and this is visible in the selection status.

A selected departure must pass task support, prespecified-shape comparison, and
training-run compatibility, as well as `min_r2`. The score penalty toward the
prespecified shape defaults to 0.01; this is a fixed discrete penalty, not a
Bayesian prior or a universally optimized strength. Insufficient voxel evidence
falls back to a supported, non-rejected parcel anchor, then the prespecified
`fallback_kernel` (default first candidate). With no common compatible library
shape, use the prespecified shape and report inconsistency/misspecification.
Parcel mode also falls back when its anchor is unsupported. Ties choose the first
candidate. `score_margin` remains descriptive. Inspect `voxel_evidence`, selected
task/shape p-values, compatible/supported shape counts, and `stability`.
Borrowed parcel shapes do not turn an unidentified voxel into a supported voxel.
The legacy mode retains only the original R-squared/count heuristics.

Inspect the maps of chosen shapes, parcel disagreements, low-support voxels, and
fallback counts. A parcel mixing distinct hemodynamics can bias either sharing
option. Define parcels independently of target-run data, or learn them inside
each training fold before supplying them. The interface accepts one fixed parcel
map; it does not learn fold-specific parcellations or choose a spatial strength.

Published precedents distinguish several forms of spatial regularization:

- [PyHRF/JDE](https://pmc.ncbi.nlm.nih.gov/articles/PMC3989728/) shares HRFs within
  parcels and combines temporal HRF regularization with spatial modeling of
  response levels. That is more extensive joint inference than this interface.
- [Chaari et al., parcel-free JDE](https://mistis.inrialpes.fr/~forbes/PAPERS/SFDS2012_CHAARI.pdf)
  describes voxel-varying HRFs tied to HRF groups with spatially regularized
  group membership. This motivates allowing local departures rather than
  assuming every voxel in an anatomical parcel has identical hemodynamics.
- [Cherkaoui et al., 2021](https://doi.org/10.1016/j.neuroimage.2021.118418) uses
  hemodynamic parcels, constrained HRF dilation, and a multivariate semi-blind
  model. Adopting that full model would be a separate research implementation.
- [GLMsingle](https://pmc.ncbi.nlm.nih.gov/articles/PMC9708069/) demonstrates a
  constrained HRF library for voxel-level fits. Its HRF selection uses fit rather
  than cross-validation; this PR adds an explicit exclusion of the target run.

## Frozen inverse, scaling, and computational cost

`deconvolve_crossrun_hrf()` requires an explicit lambda. It never reruns HRF-grid
selection on the target data. Candidate/library choices, spatial strength, and
inverse lambda should be prespecified, learned externally, or selected within
training data. With two runs, choose conservative fixed settings rather than
claiming nested validation from one training run. Existing within-run GCV is not
an independent test of the learned HRF or neural recovery.

The default inverse is an Rcpp masked GLS implementation. The covariance from
each run is frozen and used for its inverse; using a target run to estimate its
noise does not use its responses to select its own HRF. Missing observations are
integrated out of the AR(1) covariance: a gap of d acquisitions uses rho^d and
its accumulated innovation variance. `observation_mask = "fit"` is the default
for AR/calibrated fits. Masked BOLD values do not enter scaling or fitting, but
latent activity is estimated at all original acquisition times. The activity's
`observation_mask` attribute distinguishes observed and modeled samples.

Inverse inputs remain assumed to be consistently prepared: `nuisance_mode =
"ignore"` is the default. Optional `"joint"` fits unpenalized nuisance coefficients
by projecting both whitened convolution and BOLD designs. A free nuisance trend
can amplify low-frequency uncertainty in a continuous inverse; this option harmed
recovery in the development pilot and should not be enabled automatically.

R provides an independent GLS reference. Scalar Wiener remains available only for
OLS fits with `nuisance_mode = "ignore", observation_mask = "all"`; it does not
silently approximate masked GLS. For cpp/r, each distinct HRF/rho/design yields one
operator, applied using block BLAS. There is no factorization per voxel. Caching
reuses operators only when run lengths, masks, and nuisance designs match.
Operators are reused within a bounded 128-MiB FIFO cache; change
`operator_cache_mb` to control this budget. Changing run length, mask, or
nuisance design clears that cache. A current operator and factorization workspaces can
exceed the cache budget: the existing dense inverse has quadratic memory and
cubic factorization cost in run length. Training uses shared design factorizations and batched projections, including
pairwise encompassing spaces for the evidence tests and grid REML for noise.

`batch_size` controls temporary voxel-block memory, not estimator accuracy.
Calibrated training stores scores and three evidence matrices, approximately
`4 * 8 * runs * candidates * voxels` bytes, plus input BOLD, noise profiles,
manifests, shared bases, and fold workspaces. Legacy training stores scores only. `return_scores = FALSE` avoids keeping
the additional per-fold candidate-score tables in the returned fit. Native code
checks interrupts between voxel blocks. No new package dependency, extra process
pool, or nested OpenMP layer is introduced. Control BLAS threads at the job level
to avoid oversubscribing scheduled jobs.

Scaling is applied once to measured BOLD, before optional zero padding. Z-scoring
remains the default for standardized within-run responses; use a justified PSC or
demeaning convention when between-run units matter. Centers/divisors are returned
per voxel/run. Unit-sum HRFs and these transformations do not identify an absolute
physiological neural scale. `tail_padding = 0` is the finite-run default; setting it
to the kernel length reproduces the existing voxelwise wrapper's zero-tail
assumption, rather than correcting that assumption.

Inverse masks now follow the explicit observation policy above. Supplied nuisance
coefficients are not silently transferred from training to target runs. Nonfinite
or constant signals on usable rows, or unidentifiable AR covariance, remain
missing. Filling censored BOLD with zeros is not a substitute for the masked
forward model. Modeled activity in a gap does not restore a measured acquisition;
alignment must retain the acquisition mask.

The returned fit records normalized kernels/parameters, run/voxel IDs and lengths,
training run IDs for each target, design ranks, training signal statuses, selected
shapes, parcel anchors, evidence/fallback fields, and frozen settings. It also
retains the supplied event tables, nuisance matrices, masks, and parcel labels in
`specification`. Save this object with the activity/scaling output for review.
The BOLD arrays themselves are not copied into the fit; archive input provenance
alongside the saved model. The application checks identity/order/dimensions, not
an input-content checksum. Changing inputs requires a new fit.

## Validation and review artifacts

The repeated comparison against a fixed canonical HRF is documented in
[HRF_ESTIMATION_VALIDATION.md](HRF_ESTIMATION_VALIDATION.md). It uses an
independent continuous-time generator, 40 datasets per scenario, voxel-level
recovery measures, paired Monte Carlo intervals, and all six methods at three
fixed inverse penalties. It includes controls where learning can hurt. Use that
report to assess recovery benefits; the single-seed example below illustrates
the interface and is not sufficient evidence for choosing a default.

After installing this checkout, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
  inst/examples/crossrun_hrf_review.R /tmp/crossrun-hrf-review
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
  inst/benchmarks/crossrun_hrf.R /tmp/crossrun-hrf-benchmark 10000 3
```

The review script writes the library, per-fold voxel selections, parcel anchors,
design and training-signal manifests, frozen models, simulated aligned epochs,
graded-RPE recovery summaries, and a plot. It compares canonical and learned
HRFs using known generating shapes and a prespecified graded input. The example
does not select settings by maximizing its RPE association. It also writes
`sensitivity_recovery.csv` and `sensitivity_models.rds` for off-library generating
HRFs with AR(1) noise and heterogeneous parcels, comparing hard pooling and the
soft voxel penalty under frozen settings. These are controlled simulations, not
participant-data validation or an uncertainty analysis. The example's alignment
is simple grid interpolation for plotting and
does not depend on PR #44's alignment API.

The benchmark compares scores with independent per-voxel `lm.fit` models on a
subset, then times the native implementation and all-run selection/application
at the requested voxel count. Its R baseline is a conventional fitting loop, not
an optimized R projection algorithm. Compilation/loading are excluded; environment
and BLAS-thread settings are saved with results. See the PR for measured results.

The following historical measurements use the original OLS/heuristic workflow;
new covariance/evidence costs and validation are in
[NOISE_AWARE_HRF_RECOVERY.md](NOISE_AWARE_HRF_RECOVERY.md).

Measured on this host with R 4.5.1, OpenBLAS 0.3.30, one BLAS/OpenMP thread,
and median elapsed time over three repetitions (compilation/loading excluded):

| Operation | Data | Seconds |
| --- | --- | ---: |
| Conventional per-voxel `lm.fit` scoring | One run, 300 volumes, 300 voxels, nine HRFs | 0.322 |
| Native batched scoring, identical models | Same | 0.002 |
| All cross-run selections, condition amplitudes | Three runs, 300 volumes each, 50,000 voxels, nine HRFs | 1.721 |
| All cross-run selections, free event amplitudes | Same; 56 event amplitudes per run | 2.530 |
| All frozen cpp inverses, including scaling/operator construction | Three runs, 300 volumes each, 50,000 voxels | 17.072 |
| All frozen scalar Wiener inverses | Same | 33.014 |

Native scores differed from independent least-squares scores by at most
`1.45e-15` in the comparison. The speedup reflects shared fitting work as well as
native code; an optimized R projection implementation would be a different
baseline. Wiener was slower than cpp in this workload despite its lower
asymptotic application cost. Run length, event-model rank, selected-kernel count,
hardware, and BLAS configuration affect these results. No candidate pruning,
voxel subsampling, or approximate projections were used to obtain these timings.

Tests compare native designs with independent causal convolution, native scores
with masked least squares, grouped inverses with direct existing-API fits, and
the spatial objective with its explicit formula. They check target-run exclusion,
fallbacks, rank deficiencies, nuisance scaling, two-run application, identities,
missing values, unequal run lengths, voxel blocks, cache budgets, and zero-mean
signed trial responses. An isolated package installation and **43 focused test
cases / 761 assertions** passed, with no test failures, errors, warnings, or skips.
The focused checks include the existing linear, voxelwise, and native NIfTI
extraction regressions. R parsing, Rd checks, and whitespace checks passed. The
review script generated the documented artifacts. The full unrelated package
suite and participant-data validation were not run.

In the matched-library example, all 160 voxel/fold HRF choices matched their
generating candidate. Mean onset-response correlations with injected RPE were
0.841 versus 0.826 for canonical HRFs in the early parcel, and 0.809 versus 0.751
in the late parcel. Recovered onset slopes remained 0.129--0.177 for a known input
slope of 0.4: correct shape selection does not remove regularization attenuation
or finite-grid effects. In the heterogeneous-parcel/AR(1)/off-library example,
the soft voxel model reduced mean HRF-shape MSE from `5.41e-4` under hard pooling
to `1.29e-4`. These are descriptive results for the fixed example and settings,
not estimates of expected participant-data performance.

Further scientific validation should broaden these mismatch/noise/parcel tests
and include neural latency/duration variation, actual schedules, and evaluation
of any parameter-selection procedure.
A delayed neural response and delayed vascular response remain difficult to
separate. Estimating on another run does not remove that ambiguity. Compare
inference across plausible HRFs and preprocessing assumptions, and avoid claiming
EEG-like temporal resolution from interpolated BOLD inverses.
