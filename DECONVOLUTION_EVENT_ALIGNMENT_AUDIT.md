# Voxelwise deconvolution and event alignment audit

Date: November 6, 2026

Package: `fmri.pipeline`

Reviewed checkout: `c6672b2e2264829550a84a1c1e2a4392aeed5ba2`

## Purpose and overall assessment

The scientific goal is to estimate trial-level time series of latent neural activity from BOLD data, analogous to event-locked analysis in EEG, while retaining sensitivity to graded inputs such as signed reward prediction errors (RPEs). This audit covers the Bush/Cisler algorithms, regularized linear deconvolution, event alignment, overlap censoring, and computational performance.

The regularized linear method is the best development target for this goal. Its signed, unbounded activity estimates provide a natural parameterization for graded increases and decreases, and its linear algebra offers practical opportunities for acceleration. However, correctness problems in the surrounding pipeline should be repaired before interpreting trial-level results.

Basic off-grid event alignment passed a synthetic timing check. The current censoring correctly applies nearest-event boundaries to source samples, but it does not guarantee nonoverlapping trial windows or remove all contributions from neighboring events.

The original audit was a review and diagnostic exercise, with no package implementation changes at that stage. Subsequent remediation is recorded below. Synthetic results establish specific behavior under the tested conditions; they do not establish performance on participants' data.

## Scope and verification

### Remediation update: first two P1 findings

The first two P1 findings have been addressed in the working tree following this
audit. The findings below preserve the original reviewed behavior and line
references for historical context.

- Voxelwise extraction now uses an internal RNifti/Rcpp function,
  `extract_nifti_voxel_ts()`, with the exact same masked linear voxel indices as
  the metadata. The `3dmaskdump` dependency, temporary text extraction, and AFNI
  coordinate-order conversion have been removed. Atlas, mask, and input spatial
  geometry are checked, and output coordinates use the checked affine.
  `afni_dir` is accepted for compatibility and ignored.
- Constant and nonfinite signals are screened before scaling and excluded from
  fitting and tuning. Invalid signals retain their output positions with `NA`
  estimates. Voxelwise CSVs include `decon_status`; standalone linear diagnostics
  include `signal_status`. The status values are `ok`, `constant`, and
  `nonfinite`. Partial missingness excludes the entire signal rather than
  deleting or interpolating observations. All-invalid batches return missing
  estimates without invoking a solver.
- Regression coverage includes native voxel ordering, compressed and uncompressed
  input, integer scaling, nonfinite values, single-voxel output, incompatible
  grids, all-invalid batches, and isolation of valid estimates and tuning across
  the C++, R, and Wiener backends.

The package compiled and installed successfully into an isolated temporary R
library. Focused installed-package tests passed for linear deconvolution,
voxelwise deconvolution, native extraction, and the shared ROI geometry checks.
An installed-package integration check with two local PSOCK workers also passed,
including correct masked labels and preservation of invalid rows. Local worker
sockets required execution outside the filesystem/network sandbox.
The regression tests also cover the exact `lambda` lookup required to use a
`lambda_grid` without R's list partial matching treating the grid as a fixed
lambda.

Internal extraction holds the native image buffer while creating the selected
voxel matrix. It avoids a full double-precision R image copy and text
serialization; it is not a streaming reader. The audit's timing tables below
describe the original deconvolution benchmarks, not a benchmark of this change.


### Remediation update: explicit normalization and Wiener tuning

The two P2 findings below are now addressed. Their original reproductions are
retained as historical findings, rather than descriptions of the updated code.

**Recommended scaling.** Z-scoring remains the default for standardized trial
responses. The voxelwise Bush/Cisler paths require it: the fixed learning rate,
sigmoid encoding, and standardized residual resampling do not provide a calibrated
arbitrary-unit estimator. Standalone Bush 2011 retains its existing expert
`normalize = FALSE` option. Bush resampling now defaults to normalizing its
observations consistently for both the base fit and residual construction; its
new `normalize = FALSE` option accepts caller-standardized data. The voxelwise
Bush 2015 path also forwards the requested `beta`.

Linear deconvolution accepts explicit `scaling = "zscore"`, `"demean"`, `"none"`,
or `"psc"`; the voxelwise wrapper accepts the same field in `decon_settings`.
Explicit scaling overrides the legacy flags. Otherwise `normalize = FALSE`
selects demeaning by default, or no transformation when `demean = FALSE` too.
For the fixed linear inverse, multiplying the input does not change the
conditioning of the inverse matrix. Regression checks across gains from `1e-4`
to `1e6` verify proportional estimates and unchanged globally selected lambda
when the same gain is applied to all signals. Relative gains between signals can
change global tuning weights. Z-scoring is therefore a sensible standardized
analysis default, rather than a numerical requirement for this linear solver.

Use demeaning to retain a consistent input scale in residualized data; use PSC
only with a meaningful positive BOLD baseline. PSC rejects nonpositive and
near-zero means relative to signal spread. Raw (`none`) input is mainly for
already-centered data or simulations. Neither preserving BOLD units nor
standardizing identifies a physiological neural amplitude or removes penalty
attenuation. Z-scoring retains relative trial variation within a voxel/run but
changes comparisons across runs and voxels.

The wrapper applies its selected transform once to measured volumes before the
existing tail padding. It disables subsequent solver normalization and linear
solver demeaning. This intentionally removes the previous double scaling/centering
of padded data; regenerate existing files to adopt the new behavior. Per-voxel
`.scaling.csv` sidecars record mode, removed center, divisor, and signal status.
Identity sidecars additionally record algorithm, scaling, output units, and the
scaling filename. Standalone linear matrices/diagnostics record the same transform.
Original-series outputs expose `BOLD_input`; `BOLD_z` remains an alias only when
z-scoring actually occurred. These are descriptions of the BOLD transformation;
the removed BOLD mean is not an identified neural baseline.

**Wiener tuning.** Scalar Wiener GCV now uses the exact trace of its finite
observed smoother, tested against dense responses to observed unit impulses,
including short runs, different FFT lengths, and fixed AR(1) whitening settings.
If `h` is the implemented causal forward kernel and `f` the circular impulse
response of its inverse filter, the trace is
`sum((n_time - lag) * h[lag] * f[-lag])` over supported nonnegative lags shorter
than the run (zero-based indices; inverse indices wrap within the FFT).
This counts only observed diagonal entries after extension and truncation,
rather than summing padded Fourier gains. The audit reproduction now gives
approximately `40.32813`, matching the directly constructed smoother.

This trace conditions on preprocessing and the AR(1) coefficient. Adaptive PSD
filters also depend on the observations through their estimated spectra: holding
spectra fixed does not give their full Jacobian degrees of freedom. Consequently,
PSD fits require one lambda and one HRF; automatic GCV selection now fails with
an actionable error. Their diagnostics label `df_method = "conditional_fixed_psd"`
and return missing GCV values. Scalar diagnostics label
`df_method = "finite_observation_trace"`. Fixed-parameter PSD activity is preserved.

**Accuracy-preserving acceleration.** The C++ backend constructs the exact
finite-difference penalty crossproduct from its short stencil, including boundary
entries, instead of multiplying dense difference matrices. A single fixed
candidate with matrix-only output skips RSS, trace, and GCV calculations. Native
scaling is retained for ordinary C++ modes to avoid extra R matrix copies.
Wiener fitting caches input/kernel FFTs across candidate lambdas and skips
unrequested scoring for fixed-parameter fits. Fixed estimates are compared with
the original implementation; C++ tuning and diagnostics are also checked against
the independent dense R reference across penalties, HRFs, and whitening settings.
Explicit dense Cholesky reuse and combined-RHS reuse were benchmarked and rejected
because they slowed these problems; the existing band-aware solver is retained.
The larger sparse/banded redesign remains separate work.

The included `inst/benchmarks/deconvolution.R` compares the current checkout
against baseline `bfbffcb`. It excludes compilation/loading, uses 200 signals,
a canonical HRF at TR = 1, default z-scoring, fixed lambda `0.01` or ten candidate
lambdas, and reports median elapsed time. Set `OPENBLAS_NUM_THREADS=1` and
`OMP_NUM_THREADS=1` before starting R. These measurements concern computation,
not participant-data validation or improved scientific recovery. Event overlap
policies, gap preservation, artificial boundary observations, amplitude-calibration
simulations, and Bush uncertainty/convergence issues remain separate audit items.

Median seconds from five repetitions in this environment:

| Timepoints | Backend | Fixed: before / after | Tuned: before / after |
| --- | --- | --- | --- |
| 300 | C++ | 0.019 / 0.015 | 0.076 / 0.091 |
| 1,200 | C++ | 0.207 / 0.113 | 1.136 / 1.152 |
| 300 | Wiener | 0.052 / 0.035 | 0.384 / 0.327 |
| 1,200 | Wiener | 0.199 / 0.111 | 1.558 / 1.428 |

These runs show the clearest benefit for fixed-parameter fits. Tuned C++ did not
improve in this benchmark (and was slower for the smaller problem); no universal
speedup is claimed. All activity differences in these particular comparisons
were zero, including the tuned Wiener examples, which happened to select the
same candidate. Corrected Wiener GCV can select a different candidate in other
data; dense-oracle scoring tests establish the corrected selection behavior.


**Validation:** The updated package compiled and installed in an isolated library.
Focused installed-package checks passed: 56 test cases and 560 assertions, with
no failures, errors, or warnings. Coverage includes scaling and gain invariance,
finite-smoother trace and GCV selection oracles, C++/R agreement, fixed-fit fast
paths, invalid signals, voxel metadata, identity, alignment, and ROI geometry.
A local two-worker integration check matched serial estimates and scaling
sidecars for z-score/C++ and PSC/Wiener runs. Bush wrapper checks cover single
voxels, normalization, retained default settings, and `beta` forwarding; its
resampling-input check isolates the base/residual normalization contract.
Changed R/Rd files parse and Rd checks pass, as does `git diff --check`.
The full unrelated package suite and participant-data validation were not run.


### Remediation update: overlap support and fixed-operator diagnostics

See [EVENT_ALIGNMENT_OVERLAP_GUIDE.md](EVENT_ALIGNMENT_OVERLAP_GUIDE.md) for the
implementation contract, changed statistics, field dictionary, limitations, and
runnable review example. Alignment now preserves explicit missing intervals and
acquisition-clock gaps, interpolating each unit before aggregation. Per-lag
contributor counts and acquisition-time memberships make changing support visible.
Matrix interpolation on common grids is checked against per-unit interpolation;
ragged grids use that fallback, processing one output group at a time.

Overlapping epochs are retained by default when no collision settings were
specified. Explicit legacy collision settings continue to select truncation;
`allow` can override them, while `exclude` removes whole intersecting requested
epochs. Neighbor/event/edge annotations and a source-clock mapping support audits.
These policies manage data selection and shared acquisition times, not physiology.

`deconvolution_impulse_diagnostics()` computes conditional recovery and entry into
other event windows using fixed HRF, lambda, whitening, scaling and padding choices.
Missing/outside/unobserved impulses remain explicit. Off-grid impulses use linear
acquisition-bin deposition, not an assertion of sub-TR neural resolution. Directed
window area ratios can exceed one and are not attribution probabilities. The tool
supports C++, R, and scalar Wiener; no adaptive PSD, Bush, noise simulation,
automatic tuning, or observation-weighted refit is implied.

Joint event-informed estimation, observation masks inside the forward model,
removal of artificial boundary observations, and realistic full-pipeline RPE/noise
validation remain separate work. Post-hoc alignment masks cannot undo artifacts
that already entered deconvolution. Validation passed 74 focused installed-package
test cases and 708 assertions without failures/errors/warnings, plus serial versus
two-worker checks for support, overlap annotations, and source clocks. The review
script generated policy tables and centered/uncentered impulse comparisons. No
participant-data validation or full unrelated package suite is claimed.

Primary files reviewed:

- [R/voxelwise_deconvolution.R](R/voxelwise_deconvolution.R)
- [R/deconvolve_reglin.R](R/deconvolve_reglin.R)
- [src/deconvolve_reglin.cpp](src/deconvolve_reglin.cpp)
- [src/deconvolve_nlreg.cpp](src/deconvolve_nlreg.cpp)
- [R/deconvolve_nlreg_resample.R](R/deconvolve_nlreg_resample.R)
- [R/event_lock_timeseries.R](R/event_lock_timeseries.R)
- [R/fmri_ts.R](R/fmri_ts.R)
- [tests/testthat/test-deconvolve_reglin.R](tests/testthat/test-deconvolve_reglin.R)
- [tests/testthat/test-voxelwise_deconvolution.R](tests/testthat/test-voxelwise_deconvolution.R)

### Remediation update: remaining P1 findings

All five P1 findings in this audit have now been addressed in the working tree.
The original findings and line references below describe the reviewed checkout;
P2 remediation is recorded above; broader methodological recommendations remain
separate from the P1 fixes.

- **Missing event times:** Target-onset validation precedes collision searches,
  which use finite candidates only. Trials with missing onsets and unsupported
  windows retain explicit missing epochs. Original event variable mappings are
  preserved through copying, time remapping, and trial keying.
- **Interpolation grids:** Zero-width, one-sided, and sub-sample windows have
  unique times within their bounds, including zero. An exact singleton sample
  is retained without extrapolation. Duplicate reshape keys trigger an error
  instead of count aggregation.
- **File discovery and identity:** Default output names start with `decon_` and
  are visible. New files have `.metadata.csv` sidecars containing the source
  NIfTI path and supplied `id`, `run`, and optional `session` metadata. Discovery
  includes historical hidden outputs while excluding sidecars. Explicit
  `file_manifest` data override automatic identity resolution; legacy/BIDS
  filename or source-path parsing remains available as a fallback. Exact string
  matches preserve leading zeros and alphanumeric identifiers. Ambiguous
  participant/run/session matches fail rather than combining runs silently.
- **File failures:** Every input receives a status row. Alignment stops before
  writing combined data if any file fails, unless `fail_on_error = FALSE` is
  explicitly supplied. Partial output then carries a warning. The returned
  `file_status` attribute records failures and row counts; when an output path
  is supplied, the manifest is written to `paste0(out_file, ".files.csv")`,
  including on failure. The default error condition also contains the manifest.

Interpolated output grouped by trial includes `epoch_status`: `ok`, `partial`,
`insufficient_data`, `missing_onset`, `outside_run`, `censored`, or `no_data`.
These describe timing and output support. Nearest-event truncation still does
not guarantee disjoint windows or physiological isolation, and interpolation
across interior gaps has not been redesigned by these P1 fixes.

For custom filenames without embedded identity, provide `id`, `run`, and optional
`session` in `add_metadata` when deconvolving, or supply a manifest at alignment:

```r
manifest <- data.frame(
  file = decon_files[1L],
  id = "001",
  run = 1L,
  session = "B"
)
aligned <- evt_align_decon_files(
  decon_files[1L], trial_df,
  alignment = list(evt_col = "feedback_onset", time_before = -2, time_after = 4),
  tr = 1, ncpus = 1, file_manifest = manifest
)
attr(aligned, "file_status")
```

Verification: the updated package installed into an isolated temporary library.
All 45 focused tests passed, with 341 passing assertions, zero failures, and zero
test warnings. Coverage includes all five P1 findings, the existing linear
deconvolution and native extraction tests, and shared ROI geometry tests. A
two-worker installed-package alignment check passed with successful files,
missing-onset trials, and an intentionally failed file reported in the manifest.
The full package suite and participant-data validation were not run.

Original audit verification used R 4.5.1, freshly compiled copies of the current deconvolution C++ sources, synthetic signals, and an AFNI extraction test with small NIfTI images. The AFNI test reported version `AFNI_26.0.00`.

The nine existing linear-deconvolution test cases and the wrapper's retired-backend rejection test passed in the final focused harness. The harness sourced the checkout, redirected package-qualified helper references to those sourced helpers, and explicitly selected testthat edition 3. Earlier harness runs used the wrong test edition and reported two attribute-only comparison failures; the final edition-3 run passed. A separate numerical C++/R comparison differed by approximately `7.1e-15`.

The tests available at the original audit did not cover the alignment failures, mask-label mismatch, and other integration problems documented below. Subsequent regression coverage is recorded above. The full package test suite and real-data validation were not run.

## Correctness findings

Priorities: **P1** denotes a high-priority correctness issue that can mislabel, invalidate, corrupt, or omit results. **P2** denotes a consequential estimation or configuration issue.

### 1. P1: An additional mask can assign signals to incorrect voxel and parcel labels

**Location:** `R/voxelwise_deconvolution.R`, approximately lines 123–146, 185–192, and 270–275.

The metadata are constructed from the atlas after multiplication by the additional mask. The AFNI extraction command still uses the original atlas file. Extracted signals are then joined to metadata by a sequential voxel number rather than by the returned voxel coordinates.

Consequently, the extracted and metadata rows can describe different voxels.

**Reproduction:** An eight-voxel atlas labeled 1 through 8 was paired with a mask retaining only voxels 2 and 5. Extraction returned all eight voxels. The first signal received label 2, the second received label 5, and the remaining six received missing labels. The first labeled signal correlated `1.0` with voxel 1's generated signal and approximately `-0.3301` with voxel 2's signal.

**Recommendation:** Apply the same voxel selection to extraction and metadata, and join using explicit voxel coordinates after reconciling indexing conventions. Validate image dimensions and spatial transforms before extraction. Check that every extracted voxel maps to exactly one metadata row.

### 2. P1: A constant or nonfinite voxel can invalidate the entire linear batch

**Location:** `R/voxelwise_deconvolution.R`, approximately line 195; `R/deconvolve_reglin.R`, input validation and global candidate selection; `src/deconvolve_reglin.cpp`, approximately lines 207–238.

The wrapper unconditionally z-scores every voxel. Constant signals therefore become nonfinite. The linear solver does not reject or isolate these columns before global GCV selection. Nonfinite scores can prevent any candidate from being selected, leaving valid signals missing as well.

**Reproduction:** A matrix with one varying signal and one constant signal, after wrapper-style scaling, returned no finite activity estimates in either column. A valid signal paired with an entirely missing signal also invalidated the tuned batch.

**Recommendation:** Screen signals before scaling, fit valid columns, retain an explicit validity/status map, and restore invalid columns afterward. Define separate behavior for constant columns and partially missing observations. Invalid columns must not participate in global tuning.

### 3. P1: Missing event times can abort alignment

**Location:** `R/event_lock_timeseries.R`, approximately lines 60–75 and 99–103.

Collision searches execute before the missing target-onset check. Comparisons containing missing times can make `if (any(...))` fail. Missing collision times can also propagate into the selected boundary. A separate list-management error occurs when a missing final target onset is skipped and the code assigns more names than the trial list contains.

**Reproduction:** Missing target or collision onsets produced `missing value where TRUE/FALSE needed`. A missing final onset without collision settings produced a list-name length mismatch.

**Recommendation:** Validate the target onset first; restrict collision candidates to finite values; use stable, preallocated trial bookkeeping; and retain a trial-level exclusion reason instead of failing the run.

### 4. P1: Zero or sub-sample windows can replace signal values with counts

**Location:** `R/event_lock_timeseries.R`, approximately lines 280–288 and 311.

The interpolation grid uses colon expressions such as `1:floor(time_after/output_resolution)`. In R, `1:0` is a descending sequence rather than an empty sequence. This creates duplicate times and times outside the requested window when either side is zero or shorter than the output resolution.

Duplicate trial/time combinations subsequently cause `dcast()` to use counts rather than the interpolated signal values.

**Reproduction:** Both a `[0, 0]` window and a `[-0.2, 0.2]` window at one-second output resolution returned unique times `-1, 0, 1`, accompanied by duplicate-cell aggregation warnings.

**Recommendation:** Construct positive and negative grid branches with explicitly empty sequences when their lengths are zero. Assert that times are unique and lie inside the requested window. Reject duplicate reshape keys rather than accepting count aggregation.

### 5. P1: Default output naming and downstream discovery are incompatible

**Location:** `R/voxelwise_deconvolution.R`, approximately lines 154–156; `R/event_lock_timeseries.R`, approximately lines 532 and 613–617.

The default naming expression replaces path separators with periods. An absolute input path consequently produces a hidden filename beginning with a period. The alignment stage's default `list.files()` excludes such files.

Even when files are discovered, alignment derives identity from filenames matching numeric `sub..._run...` patterns, despite the upstream function allowing arbitrary output expressions. The default path-derived names do not satisfy that parser. This also limits support for alphanumeric participant identifiers and session information.

**Reproduction:** A default-style filename was absent from ordinary file discovery but present with `all.files = TRUE`; its parsed subject and run were both missing.

The alignment loop additionally uses `.errorhandling = "remove"`, allowing worker errors to omit files without a complete returned failure manifest.

**Recommendation:** Use an explicit file manifest containing subject, session, run, and output path. Choose visible default filenames. Return successes, exclusions, and errors together, and verify expected versus observed files before combining outputs.

### 6. P2: The voxelwise wrapper overrides amplitude-preservation settings

**Status:** Addressed in the normalization/Wiener remediation update above. The following describes the original finding.

**Location:** `R/voxelwise_deconvolution.R`, approximately lines 195 and 235.

The wrapper standardizes input before passing `decon_settings$normalize` to the linear solver. Setting that option to `FALSE` cannot preserve the original amplitude scale.

**Reproduction:** Standalone linear deconvolution with normalization disabled preserved a twofold difference between two signals. Wrapper-style z-scoring reduced the output scale ratio to one.

Within a run and voxel, this standardization does not itself remove graded trial variation: it applies a common scaling factor. It does remove between-run and between-voxel scale differences and changes the interpretation of parcel averages. Preserving BOLD units also does not, by itself, identify a common physiological neural scale across voxels.

**Recommendation:** Make preprocessing algorithm-specific and explicit. Support raw/demeaned, percent-signal-change, and standardized inputs as appropriate. Save scaling factors and estimation settings with the outputs.

### 7. P2: Wiener GCV uses incorrect effective degrees of freedom

**Status:** Addressed in the normalization/Wiener remediation update above. The following describes the original finding.

**Location:** `R/deconvolve_reglin.R`, approximately lines 643–662.

The degrees-of-freedom calculation sums frequency-domain gains over the padded FFT without correctly accounting for the actual observed and truncated smoothing operator. Automatic regularization selection is therefore affected.

**Reproduction:** For 80 observations, kernel `c(0.05, 0.25, 0.45, 0.2, 0.05)`, and lambda `0.1`, the reported degrees of freedom were `66.11319`. The trace of the actual smoother, obtained by applying deconvolution and reconvolution to an identity matrix, was `40.32813`.

**Recommendation:** Derive the trace for the implemented finite-observation operator and verify it numerically on small problems. Treat scalar and data-adaptive PSD shrinkage separately; the latter requires additional care because its filter depends on the data. Keep automatic Wiener tuning experimental until validated.

## Event alignment and overlap censoring

### What passed

The basic time convention is internally consistent: volume 1 is assigned time zero, `time_offset` shifts that clock, and event onset is subtracted before interpolation. A synthetic ramp aligned to noninteger onsets returned exactly the expected absolute-time values, with maximum error zero.

The nearest-event logic excludes retained source samples at or beyond its selected boundaries. Linear interpolation does not extrapolate outside the retained support. Removing source samples before interpolation is conservative and can lose usable output near a boundary when no valid bracketing samples remain.

For actual acquisitions, independently verify that behavioral onsets, discarded initial volumes, and the slice-timing reference use the same clock. These conventions cannot be established from a synthetic timing check alone.

### What the current censoring does not guarantee

1. **Exclusive samples across trials.** Nearest-event truncation is not a disjoint-window policy. With events at 0 and 6 seconds and requested windows of plus/minus 4 seconds, samples between 2 and 4 seconds can occur in both epochs. The diagnostic with onsets 4.5, 10.5, and 16.5 seconds reproduced shared samples at times 7, 8, 13, and 14 seconds.
2. **Removal of neighboring-event contributions.** A regularized inverse filter spreads an event's influence in both temporal directions. Truncation at another event's onset cannot remove all such leakage.
3. **Preservation of missing intervals.** Interpolation silently bridges interior gaps. A synthetic three-sample missing interval was replaced with finite values. This behavior is unsuitable for censoring that is intended to retain missingness.
4. **A complete trial manifest.** Trials with insufficient retained data disappear rather than remaining as explicit missing epochs with a status.
5. **A constant analysis sample across epoch time.** Censoring can leave different trials contributing at different lags. If RPE, response time, or another covariate relates to timing or exclusions, estimates at different lags can reflect different trial populations.

### Recommended censoring contract

Expose and document separate policies for nearest-event truncation, disjoint windows, and whole-epoch exclusion. If exclusive sample ownership is required, use an explicit allocation rule, such as specified boundaries or priorities, and test it directly. Such a rule still does not establish physiological isolation.

Retain every requested trial and report onset validity, collision bounds, interpolation support, retained sample counts, missing intervals, and exclusion reasons. Do not interpolate across designated censored gaps. Define behavior for simultaneous events and boundary equality explicitly.

For event isolation, model overlapping contributions jointly where the design permits identification. Regression-based EEG overlap correction uses event timings and covariates to estimate contributions together; it is a useful conceptual precedent rather than evidence that truncation alone separates events. See [Ehinger and Dimigen, 2019](https://pubmed.ncbi.nlm.nih.gov/31660265/).

## Improving neural activity and parametric-response estimation

### Validate gain and temporal resolution, not only correlation

A noiseless canonical-HRF simulation with signed impulses of graded amplitude produced a correlation above `0.999` between true amplitudes and linear-deconvolution peaks. Yet the peak regression slope was only approximately `0.34` with lambda `0.01` and the second-difference penalty.

A separate isolated-impulse example illustrates the difference between peaks and area:

| Method | Penalty | Recovered unit-impulse peak | Recovered area over plus/minus 10 bins |
| --- | --- | ---: | ---: |
| Linear C++ | Ridge | 0.3620 | 0.9899 |
| Linear C++ | First difference | 0.3581 | 0.9981 |
| Linear C++ | Second difference | 0.3542 | 0.9976 |
| Scalar Wiener | Ridge-style | 0.3620 | 0.9899 |

The peak occurred at the correct bin in these examples, but the estimates showed temporal spreading and negative sidelobes. These are favorable, noiseless tests using a matched HRF; the numerical gains are not general calibration constants.

A strong RPE correlation can therefore coexist with amplitude attenuation and temporal smearing. Integrated activity can be more stable for isolated responses, but its suitability depends on overlap, chosen support, and the scientific estimand. It should not be assumed to solve the problem for rapid event sequences.

Pass isolated impulses through the complete pipeline to quantify peak attenuation, temporal width, pre-event ringing, and neighboring-trial leakage. Interpolating onto a finer output grid improves alignment convenience but does not establish finer neural resolution.

### Constrain HRF estimation with independent information

Use known event timing, independent runs, or a localizer to constrain HRF shape before estimating unrestricted latent activity. Consider regional or hierarchical HRF models when voxelwise estimation is unstable.

There is an identifiability limitation in the continuous model: shifting the HRF can be offset by shifting the latent activity. Consequently, blind selection among HRF lags cannot uniquely establish neural latency. A synthetic test that selects the generating kernel does not establish general identification; boundary conditions and regularization also influence the choice.

### Incorporate nuisance effects and observation validity in the forward model

A useful extension is a penalized model of the form:

```text
minimize over x and gamma:
    || W (y - H x - C gamma) ||^2 + lambda || L x ||^2
```

Here `H` models HRF convolution, `C` represents nuisance regressors, `L` defines the activity penalty, and `W` represents the chosen observation weighting and noise treatment. The handling of missing observations and whitening must be specified consistently.

Preserve the regular acquisition clock when volumes are invalid. If preprocessing filters or residualizes the data, apply the corresponding transformation consistently to the forward model. Estimate temporal noise from suitable residuals and check whether pooled AR(1) treatment is adequate.

The wrapper currently appends zero-valued BOLD samples and treats them as observations. This imposes an artificial future-data constraint. In one late-event example, zero padding changed the estimated peak from `0.3410` to `0.3584`, while the middle-of-run estimate barely changed. Neither value establishes truth at an incompletely observed boundary. Use explicit boundary assumptions and report reduced support near run edges rather than treating invented samples as measured data.

### Add an event-informed option that still estimates trial time series

Fit a small temporal basis per trial, with shrinkage and shared event components, while modeling cue, response, feedback, and other relevant events together. Coarse splines or another low-dimensional basis can reduce the number of unconstrained latent samples. Trial deviations can be partially pooled around a shared response.

The resulting estimates remain dependent on the basis, penalty, HRF assumptions, and design identifiability. Closely spaced or simultaneous components cannot be separated merely by adding parameters.

GLMsingle offers relevant precedents for HRF selection, nuisance modeling, and regularization, but its scalar trial estimates do not replace the requested trial time series. Its cross-validated procedures also require some repeated conditions; suitability for an RPE task must be assessed rather than assumed. See [Prince et al., 2022](https://elifesciences.org/articles/77599).

Sparse activity or innovation penalties are additional candidates to compare with smooth quadratic penalties. They impose different assumptions and can introduce their own amplitude bias. Paradigm Free Mapping and Total Activation provide useful methodological context: [Uruñuela et al., Hemodynamic Deconvolution Demystified](https://arxiv.org/abs/2107.12026).

### Use realistic, independent validation

Simulate the actual event schedules, including off-grid onsets, durations, signed RPE amplitudes, neighboring event types, HRF variability, drift, motion, temporally correlated noise, and invalid volumes. Include relationships between timing, nuisance effects, and RPE.

Evaluate amplitude slope and bias, sign recovery, peak timing, temporal width, neighboring-trial leakage, interval coverage, exclusion patterns, and false-positive behavior under null inputs. Compare continuous deconvolution against event-informed alternatives and a conventional parametric BOLD model.

Use held-out runs or appropriately separated temporal blocks for validation. Do not select a deconvolution configuration by maximizing its observed RPE association and then treat that same association as independent validation. If RPE is explicitly included in a supervised event model, make that modeling choice and its inferential target clear.

## Bush/Cisler-specific assessment

The Bush/Cisler method can encode graded differences, but its sigmoid-bounded representation complicates linear amplitude interpretation. The original method defines continuous encodings between zero and one; it is not simply a binary event detector. Relative increases and decreases can be represented around a baseline, but the encoding is not automatically a calibrated signed neural-amplitude scale. See [Bush and Cisler, 2013](https://pubmed.ncbi.nlm.nih.gov/23602664/).

Additional implementation concerns:

- **Infinite initialization:** `src/deconvolve_nlreg.cpp`, approximately lines 109–114, maps exact minimum and maximum pre-encodings through a logit, creating infinite activations. The associated sigmoid derivatives vanish. Clip initialization away from zero and one and verify finite optimization states.
- **Termination:** The optimization has no iteration cap, and its stopping comparison uses absolute error change despite documentation describing a relative criterion. Add a maximum iteration count, finite-state checks, and explicit convergence diagnostics.
- **Ignored beta setting:** The `bush2015` branch in `R/voxelwise_deconvolution.R`, approximately line 215, does not forward `decon_settings$beta` to the resampling function, which therefore uses its own default.
- **Temporal dependence in resampling:** `R/deconvolve_nlreg_resample.R`, approximately lines 65–76, independently shuffles residuals. This destroys residual autocorrelation. Consider a validated block-resampling or fitted-noise procedure.
- **Interval interpretation:** Approximately lines 103–106 and 123–126 apply a t interval to the surrogate estimates. This describes uncertainty in the resampling mean rather than the distribution of plausible underlying activity estimates. In a diagnostic, the ratio of interval half-width to surrogate standard deviation fell from `0.8360` with 8 samples to `0.3605` with 32 samples. More computation should improve Monte Carlo precision, not by itself establish greater certainty about the participant's latent activity. Define the intended interval target and validate its coverage.

The 2015 resampling method was developed to estimate confidence and improve precision of event estimates; this does not establish amplitude calibration for the present trial-level goal. See [Bush et al., 2015](https://pubmed.ncbi.nlm.nih.gov/26226647/).

## Computational performance

### Demonstrated sparse-solver improvement

The convolution and finite-difference penalties have banded structure, but the current linear implementation constructs dense matrices. A diagnostic sparse prototype used sparse convolution and second-difference matrices, sparse Cholesky factorization, and a fixed lambda. It included matrix construction in the timed operation.

| Fixed-lambda problem | Signals | Current dense C++ | Sparse prototype | Speedup | Maximum absolute difference |
| --- | ---: | ---: | ---: | ---: | ---: |
| 300 timepoints | 200 | 0.062 s | 0.027 s | 2.30x | 1.87e-10 |
| 1,200 timepoints | 200 | 0.246 s | 0.084 s | 2.93x | 1.99e-10 |

These are local synthetic benchmarks, median of three runs, using a canonical HRF sampled at TR = 1 second, lambda `0.01`, a second-difference penalty, and no normalization or demeaning. The benchmark commands set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`. The prototype is not a production implementation and does not implement the full tuning and diagnostics interface.

### Wiener is not consistently faster in this implementation

A separate comparison used the current backends with a single canonical HRF and either fixed lambda `0.01` or the default 20-point tuning grid:

| Timepoints | Signals | Backend | Fixed lambda | 20-candidate tuning |
| --- | ---: | --- | ---: | ---: |
| 300 | 200 | C++ | 0.015 s | 0.145 s |
| 300 | 200 | Wiener | 0.044 s | 0.623 s |
| 1,200 | 200 | C++ | 0.240 s | 2.151 s |
| 1,200 | 200 | Wiener | 0.193 s | 2.977 s |

These are also medians of three local runs. Timings vary across runs and configurations. FFT-based asymptotic advantages do not establish an end-to-end speed advantage for the current R implementation, especially with tuning. The Wiener GCV issue must also be resolved before comparing automatically selected estimates scientifically.

### Further opportunities

1. **Reuse factorizations.** The C++ path solves the same left-hand-side matrix separately for activity and degrees of freedom. Factor it once per candidate and reuse it.
2. **Skip unnecessary diagnostics.** With fixed HRF and lambda and no requested diagnostics, omit tuning-only residual and trace calculations.
3. **Exploit banded structure throughout.** Avoid constructing dense convolution and penalty matrices. Preserve the intended boundary and whitening operators when optimizing.
4. **Separate tuning from application.** For global tuning, evaluate representative valid voxels, freeze the selected settings, and apply the operator in memory-bounded batches. Verify that the tuning sample represents relevant regions and noise levels.
5. **Avoid repeated deconvolution across overlapping atlases.** Extract and deconvolve a union of required voxels where settings permit, then map results into atlases.
6. **Aggregate before repeated alignment scans.** The alignment code repeatedly scans voxel data by trial. Preaggregation and indexed time-window extraction can reduce work substantially. Linear interpolation of a mean commutes with averaging under common support; median and standard deviation do not share that identity.
7. **Average before deconvolution when mathematically appropriate.** Parcel averaging commutes with a shared, fixed linear deconvolution operator. This does not hold generally for independently tuned voxel operators or nonlinear Bush estimation, and the intended voxel scaling must be respected.
8. **Reduce serialization cost.** Long voxel-by-time CSV outputs repeat extensive metadata. Consider matrix or image storage with separate metadata and validity information. Control worker counts and numerical-library threads together to avoid oversubscription.

## Recommended implementation order

1. Repair voxel selection and coordinate-based metadata mapping.
2. Isolate invalid signals and make normalization settings effective.
3. Repair missing-onset handling, interpolation-grid construction, and gap preservation.
4. Define censoring policies explicitly and return a complete trial/file status manifest.
5. Correct filename/identity integration and stop silently removing failed workers.
6. Correct and verify Wiener degrees of freedom; address Bush configuration and numerical issues if that path remains supported.
7. Implement sparse or banded linear fitting, factorization reuse, and a fast fixed-parameter application path.
8. Validate amplitude gain, timing, leakage, and uncertainty using realistic event schedules before selecting scientific defaults.
9. Develop and compare an event-informed trial-time-series model against the continuous method.

## Diagnostic artifacts

At the time of the audit, diagnostic scripts and output were stored in `/tmp/fmri-decon-review/`:

| Artifact | Contents |
| --- | --- |
| `review.R` / `results.txt` | Initial synthetic alignment, censoring, missingness, amplitude, and Wiener-trace checks |
| `followup.R` / `followup_results.txt` | Additional impulse, boundary, resampling, and backend timing checks |
| `integration.R` / `integration_results.txt` | Successful AFNI mask-label reproduction, filename discovery, and larger timing checks |
| `sparse_benchmark.R` / `sparse_results.txt` | Sparse/dense numerical comparison and final successful edition-3 focused tests |
| `images/` | Synthetic NIfTI inputs and outputs for the mask-label check |

These temporary artifacts are not committed and may be removed by temporary-directory cleanup. The follow-up script's first integration attempt lacked an exposed internal `mat2df` helper; `integration.R` supplied that helper and completed the actual extraction check. That harness issue is distinct from the reproduced package mask-label bug.

## References

- Bush and Cisler (2013). [Decoding neural events from fMRI BOLD signal: a comparison of existing approaches and development of a new algorithm](https://pubmed.ncbi.nlm.nih.gov/23602664/).
- Bush et al. (2015). [Improving the precision of fMRI BOLD signal deconvolution with implications for connectivity analysis](https://pubmed.ncbi.nlm.nih.gov/26226647/).
- Ehinger and Dimigen (2019). [Unfold: an integrated toolbox for overlap correction, non-linear modeling, and regression-based EEG analysis](https://pubmed.ncbi.nlm.nih.gov/31660265/).
- Prince et al. (2022). [Improving the accuracy of single-trial fMRI response estimates using GLMsingle](https://elifesciences.org/articles/77599).
- Uruñuela et al. [Hemodynamic Deconvolution Demystified: Sparsity-Driven Regularization at Work](https://arxiv.org/abs/2107.12026).
