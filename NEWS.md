# fmri.pipeline 0.4-1

## Banded linear deconvolution solver and exact Takahashi trace

* `deconvolve_reglin_cpp()` now exploits the strict banded structure of
  $H_k^\top H_k + \lambda P + \delta I$ (half-bandwidth $b = \max(\text{kernel\_len} - 1 + \mathbb{I}(\text{prewhiten}), \; b_P)$),
  replacing dense $O(n_{\text{latent}}^3)$ normal-equation solves and dense
  matrix-RHS trace evaluations with an in-place $O(n_{\text{latent}} b^2)$
  lower-banded Cholesky factorization, contiguous multi-signal forward/backward
  substitution, and the exact $O(n_{\text{latent}} b^2)$ Takahashi banded-inverse
  recursion for effective degrees of freedom $\operatorname{tr}(A(\lambda)^{-1} H_k^\top H_k)$.
* `hrf_inverse_operator_cpp()`, `hrf_score_run_cpp()`, and `hrf_noise_profile_cpp()`
  now eliminate dense $m \times m$ AR(1) whitening matrix construction and
  dense $O(m^2 L)$ / $O(L m m_0)$ matrix multiplications in favor of
  in-place column-wise bidiagonal whitening, tridiagonal $W^\top W$ stencil
  accumulation, L1-tiled banded Cholesky solves (for `nuisance_mode = "ignore"`),
  and single-pass batch normalization across the REML noise grid.

## Coefficient recovery validation

* Test fixed window-sum and calibrated-template amplitude readouts, and
  regression through the frozen inverse's transformed event design. Include
  active/confounded nulls, off-library HRFs, matched references and a direct
  forward GLS benchmark. See `COEFFICIENT_RECOVERY_VALIDATION.md`. These are
  validation candidates; package algorithms and defaults are unchanged.

## Parametric recovery validation

* Add a separate simulation study with genuine neural trial variability,
  task-active RPE nulls, correlated task confounds, and paired true-HRF,
  true-covariance, noiseless and matched-grid references. Audit coefficient bias,
  coverage, false positives and power across runs and independent datasets at
  eight fixed lambdas. Native inverse operators and trial regression
  factorizations are reused. See `PARAMETRIC_RECOVERY_VALIDATION.md`.

## Noise-aware HRF evidence and event recovery diagnostics

* Cross-run HRF fitting now defaults to pooled grid-REML AR(1) covariance,
  conditional GLS task/shape evidence with library corrections, training-run
  compatibility, and a fixed penalty toward the prespecified shape. Unsupported
  adaptation has an explicit parcel/prespecified fallback. Tests are conditional
  on the noise/mean model and are not whole-brain significance or posterior
  confidence. `noise_model = "ols", evidence = "legacy"` retains the prior fit.
* The native inverse shares operators across HRF/noise strata and uses original
  acquisition times across censor gaps. AR/calibrated fits use observation masks
  during scaling and fitting by default. Nuisance adjustment remains explicit;
  optional joint projection applies to both convolution and observations.
  Unidentifiable covariance is reported rather than silently assumed white.
* `diagnose_crossrun_recovery()` reports conditional event amplitude gain and
  excess leakage relative to the ideal readout, accounting for natural neural
  window overlap, fractional onsets, masks, and the frozen inverse. It does not
  automatically rescale amplitudes. See `NOISE_AWARE_HRF_RECOVERY.md` for review,
  independent reference checks, repeated validation, and throughput benchmarks.

## Cross-run HRF estimation and spatial regularization

* `hrf_library()` constructs a small unit-sum HRF library.
  `estimate_crossrun_hrf()` selects HRFs using other runs only, with batched
  Rcpp/BLAS training fits, nuisance regressors, observation masks, explicit
  support gates, and a recorded parcel/prespecified fallback hierarchy.
  Training can use condition means or free event amplitudes for signed trial
  responses that would otherwise average away, without supplying RPE values.
* Explicit parcel pooling and optional voxel shape penalties toward
  training-derived parcel HRFs provide spatial regularization. These constrain
  HRF selection without spatially smoothing recovered neural activity. Strength
  and effect-size thresholds are user-controlled. The original heuristic mode
  does not supply significance tests.
* `deconvolve_crossrun_hrf()` applies frozen HRFs with required fixed lambda,
  grouped inverse operators, bounded caching, block multiplication, and recorded
  per-voxel scaling. Existing NIfTI/file deconvolution defaults are unchanged.
  See `CROSSRUN_HRF_GUIDE.md` and the reproducible review/benchmark scripts.
* `HRF_ESTIMATION_VALIDATION.md` reports 480 independent simulated datasets
  comparing learned, canonical, and true HRFs at three fixed inverse penalties.
  Learning helps with stable shape mismatch but can hurt with low SNR,
  zero-mean signed responses, run changes, or neural-duration mismatch. The
  original support gate accepted all null free-event fits under colored noise;
  support labels are not confidence measures. Estimator settings were not tuned
  to these results. Reproducible scripts, compact tables, and figures are included.

## Auditable event overlap and impulse recovery

* Alignment now interpolates within each unit's valid acquisition segments before
  aggregation. Explicit invalid values/masks block interpolation, and brackets
  wider than the acquisition TR are missing by default. Off-grid median/SD now
  summarize interpolated units; means with changing contributors also change.
* `overlap_policy` supports `allow`, legacy `truncate`, and whole-window `exclude`.
  NULL retains collision-column truncation when explicitly configured, otherwise
  allows overlap. `valid_col`, `neighbor_events`, and `max_gap` pass through file
  and batch alignment settings. No policy claims physiological isolation.
* Outputs report per-lag contributor/reason counts, acquisition-time IDs and epoch
  reuse, neighboring events, simultaneous events, window overlap, and run edges.
  Source clocks are returned as attributes; file alignment adds `alignment_file`
  and writes a `.sources.csv` clock sidecar beside requested outputs.
* `deconvolution_impulse_diagnostics()` evaluates isolated impulses on actual
  event schedules under frozen linear settings, including specified HRF mismatch,
  centering, whitening and padding. Recovery and source-to-target window metrics
  are conditional diagnostics, not measured neural-attribution fractions.
* `EVENT_ALIGNMENT_OVERLAP_GUIDE.md` documents changed semantics, limitations,
  review checks, and the runnable `inst/examples/event_overlap_review.R` example.
  Regenerate existing alignment outputs to adopt the new behavior.

## Deconvolution scaling, Wiener tuning, and performance

* Linear deconvolution supports explicit `scaling = "zscore"`, `"demean"`,
  `"none"`, or `"psc"`; z-scoring remains the default. Legacy normalization
  flags still work. Voxelwise settings now control the actual preprocessing,
  applied once before padding. Bush voxelwise fits require z-scored input.
* Per-voxel `.scaling.csv` sidecars record the applied centers and divisors;
  file metadata records scaling, algorithm, and units. Original-series files
  expose `BOLD_input`, with `BOLD_z` retained only for z-scored observations.
  Regenerate existing outputs to adopt these changes.
* Bush resampling uses the same normalized observations for its base fit and
  residual construction. Its new `normalize` argument defaults to `TRUE`;
  voxelwise callers disable redundant normalization. The voxelwise Bush 2015
  path now forwards `beta` and retains matrix dimensions for single-voxel fits.
* Scalar Wiener GCV uses the exact finite-observation smoother trace rather
  than summing padded Fourier gains. Adaptive PSD Wiener requires fixed lambda
  and HRF; its conditional filter trace is labeled and GCV is missing, avoiding
  invalid automatic selection from a data-dependent filter.
* The C++ linear solver builds difference-penalty crossproducts directly
  and skips scoring for fixed-parameter matrix-only fits. Wiener fitting reuses input/kernel FFTs
  across candidates and skips unrequested fixed-parameter scoring. The R
  reference implementation remains available for numerical comparisons.

## Event alignment and file identity

* Missing/nonfinite event onsets and unsupported windows retain explicit missing
  epochs with `epoch_status`; nonfinite collision times no longer abort runs.
  Original event variable mappings are preserved during alignment.
* Zero-width, one-sided, and sub-sample windows use a unique grid within the
  requested bounds. An exact singleton observation can be retained without
  interpolation; duplicate reshape keys cannot silently turn signals into counts.
* Default deconvolution names are visible, and each new output has an identity
  metadata sidecar. Alignment supports explicit file manifests, alphanumeric IDs,
  sessions, and leading-zero labels, with legacy filenames as a fallback.
  Discovery includes old hidden outputs and excludes sidecars.
* Alignment records a status for every input file and stops on failures by
  default. Explicit `fail_on_error = FALSE` permits partial output with a warning.
  Status is returned in `file_status` and saved beside requested CSV outputs.

## Voxelwise deconvolution correctness

* Replaced `3dmaskdump` with internal RNifti/Rcpp extraction using the same
  masked voxel indices as the output metadata. Atlas, additional mask, and
  input spatial geometry are validated before extraction. `afni_dir` remains
  accepted for compatibility but is ignored.
* Constant and nonfinite signals are excluded before scaling and fitting, so
  they cannot invalidate other voxels or influence global tuning. Output rows
  and columns are preserved with missing estimates and explicit status values
  (`decon_status` in voxelwise output, `signal_status` in linear diagnostics).

## Confound cache invalidation and batch exit status

* Confound caches now record settings, demeaning, run inputs, and source-file
  fingerprints. Changes rebuild confounds, column metadata, exclusion decisions,
  and truncation data together, including cached raw confound and motion inputs.
  Original image provenance prevents repeated volume dropping when finalized run
  data are reused. Legacy caches rebuild from original inputs; a previously
  truncated run without provenance requires its original image and volume count.
  Existing fitted models must still be explicitly re-estimated after changes.

* Generated batch scripts now retain nonzero task exit statuses in both serial
  and background execution, wait for every background task, and preserve failures
  through cleanup. Each task runs in its own subshell with `set -e`; setup and
  cleanup failures also produce a nonzero status. Standard shell rules for
  conditionals and pipelines still apply.
## L2 model refitting and concatenated confounds

* L2 refits now preserve configured covariate transformations and factor reference
  levels. Transformation parameters are recomputed from complete retained runs
  within each subject/session or pooled subject, according to the L2 scope.

* Concatenated L1 confounds now keep each session/run's spike regressors separate,
  with zeros outside that run. Single-spike regressors retain their matrix shape
  and column names so they are correctly identified during concatenation.

* Existing analyses need regenerated per-run confounds (including column metadata),
  concatenated confounds, and L2 designs, followed by re-estimation of affected models.

## Declarative FWE correction workflow

* Added a serializable `fwe_spec()` API for describing level-3 FWE corrections
  with semantic model/contrast selectors, method options, correction masks, and
  scheduler settings. Specifications can be read from and written to YAML.

* Added pTFCE plan and execution support through `plan_fwe_correction()`,
  `run_fwe_plan()`, and `refresh_fwe_plan()`. A single specification can fan
  out across all selected level-3 z-statistic maps, with explicit preflight
  status for required inputs and optional dry-run command inspection.

* Added persistent FWE result snapshots and artifact manifests via
  `collect_fwe_results()` and `fwe_result_artifacts()`, providing stable,
  semantically selectable output paths for downstream analyses.

## ROI extraction and scheduler error handling

* ROI extraction now validates every statistic image against the mask's spatial
  dimensions, voxel sizes, units, and transform before indexing. Misaligned
  inputs stop extraction with a resampling/registration diagnostic. Reported
  coordinates use the same sform-first transform as the geometry checks.
* Completed L3 ROI extraction using the shared group-output lookup, preserving
  lower-level contrast identities and pooled/per-session metadata without
  duplicating maps for each subject. Group-only extraction no longer requires
  lower-level setup tables, and L1-to-L3 analyses are supported.
* Slurm cancellations, including `CANCELLED by <uid>` responses, are terminal
  failures in `wait_for_job()`. Failed submissions with `fail_on_error = FALSE`
  warn and return `NULL` without attempting to track a nonexistent job.
* Removed the redundant package reload from test setup and added package-local
  command guards so unit tests fail on unmocked scheduler submissions.

## Bug fixes in `build_design_matrix` / `fmri.stimulus` convolution

* **Fix: overlapping events silently overwritten during stimulus construction.** `construct_stimulus`
  (inside `fmri.stimulus`) used assignment (`stim[idx] <- values[i]`) when building the microtime
  stimulus vector. When two events overlapped in time (e.g., long-duration events with short ITIs),
  later events overwrote earlier ones, silently producing incorrect regressors. Fixed to use additive
  accumulation (`stim[idx] <- stim[idx] + values[i]`), which is correct by linearity of convolution.
  This affected all normalization modes (`none`, `durmax_1`, `evtmax_1`).

* **Fix: microtime bin truncation causing position-dependent convolution peaks.** Floating-point
  arithmetic in `fmri.stimulus` when converting onset times to microtime bins could produce values
  like `2013.999...` instead of `2014`, which R then truncated to the wrong bin. This shifted events
  by one microtime step, causing identical events at different run positions to have ~0.4% different
  peak amplitudes. Fixed by applying `round()` to microtime onset and duration indices.

* **Optimized `evtmax_1` normalization.** Replaced the previous approach of N separate full-run FFT
  convolutions (one per event) with precomputed HRF peaks via small padded windows, followed by a
  single combined convolution. This also eliminates the previous fragile workaround for end-of-run
  events (re-convolving at the center of the run to estimate the "true" peak). Speedups of 2-9x
  depending on number of events, number of unique durations, and run length.

* Thread `microtime_resolution` parameter through `convolve_regressor` to ensure peak lookups and
  grid-phase computation always use the same resolution as the main convolution.

* **Fix: insufficient convolution padding with sub-second TRs.** `convolve_stimulus` used a fixed
  20-TR padding for edge-effect suppression, which at TR=0.5s was only 10 seconds — not enough for
  the double-gamma HRF undershoot to decay (still ~-0.02 at 20s). Changed to 30 seconds of absolute
  time padding regardless of TR, preventing circular convolution artifacts near run boundaries in
  multiband acquisitions.

## Repeated timing rows and within-trial occurrences

* Fixed L1 setup for occurrence-keyed signals when subject events are aggregated as a `data.table`.
  The occurrence validation added in July used data-frame column-selection syntax that `data.table`
  interpreted as a request for a literal `join_cols` column, causing every affected design build to
  fail and leaving an empty FSL setup table.

* `run_feat_sepjobs()` now validates its FSL setup table before model filtering and reports a clear
  upstream-setup diagnostic instead of the misleading `object 'l1_model' not found` error.

* Pipeline inputs and time-series multipliers now copy incoming `data.table`s to ordinary
  `data.frame`s only when necessary. Fresh `rbindlist()` aggregations are likewise normalized at
  data-frame API boundaries without copying caller-owned objects or converting existing data.frames.

* `mixed_by()` and `fill_atlas_with_stats()` now isolate their internal `data.table` operations from
  caller-owned inputs. Keying, sorting, and conversion inside these functions no longer change an
  input object's class, key, or row order by reference.

* `trial_data` now supports ragged timing rows: a row may contain trial-level timings alongside a
  within-trial occurrence, and sparse onset columns use `NA` to indicate no occurrence on that row.
  Trial identifiers may be `NA` for free-running occurrences.

* Exact duplicate timing/value records are collapsed before convolution to prevent copied trial rows
  from inflating regressors. Set `keep_duplicate_occurrences: true` on a YAML signal to retain them.

* Parametric modulators and within-subject factors now align to an internal timing-occurrence key,
  rather than `run_number + trial`, so repeated trials do not create many-to-many joins. Modulator
  values must be present on each modeled timing row; missing values omit that occurrence.

* Beta-series signals now require exactly one distinct timing occurrence per trial after duplicate
  collapse and reject missing trial identifiers or ambiguous multi-occurrence trials.

* Removed unused `subtrial` classes, mappings, arguments, and documentation.

## Backend system overhaul

* Per-level and per-model backend selection via `gpa$level_backends`, `execution_backend`, and `producer_backend` fields
* Artifact-based dependency resolution: backends now declare which artifacts they produce and require, and the pipeline validates these at submission time
* Preflight report logged at pipeline start showing resolved backends, producer mappings, and artifact satisfaction for each analysis level
* `harvest_l3_inputs()` replaces the former `harvest_fsl_copes()` with a dispatcher that routes to backend-specific harvesters based on the declared `producer_backend`
* Run-time backend overrides via `level_backends` and `backend_overrides` arguments to `run_glm_pipeline()`

## FSL/FLAME backend reliability

* **Fix: FSL L3 jobs now wait on the backend-specific L2 producer job.** FSL L3 setup/run jobs now
  depend on `setup_run_l2_fsl` when they consume FSL L2 outputs, avoiding races where L3 setup could
  start before the FSL L2 cache and FEAT outputs were ready.

* **Improved slice-parallel FLAME1+2 recovery.** The local `flame_runner` keeps FSL's default
  z-slice FLAME parallelization, but no longer fails the full model immediately when one FLAME12
  slice crashes. Failed `--runmode=flame12` slice commands are retried as `--runmode=flame1`, with
  FLAME12-only MCMC options removed from the retry command.

* **Surfaced FLAME12 slice fallbacks in logs.** When a slice is recovered by FLAME1, the failed
  FLAME12 slice log directory is moved aside as `statsNNNN.flame12_failed.*`, and the FEAT directory
  receives a `flame_runner_fallbacks.tsv` audit file. These records are also promoted into the
  level-specific `lgr` model-estimation logs (`logs/l1_estimation.*`, `logs/l2_estimation.*`, or
  `logs/l3_estimation.*`) as warnings so users can see that estimation recovered by falling back to
  FLAME1 for one or more slices.

## Level 1 setup refactor

* Refactor `setup_l1_models` into smaller helper functions (`setup_l1_subject`, `load_cached_l1_bdm`, `drop_runs_without_events`)
* Improved `build_design_matrix` error messages for missing onset/duration columns and empty event tables

## Other improvements

* Support for SPM-style concatenated run handling at Level 1
* AFNI `3dLMEr` Level 3 execution backend with automatic `-qVars` detection and `dataTable.txt` generation
* Backend specification defaults (`default_glm_backend_specs`) for FSL, SPM, and AFNI
* `lookup_feat_outputs()` now returns stable, level-specific `$l1`, `$l2`, and
  `$l3` tables by default, with an explicit `format = "long"` option. L2
  passthroughs are represented as logical L2 rows and identify their L1
  materialization and source image.

# fmri.pipeline 0.2-1

* Fix to use `run_data$tr` inside `build_design_matrix` call
* Always include evt_time 0 in event-aligned outputs for MEDuSA
* Provide `time_audit` argument in event alignment to output additional columns pertaining to alignment timing
* include `volume` as column in outputs of `voxelwise_deconvolution` to support time audits
