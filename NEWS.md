# fmri.pipeline 0.4-1

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
