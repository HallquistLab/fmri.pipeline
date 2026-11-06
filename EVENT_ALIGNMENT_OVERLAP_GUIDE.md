# Reviewing event overlap and impulse recovery

This change makes overlapping event windows auditable and preserves missing data.
It does not separate neighboring neural responses or introduce a joint event model.
For the primary analysis, retain overlapping epochs, inspect their support, and
assess the fitted inverse's temporal spread using a specified generating model.

## What changed

1. Alignment interpolates each unit (usually a voxel) within its valid source
   segments **before** computing the requested group's mean, median, and SD.
   Adjacent source samples must both be finite and uncensored. An explicit
   missing row blocks interpolation; deleted acquisitions create a gap when the
   bracket exceeds `max_gap`, which defaults to the acquisition TR.
2. Every output value includes contributor counts, reasons for unavailable
   contributors, acquisition-time IDs, and the maximum distinct-epoch reuse of
   those IDs within the returned grouping.
3. Trial-grouped output includes neighboring event times, simultaneous-event
   counts, requested-window overlap counts, distances to run edges, and the policy.
4. `deconvolution_impulse_diagnostics()` measures recovery and spread into other
   windows under **fixed** linear deconvolution settings. It never tunes settings
   or z-scores each impulse.

Unmissing mean estimates remain equivalent to the previous interpolate-after-
aggregation approach on common grids, up to floating-point precision. With
missing contributors they can change intentionally: interpolation between two
means based on different voxels no longer supplies a nonexistent voxel response.
Off-grid median and SD also change intentionally: they now summarize interpolated
unit values, rather than interpolating the medians/SDs of source distributions.

Existing output files are not automatically invalidated. Regenerate alignment
outputs (`overwrite = TRUE` for `run_decon_alignment`) to adopt this behavior.

## Run a complete review example

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
  inst/examples/event_overlap_review.R /tmp/overlap-review
```

The script loads this checkout when `pkgload` is available, otherwise the installed
package. It writes:

- `alignment_allow.csv`, `alignment_truncate.csv`, and `alignment_exclude.csv`;
- `source_clock.csv`, mapping source IDs to absolute acquisition times;
- `impulse_recovery.csv`, with one status and recovery summary per injected event;
- `source_target_windows.csv`, with directed source-to-target window metrics;
- `impulse_recovery_uncentered.csv` and `source_target_windows_uncentered.csv`,
  providing a comparison with measured-run centering disabled;
- `impulse_timeseries.csv`, `frozen_settings.rds`, and `overlap_review.png`.

Inspect the middle event's alignment: the deliberately censored interval remains
missing even with half-TR output. In the recovery experiment, only the selected
source is injected in each simulation. A nonzero estimate in another trial's
window therefore describes spread of that injected event under the specified
model, without any actual activity from the neighboring trial.

## Alignment policies

| Policy | Behavior | Intended use |
| --- | --- | --- |
| `allow` | Retain the full requested window and available interpolation padding. | Main analysis with overlap annotations. |
| `truncate` | Retain source samples strictly after the closest supplied earlier collision and before the closest later collision. | Compatibility with previous collision settings; sensitivity analysis. |
| `exclude` | Return a missing epoch if its closed requested window intersects any other target epoch's requested window. | Conservative whole-epoch sensitivity analysis. |

`overlap_policy = NULL` preserves explicit legacy intent: collision columns select
`truncate`, otherwise the policy is `allow`. Set `allow` explicitly to retain
windows even when collision columns are supplied. `truncate` with no collision
columns has no truncation bounds.

`exclude` compares **requested target windows**, not interpolation padding and
not other event types. Touching endpoints count as an intersection; both affected
epochs are excluded. Simultaneous target onsets remain distinct trial rows under
`allow`/`truncate`, and are excluded under `exclude`. Strict truncation considers
only strictly earlier/later collisions, matching the previous implementation.

`neighbor_events` controls timing annotations independently of censorship. Supply
columns such as `c("cue", "response", "feedback")`. The exact target event record
is omitted from its own neighbor search; other event records, including another
event type on the same trial, are retained. Simultaneous counts count records.
Window-overlap counts always compare target epochs, regardless of this setting.

Neither truncation nor exclusion establishes physiological isolation. Even
nonintersecting requested windows can use the same interpolation brackets.
Midpoint allocation is not implemented as a presumed decontamination procedure.

## Supplying validity and inspecting support

```r
aligned <- get_medusa_interpolated_ts(
  object, event = "feedback", time_before = -3, time_after = 5,
  group_by = c("trial", "region"), output_resolution = 0.5,
  overlap_policy = "allow", neighbor_events = c("cue", "response", "feedback"),
  valid_col = "usable", max_gap = object$tr, time_audit = TRUE
)
```

The optional `usable` column must be logical in the source time-series table.
`FALSE` and `NA` make that row unavailable for all value columns while retaining
its time. Without a mask, nonfinite values remain unavailable independently for
each value column. Retain acquisition rows and mark them invalid whenever possible.
For deleted rows, `max_gap` is the longest allowed source bracket in **seconds**;
it is independent of the requested output resolution. Increasing it can bridge
deleted rows, but cannot bridge explicitly invalid rows. There is no extrapolation.

For a value named `decon`, these columns describe each output lag:

| Field | Meaning |
| --- | --- |
| `decon_mean`, `decon_median`, `decon_sd` | Statistics of finite interpolated unit estimates; SD requires two units. |
| `decon_n`, `decon_n_total` | Available and represented contributors. With trial grouping, contributors are units; when pooling trials, they are unit/epoch combinations. |
| `decon_n_censored` | Contributors blocked by an explicit validity mask at the requested sample or bracket. |
| `decon_n_missing` | Contributors blocked by nonfinite values, including placeholder units in missing/excluded epochs. Check `epoch_status` to distinguish the latter. |
| `decon_n_gap` | Contributors blocked by a deleted-row bracket wider than `max_gap`. |
| `decon_n_outside` | Contributors lacking a source sample or bracketing support, including run/collision boundaries. |
| `decon_source_ids` | Semicolon-separated sorted acquisition-time IDs used by finite contributors; empty when none contribute. |
| `decon_source_max_epochs` | Largest number of distinct epochs using any of those acquisition-time IDs within the same value/output group; zero with no contributors. Repeated lags within one epoch do not inflate it. |

The available count plus the four unavailable-reason counts equals the total.
`epoch_status` summarizes lag availability of any value's finite mean, not a
complete contributor cohort; inspect per-value counts for changing support.
When several reasons could apply to a bracket, explicit censoring takes precedence,
then nonfinite values, then excessive bracket duration. Exact valid samples remain
available even next to a gap. Duplicate times within a unit/epoch are errors.

IDs index sorted observed acquisition times within the input run, **not voxels**.
Reuse indicates a shared acquisition time, and does not assert identical voxel
subsets contributed to each epoch. The `source_clock` attribute maps IDs to times.
It is present on raw `event_lock_ts()` results and interpolated data. For batch
alignment, `alignment_file` disambiguates source files; the combined clock includes
file/identity columns and is written to `<out_file>.sources.csv`. The existing
`file_status` attribute and `<out_file>.files.csv` continue to record file failures.
IDs are not globally unique across files. No neural independence is implied by a
reuse count of one.

Per-trial fields are `align_onset`, `align_previous_event`, `align_next_event`,
`align_simultaneous_n`, `align_window_overlap_n`, `align_run_start_distance`,
`align_run_end_distance`, and `align_overlap_policy`. Run distances are signed;
an onset outside the run has a negative distance to one edge. Missing onsets have
missing timing annotations. `overlap_excluded` is the new epoch status; other
status meanings are retained. If trial is omitted from `group_by`, statistics can
pool epochs, but trial-level annotations and epoch statuses are not attached.
The `align_*` names and internal `.align_*` fields are reserved in source data.

The same four settings (`overlap_policy`, `valid_col`, `neighbor_events`, `max_gap`)
can be included in an alignment list passed to `evt_align_decon_files()` or
`run_decon_alignment()`. Validity masks apply to **post-deconvolution alignment**.
They do not refit deconvolution with masked BOLD observations, nor repair effects
of invalid BOLD samples that already entered an inverse fit. Observation-weighted
forward modeling remains separate work.

## Measuring temporal spread with frozen settings

```r
recovery <- deconvolution_impulse_diagnostics(
  event_data = trials, onset = "feedback", trial = "trial",
  n_time = n_volumes, tr = TR, time_offset = acquisition_offset,
  kernel = selected_hrf, lambda = selected_lambda,
  backend = "cpp", penalty = "diff2", prewhiten_rho = selected_rho,
  input_scale = fitted_bold_sd, input_demean = TRUE,
  tail_padding = length(selected_hrf),
  time_before = -3, time_after = 5,
  source_trials = trials$trial[unique(c(1, ceiling(nrow(trials) / 2), nrow(trials)))],
  return_timeseries = TRUE
)
```

Use the actual fitted settings from the run/voxel or a documented representative
operator. `lambda` is required. `kernel` is the selected HRF, not a candidate grid;
`prewhiten_rho` is fixed rather than re-estimated. The supported backends are
`cpp`, `r`, and **scalar** `wiener`. Adaptive PSD and Bush estimates are not fixed
linear operators and are not supported by this diagnostic.

For a fixed operator, the relevant response is `D H x`, including the specified
centering and boundary treatment. `input_scale` uses a frozen positive divisor;
outputs are multiplied back into injected-impulse units, so changing this divisor
alone does not manufacture a different recovery gain. `input_demean = TRUE`
applies the same linear measured-run mean-removal operation to each generated
column before padding. It does not divide each impulse by its own SD.

`tail_padding = 0` is the standalone finite-run case. Use the original kernel
length to reproduce the wrapper's appended zero BOLD observations. The diagnostic
retains that boundary assumption so its effects can be assessed; it does not fix
it. A different `forward_kernel` tests a specified HRF mismatch; both generating
and inverse HRFs are normalized to unit sum, consistent with the linear API.

Each selected source event is injected separately with zero pre-run history.
Off-grid onsets distribute unit mass between their two surrounding acquisition
bins. Thus the injected samples sum to one, and area is `tr`; this approximation
does not claim sub-TR neural resolution. `truth_at_onset` records the interpolated
grid truth, which can be less than one off-grid. `onset_gain` divides recovery by
that grid truth, while `mass_gain` compares summed recovered samples with unit
injected mass. Peak timing is reported relative to the specified onset and is
limited by the acquisition grid. With `input_demean = TRUE`, mean removal adds a
run-wide baseline component to the impulse response. Signed recovered mass is no
longer expected to equal one; pre-onset and distant-window activity can include
this centering effect as well as inverse-filter spreading. The example compares
centered and uncentered inputs so these effects can be inspected separately.
Missing/outside onsets and impulses producing
zero or constant observed BOLD retain statuses and missing recovery metrics.
Events outside the measured acquisition span are not injected, even though a
pre-run event can physiologically affect BOLD acquired later.

The directed `windows` table compares each selected source with every target
window. It reports finite lag counts, absolute-area and signed-area sums, maximum
absolute response, and `window_l1_ratio` relative to that source's own window.
Area summaries are rectangular grid sums times `output_resolution` over available
lags. Compare coverage before comparing areas near run edges. Ratios can exceed
one, and overlapping windows count some response repeatedly: **they are not
fractions of measured neural attribution** and need not sum to one.

No noise, RPE values, automatic tuning, observation-mask refitting, or uncertainty
intervals are included. These are conditional operator diagnostics. Validating
scientific recovery still requires realistic HRFs, noise, actual schedules and
RPE sequences, and a separate evaluation of the full parameter-selection process.
For large schedules, begin with selected early/middle/late source events: window
output grows as selected sources times all target events. Full time series are
optional and scale as selected sources times measured volumes. Impulse fits are
batched; alignment processes one output group at a time and uses shared matrix
brackets for common rectangular grids, with a per-unit fallback for ragged grids.

## Review checklist

- Verify behavioral onsets, `time_offset`, discarded volumes, and acquisition
  timing use the same clock.
- Inspect the example's missing interval and compare counts with the original mask.
- Confirm source IDs resolve through `source_clock` and repeated lags do not
  increase distinct-epoch reuse counts.
- Compare policies and report changes in included trials and contributors by lag.
- Inspect whether exclusions/support relate to RPE, reaction time, or condition.
- Compare impulse recovery across run positions, lambda/HRF choices, and padding
  assumptions. Do not select a pipeline by maximizing its observed RPE association.
- Account for repeated acquisitions and temporal dependence in downstream
  inference; avoid treating overlapping epoch rows as independent observations.

## Validation performed

The package installed successfully in an isolated library. Focused installed-
package checks passed **74 test cases and 708 assertions**, with no failures,
errors, or warnings. They cover explicit/deleted gaps, masks, contributor changes,
simultaneous events, policy boundaries, source clocks, empty clocks, legacy/custom
identity mappings, multivariate output, native extraction/geometry, and the
existing linear/voxelwise regressions. Matrix interpolation agrees with the
per-unit oracle on masked grids. Impulse recovery agrees with a dense finite
inverse, and directed window metrics agree with independently sampled responses.
Tests also cover fixed centering/scaling/whitening/padding, batches, off-grid
onsets, missing/outside/unobserved impulses, and specified HRF mismatch.

A two-worker installed-package check matched serial estimates, contributor
counts, overlap annotations, and source clocks for both allow and truncate
policies. R/Rd parsing, Rd checks, and whitespace checks passed. The review script
ran and generated the documented tables and plot. The full unrelated package
suite and participant-data scientific validation were not run.
