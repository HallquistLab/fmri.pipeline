# Align a set of deconvolved time series files from `voxelwise_deconvolution` to a set of events in a task-related design

Align a set of deconvolved time series files from
`voxelwise_deconvolution` to a set of events in a task-related design

## Usage

``` r
run_decon_alignment(
  atlas_files,
  decon_dir,
  trial_df,
  alignments = list(),
  nbins = 12,
  aggregate_by = "atlas_value",
  atlas_subset = NULL,
  overwrite = FALSE,
  tr = NULL,
  ncpus = 8,
  mem_per_cpu = "4g",
  walltime = "1:00:00",
  scheduler = "slurm",
  time_audit = FALSE,
  file_manifest = NULL,
  fail_on_error = TRUE
)
```

## Arguments

- atlas_files:

  A character vector containing filenames of atlases that used for
  aggregating voxel-level deconvolved time series

- decon_dir:

  The directory containing the deconvolved data from
  `voxelwise_deconvolution`. Should be the same as `out_dir` from
  `voxelwise_deconvolution`.

- trial_df:

  Trial-level data used for event alignment.

- alignments:

  A list of alignment specifications, each of which will be computed by
  this function. See `get_medusa_interpolated_ts` for help with
  arguments, and see Details. Settings include overlap_policy,
  valid_col, neighbor_events, and max_gap.

- nbins:

  For atlases with continuous values, how many bins should be used to
  discretize the mask values, leading to aggregation by bin.

- aggregate_by:

  The column name in the individual deconvolved files used for averaging
  repeated units (e.g., voxels) into single event-aligned time series.
  Most commonly, this is "atlas_value", which will lead to averaging of
  voxels within each parcel in the mask.

- atlas_subset:

  Optional subset of atlas values to retain.

- overwrite:

  if TRUE, overwrite existing output files

- tr:

  The repetition time of the scan in seconds

- ncpus:

  The number of cores to use for each parallel job

- mem_per_cpu:

  The amount of memory (RAM) allocated for each CPU job. Default is "4g"
  (4 gigabytes)

- walltime:

  The time requested for each event alignment job. Default: 1:00:00 (1
  hour).

- scheduler:

  The job scheduler to be used in alignment. Options are "slurm"
  (default), "torque", or "local"

- time_audit:

  If TRUE, additional columns will be added to the output showing how
  alignment is calculated vis-a-vis event timing

- file_manifest:

  Optional data.frame with `file`, `id`, `run`, and optional `session`
  columns. Paths must identify the deconvolution files. Overrides
  per-file metadata and legacy filename parsing.

- fail_on_error:

  Stop if any file fails alignment (default TRUE). Set FALSE explicitly
  to allow partial output with a warning and a file status manifest.

## Details

This function will create a separate R batch job on the scheduler for
each combination of atlas files and alignments. Each of these jobs can
be allocated a number of cores, which will be used to loop over the
deconvolved time series files for event alignment. The number of cores
requested for each alignment job is specified using `ncpus`, which
defaults to 8.

The alignments argument is a list containing information about how to
event-align the deconvolved time series. Elements include:

- `$evt_col`: (Required) The name of the column in `trial_df` containing
  the event of interest.

- `$time_before`: (Required) The number of seconds before the event of
  interest to include in the interpolation window.

- `$time_after`: (Required) The number of seconds after the event of
  interest to include in the interpolation window.

- `$pad_before`: (Optional) The number of seconds before the earliest
  event to include in the interpolation window to avoid missing values
  at the edge of the interpolation window. Default: -1.5

- `$pad_after`: (Optional) The number of seconds after the latest event
  to include in the interpolation window to avoid missing values at the
  edge of the interpolation window. Default: 1.5

## Examples

``` r

if (FALSE) { # \dontrun{
  atlas_dir <- file.path(
    "/proj/mnhallqlab/projects/clock_analysis/fmri",
    "ph_da_striatum/masks"
  )
  atlas_files <- file.path(
    atlas_dir,
    c(
      "bilateral_striatum_tight_7Networks_2.3mm.nii.gz",
      "pauli_combined_integermask_2.3mm.nii.gz"
    )
  )

  # Outputs of voxelwise_deconvolution() for these atlases.
  decon_dir <- "/proj/mnhallqlab/users/michael/sceptic_decon"
  trial_df <- get_trial_data(
    repo_directory = "/proj/mnhallqlab/projects/clock_analysis",
    dataset = "mmclock_fmri",
    groupfixed = TRUE
  )

  alignments <- list(
    clock_long = list(
      evt_col = "clock_onset",
      time_before = -5,
      time_after = 10,
      collide_before = "iti_onset", # censor data if we bump into the end of the prior trial
      collide_after = "clock_onset" # censor data it we hit the next trial
    ),
    rt_long = list(
      evt_col = "rt_time",
      time_before = -4,
      time_after = 7,
      collide_before = "iti_onset",
      collide_after = "clock_onset"
    )
  )

  # run all atlases and alignments as separate slurm jobs (2 hours each, 8 cpus) 
  run_decon_alignment(atlas_files, decon_dir, trial_df, alignments, 
    overwrite = TRUE, tr = 1.0, ncpus = 8, walltime = "2:00:00", scheduler = "slurm")

} # }
```
