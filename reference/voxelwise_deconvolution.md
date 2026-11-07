# Function to perform voxelwise deconvolution on an fMRI dataset using the fMRI model arguments object

Function to perform voxelwise deconvolution on an fMRI dataset using the
fMRI model arguments object

## Usage

``` r
voxelwise_deconvolution(
  niftis,
  add_metadata = NULL,
  out_dir = getwd(),
  out_file_expression = NULL,
  log_file = file.path(out_dir, "deconvolve_errors"),
  TR = NULL,
  time_offset = 0,
  atlas_files = NULL,
  mask = NULL,
  nprocs = 20,
  save_original_ts = TRUE,
  algorithm = "bush2011",
  decon_settings = list(nev_lr = 0.01, epsilon = 0.005, beta = 60, kernel =
    spm_hrf(TR)$hrf, n_resample = 25),
  afni_dir = NULL
)
```

## Arguments

- niftis:

  A vector of processed fMRI timeseries images (4D files) to be
  deconvolved

- add_metadata:

  A data.frame with one row per value of `niftis`, saved as atlas-level
  metadata. The `id`, `run`, and optional `session` columns are also
  saved in per-file sidecars for event alignment.

- out_dir:

  Base output directory for deconvolved time series files. Default is
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- out_file_expression:

  Expression evaluated to resolve the filename for the deconvolved csv
  files. Default is to prefix `"decon_"`, replace slashes in `niftis`
  with periods, and add the atlas image name. Note that the suffix
  `_deconvolved.csv.gz` or `_original.csv.gz` will be added to the
  expression, so don't pass this piece.

- log_file:

  Name (and path) of log file for any deconvolution errors or messages

- TR:

  the repetition time of the sequence in seconds. Required

- time_offset:

  The number of seconds that will be subtracted or added to the time
  field. Default: 0. Useful if some number of volumes have been dropped
  from the NIfTI data prior to deconvolution.

- atlas_files:

  optional atlas file vector specifying voxels used in deconvolution. If
  omitted, perform whole-brain deconvolution

- mask:

  an optional character string specifying a mask that should be used to
  constrain bounds of deconvolution.

- nprocs:

  The number of processors to use simultaneously for deconvolution

- save_original_ts:

  Whether to save the voxelwise BOLD data prior to deconvolution (for
  comparison/diagnosis). Default: TRUE

- algorithm:

  Which deconvolution algorithm to use for deconvolving voxelwise time
  series. Default: "bush2011". Alternative is "bush2015", which
  implements a resampling approach as well, or
  "reglin"/"regularized_linear", which estimates a continuous latent
  activity time series using regularized linear deconvolution with
  optional HRF tuning.

- decon_settings:

  A list of settings passed to the deconvolution algorithm.

- afni_dir:

  Retained for compatibility; ignored. NIfTI extraction is internal and
  does not require AFNI.

## Value

Nothing (invisible NULL).

## Details

The Bush 2011 algorithm uses the package's internal RcppArmadillo
implementation and does not require an external executable. Atlas, mask,
and input images must share a spatial grid. Constant signals and signals
containing nonfinite observations are excluded before fitting; their
output rows are retained with missing estimates. The `decon_status`
column identifies `"ok"`, `"constant"`, or `"nonfinite"` signals. Each
deconvolution file has a `.metadata.csv` sidecar containing its source
NIfTI path and any `id`, `run`, and `session` columns in `add_metadata`.
These identify files during event alignment, including custom filenames.
If absent, alignment falls back to the source path or requires an
explicit `file_manifest`.

If you want to use subject metadata to name the output file, use
`this_subj` in your `out_file_expression`, which will give you access to
a one-row data.frame containing the metadata for the current subject in
the loop.

Set `decon_settings$scaling` explicitly to `"zscore"`, `"demean"`,
`"none"`, or `"psc"`. The default is z-scoring. An explicit mode
overrides `normalize`/`demean`; otherwise `normalize = FALSE` selects
demeaning unless `demean = FALSE` too. The voxelwise Bush paths require
z-scoring because their learning-rate and bounded-encoding/resampling
behavior is sensitive to scale. Their output remains a bounded encoding.
Use the linear algorithm for estimates in alternative units. For linear
methods, z-scoring is recommended for standardized trial responses;
demeaning retains input units, and percent signal change requires an
appropriate positive baseline rather than residualized/zero-mean data.
Raw scale does not by itself degrade the conditioning of the fixed
linear inverse, but it changes signal weights during global tuning. See
[`deconvolve_reglin`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_reglin.md).

Scaling uses the observed run once, before any zero tail is appended
(`decon_settings$tail_padding`, defaulting to `0L` for
`reglin`/`regularized_linear` and `length(decon_settings$kernel)` for
`bush2011`/`bush2015`); solver normalization is disabled, as is
linear-solver demeaning of padded data. Each output has a `.scaling.csv`
sidecar containing voxel number, mode, subtracted center, divisor, and
validity status. Its `.metadata.csv` records algorithm, scaling, output
units, and the scaling filename. Original-series files contain
`BOLD_input`, the transformed observations supplied before padding, plus
the legacy `BOLD_z` alias only for z-scored inputs. Existing output
files are skipped; regenerate them to adopt the new scaling behavior and
sidecars.

## Examples

``` r
  if (FALSE) { # \dontrun{

    #name outputs according to subject metadata
    xx <- voxelwise_deconvolution(
      niftis="/proj/mnhallqlab/user/michael/test_nifti.nii.gz",
      out_dir="/proj/mnhallqlab/user/michael/decon_outputs",
      out_file_expression = expression(
        paste0(
          this_subj$subid,
          "_run",
          this_subj$run_num,
          "_",
          atlas_img_name
        )
      )
    )
  } # }
```
