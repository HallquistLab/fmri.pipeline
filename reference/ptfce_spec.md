# R6 class for pTFCE specification for one or more z-statistic images

R6 class for pTFCE specification for one or more z-statistic images

## Active bindings

- `fwe_p`:

  a vector of p-values used for familywise error (FWE) z-statistic
  threshold calculations in pTFCE.

- `two_sided`:

  If `TRUE`, use native two-sided pTFCE on a shared threshold grid. If
  `FALSE`, only run pTFCE on the positive tail (z \> 0).

## Methods

### Public methods

- [`ptfce_spec$new()`](#method-ptfce_spec-initialize)

- [`ptfce_spec$get_ptfce_calls()`](#method-ptfce_spec-get_ptfce_calls)

- [`ptfce_spec$get_expected_files()`](#method-ptfce_spec-get_expected_files)

- [`ptfce_spec$run()`](#method-ptfce_spec-run)

- [`ptfce_spec$submit()`](#method-ptfce_spec-submit)

- [`ptfce_spec$is_complete()`](#method-ptfce_spec-is_complete)

- [`ptfce_spec$get_clusters()`](#method-ptfce_spec-get_clusters)

- [`ptfce_spec$clone()`](#method-ptfce_spec-clone)

------------------------------------------------------------------------

### `ptfce_spec$new()`

Create a new ptfce_spec object

#### Usage

    ptfce_spec$new(
      gfeat_dir = NULL,
      zstat_numbers = NULL,
      fsl_smoothest_file = NULL,
      dof = NULL,
      residuals_file = NULL,
      z_files = NULL,
      mask_files = NULL,
      fwe_p = 0.05,
      two_sided = TRUE,
      write_thresh_imgs = TRUE,
      scheduler = NULL,
      time_per_zstat = NULL,
      memgb_per_command = NULL
    )

#### Arguments

- `gfeat_dir`:

  One or more .gfeat folders containing a higher-level FSL analysis.
  These will be used for zstat images, mask files, and fsl residual
  smoothness estimates.

- `zstat_numbers`:

  if a `gfeat_dir` is used, a vector of zstat numbers can also be
  provided to subset the zstat images that are used in pTFCE correction.
  Ignored if gfeat_dir is not provided.

- `fsl_smoothest_file`:

  The smoothness file created by the smoothest command in FSL. Created
  by FEAT automatically or can be run manually with smoothest on the
  res4d file. Used by ptfce to compute resels.

- `dof`:

  The degrees of freedom for the test of interest. Found in the dof file
  created by FSL and used by ptfce for FWE correction.

- `residuals_file`:

  The residuals file from the group analysis used by ptfce for FWE
  correction. In FSL, this is the res4d.nii.gz file created by FEAT.

- `z_files`:

  A vector of z-statistic filenames that should be corrected using pTFCE

- `mask_files`:

  A vector of mask filenames that correspond to `z_files`. If this is of
  length 1, then the mask file will be recycled for all zstat images.

- `fwe_p`:

  A vector of p-values for which z-statistic thresholds will be
  calculated.

- `two_sided`:

  If TRUE, p-values for `fwe_p` are treated as two-tailed (i.e.,
  p-values are divided by 2 in the TFCE z-threshold calculation.)

- `write_thresh_imgs`:

  If TRUE, then pTFCE thresholds for each fwe_p will be applied to the
  TFCE image and saved to the same folder as the z-statistic. These
  thresholded files let you look at the map at a given FWE threshold.

- `scheduler`:

  Which scheduler to use for submitting jobs. Options are 'local',
  'slurm', and 'torque'.

- `time_per_zstat`:

  The amount of time to budget for each zstat to run through pTFCE in
  dd-hh:mm:ss format. Default is 00:00:45 (45 seconds).

- `memgb_per_command`:

  How many GB of memory should be requested for each pTFCE command/run.
  If not provided, defaults to 1.

------------------------------------------------------------------------

### `ptfce_spec$get_ptfce_calls()`

method to return calls to external ptfce_zstat.R script for each zstat

#### Usage

    ptfce_spec$get_ptfce_calls(include_complete = FALSE)

#### Arguments

- `include_complete`:

  if TRUE, return calls for zstats that already appear to have
  pTFCE-corrected images in place. Default: FALSE.

------------------------------------------------------------------------

### `ptfce_spec$get_expected_files()`

return the vector of expected output files

#### Usage

    ptfce_spec$get_expected_files()

------------------------------------------------------------------------

### `ptfce_spec$run()`

Run each pending map sequentially in the current R session.

#### Usage

    ptfce_spec$run(force = FALSE)

#### Arguments

- `force`:

  if TRUE, re-run pTFCE on an existing output method to submit all
  ptfce_zstat.R calls to a cluster based on the scheduler specified

------------------------------------------------------------------------

### `ptfce_spec$submit()`

#### Usage

    ptfce_spec$submit(force = FALSE)

#### Arguments

- `force`:

  if TRUE, re-run pTFCE for zstat images that already appear to have
  pTFCE-corrected outputs in place

------------------------------------------------------------------------

### `ptfce_spec$is_complete()`

returns `TRUE` when all expected outputs exist and their threshold
tables match the requested settings

#### Usage

    ptfce_spec$is_complete()

------------------------------------------------------------------------

### `ptfce_spec$get_clusters()`

for each input file, obtain 3dClusterize objects that reflect the
pTFCE-corrected clusters

#### Usage

    ptfce_spec$get_clusters(
      fwep = 0.05,
      clust_nvox = 10,
      NN = 1L,
      add_whereami = TRUE,
      whereami_atlases = NULL
    )

#### Arguments

- `fwep`:

  the whole-brain familywise error rate to use (often .05)

- `clust_nvox`:

  The minimum number of voxels to allow in a given cluster

- `NN`:

  The cluster definition in AFNI terms. 1 = faces touch, 2 = edges
  touch, 3 = corners touch

- `add_whereami`:

  if TRUE, lookup labels for each cluster

- `whereami_atlases`:

  The atlases to request in the whereami lookup. If NULL, it uses the
  defaults

#### Details

Note that even though pTFCE enhances clusters, you still see some very
small clusters in some cases. Hence, there is no requirement to have
clust_nvox \> 1, but it may be a good idea for your sanity.

------------------------------------------------------------------------

### `ptfce_spec$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ptfce_spec$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
