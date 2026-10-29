# Extract FSL GLM statistics within spatially aligned masks

Extract FSL GLM statistics within spatially aligned masks

## Usage

``` r
extract_glm_betas_in_mask(
  gpa,
  mask_files,
  what = c("cope", "zstat"),
  out_dir = getwd(),
  extract_l1 = "all",
  extract_l2 = "all",
  extract_l3 = "all",
  ncores = NULL,
  scheduler = "local",
  aggregate = TRUE,
  aggFUN = mean,
  return_data = TRUE,
  write_data = TRUE
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object for which results of GLM analyses
  are already available

- mask_files:

  vector of filenames for NIfTI-based atlases. The function will loop
  over each and extract coefficients from each unique value in each
  mask.

- what:

  which statistics to extract from each parcel. Default is 'cope' (aka
  'beta') and 'zstat'.

- out_dir:

  the directory to which statistics are written as .csv.gz files.

- extract_l1:

  a character vector of l1 models from which to extract run-level
  (level 1) statistics. If "none", l1 statistic extraction will be
  skipped. If "all", then all l1 models will be extracted.

- extract_l2:

  a character vector of l2 models from which to extract subject-level
  (level 2) statistics. If "none", l2 statistic extraction will be
  skipped. If "all", then all l2 models will be extracted (you will
  still be prompted for which l1 models you want to extract from each l2
  model).

- extract_l3:

  a character vector of l3 models from which to extract group-level
  (level 3) statistics. If "none", l3 statistic extraction will be
  skipped. If "all", then all l3 models will be extracted (you will
  still be prompted for which l1 and l2 models you want to extract from
  each l3 model).

- ncores:

  The number of cores to use for extracting statistics from lower-level
  imgages. If NULL, lookup value in
  gpa\$parallel\$extract_glm_betas_ncores. If this value is \> 1 and
  scheduler = "local", then mclapply will be used locally to extract.

- scheduler:

  The scheduler to use for extracting statistics. If "local", use
  lapply/mclapply within the current compute session. If 'slurm', use
  doFuture with multiple slurm jobs to extract.

- aggregate:

  whether to take the average (or other central tendency measure) of
  voxels within a given mask value. This only pertains to integer-valued
  masks, not continuous ones.

- aggFUN:

  the function used to aggregate statistics for voxels within a given
  mask value. Default is mean.

- return_data:

  if TRUE, then extracted statistics will be returned as a list of
  data.frames with elements l1, l2, and l3.

- write_data:

  if TRUE, then extracted statistics will be written as .csv.gz files to
  `out_dir`.

## Details

Every statistic image must match the mask's spatial dimensions, voxel
sizes, spatial units, and voxel-to-world transform (absolute tolerance
1e-4). The sform is preferred, with qform and voxel-size fallbacks, and
is also used for output coordinates. A mismatch stops extraction; masks
are not automatically resampled. Group-level rows have a missing subject
id and retain their setup session, including 0 for pooled sessions.
