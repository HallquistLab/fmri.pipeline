# Lookup FSL FEAT output images by model and contrast

Build a tidy lookup table that maps the human-readable pipeline model
and contrast names onto the FEAT folder structure for level 1, 2, and/or
3 FSL outputs. The default result separates the structurally different
FEAT levels into `$l1`, `$l2`, and `$l3` tables with stable
level-specific schemas; levels that were not requested are present as
empty tables.

## Usage

``` r
lookup_feat_outputs(
  gpa,
  level = c(1L, 2L, 3L),
  what = c("cope", "varcope", "zstat", "tstat"),
  include_missing = TRUE,
  include_internal = FALSE,
  source = c("auto", "setup", "cache", "filesystem"),
  cache_dir = NULL,
  refresh_status = FALSE,
  lg = NULL,
  format = c("by_level", "long")
)

# S3 method for class 'feat_output_lookup'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'feat_output_lookup'
print(x, ...)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object. The function is most informative
  when populated FSL setup tables are present, but `source = "auto"` can
  also search scheduler caches and crawl FEAT folders.

- level:

  integer vector selecting FEAT levels to lookup. Valid values are `1`,
  `2`, and `3`. Defaults to all levels.

- what:

  statistic image types to include. Valid values are `"cope"`,
  `"varcope"`, `"zstat"`, and `"tstat"`.

- include_missing:

  if `TRUE`, include expected output paths even when the image does not
  exist. If `FALSE`, return only existing images.

- include_internal:

  if `TRUE`, retain setup/debug columns such as FEAT execution
  timestamps and raw FEAT status fields. The default `FALSE` returns
  compact, level-specific user-facing tables.

- source:

  where to look for output metadata. `"setup"` uses only the
  `gpa$l*_model_setup$fsl` tables, `"cache"` searches scheduler
  `run_pipeline_cache*.RData` files, `"filesystem"` crawls FEAT folders,
  and `"auto"` tries these in that order.

- cache_dir:

  optional directory containing scheduler batch caches. If `NULL`,
  caches are searched below `gpa$output_directory` and
  `gpa$output_locations$scheduler_scripts`.

- refresh_status:

  if `TRUE`, refresh the FEAT status columns in the setup table before
  building the lookup.

- lg:

  optional
  [`lgr::Logger`](https://s-fleck.github.io/lgr/reference/Logger.html)
  object.

- format:

  return format. `"by_level"` returns a `feat_output_lookup` object with
  `$l1`, `$l2`, and `$l3` tables. `"long"` returns their combined
  data.frame form.

- x:

  a `feat_output_lookup` object.

- row.names, optional:

  standard `as.data.frame` arguments.

- ...:

  unused.

## Value

With `format = "by_level"`, a `feat_output_lookup` object containing
level-specific `$l1`, `$l2`, and `$l3` data.frames. With
`format = "long"`, a data.frame with one row per analysis instance,
contrast, and statistic image. The long form adds `output_model`,
`output_cope_number`, and `output_cope_name` as level-independent
convenience columns.

## Details

L2 passthroughs remain logical L2 rows (`output_level = 2`) but are
marked with `is_passthrough = TRUE`, `materialized_level = 1`, and
`analysis_status = "passthrough"`. Because no L2 FEAT analysis is
materialized, `analysis_fsf` and `analysis_dir` are missing and
`image_feat_dir` identifies the originating L1 FEAT directory. A
passthrough has only a `"cope"` image row.

Higher-level rows report `feat_requested_estimator` and
`feat_actual_estimator` (`ols`, `fe`, `flame1`, `flame12`, or
`randomise`). A partial FLAME1 fallback is `flame12+flame1`;
`feat_flame1_fallback_count` counts affected slices.
`feat_requested_outlier_deweighting` and
`feat_actual_outlier_deweighting` report the FSF request and the
retained fit's outlier option. `feat_auto_retried` identifies a
whole-model outlier retry. Actual settings are missing for unvalidated
fits or when execution evidence is unavailable. Use
`refresh_status = TRUE` to update a setup table saved before estimation.
