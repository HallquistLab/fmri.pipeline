test_that("legacy pre-PR-42 MEDuSA exports emit deprecation warnings pointing to medusa", {
  legacy_fns <- c(
    "voxelwise_deconvolution",
    "deconvolve_nlreg",
    "deconvolve_nlreg_resample",
    "deconvolve_reglin",
    "event_lock_ts",
    "get_medusa_interpolated_ts",
    "get_medusa_compression_score",
    "compress_mts_pca",
    "evt_align_decon_files",
    "run_decon_alignment",
    "mixed_by",
    "meta_mixed_by",
    "jlmer"
  )

  # Verify all 14 legacy symbols are exported from fmri.pipeline
  pkg_exports <- getNamespaceExports("fmri.pipeline")
  expect_true(all(c(legacy_fns, "fmri_ts") %in% pkg_exports))

  # Verify post-PR-42 MEDuSA symbols are NOT exported from fmri.pipeline
  post_pr42_fns <- c(
    "hrf_library",
    "estimate_crossrun_hrf",
    "deconvolve_crossrun_hrf",
    "decon_epoch_regression",
    "diagnose_crossrun_recovery",
    "deconvolution_impulse_diagnostics",
    "decon_connectivity_repair",
    "medusa_between_subject_repair",
    "medusa_fit_mmom",
    "medusa_multilag_infer"
  )
  expect_false(any(post_pr42_fns %in% pkg_exports))

  # Mock internal lookup helpers to test forwarding of all 13 legacy functions
  called_with <- NULL
  local_mocked_bindings(
    medusa_is_installed = function() TRUE,
    medusa_get_export = function(name) {
      force(name)
      function(...) {
        called_with <<- list(name = name, args = list(...))
        "forwarded_ok"
      }
    },
    .package = "fmri.pipeline"
  )

  for (fn_name in legacy_fns) {
    fn <- getExportedValue("fmri.pipeline", fn_name)
    expect_warning(
      res <- fn(x = 1, y = "two"),
      regexp = paste0("medusa::", fn_name),
      class = "deprecatedWarning"
    )
    expect_identical(res, "forwarded_ok")
    expect_identical(called_with$name, fn_name)
    expect_identical(called_with$args, list(x = 1, y = "two"))
  }

  # Also verify fmri_ts$new(...)
  local_mocked_bindings(
    medusa_is_installed = function() TRUE,
    medusa_get_export = function(name) {
      expect_identical(name, "fmri_ts")
      list(new = function(...) list(constructed = TRUE, args = list(...)))
    },
    .package = "fmri.pipeline"
  )
  expect_warning(
    ts_obj <- fmri.pipeline::fmri_ts$new(tr = 1.0),
    regexp = "medusa::fmri_ts",
    class = "deprecatedWarning"
  )
  expect_true(isTRUE(ts_obj$constructed))
  expect_identical(ts_obj$args, list(tr = 1.0))

  # Verify informative error when medusa is not installed
  local_mocked_bindings(
    medusa_is_installed = function() FALSE,
    .package = "fmri.pipeline"
  )
  expect_warning(
    expect_error(
      fmri.pipeline::mixed_by(data.frame(y = 1)),
      regexp = "not installed"
    ),
    class = "deprecatedWarning"
  )
  expect_warning(
    expect_error(
      fmri.pipeline::fmri_ts$new(tr = 1.0),
      regexp = "not installed"
    ),
    class = "deprecatedWarning"
  )
})

test_that("deprecated fmri.pipeline wrappers forward live to installed medusa", {
  skip_if_not_installed("medusa")

  ts_df <- data.frame(id = 1L, run = 1L, time = seq(0, 4, by = 1), v1 = rnorm(5))
  expect_warning(
    obj <- fmri.pipeline::fmri_ts$new(ts_data = ts_df, vm = list(value = "v1"), tr = 1.0),
    regexp = "medusa::fmri_ts",
    class = "deprecatedWarning"
  )
  expect_s3_class(obj, "fmri_ts")
})
