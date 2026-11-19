test_that("ClustSim requires exact thresholds and applies its simulation mask", {
  root <- tempfile("clustsim_mask_")
  dir.create(root)
  mask <- write_fwe_test_image(file.path(root, "mask.nii.gz"), 1)
  statistic <- write_fwe_test_image(file.path(root, "z.nii.gz"))
  sim <- afni_3dclustsim$new(out_dir = root, clustsim_mask = mask, acf_params = c(0.5, 2, 4))
  private <- sim$.__enclos_env__$private
  private$clustsim_df <- data.frame(nn = 1L, sided = "bi", pthr = c(0.001, 0.01), athr = 0.05, nvoxels = c(10.2, 40))
  private$clustsim_complete <- TRUE
  args <- NULL
  testthat::local_mocked_bindings(afni_3dclusterize = list(new = function(...) {
    args <<- list(...)
    list(run = function(...) invisible(NULL))
  }), .package = "fmri.pipeline")
  sim$apply_clustsim(statistic, add_whereami = FALSE)
  expect_identical(normalizePath(args$mask), normalizePath(mask))
  expect_equal(args$clust_nvox, 11)
  expect_equal(args$upper_thresh, qnorm(0.001 / 2, lower.tail = FALSE))
  expect_error(sim$apply_clustsim(statistic, pthr = 0.0005, add_whereami = FALSE), "Specified pthr")
  expect_error(sim$apply_clustsim(statistic, athr = 0.01, add_whereami = FALSE), "Specified athr")
})
