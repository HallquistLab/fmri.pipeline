test_that("pTFCE threshold lookup preserves sidedness for new alpha values", {
  skip_if_not_installed("pTFCE", "0.4.0")
  for (sided in c(TRUE, FALSE)) {
    tab <- data.frame(p_value = 0.05, two_sided = sided, z_ptfce = 4, number_of_resels = 100)
    for (alpha in c(0.01, 0.1)) {
      result <- fmri.pipeline:::ptfce_threshold_from_table(tab, alpha)
      expect_equal(result$z, pTFCE::fwe.p2z(100, FWEP = if (sided) alpha / 2 else alpha))
      expect_identical(result$two_sided, sided)
    }
    expect_equal(fmri.pipeline:::ptfce_threshold_from_table(tab, 0.05)$z, 4)
  }
})

test_that("pTFCE uses correction volume and rejects incompatible masks", {
  root <- tempfile("ptfce_mask_")
  dir.create(root)
  z <- write_fwe_test_image(file.path(root, "z.nii.gz"))
  mask <- file.path(root, "roi.nii.gz")
  region <- array(0, c(5, 5, 5))
  region[2:4, 2:4, 2:4] <- 1
  RNifti::writeNifti(region, mask)
  smoothness <- file.path(root, "smoothness")
  writeLines(c("DLH 0.1", "VOLUME 125", "RESELS 10"), smoothness)
  inside <- fmri.pipeline:::ptfce_read_mask(mask, z)
  expect_equal(fmri.pipeline:::ptfce_smoothness(smoothness, inside), list(V = 27L, Rd = 2.7, resels = 10))
  write_fwe_test_image(mask, 1, c(4, 5, 5))
  expect_error(fmri.pipeline:::ptfce_read_mask(mask, z), "geometry")
  write_fwe_test_image(mask, 0)
  expect_error(fmri.pipeline:::ptfce_read_mask(mask, z), "nonempty")
  write_fwe_test_image(mask, 0.5)
  expect_error(fmri.pipeline:::ptfce_read_mask(mask, z), "binary")
})

test_that("real pTFCE matches direct calculation and writes distinct small-alpha maps", {
  skip_if_not_installed("pTFCE", "0.4.0")
  root <- tempfile("ptfce_real_")
  dir.create(root)
  values <- array(0, c(7, 7, 7))
  values[2:3, 2:3, 2:3] <- 5
  values[5:6, 5:6, 5:6] <- -4
  zfile <- file.path(root, "z.nii.gz")
  mask <- write_fwe_test_image(file.path(root, "mask.nii.gz"), 1, c(7, 7, 7))
  RNifti::writeNifti(values, zfile)
  smoothness <- file.path(root, "smoothness")
  writeLines(c("DLH 0.1", "VOLUME 343", "RESELS 10"), smoothness)
  out <- file.path(root, "out")
  alphas <- c(0.05, 0.0005, 0.0002, 0.0001)
  thresholds <- fmri.pipeline:::run_ptfce_image(zfile, mask, out, fwe_alpha = alphas, fsl_smoothest = smoothness)
  reference <- pTFCE::ptfce(zfile, mask, V = 343, Rd = 34.3, resels = 10,
    sided = "two.sided", fwe_p = alphas, verbose = FALSE)
  expected <- reference$Z
  actual <- RNifti::readNifti(file.path(out, "z_ptfce.nii.gz"))
  expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-6)
  expect_equal(thresholds$z_ptfce, vapply(alphas, function(a) pTFCE::fwe.p2z(34.3, a / 2), numeric(1)))
  for (i in seq_along(alphas)) {
    name <- paste0("z_ptfce_fwep_", fmri.pipeline:::format_fwe_probability(alphas[i]), ".nii.gz")
    expected_thresholded <- expected
    expected_thresholded[abs(expected) < thresholds$z_ptfce[i]] <- 0
    expect_equal(as.numeric(RNifti::readNifti(file.path(out, name))), as.numeric(expected_thresholded), tolerance = 1e-6)
  }
  expect_equal(length(unique(list.files(out, pattern = "fwep"))), length(alphas))
  # A completed map with no discoveries is a valid result.
  empty <- file.path(root, "empty.nii.gz")
  write_fwe_test_image(empty, 0, c(7, 7, 7))
  zero <- fmri.pipeline:::run_ptfce_image(empty, mask, file.path(root, "zero"), fsl_smoothest = smoothness)
  expect_equal(zero$significant_voxels, 0L)
})

test_that("the built-in worker succeeds in a clean subprocess with no discoveries", {
  skip_if_not_installed("pTFCE", "0.4.0")
  gpa <- make_fwe_execution_test_gpa(tempfile("ptfce_process_"))
  plan <- plan_fwe_correction(gpa, fwe_spec("builtin", list(l3_cope_name = "age")), source = "setup")
  run <- run_fwe_plan(plan, scheduler = "local")
  expect_identical(run$execution$exit_status, 0L)
  expect_true(run$result$complete)
  table <- read.csv(run$result$manifest$artifact_file[run$result$manifest$artifact_role == "threshold_table"])
  expect_equal(table$significant_voxels, 0L)
  expect_equal(table$search_voxels, 125L)
})

test_that("the worker rejects masks outside model support and nonfinite input", {
  skip_if_not_installed("pTFCE", "0.4.0")
  root <- tempfile("ptfce_bad_support_")
  dir.create(root)
  z <- write_fwe_test_image(file.path(root, "z.nii.gz"))
  roi <- write_fwe_test_image(file.path(root, "roi.nii.gz"), 1)
  values <- array(1, c(5, 5, 5))
  values[1, 1, 1] <- 0
  model <- file.path(root, "model.nii.gz")
  RNifti::writeNifti(values, model)
  expect_error(fmri.pipeline:::run_ptfce_image(z, roi, model_mask = model), "outside the model mask")
  values[1, 1, 1] <- NaN
  RNifti::writeNifti(values, z)
  expect_error(fmri.pipeline:::run_ptfce_image(z, roi), "Non-finite Z")
})

test_that("the directory CLI runs all contrasts locally and honors positive-only inference", {
  skip_if_not_installed("pTFCE", "0.4.0")
  group <- tempfile("group with spaces", fileext = ".gfeat")
  for (i in 1:2) {
    feat <- file.path(group, paste0("cope", i, ".feat"))
    stats <- file.path(feat, "stats")
    write_fwe_test_image(file.path(feat, "mask.nii.gz"), 1)
    z <- array(-8, c(5, 5, 5))
    z[2:3, 2:3, 2:3] <- 5
    dir.create(stats)
    RNifti::writeNifti(z, file.path(stats, "zstat1.nii.gz"))
    writeLines(c("DLH 0.1", "VOLUME 125", "RESELS 10"), file.path(stats, "smoothness"))
  }
  fmri.pipeline:::ptfce_cli(c("--gfeat", group, "--onesided", "--fwep", "0.05"))
  spec <- ptfce_spec$new(gfeat_dir = group, two_sided = FALSE, write_thresh_imgs = FALSE)
  expect_true(spec$is_complete())
  expect_length(spec$get_ptfce_calls(), 0L)
  signed <- ptfce_spec$new(gfeat_dir = group, two_sided = TRUE, write_thresh_imgs = FALSE)
  expect_false(signed$is_complete())
  expect_length(signed$get_ptfce_calls(), 2L)
  for (map in spec$get_expected_files()) {
    value <- RNifti::readNifti(map$z_ptfce)
    expect_equal(as.numeric(value[1, 1, 1]), 0)
    expect_gt(as.numeric(value[2, 2, 2]), 0)
    thresholds <- read.csv(map$z_csv)
    expect_false(thresholds$two_sided)
    expect_equal(thresholds$z_ptfce, pTFCE::fwe.p2z(12.5, 0.05))
  }
})

test_that("native smoothness estimation uses residuals and degrees of freedom", {
  skip_if_not_installed("pTFCE", "0.4.0")
  withr::local_seed(2026)
  root <- tempfile("ptfce_residual_")
  dims <- c(9L, 9L, 9L)
  mask <- write_fwe_test_image(file.path(root, "mask.nii"), 1, dims)
  z <- file.path(root, "z.nii")
  residual <- file.path(root, "residual.nii")
  # Adjacent slices share signal so the smoothness estimator has a positive spatial correlation.
  smooth_noise <- function() {
    x <- array(rnorm(prod(dims)), dims)
    ix <- c(1L, seq_len(8L))
    x + x[ix,,] + x[,ix,] + x[,,ix]
  }
  RNifti::writeNifti(smooth_noise(), z)
  RNifti::writeNifti(array(replicate(12L, smooth_noise()), c(dims, 12L)), residual)
  for (use_residual in c(FALSE, TRUE)) {
    resid <- if (use_residual) residual
    out <- file.path(root, as.character(use_residual))
    table <- fmri.pipeline:::run_ptfce_image(z, mask, out, residuals = resid, dof = 10)
    reference <- pTFCE::ptfce(z, mask, sided = "two.sided", residual = resid, dof = 10,
      verbose = FALSE)
    expect_equal(table$number_of_resels, reference$number_of_resels)
    expect_equal(as.numeric(RNifti::readNifti(file.path(out, "z_ptfce.nii"))),
      as.numeric(reference$Z), tolerance = 1e-6)
  }
})
