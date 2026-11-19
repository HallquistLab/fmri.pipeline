test_that("pTFCE threshold lookup preserves sidedness for new alpha values", {
  skip_if_not_installed("pTFCE")
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
  skip_if_not_installed("pTFCE")
  root <- tempfile("ptfce_real_")
  dir.create(root)
  values <- array(0, c(7, 7, 7))
  values[2:3, 2:3, 2:3] <- 5
  values[5:6, 5:6, 5:6] <- -5
  zfile <- file.path(root, "z.nii.gz")
  mask <- write_fwe_test_image(file.path(root, "mask.nii.gz"), 1, c(7, 7, 7))
  RNifti::writeNifti(values, zfile)
  smoothness <- file.path(root, "smoothness")
  writeLines(c("DLH 0.1", "VOLUME 343", "RESELS 10"), smoothness)
  out <- file.path(root, "out")
  alphas <- c(0.05, 0.0005, 0.0002, 0.0001)
  thresholds <- fmri.pipeline:::run_ptfce_image(zfile, mask, out, fwe_alpha = alphas, fsl_smoothest = smoothness)
  Z <- oro.nifti::readNIfTI(zfile, reorient = FALSE)
  M <- oro.nifti::readNIfTI(mask, reorient = FALSE)
  reference <- function(z) {
    withr::local_package("oro.nifti")
    result <- pTFCE::ptfce(z, M, V = 343, Rd = 34.3, resels = 10, verbose = FALSE,
      logpmax = -pnorm(max(z), lower.tail = FALSE, log.p = TRUE))$Z
    result[result < 0] <- 0
    result
  }
  expected <- reference(Z) - reference(-Z)
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
  skip_if_not_installed("pTFCE")
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
  skip_if_not_installed("pTFCE")
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
