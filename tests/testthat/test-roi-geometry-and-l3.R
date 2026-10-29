roi_group_fixture <- function() {
  root <- tempfile("roi_group_")
  dir.create(root)
  mask <- file.path(root, "mask.nii.gz")
  RNifti::writeNifti(array(1, c(2, 2, 2)), mask)
  setup <- data.frame(l1_model = "task", l1_cope_name = c("a", "b"),
    l2_model = "subject", l2_cope_name = "mean", l3_model = "group", session = 0L,
    feat_dir = file.path(root, c("a.gfeat", "b.gfeat")),
    feat_dir_exists = TRUE, feat_complete = TRUE)
  metadata <- tidyr::crossing(l1_cope_name = c("a", "b"), l3_cope_number = 1:2)
  metadata$l1_model <- "task"
  metadata$l1_cope_number <- match(metadata$l1_cope_name, c("a", "b"))
  metadata$l2_model <- "subject"
  metadata$l2_cope_name <- "mean"
  metadata$l2_cope_number <- 1L
  metadata$l3_model <- "group"
  metadata$l3_cope_name <- c("overall", "slope")[metadata$l3_cope_number]
  # Repeated subject metadata must not duplicate group maps in the results.
  metadata <- dplyr::bind_rows(transform(metadata, id = "s1"), transform(metadata, id = "s2"))
  for (i in 1:2) {
    directory <- file.path(setup$feat_dir[i], "cope1.feat", "stats")
    dir.create(directory, recursive = TRUE)
    for (stat in c("cope", "varcope", "zstat", "tstat")) {
      for (j in 1:2) RNifti::writeNifti(array(i * 10 + j, c(2, 2, 2)),
        file.path(directory, paste0(stat, j, ".nii.gz")))
    }
  }
  gpa <- structure(list(lgr_threshold = "fatal", multi_run = TRUE, glm_software = "fsl",
    subject_data = data.frame(id = c("s1", "s2"), session = 1L, exclude_subject = FALSE),
    l3_model_setup = structure(list(fsl = setup, metadata = metadata), class = c("l3_setup", "list"))),
    class = c("glm_pipeline_arguments", "list"))
  list(root = root, mask = mask, gpa = gpa,
    extract = data.frame(l1_model = "task", l2_model = "subject", l3_model = "group"))
}

extract_roi_group <- function(fixture, ...) {
  fmri.pipeline:::extract_fsl_betas(fixture$gpa, extract = fixture$extract,
    level = 3L, mask_file = fixture$mask, lg = lgr::get_logger("test/roi"), ...)
}

test_that("L3 extraction returns one row per group contrast and statistic", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  result <- extract_roi_group(f, what = c("cope", "varcope", "zstat", "tstat"))
  expect_equal(nrow(result), 16L)
  expect_true(all(is.na(result$id)))
  expect_true(all(result$session == 0L))
  expect_equal(sort(unique(result$value)), c(11, 12, 21, 22))
  expect_equal(sort(unique(result$l3_cope_name)), c("overall", "slope"))
  f$extract$l1_cope_name <- "b"
  f$extract$l3_cope_number <- 2L
  selected <- extract_roi_group(f, what = "cope")
  expect_equal(nrow(selected), 1L)
  expect_equal(selected$value, 22)
})

test_that("L3 extraction supports L1-to-L3 and group-only public extraction", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$gpa$multi_run <- FALSE
  f$gpa$l3_model_setup$fsl$l2_model <- NULL
  f$gpa$l3_model_setup$fsl$l2_cope_name <- NULL
  f$gpa$l3_model_setup$metadata <- f$gpa$l3_model_setup$metadata[,
    !grepl("^l2_", names(f$gpa$l3_model_setup$metadata))]
  f$extract$l2_model <- NULL
  expect_equal(nrow(extract_roi_group(f, what = "cope")), 4L)
  local_mocked_bindings(build_beta_extraction = function(...) list(l1 = NULL, l2 = NULL, l3 = f$extract),
    .package = "fmri.pipeline")
  result <- extract_glm_betas_in_mask(f$gpa, mask_files = f$mask,
    extract_l1 = "none", extract_l2 = "none", extract_l3 = "all", what = "cope",
    write_data = TRUE, out_dir = file.path(f$root, "extracted"))
  expect_equal(nrow(result$l3), 4L)
  expect_true(file.exists(file.path(f$root, "extracted", "L1m-task", "mask_cope_l3.csv.gz")))
})

test_that("ROI extraction rejects translation, orientation, spacing, and unit mismatches", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  original <- RNifti::readNifti(f$mask)
  for (kind in c("translation", "orientation", "spacing", "units", "qform")) {
    header <- RNifti::niftiHeader(original)
    header$sform_code <- 1L
    header$srow_x <- c(1, 0, 0, if (kind == "translation") 100 else 0)
    header$srow_y <- c(0, 1, 0, 0)
    header$srow_z <- c(0, 0, 1, 0)
    if (kind == "orientation") header$srow_x[1] <- -1
    if (kind == "spacing") header$pixdim[2:4] <- 2
    if (kind == "units") header$xyzt_units <- 2L
    if (kind == "qform") {
      header$sform_code <- 0L
      header$qform_code <- 1L
      header$qoffset_x <- 100
    }
    RNifti::writeNifti(RNifti::updateNifti(original, header), f$mask)
    expect_error(extract_roi_group(f, what = "cope"), "incompatible spatial geometry", info = kind)
  }
})

test_that("ROI checks all images and uses the checked affine for coordinates", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  # With four statistics there are sixteen images; this mismatch is past five.
  last <- file.path(f$gpa$l3_model_setup$fsl$feat_dir[2], "cope1.feat", "stats", "tstat2.nii.gz")
  img <- RNifti::readNifti(last)
  affine <- diag(4)
  affine[1, 4] <- 100
  attr(affine, "code") <- 1L
  RNifti::sform(img) <- affine
  RNifti::writeNifti(img, last)
  expect_error(extract_roi_group(f, what = c("cope", "varcope", "zstat", "tstat")),
    "incompatible spatial geometry")
  # Tiny header rounding differences are harmless; matching translated grids
  # must return coordinates derived from that same transform.
  mask <- RNifti::readNifti(f$mask)
  RNifti::sform(mask) <- affine
  RNifti::writeNifti(mask, f$mask)
  f$extract$l1_cope_name <- "b"
  f$extract$l3_cope_number <- 2L
  out <- extract_roi_group(f, what = "tstat", aggregate = FALSE)
  expect_equal(out$x, out$i - 1 + 100)
  expect_equal(out$y, out$j - 1)
  affine[1, 4] <- 100 + 1e-5
  RNifti::sform(mask) <- affine
  RNifti::writeNifti(mask, f$mask)
  expect_equal(nrow(extract_roi_group(f, what = "tstat")), 1L)
})

test_that("L3 extraction preserves session-specific group rows", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$gpa$l3_model_setup$fsl <- dplyr::bind_rows(
    transform(f$gpa$l3_model_setup$fsl, session = 1L),
    transform(f$gpa$l3_model_setup$fsl, session = 2L))
  result <- extract_roi_group(f, what = "cope")
  expect_equal(nrow(result), 8L)
  expect_equal(as.integer(table(result$session)), c(4L, 4L))
  f$extract$session <- 2L
  expect_true(all(extract_roi_group(f, what = "cope")$session == 2L))
})

test_that("L3 extraction can reconstruct contrast metadata from model definitions", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  f$gpa$l3_model_setup$metadata <- NULL
  f$gpa$l1_cope_names <- list(task = c("a", "b"))
  f$gpa$l3_models <- list(models = list(group = list(
    contrasts = matrix(1, 2, 1, dimnames = list(c("overall", "slope"), "Intercept")))))
  out <- extract_roi_group(f, what = "cope")
  expect_equal(nrow(out), 4L)
  expect_equal(sort(out$value), c(11, 12, 21, 22))
})

test_that("ROI extraction rejects mismatched dimensions and empty masks", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  RNifti::writeNifti(array(1, c(3, 2, 2)), f$mask)
  expect_error(extract_roi_group(f, what = "cope"), "incompatible spatial geometry")
  RNifti::writeNifti(array(0, c(2, 2, 2)), f$mask)
  expect_error(extract_roi_group(f, what = "cope"), "no nonzero voxels")
})

test_that("4D masks use their first volume but statistics must be single-volume", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  mask <- array(0, c(2, 2, 2, 2))
  mask[, , , 1] <- 1
  RNifti::writeNifti(mask, f$mask)
  expect_equal(nrow(extract_roi_group(f, what = "cope")), 4L)
  statistic <- file.path(f$gpa$l3_model_setup$fsl$feat_dir[1], "cope1.feat", "stats", "cope1.nii.gz")
  RNifti::writeNifti(array(1, c(2, 2, 2, 2)), statistic)
  expect_error(extract_roi_group(f, what = "cope"), "incompatible spatial geometry")
})

test_that("missing group statistic images yield empty public results", {
  f <- roi_group_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  unlink(list.files(f$root, pattern = "^cope[12]\\.nii\\.gz$", recursive = TRUE, full.names = TRUE))
  local_mocked_bindings(build_beta_extraction = function(...) list(l1 = NULL, l2 = NULL, l3 = f$extract),
    .package = "fmri.pipeline")
  out <- extract_glm_betas_in_mask(f$gpa, mask_files = f$mask,
    extract_l1 = "none", extract_l2 = "none", extract_l3 = "all", what = "cope",
    write_data = TRUE, out_dir = file.path(f$root, "extracted"))
  expect_equal(nrow(out$l3), 0L)
  expect_false(dir.exists(file.path(f$root, "extracted")))
})
