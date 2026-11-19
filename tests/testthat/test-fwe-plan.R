make_fwe_plan_test_gpa <- function(root, create_mask = TRUE,
                                   create_smoothness = TRUE,
                                   scheduler = "slurm") {
  l3_feat <- file.path(
    root, "feat_l3", "L3m-l3_group", "L1m-facehouse",
    "L2m-l2_session", "l2c-overall", "FEAT_l1c-EV_face.gfeat"
  )
  image_feat_dir <- file.path(l3_feat, "cope1.feat")
  stats_dir <- file.path(image_feat_dir, "stats")
  dir.create(stats_dir, recursive = TRUE)
  for (file in file.path(stats_dir, c("zstat1.nii.gz", "zstat2.nii.gz"))) write_fwe_test_image(file)
  if (isTRUE(create_mask)) write_fwe_test_image(file.path(image_feat_dir, "mask.nii.gz"), 1)
  if (isTRUE(create_smoothness)) {
    writeLines(c("DLH 0.1", "VOLUME 1000", "RESELS 100"), file.path(stats_dir, "smoothness"))
  }

  l3_setup <- data.frame(
    l1_model = "facehouse",
    l1_cope_name = "EV_face",
    l2_model = "l2_session",
    l2_cope_name = "overall",
    l3_model = "l3_group",
    session = 0L,
    l3_input_mode = "all_sessions",
    feat_fsf = sub("\\.gfeat$", ".fsf", l3_feat),
    feat_dir = l3_feat,
    feat_complete = TRUE,
    feat_failed = FALSE,
    feat_dir_exists = TRUE,
    feat_fsf_exists = FALSE,
    stringsAsFactors = FALSE
  )
  l3_metadata <- data.frame(
    id = c("01", "01"),
    session = c(1L, 1L),
    l1_model = c("facehouse", "facehouse"),
    l1_cope_number = c(1L, 1L),
    l1_cope_name = c("EV_face", "EV_face"),
    l2_model = c("l2_session", "l2_session"),
    l2_cope_number = c(1L, 1L),
    l2_cope_name = c("overall", "overall"),
    l3_model = c("l3_group", "l3_group"),
    l3_cope_number = c(1L, 2L),
    l3_cope_name = c("group", "age"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      output_directory = root,
      scheduler = scheduler,
      lgr_threshold = "warn",
      multi_run = TRUE,
      l3_model_setup = structure(
        list(metadata = l3_metadata, fsl = l3_setup),
        class = c("l3_setup", "list")
      )
    ),
    class = c("glm_pipeline_arguments", "list")
  )
}

test_that("pTFCE planning fans out over all selected L3 contrasts", {
  root <- tempfile("fwe_plan_")
  gpa <- make_fwe_plan_test_gpa(root)
  spec <- fwe_spec(
    "all_ptfce",
    targets = list(l3_model = "l3_group"),
    fwe_alpha = c(0.05, 0.01)
  )

  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  expect_s3_class(plan, "fwe_plan")
  expect_invisible(fmri.pipeline:::validate_fwe_plan(plan, require_inputs = TRUE))
  expect_identical(plan$method, "ptfce")
  expect_identical(plan$scheduler, "slurm")
  expect_equal(nrow(plan$analyses), 1L)
  expect_equal(nrow(plan$tasks), 2L)
  expect_setequal(plan$tasks$l3_cope_name, c("group", "age"))
  expect_true(all(plan$tasks$status == "ready"))
  expect_true(all(vapply(plan$tasks$fwe_alpha, identical, logical(1L), c(0.05, 0.01))))
  expect_equal(nrow(plan$artifacts), 8L)
  expect_equal(
    table(plan$artifacts$artifact_role),
    table(c(rep("ptfce_statistic", 2L), rep("threshold_table", 2L),
            rep("thresholded_statistic", 4L)))
  )
  expect_false(dir.exists(plan$output_directory))
  expect_false(any(dir.exists(plan$tasks$output_directory)))
})

test_that("pTFCE planning applies a custom correction mask to every task", {
  root <- tempfile("fwe_plan_mask_")
  gpa <- make_fwe_plan_test_gpa(root)
  custom_mask <- file.path(root, "custom_mask.nii.gz")
  write_fwe_test_image(custom_mask, 1)
  spec <- fwe_spec(
    "custom_mask",
    targets = list(l3_model = "l3_group"),
    correction_mask = custom_mask
  )

  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  expect_identical(
    unique(plan$tasks$correction_mask_file),
    as.character(fs::path_abs(custom_mask))
  )
  expect_true(all(plan$tasks$status == "ready"))
  expect_true(all(plan$analyses$model_mask_exists))
})

test_that("pTFCE planning reports missing inputs or returns a blocked preview", {
  root <- tempfile("fwe_plan_blocked_")
  gpa <- make_fwe_plan_test_gpa(root, create_smoothness = FALSE)
  spec <- fwe_spec("blocked", targets = list(l3_model = "l3_group"))

  expect_error(
    plan_fwe_correction(gpa, spec, source = "setup"),
    "Missing required pTFCE input"
  )

  plan <- plan_fwe_correction(
    gpa, spec, source = "setup", require_inputs = FALSE
  )
  expect_true(all(plan$tasks$status == "blocked"))
  expect_true(all(grepl("smoothness_file", plan$tasks$missing_inputs, fixed = TRUE)))
  expect_invisible(fmri.pipeline:::validate_fwe_plan(plan))
  expect_error(
    fmri.pipeline:::validate_fwe_plan(plan, require_inputs = TRUE),
    "blocked tasks"
  )
})

test_that("compiled FWE plans detect completed artifacts", {
  root <- tempfile("fwe_plan_complete_")
  gpa <- make_fwe_plan_test_gpa(root)
  spec <- fwe_spec(
    "no_threshold_images",
    targets = list(l3_cope_name = "age"),
    write_threshold_images = FALSE
  )
  preview <- plan_fwe_correction(gpa, spec, source = "setup")
  expect_equal(nrow(preview$artifacts), 2L)

  dir.create(preview$tasks$output_directory, recursive = TRUE)
  write_fwe_test_result(preview)
  complete <- plan_fwe_correction(gpa, spec, source = "setup")

  expect_identical(complete$tasks$status, "complete")
  expect_true(all(complete$artifacts$exists))
})

test_that("compiled FWE plans round trip through RDS", {
  root <- tempfile("fwe_plan_rds_")
  gpa <- make_fwe_plan_test_gpa(root)
  spec <- fwe_spec("roundtrip", targets = list(l3_cope_name = "group"))
  plan <- plan_fwe_correction(gpa, spec, source = "setup")
  path <- file.path(root, "plan.rds")

  written <- write_fwe_plan(plan, path)
  restored <- read_fwe_plan(written)

  expect_s3_class(restored, "fwe_plan")
  expect_equal(restored, plan)
  expect_error(write_fwe_plan(plan, path), "Refusing to overwrite")
  expect_silent(write_fwe_plan(plan, path, overwrite = TRUE))
})

test_that("planning fails clearly when a method adapter is unavailable", {
  root <- tempfile("fwe_plan_adapter_")
  gpa <- make_fwe_plan_test_gpa(root)
  spec <- fwe_spec(
    "clustsim",
    targets = list(l3_cope_name = "group"),
    method = "afni_3dclustsim_permutation"
  )

  expect_error(
    plan_fwe_correction(gpa, spec, source = "setup"),
    "No FWE planning adapter"
  )
})
