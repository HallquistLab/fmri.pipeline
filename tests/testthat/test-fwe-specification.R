make_fwe_spec_test_gpa <- function(root, create_images = c(TRUE, TRUE),
                                   complete = TRUE) {
  l3_feat <- file.path(
    root, "feat_l3", "L3m-l3_group", "L1m-facehouse",
    "L2m-l2_session", "l2c-overall", "FEAT_l1c-EV_face.gfeat"
  )
  stats_dir <- file.path(l3_feat, "cope1.feat", "stats")
  dir.create(stats_dir, recursive = TRUE)
  for (ii in seq_along(create_images)) {
    if (isTRUE(create_images[[ii]])) {
      file.create(file.path(stats_dir, paste0("zstat", ii, ".nii.gz")))
    }
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
    feat_complete = complete,
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

test_that("fwe_spec constructs and validates serializable method settings", {
  spec <- fwe_spec(
    name = "age_ptfce",
    targets = list(l3_model = "l3_group", l3_cope_name = "age"),
    method = "ptfce",
    fwe_alpha = c(0.05, 0.01)
  )

  expect_s3_class(spec, "fwe_spec")
  expect_identical(spec$schema_version, 1L)
  expect_identical(spec$level, 3L)
  expect_equal(spec$method$fwe_alpha, c(0.05, 0.01))
  expect_identical(spec$correction_mask, list(source = "model"))
  expect_invisible(validate_fwe_spec(spec))

  expect_error(
    fwe_spec("bad", list(unknown_model = "x")),
    "Unknown FWE target selector"
  )
  expect_error(
    fwe_spec("bad", list(l3_model = "x"), method = "ptfce", fwe_alpa = 0.05),
    "Unknown option"
  )
  expect_error(
    fwe_spec("bad", list(l3_model = "x"), level = 2L),
    "level = 3 only"
  )
})

test_that("3dClustSim specifications validate decision thresholds", {
  spec <- fwe_spec(
    "group_clustsim",
    targets = list(l3_model = "l3_group"),
    method = "3dclustsim_permutation",
    voxel_p = 0.005,
    cluster_alpha = 0.02
  )

  expect_identical(spec$method$name, "afni_3dclustsim_permutation")
  expect_equal(spec$method$voxel_p, 0.005)
  expect_error(
    fwe_spec(
      "bad_clustsim",
      targets = list(l3_model = "l3_group"),
      method = "afni_3dclustsim_permutation",
      voxel_p = 0.003
    ),
    "voxel_p must be included"
  )
})

test_that("FWE specifications round trip through YAML", {
  path <- file.path(tempdir(), paste0("fwe_spec_", basename(tempfile()), ".yaml"))
  spec <- fwe_spec(
    name = "age_ptfce",
    targets = list(
      l1_model = "facehouse",
      l3_model = "l3_group",
      l3_cope_name = c("group", "age")
    ),
    method = "ptfce",
    correction_mask = "/project/masks/brain.nii.gz",
    compute = list(scheduler = "slurm"),
    fwe_alpha = c(0.05, 0.01)
  )

  written <- write_fwe_spec(spec, path)
  restored <- read_fwe_spec(path)

  expect_true(file.exists(written))
  expect_s3_class(restored, "fwe_spec")
  expect_equal(unclass(restored), unclass(spec))
  expect_error(write_fwe_spec(spec, path), "Refusing to overwrite")
  expect_silent(write_fwe_spec(spec, path, overwrite = TRUE))
})

test_that("resolve_fwe_targets selects exact group contrasts", {
  root <- tempfile("fwe_targets_")
  gpa <- make_fwe_spec_test_gpa(root)

  age_spec <- fwe_spec(
    "age_ptfce",
    targets = list(
      l1_model = "facehouse",
      l2_cope_name = "overall",
      l3_model = "l3_group",
      l3_cope_name = "age"
    )
  )
  age <- resolve_fwe_targets(gpa, age_spec, source = "setup")

  expect_s3_class(age, "fwe_target_set")
  expect_equal(nrow(age), 1L)
  expect_identical(age$l3_cope_name, "age")
  expect_identical(age$l3_cope_number, 2L)
  expect_true(age$target_ready)
  expect_match(age$target_key, "l3_cope_name=age", fixed = TRUE)

  both_spec <- fwe_spec(
    "both_ptfce",
    targets = list(l3_model = "l3_group", l3_cope_name = c("group", "age"))
  )
  both <- resolve_fwe_targets(gpa, both_spec, source = "setup")
  expect_equal(nrow(both), 2L)
  expect_setequal(both$l3_cope_name, c("group", "age"))
  expect_equal(length(unique(both$target_id)), 2L)

  no_match <- fwe_spec("none", targets = list(l3_cope_name = "missing"))
  expect_error(
    resolve_fwe_targets(gpa, no_match, source = "setup"),
    "No level-3 outputs match"
  )
})

test_that("resolve_fwe_targets can preview missing and incomplete outputs", {
  missing_root <- tempfile("fwe_missing_")
  missing_gpa <- make_fwe_spec_test_gpa(missing_root, create_images = c(TRUE, FALSE))
  age_spec <- fwe_spec("age", targets = list(l3_cope_name = "age"))

  expect_error(
    resolve_fwe_targets(missing_gpa, age_spec, source = "setup"),
    "do not exist"
  )
  preview <- resolve_fwe_targets(
    missing_gpa, age_spec, source = "setup", require_existing = FALSE
  )
  expect_false(preview$target_ready)
  expect_false(preview$image_exists)

  incomplete_root <- tempfile("fwe_incomplete_")
  incomplete_gpa <- make_fwe_spec_test_gpa(incomplete_root, complete = FALSE)
  expect_error(
    resolve_fwe_targets(incomplete_gpa, age_spec, source = "setup"),
    "not complete"
  )
  incomplete_preview <- resolve_fwe_targets(
    incomplete_gpa, age_spec, source = "setup", require_complete = FALSE
  )
  expect_false(incomplete_preview$target_ready)
})
