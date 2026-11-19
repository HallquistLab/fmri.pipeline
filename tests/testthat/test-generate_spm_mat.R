test_that("generate_spm_mat handles concatenated ts_files and non-concat mismatch", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_mat_")
  dir.create(tmp_dir, recursive = TRUE)

  make_nifti <- function(path, nvol = 5L) {
    dims <- c(2L, 2L, 2L, nvol)
    img <- RNifti::asNifti(array(0, dim = dims))
    RNifti::writeNifti(img, path)
  }

  nifti1 <- file.path(tmp_dir, "run1.nii")
  nifti2 <- file.path(tmp_dir, "run2.nii")
  make_nifti(nifti1, 5L)
  make_nifti(nifti2, 6L)

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg1 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "feedback")
  design <- array(list(reg1, reg1), dim = c(2L, 1L), dimnames = list(c("run1", "run2"), "reg1"))
  design_concat <- list(reg1)
  names(design_concat) <- "reg1"

  bdm <- list(
    design = design,
    design_concat = design_concat,
    run_niftis = c(nifti1, nifti2),
    tr = 1.0
  )
  class(bdm) <- c("bdm", "list")

  ts_file <- file.path(tmp_dir, "confounds.txt")
  writeLines(rep("1 2 3", 11), ts_file)

  expect_error(
    generate_spm_mat(
      bdm = bdm,
      ts_files = ts_file,
      output_dir = file.path(tmp_dir, "non_concat"),
      concatenate_runs = FALSE,
      spm_path = tmp_dir
    ),
    "Length of ts_files argument is not equal to length of run_niftis"
  )

  expect_no_error(
    generate_spm_mat(
      bdm = bdm,
      ts_files = ts_file,
      output_dir = file.path(tmp_dir, "concat"),
      concatenate_runs = TRUE,
      spm_path = tmp_dir
    )
  )
})

test_that("generate_spm_mat writes projected L2 regressors for concatenated runs", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_mat_proj_")
  dir.create(tmp_dir, recursive = TRUE)

  make_nifti <- function(path, nvol = 5L) {
    dims <- c(2L, 2L, 2L, nvol)
    img <- RNifti::asNifti(array(0, dim = dims))
    RNifti::writeNifti(img, path)
  }

  nifti1 <- file.path(tmp_dir, "run1.nii")
  nifti2 <- file.path(tmp_dir, "run2.nii")
  make_nifti(nifti1, 5L)
  make_nifti(nifti2, 6L)

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg1 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "feedback")
  design <- array(list(reg1, reg1), dim = c(2L, 1L), dimnames = list(c("run1", "run2"), "reg1"))
  design_concat <- list(reg1)
  names(design_concat) <- "reg1"

  bdm <- list(
    design = design,
    design_concat = design_concat,
    run_niftis = c(nifti1, nifti2),
    runs_to_output = c(1L, 2L),
    tr = 1.0,
    spm_l2_projection = data.frame(
      run_number = c(1L, 2L),
      drug = c(1, 0),
      stringsAsFactors = FALSE
    )
  )
  class(bdm) <- c("bdm", "list")

  out_dir <- file.path(tmp_dir, "concat_proj")
  expect_no_error(
    generate_spm_mat(
      bdm = bdm,
      output_dir = out_dir,
      concatenate_runs = TRUE,
      spm_path = tmp_dir
    )
  )

  mfile <- file.path(out_dir, "glm_design_batch.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_true(any(grepl("sess\\(1\\)\\.regress\\(1\\)\\.name = 'drug'", lines)))
  expect_true(any(grepl("1 1 1 1 1 0 0 0 0 0 0", lines, fixed = TRUE)))
})

test_that("generate_spm_mat writes selected L2-by-L1 interaction modulators", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_mat_projint_")
  dir.create(tmp_dir, recursive = TRUE)

  make_nifti <- function(path, nvol = 5L) {
    dims <- c(2L, 2L, 2L, nvol)
    img <- RNifti::asNifti(array(0, dim = dims))
    RNifti::writeNifti(img, path)
  }

  nifti1 <- file.path(tmp_dir, "run1.nii")
  nifti2 <- file.path(tmp_dir, "run2.nii")
  make_nifti(nifti1, 5L)
  make_nifti(nifti2, 5L)

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg_run1 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "face_onset")
  reg_run2 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "face_onset")
  reg_concat <- mk_reg(c(1, 3, 6, 8), c(1, 1, 1, 1), c(1, 1, 1, 1), "face_onset")

  design <- array(
    list(reg_run1, reg_run2),
    dim = c(2L, 1L),
    dimnames = list(c("run1", "run2"), "face_reg")
  )
  design_concat <- list(reg_concat)
  names(design_concat) <- "face_reg"

  bdm <- list(
    design = design,
    design_concat = design_concat,
    run_niftis = c(nifti1, nifti2),
    runs_to_output = c(1L, 2L),
    run_volumes = c(5L, 5L),
    tr = 1.0,
    spm_l2_projection = data.frame(
      run_number = c(1L, 2L),
      face_emotion = c(1, -1),
      stringsAsFactors = FALSE
    ),
    spm_l2_projection_interactions = "face_emotion"
  )
  class(bdm) <- c("bdm", "list")

  out_dir <- file.path(tmp_dir, "concat_proj_interactions")
  syntax <- generate_spm_mat(
    bdm = bdm,
    output_dir = out_dir,
    concatenate_runs = TRUE,
    spm_path = tmp_dir
  )

  expect_true("face_onset_x_l2_face_emotion" %in% syntax$projection_interaction_terms)
  mfile <- file.path(out_dir, "glm_design_batch.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_true(any(grepl("face_onset_x_l2_face_emotion-pmod", lines)))
  expect_true(any(grepl("1, 1, -1, -1", lines, fixed = TRUE)))
})

test_that("generate_spm_mat drops projected L2 terms without within-run variation in non-concatenated mode", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_mat_projint_nonconcat_")
  dir.create(tmp_dir, recursive = TRUE)

  make_nifti <- function(path, nvol = 5L) {
    dims <- c(2L, 2L, 2L, nvol)
    img <- RNifti::asNifti(array(0, dim = dims))
    RNifti::writeNifti(img, path)
  }

  nifti1 <- file.path(tmp_dir, "run1.nii")
  nifti2 <- file.path(tmp_dir, "run2.nii")
  make_nifti(nifti1, 5L)
  make_nifti(nifti2, 5L)

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg_run1 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "face_onset")
  reg_run2 <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "face_onset")

  bdm <- list(
    design = array(
      list(reg_run1, reg_run2),
      dim = c(2L, 1L),
      dimnames = list(c("run1", "run2"), "face_reg")
    ),
    run_niftis = c(nifti1, nifti2),
    runs_to_output = c(1L, 2L),
    run_volumes = c(5L, 5L),
    tr = 1.0,
    spm_l2_projection = data.frame(
      run_number = c(1L, 2L),
      face_emotion = c(1, -1),
      stringsAsFactors = FALSE
    ),
    spm_l2_projection_interactions = "face_emotion"
  )
  class(bdm) <- c("bdm", "list")

  out_dir <- file.path(tmp_dir, "nonconcat_proj_interactions")
  syntax <- NULL
  expect_warning(
    syntax <- generate_spm_mat(
      bdm = bdm,
      output_dir = out_dir,
      concatenate_runs = FALSE,
      spm_path = tmp_dir
    ),
    "Skipping projected L2 main-effect regressors in non-concatenated SPM design"
  )

  expect_equal(length(syntax$projection_interaction_terms), 0L)
  mfile <- file.path(out_dir, "glm_design_batch.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_false(any(grepl("face_onset_x_l2_face_emotion-pmod", lines, fixed = TRUE)))
  expect_false(any(grepl("\\.regress\\(1\\)\\.name = 'face_emotion'", lines)))
})

test_that("generate_spm_mat handles multiple unit-height regressors and empty levels", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_mat_multi_")
  dir.create(tmp_dir, recursive = TRUE)

  dims <- c(2L, 2L, 2L, 5L)
  img <- RNifti::asNifti(array(0, dim = dims))
  nifti <- file.path(tmp_dir, "run1.nii")
  RNifti::writeNifti(img, nifti)

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg_a <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "feedback")
  reg_b <- mk_reg(c(2, 4), c(1, 1), c(1, 1), "feedback")
  design_multi <- array(
    list(reg_a, reg_b),
    dim = c(1L, 2L),
    dimnames = list("run1", c("feedback.A", "feedback.B"))
  )

  bdm_multi <- list(
    design = design_multi,
    design_concat = NULL,
    run_niftis = c(nifti),
    tr = 1.0
  )
  class(bdm_multi) <- c("bdm", "list")

  out_dir_multi <- file.path(tmp_dir, "multi_unit")
  expect_no_error(
    generate_spm_mat(
      bdm = bdm_multi,
      ts_files = NULL,
      output_dir = out_dir_multi,
      concatenate_runs = FALSE,
      spm_path = tmp_dir
    )
  )
  mfile <- file.path(out_dir_multi, "glm_design_batch.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_true(any(grepl("feedback.A", lines)))
  expect_true(any(grepl("feedback.B", lines)))

  empty_reg <- matrix(numeric(0), nrow = 0L, ncol = 4L)
  colnames(empty_reg) <- c("trial", "onset", "duration", "value")
  attr(empty_reg, "event") <- "feedback"

  design_empty <- array(
    list(reg_a, empty_reg),
    dim = c(1L, 2L),
    dimnames = list("run1", c("feedback.A", "feedback.EMPTY"))
  )
  bdm_empty <- list(
    design = design_empty,
    design_concat = NULL,
    run_niftis = c(nifti),
    tr = 1.0
  )
  class(bdm_empty) <- c("bdm", "list")

  out_dir_empty <- file.path(tmp_dir, "empty_level")
  expect_warning(
    generate_spm_mat(
      bdm = bdm_empty,
      ts_files = NULL,
      output_dir = out_dir_empty,
      concatenate_runs = FALSE,
      spm_path = tmp_dir
    ),
    "Dropping empty regressors"
  )
  mfile_empty <- file.path(out_dir_empty, "glm_design_batch.m")
  expect_true(file.exists(mfile_empty))
  lines_empty <- readLines(mfile_empty)
  expect_true(any(grepl("feedback", lines_empty)))
  expect_false(any(grepl("feedback.EMPTY", lines_empty)))
})

test_that("generate_spm_mat expands wi regressors with all levels present", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_mat_wi_")
  dir.create(tmp_dir, recursive = TRUE)

  dims <- c(2L, 2L, 2L, 5L)
  img <- RNifti::asNifti(array(0, dim = dims))
  nifti <- file.path(tmp_dir, "run1.nii")
  RNifti::writeNifti(img, nifti)

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  # Two within-subject levels, both present with unit height
  reg_a <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "feedback")
  reg_b <- mk_reg(c(2, 4), c(1, 1), c(1, 1), "feedback")
  design <- array(
    list(reg_a, reg_b),
    dim = c(1L, 2L),
    dimnames = list("run1", c("feedback.A", "feedback.B"))
  )

  bdm <- list(
    design = design,
    design_concat = NULL,
    run_niftis = c(nifti),
    tr = 1.0
  )
  class(bdm) <- c("bdm", "list")

  out_dir <- file.path(tmp_dir, "wi_full")
  expect_no_warning(
    generate_spm_mat(
      bdm = bdm,
      ts_files = NULL,
      output_dir = out_dir,
      concatenate_runs = FALSE,
      spm_path = tmp_dir
    )
  )

  mfile <- file.path(out_dir, "glm_design_batch.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_true(any(grepl("feedback.A", lines)))
  expect_true(any(grepl("feedback.B", lines)))
  expect_false(any(grepl("\\.pmod\\(", lines)))
})

test_that("generate_spm_mat writes estimation method and residual flags", {
  skip_if_not_installed("RNifti")

  tmp_dir <- tempfile("spm_mat_est_")
  dir.create(tmp_dir, recursive = TRUE)

  dims <- c(2L, 2L, 2L, 5L)
  img <- RNifti::asNifti(array(0, dim = dims))
  nifti <- file.path(tmp_dir, "run1.nii")
  RNifti::writeNifti(img, nifti)

  mk_reg <- function(onsets, durations, values, event_name) {
    m <- cbind(
      trial = seq_along(onsets),
      onset = onsets,
      duration = durations,
      value = values
    )
    colnames(m) <- c("trial", "onset", "duration", "value")
    attr(m, "event") <- event_name
    m
  }

  reg <- mk_reg(c(1, 3), c(1, 1), c(1, 1), "feedback")
  design <- array(list(reg), dim = c(1L, 1L), dimnames = list("run1", "reg1"))

  bdm <- list(
    design = design,
    design_concat = NULL,
    run_niftis = c(nifti),
    tr = 1.0
  )
  class(bdm) <- c("bdm", "list")

  out_dir_bayes <- file.path(tmp_dir, "bayes2")
  expect_warning(
    generate_spm_mat(
      bdm = bdm,
      output_dir = out_dir_bayes,
      concatenate_runs = FALSE,
      spm_path = tmp_dir,
      estimation_method = "Bayesian2",
      write_residuals = TRUE
    ),
    "write_residuals is only supported for Classical"
  )
  run_glm_bayes <- readLines(file.path(out_dir_bayes, "run_glm.m"))
  expect_true(any(grepl("method\\.Bayesian2 = 1", run_glm_bayes)))
  expect_true(any(grepl("write_residuals = 0", run_glm_bayes)))

  out_dir_classical <- file.path(tmp_dir, "classical")
  expect_no_warning(
    generate_spm_mat(
      bdm = bdm,
      output_dir = out_dir_classical,
      concatenate_runs = FALSE,
      spm_path = tmp_dir,
      estimation_method = "reml",
      write_residuals = TRUE
    )
  )
  run_glm_classical <- readLines(file.path(out_dir_classical, "run_glm.m"))
  expect_true(any(grepl("method\\.Classical = 1", run_glm_classical)))
  expect_true(any(grepl("write_residuals = 1", run_glm_classical)))
})
