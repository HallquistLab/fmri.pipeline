test_that("fsl_l3_model rewrites existing FSF when model is queued to run", {
  outdir <- tempfile("l3_fsf_refresh_")
  dir.create(outdir, recursive = TRUE)

  l3_df <- data.frame(
    id = "sub1",
    session = 1L,
    l1_model = "l1m",
    l2_model = "l2m",
    l3_model = "l3m",
    l1_cope_name = "EV_choice",
    l2_cope_name = "EV_Intercept",
    cope_file = "dummy.nii.gz",
    stringsAsFactors = FALSE
  )

  gpa <- list(
    lgr_threshold = "warn",
    multi_run = TRUE,
    l3_models = list(models = list(l3m = list())),
    output_locations = list(feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf"),
    additional = list(feat_l3_args = NULL),
    glm_settings = list(fsl = list(force_l3_creation = FALSE))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  l3_fsf <- file.path(outdir, "FEAT_l1c-EV_choice.fsf")
  writeLines(c("# SENTINEL_OLD_FILE", "set fmri(robust_yn) 1"), con = l3_fsf)

  mock_feat_status <- function(feat_dir, feat_fsf, feat_complete) {
    data.frame(
      feat_fsf = feat_fsf,
      feat_fsf_modified_date = Sys.time(),
      feat_fsf_exists = file.exists(feat_fsf),
      feat_dir = feat_dir,
      feat_dir_exists = dir.exists(feat_dir),
      feat_execution_start = as.POSIXct(NA),
      feat_execution_end = as.POSIXct(NA),
      feat_execution_min = NA_real_,
      feat_complete = feat_complete,
      feat_failed = !feat_complete,
      stringsAsFactors = FALSE
    )
  }

  testthat::local_mocked_bindings(
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = new_data[, c("id", "session"), drop = FALSE],
        model_matrix = matrix(1, nrow = nrow(new_data), ncol = 1),
        contrasts = matrix(1, nrow = 1, ncol = 1),
        fsl_outlier_deweighting = FALSE
      )
    },
    fsl_generate_fsf_ev_syntax = function(inputs, dmat) character(0),
    fsl_generate_fsf_contrast_syntax = function(contrasts) character(0),
    get_output_directory = function(...) outdir,
    get_feat_status = function(feat_dir, feat_fsf, lg = NULL) {
      mock_feat_status(feat_dir, feat_fsf, feat_complete = FALSE)
    },
    .package = "fmri.pipeline"
  )

  result <- fmri.pipeline:::fsl_l3_model(l3_df = l3_df, gpa = gpa)
  expect_true(isTRUE(result$to_run[1L]))

  rewritten <- readLines(l3_fsf)
  expect_true(any(grepl("^set fmri\\(robust_yn\\) 0$", rewritten)))
  expect_false(any(grepl("SENTINEL_OLD_FILE", rewritten, fixed = TRUE)))
})


test_that("fsl_l3_model refreshes changed specifications even when cached status says complete", {
  outdir <- tempfile("l3_fsf_skip_")
  dir.create(outdir, recursive = TRUE)

  l3_df <- data.frame(
    id = "sub1",
    session = 1L,
    l1_model = "l1m",
    l2_model = "l2m",
    l3_model = "l3m",
    l1_cope_name = "EV_choice",
    l2_cope_name = "EV_Intercept",
    cope_file = "dummy.nii.gz",
    stringsAsFactors = FALSE
  )

  gpa <- list(
    lgr_threshold = "warn",
    multi_run = TRUE,
    l3_models = list(models = list(l3m = list())),
    output_locations = list(feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf"),
    additional = list(feat_l3_args = NULL),
    glm_settings = list(fsl = list(force_l3_creation = FALSE))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  l3_fsf <- file.path(outdir, "FEAT_l1c-EV_choice.fsf")
  writeLines(c("# SENTINEL_OLD_FILE", "set fmri(robust_yn) 1"), con = l3_fsf)

  mock_feat_status <- function(feat_dir, feat_fsf, feat_complete) {
    data.frame(
      feat_fsf = feat_fsf,
      feat_fsf_modified_date = Sys.time(),
      feat_fsf_exists = file.exists(feat_fsf),
      feat_dir = feat_dir,
      feat_dir_exists = dir.exists(feat_dir),
      feat_execution_start = as.POSIXct(NA),
      feat_execution_end = as.POSIXct(NA),
      feat_execution_min = NA_real_,
      feat_complete = feat_complete,
      feat_failed = !feat_complete,
      stringsAsFactors = FALSE
    )
  }

  testthat::local_mocked_bindings(
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = new_data[, c("id", "session"), drop = FALSE],
        model_matrix = matrix(1, nrow = nrow(new_data), ncol = 1),
        contrasts = matrix(1, nrow = 1, ncol = 1),
        fsl_outlier_deweighting = FALSE
      )
    },
    fsl_generate_fsf_ev_syntax = function(inputs, dmat) character(0),
    fsl_generate_fsf_contrast_syntax = function(contrasts) character(0),
    get_output_directory = function(...) outdir,
    get_feat_status = function(feat_dir, feat_fsf, lg = NULL) {
      mock_feat_status(feat_dir, feat_fsf, feat_complete = TRUE)
    },
    .package = "fmri.pipeline"
  )

  result <- fmri.pipeline:::fsl_l3_model(l3_df = l3_df, gpa = gpa)
  expect_false(isTRUE(result$to_run[1L]))

  existing <- readLines(l3_fsf)
  expect_false(any(grepl("SENTINEL_OLD_FILE", existing, fixed = TRUE)))
  expect_true(any(grepl("^set fmri\\(robust_yn\\) 0$", existing)))
})

test_that("fsl_l3_model allows pooled multi-session input for pooled_sessions_subject_ev", {
  outdir <- tempfile("l3_pooled_ok_")
  dir.create(outdir, recursive = TRUE)

  l3_df <- data.frame(
    id = c("sub1", "sub1", "sub2", "sub2"),
    session = c(1L, 2L, 1L, 2L),
    l1_model = "l1m",
    l2_model = "l2m",
    l3_model = "l3m",
    l3_input_mode = "pooled_sessions_subject_ev",
    l1_cope_name = "EV_choice",
    l2_cope_name = "EV_Intercept",
    cope_file = paste0("dummy_", seq_len(4), ".nii.gz"),
    stringsAsFactors = FALSE
  )

  gpa <- list(
    lgr_threshold = "warn",
    multi_run = TRUE,
    l3_models = list(models = list(l3m = list(l3_input_mode = "pooled_sessions_subject_ev"))),
    output_locations = list(feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf"),
    additional = list(feat_l3_args = NULL),
    glm_settings = list(fsl = list(force_l3_creation = TRUE)),
    run_data = data.frame(
      id = c("sub1", "sub1", "sub2", "sub2"),
      session = c(1L, 2L, 1L, 2L),
      run_number = c(1L, 1L, 1L, 1L),
      stringsAsFactors = FALSE
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  testthat::local_mocked_bindings(
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = unique(new_data[, c("id", "session"), drop = FALSE]),
        model_matrix = matrix(1, nrow = nrow(unique(new_data[, c("id", "session"), drop = FALSE])), ncol = 1),
        contrasts = matrix(1, nrow = 1, ncol = 1),
        fsl_outlier_deweighting = FALSE
      )
    },
    fsl_generate_fsf_ev_syntax = function(inputs, dmat) character(0),
    fsl_generate_fsf_contrast_syntax = function(contrasts) character(0),
    get_output_directory = function(...) outdir,
    get_feat_status = function(feat_dir, feat_fsf, lg = NULL) {
      data.frame(
        feat_fsf = feat_fsf,
        feat_fsf_modified_date = Sys.time(),
        feat_fsf_exists = file.exists(feat_fsf),
        feat_dir = feat_dir,
        feat_dir_exists = dir.exists(feat_dir),
        feat_execution_start = as.POSIXct(NA),
        feat_execution_end = as.POSIXct(NA),
        feat_execution_min = NA_real_,
        feat_complete = FALSE,
        feat_failed = FALSE,
        stringsAsFactors = FALSE
      )
    },
    .package = "fmri.pipeline"
  )

  expect_no_error(
    result <- fmri.pipeline:::fsl_l3_model(l3_df = l3_df, gpa = gpa)
  )
  expect_equal(result$session[1L], 0L)
  expect_equal(result$l3_input_mode[1L], "pooled_sessions_subject_ev")
})

test_that("fsl_l3_model adds session suffix for per_session in multi-session projects", {
  outdir <- tempfile("l3_session_suffix_")
  dir.create(outdir, recursive = TRUE)

  l3_df <- data.frame(
    id = c("sub1", "sub2"),
    session = c(2L, 2L),
    l1_model = "l1m",
    l2_model = "l2m",
    l3_model = "l3m",
    l3_input_mode = "per_session",
    l1_cope_name = "EV_choice",
    l2_cope_name = "EV_Intercept",
    cope_file = c("dummy1.nii.gz", "dummy2.nii.gz"),
    stringsAsFactors = FALSE
  )

  gpa <- list(
    lgr_threshold = "warn",
    multi_run = TRUE,
    l3_models = list(models = list(l3m = list(l3_input_mode = "per_session"))),
    output_locations = list(feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf"),
    additional = list(feat_l3_args = NULL),
    glm_settings = list(fsl = list(force_l3_creation = TRUE)),
    run_data = data.frame(
      id = c("sub1", "sub1", "sub2", "sub2"),
      session = c(1L, 2L, 1L, 2L),
      run_number = c(1L, 1L, 1L, 1L),
      stringsAsFactors = FALSE
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  testthat::local_mocked_bindings(
    respecify_l3_model = function(mobj, new_data) {
      list(
        metadata = unique(new_data[, c("id", "session"), drop = FALSE]),
        model_matrix = matrix(1, nrow = nrow(unique(new_data[, c("id", "session"), drop = FALSE])), ncol = 1),
        contrasts = matrix(1, nrow = 1, ncol = 1),
        fsl_outlier_deweighting = FALSE
      )
    },
    fsl_generate_fsf_ev_syntax = function(inputs, dmat) character(0),
    fsl_generate_fsf_contrast_syntax = function(contrasts) character(0),
    get_output_directory = function(...) outdir,
    get_feat_status = function(feat_dir, feat_fsf, lg = NULL) {
      data.frame(
        feat_fsf = feat_fsf,
        feat_fsf_modified_date = Sys.time(),
        feat_fsf_exists = file.exists(feat_fsf),
        feat_dir = feat_dir,
        feat_dir_exists = dir.exists(feat_dir),
        feat_execution_start = as.POSIXct(NA),
        feat_execution_end = as.POSIXct(NA),
        feat_execution_min = NA_real_,
        feat_complete = FALSE,
        feat_failed = FALSE,
        stringsAsFactors = FALSE
      )
    },
    .package = "fmri.pipeline"
  )

  result <- fmri.pipeline:::fsl_l3_model(l3_df = l3_df, gpa = gpa)
  expect_true(grepl("_ses-2\\.fsf$", result$feat_fsf[1L]))
})
