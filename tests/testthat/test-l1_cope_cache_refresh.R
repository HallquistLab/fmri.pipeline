library(testthat)

test_that("refresh_l1_cope_names populates cache from L1 contrasts", {
  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 1, include_l1_models = FALSE)

  l1_cmat <- matrix(
    c(1, 0),
    nrow = 2,
    ncol = 1,
    dimnames = list(c("cope_face", "cope_house"), "Intercept")
  )
  gpa$l1_models <- list(
    models = list(
      model1 = list(
        name = "model1",
        contrasts = l1_cmat
      )
    )
  )
  class(gpa$l1_models) <- c("l1_model_set", "list")
  gpa$l1_cope_names <- NULL

  out <- fmri.pipeline:::refresh_l1_cope_names(gpa)
  expect_equal(out$l1_cope_names$model1, c("cope_face", "cope_house"))
})

test_that("setup_l3_models refreshes missing L1 cope cache before mapping copes", {
  spm_dir <- file.path(tempdir(), paste0("spm_l1_cache_refresh_", as.integer(Sys.time())))
  dir.create(spm_dir, recursive = TRUE, showWarnings = FALSE)
  make_spm_completed_fixture(spm_dir)

  subj_ids <- paste0("sub", 1:4)
  l1_models <- list(
    models = list(
      model1 = list(
        name = "model1",
        contrasts = matrix(1, nrow = 1, ncol = 1, dimnames = list("cope1", "Intercept"))
      )
    )
  )
  class(l1_models) <- c("l1_model_set", "list")

  l3_models <- list(
    models = list(
      l3_model1 = list(
        contrasts = matrix(1, nrow = 1, ncol = 1, dimnames = list("cope1", "Intercept"))
      )
    )
  )
  class(l3_models) <- c("hi_model_set", "list")

  l1_setup <- list(
    spm = data.frame(
      id = subj_ids,
      session = 1L,
      l1_model = "model1",
      spm_dir = rep(spm_dir, length(subj_ids)),
      spm_complete = TRUE,
      spm_contrast_exists = TRUE,
      stringsAsFactors = FALSE
    )
  )
  class(l1_setup) <- c("l1_setup", "list")

  fake_spm_status <- function(spm_dir, lg = NULL) {
    data.frame(spm_dir = spm_dir, spm_complete = TRUE, stringsAsFactors = FALSE)
  }
  fake_spm_l3_setup <- function(l3_df, gpa) {
    data.frame(
      l1_model = l3_df$l1_model[1L],
      l1_cope_name = l3_df$l1_cope_name[1L],
      l3_model = l3_df$l3_model[1L],
      spm_dir = l3_df$spm_dir[1L],
      stringsAsFactors = FALSE
    )
  }

  gpa <- list(
    glm_software = "spm",
    multi_run = FALSE,
    lgr_threshold = "warn",
    log_txt = FALSE,
    log_json = FALSE,
    output_locations = list(
      setup_l3_log_txt = file.path(tempdir(), "setup_l3_models.txt"),
      setup_l3_log_json = file.path(tempdir(), "setup_l3_models.json")
    ),
    subject_data = data.frame(id = subj_ids, session = 1L, exclude_subject = FALSE),
    run_data = data.frame(id = subj_ids, session = 1L, run_number = 1L, exclude_run = FALSE, exclude_subject = FALSE),
    l1_models = l1_models,
    l1_cope_names = NULL,
    l3_models = l3_models,
    l1_model_setup = l1_setup,
    glm_backend_specs = list(
      spm = list(
        name = "spm",
        runs_l1 = TRUE, runs_l2 = FALSE, runs_l3 = TRUE,
        multi_run_strategy = "concat_in_l1",
        produced_artifacts = c("run_level_contrasts", "subject_session_contrasts", "group_level_stats"),
        l1_status = fake_spm_status,
        l3_status = fake_spm_status,
        l1_status_inputs = c("spm_dir"),
        l3_status_inputs = c("spm_dir"),
        l3_setup = fake_spm_l3_setup
      )
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  foreach::registerDoSEQ()
  res <- setup_l3_models(
    gpa,
    l3_model_names = "l3_model1",
    l1_model_names = "model1",
    backend = "spm"
  )

  expect_equal(res$l1_cope_names$model1, "cope1")
  expect_true(inherits(res$l3_model_setup, "l3_setup"))
  expect_true("spm" %in% names(res$l3_model_setup))
  expect_equal(nrow(res$l3_model_setup$spm), 1L)
})

test_that("setup_l3_models errors clearly when L1 contrasts are missing", {
  spm_dir <- file.path(tempdir(), paste0("spm_l1_no_contrasts_", as.integer(Sys.time())))
  dir.create(spm_dir, recursive = TRUE, showWarnings = FALSE)
  make_spm_completed_fixture(spm_dir)

  subj_ids <- paste0("sub", 1:4)
  l1_models <- list(
    models = list(
      model1 = list(
        name = "model1",
        contrasts = matrix(numeric(0), nrow = 0, ncol = 1, dimnames = list(NULL, "Intercept"))
      )
    )
  )
  class(l1_models) <- c("l1_model_set", "list")

  l3_models <- list(
    models = list(
      l3_model1 = list(
        contrasts = matrix(1, nrow = 1, ncol = 1, dimnames = list("cope1", "Intercept"))
      )
    )
  )
  class(l3_models) <- c("hi_model_set", "list")

  l1_setup <- list(
    spm = data.frame(
      id = subj_ids,
      session = 1L,
      l1_model = "model1",
      spm_dir = rep(spm_dir, length(subj_ids)),
      spm_complete = TRUE,
      spm_contrast_exists = TRUE,
      stringsAsFactors = FALSE
    )
  )
  class(l1_setup) <- c("l1_setup", "list")

  fake_spm_status <- function(spm_dir, lg = NULL) {
    data.frame(spm_dir = spm_dir, spm_complete = TRUE, stringsAsFactors = FALSE)
  }
  fake_spm_l3_setup <- function(l3_df, gpa) {
    data.frame(
      l1_model = l3_df$l1_model[1L],
      l1_cope_name = l3_df$l1_cope_name[1L],
      l3_model = l3_df$l3_model[1L],
      spm_dir = l3_df$spm_dir[1L],
      stringsAsFactors = FALSE
    )
  }

  gpa <- list(
    glm_software = "spm",
    multi_run = FALSE,
    lgr_threshold = "warn",
    log_txt = FALSE,
    log_json = FALSE,
    output_locations = list(
      setup_l3_log_txt = file.path(tempdir(), "setup_l3_models.txt"),
      setup_l3_log_json = file.path(tempdir(), "setup_l3_models.json")
    ),
    subject_data = data.frame(id = subj_ids, session = 1L, exclude_subject = FALSE),
    run_data = data.frame(id = subj_ids, session = 1L, run_number = 1L, exclude_run = FALSE, exclude_subject = FALSE),
    l1_models = l1_models,
    l1_cope_names = NULL,
    l3_models = l3_models,
    l1_model_setup = l1_setup,
    glm_backend_specs = list(
      spm = list(
        name = "spm",
        runs_l1 = TRUE, runs_l2 = FALSE, runs_l3 = TRUE,
        multi_run_strategy = "concat_in_l1",
        produced_artifacts = c("run_level_contrasts", "subject_session_contrasts", "group_level_stats"),
        l1_status = fake_spm_status,
        l3_status = fake_spm_status,
        l1_status_inputs = c("spm_dir"),
        l3_status_inputs = c("spm_dir"),
        l3_setup = fake_spm_l3_setup
      )
    )
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  foreach::registerDoSEQ()
  expect_error(
    setup_l3_models(
      gpa,
      l3_model_names = "l3_model1",
      l1_model_names = "model1",
      backend = "spm"
    ),
    "No L1 contrasts available for model 'model1'"
  )
})
