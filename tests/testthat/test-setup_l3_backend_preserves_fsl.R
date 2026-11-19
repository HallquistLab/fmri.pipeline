# backend-scoped l3 setup should not overwrite other backends

test_that("setup_l3_models(backend='spm') preserves existing FSL l3 tables", {
  spm_dir <- file.path(tempdir(), "spm_l1_test")
  dir.create(spm_dir, recursive = TRUE, showWarnings = FALSE)
  make_spm_completed_fixture(spm_dir)
  subj_ids <- paste0("sub", 1:4)

  l1_models <- list(models = list(model1 = list()))
  class(l1_models) <- c("l1_model_set", "list")

  l3_models <- list(
    models = list(
      l3_model1 = list(
        contrasts = matrix(1, nrow = 1, ncol = 1, dimnames = list("cope1", "intercept"))
      )
    )
  )
  class(l3_models) <- c("hi_model_set", "list")

  l1_setup <- list(
    spm = data.frame(
      id = subj_ids, session = 1L, l1_model = "model1",
      spm_dir = rep(spm_dir, length(subj_ids)),
      spm_complete = TRUE, spm_contrast_exists = TRUE,
      stringsAsFactors = FALSE
    )
  )
  class(l1_setup) <- c("l1_setup", "list")

  fsl_l3 <- data.frame(
    l1_model = "model1", l1_cope_name = "cope1", l3_model = "l3_model1",
    feat_fsf = "old.fsf",
    stringsAsFactors = FALSE
  )
  l3_setup <- list(fsl = fsl_l3)
  class(l3_setup) <- c("l3_setup", "list")

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
    lgr_threshold = "info",
    log_txt = FALSE,
    log_json = FALSE,
    output_locations = list(
      setup_l3_log_txt = file.path(tempdir(), "setup_l3_models.txt"),
      setup_l3_log_json = file.path(tempdir(), "setup_l3_models.json")
    ),
    subject_data = data.frame(id = subj_ids, session = 1L, exclude_subject = FALSE),
    run_data = data.frame(id = subj_ids, session = 1L, run_number = 1L, exclude_run = FALSE, exclude_subject = FALSE),
    l1_models = l1_models,
    l1_cope_names = list(model1 = c("cope1")),
    l3_models = l3_models,
    l1_model_setup = l1_setup,
    l3_model_setup = l3_setup,
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

  expect_true(inherits(res$l3_model_setup, "l3_setup"))
  expect_equal(res$l3_model_setup$fsl, fsl_l3)
  expect_true("spm" %in% names(res$l3_model_setup))
})
