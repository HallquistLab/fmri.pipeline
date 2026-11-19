test_that("run_feat_sepjobs returns an empty character vector when nothing needs execution", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel$fsl <- list(
    l1_feat_time = "0:10:00",
    l1_feat_memgb = 1,
    l1_feat_cpus_per_job = 1L,
    l1_feat_runs_per_cpu = 1L
  )
  gpa$glm_settings <- list(
    fsl = list(
      failed_l1_folder_action = "keep",
      incomplete_l1_folder_action = "keep"
    )
  )
  gpa$l1_model_setup <- structure(
    list(
      fsl = data.frame(
        l1_model = "facehouse",
        feat_fsf = "facehouse.fsf",
        feat_dir = "facehouse.feat",
        feat_complete = TRUE,
        feat_failed = FALSE,
        to_run = FALSE,
        stringsAsFactors = FALSE
      )
    ),
    class = "l1_setup"
  )

  fit <- make_valid_feat(withr::local_tempdir())
  gpa$l1_model_setup$fsl$feat_fsf <- fit$fsf
  gpa$l1_model_setup$fsl$feat_dir <- fit$out
  res <- fmri.pipeline::run_feat_sepjobs(gpa, level = 1L, model_names = "facehouse")
  expect_type(res, "character")
  expect_length(res, 0L)
})

test_that("run_feat_sepjobs diagnoses an invalid empty setup table", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$l1_model_setup <- structure(
    list(fsl = data.table::data.table()),
    class = c("l1_setup", "list")
  )

  expect_error(
    fmri.pipeline::run_feat_sepjobs(gpa, level = 1L, model_names = "facehouse"),
    "setup step did not produce a valid FSL job table"
  )
})

test_that("R_batch_job skips waiting when child_job_ids is empty", {
  batch_dir <- tempfile("batch_empty_children_")
  dir.create(batch_dir)
  sqlite_db <- file.path(batch_dir, "jobs.sqlite")

  job <- fmri.pipeline::R_batch_job$new(
    batch_directory = batch_dir,
    job_name = "empty_children",
    r_code = "child_job_ids <- character(0)",
    wait_for_children = TRUE,
    scheduler = "local",
    sqlite_db = sqlite_db,
    print_session_info = FALSE
  )
  job$generate()

  compute_file <- list.files(batch_dir, pattern = "^batch_run.*\\.R$", full.names = TRUE)
  expect_length(compute_file, 1L)

  lines <- readLines(compute_file, warn = FALSE)
  expect_true(any(grepl("length\\(child_job_ids\\) > 0L", lines)))
  expect_true(any(grepl("warning\\('Attempt to wait for child jobs failed due to non-existent or improper child_job_ids variable\\.'\\)", lines)))
})
