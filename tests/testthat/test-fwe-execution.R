test_that("pTFCE execution dry runs expose commands without side effects", {
  root <- tempfile("fwe_execute_dry_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root)
  spec <- fwe_spec(
    "all_maps", targets = list(l3_model = "l3_group"),
    fwe_alpha = c(0.05, 0.01)
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  run <- run_fwe_plan(
    plan, scheduler = "local", dry_run = TRUE, worker_script = worker
  )

  expect_s3_class(run, "fwe_run")
  expect_true(run$dry_run)
  expect_true(all(run$execution$execution_status == "planned"))
  expect_true(all(grepl(
    rscript_executable(),
    run$execution$command,
    fixed = TRUE
  )))
  expect_true(all(grepl("--out_dir", run$execution$command, fixed = TRUE)))
  expect_true(all(grepl("--fwep", run$execution$command, fixed = TRUE)))
  expect_false(dir.exists(plan$output_directory))
})

test_that("pTFCE commands retain structured arguments for local execution", {
  root <- tempfile("fwe_execute_arguments_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root)
  spec <- fwe_spec(
    "structured", targets = list(l3_cope_name = "group")
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  commands <- build_fwe_commands(
    plan$spec$method, plan, plan$tasks, worker_script = worker
  )

  expect_identical(commands$executable, rscript_executable())
  expect_length(commands$arguments, 1L)
  expect_identical(
    commands$arguments[[1L]][5L],
    normalizePath(worker, mustWork = TRUE)
  )
  expect_true("--zstat" %in% commands$arguments[[1L]])
})

test_that("local pTFCE execution refreshes artifacts and task completion", {
  root <- tempfile("fwe_execute_local_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root)
  spec <- fwe_spec(
    "local_maps", targets = list(l3_model = "l3_group"),
    fwe_alpha = c(0.05, 0.01)
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  r_tests <- file.path(root, "child_r_tests.R")
  writeLines("stop('R_TESTS leaked into the pTFCE worker')", r_tests)
  withr::local_envvar(c(R_TESTS = r_tests))
  run <- run_fwe_plan(plan, scheduler = "local", worker_script = worker)

  expect_identical(Sys.getenv("R_TESTS"), r_tests)
  expect_true(all(run$execution$execution_status == "completed"))
  expect_true(all(run$execution$exit_status == 0L))
  expect_true(all(file.exists(run$execution$script_file)))
  expect_true(all(file.exists(run$execution$log_file)))
  expect_true(all(run$plan$tasks$status == "complete"))
  expect_true(all(run$plan$artifacts$exists))
  expect_s3_class(run$result, "fwe_result")
  expect_true(run$result$complete)
  expect_true(file.exists(run$result$manifest_file))
  expect_true(file.exists(run$result$result_file))

  restored <- read_fwe_result(plan$output_directory)
  manifest <- read_fwe_artifact_manifest(plan$output_directory)
  expect_s3_class(restored, "fwe_result")
  expect_true(restored$complete)
  expect_equal(nrow(manifest), 8L)
  expect_identical(anyDuplicated(manifest$artifact_id), 0L)
  expect_true(all(manifest$artifact_status == "available"))

  masks <- fwe_result_artifacts(
    restored,
    artifact_role = "thresholded_statistic",
    fwe_alpha = 0.05
  )
  expect_equal(nrow(masks), 2L)
  expect_true(all(masks$can_use_as_mask))
  expect_true(all(file.exists(masks$artifact_file)))

  skipped <- run_fwe_plan(run$plan, scheduler = "local")
  expect_true(all(skipped$execution$execution_status == "skipped_complete"))
  expect_length(skipped$job_ids, 0L)

  forced <- run_fwe_plan(
    run$plan, scheduler = "local", force = TRUE, dry_run = TRUE,
    worker_script = worker
  )
  expect_true(all(forced$execution$execution_status == "planned"))
})

test_that("blocked FWE tasks can be rejected or explicitly skipped", {
  root <- tempfile("fwe_execute_blocked_")
  gpa <- make_fwe_execution_test_gpa(root, create_smoothness = FALSE)
  spec <- fwe_spec("blocked", targets = list(l3_model = "l3_group"))
  plan <- plan_fwe_correction(
    gpa, spec, source = "setup", require_inputs = FALSE
  )

  expect_error(run_fwe_plan(plan, dry_run = TRUE), "blocked tasks")
  run <- run_fwe_plan(plan, dry_run = TRUE, allow_blocked = TRUE)
  expect_true(all(run$execution$execution_status == "skipped_blocked"))
  expect_true(all(is.na(run$execution$command)))
})

test_that("local execution failures return an inspectable run report", {
  root <- tempfile("fwe_execute_failure_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root, fail = TRUE)
  spec <- fwe_spec("failure", targets = list(l3_cope_name = "age"))
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  run <- run_fwe_plan(
    plan, scheduler = "local", worker_script = worker,
    stop_on_error = FALSE
  )
  expect_identical(run$execution$execution_status, "failed")
  expect_identical(run$execution$exit_status, 3L)

  condition <- tryCatch(
    run_fwe_plan(plan, scheduler = "local", worker_script = worker),
    fwe_execution_error = identity
  )
  expect_s3_class(condition, "fwe_execution_error")
  expect_s3_class(condition$run, "fwe_run")
})

test_that("refresh_fwe_plan detects artifacts created outside the runner", {
  root <- tempfile("fwe_execute_refresh_")
  gpa <- make_fwe_execution_test_gpa(root)
  spec <- fwe_spec(
    "refresh", targets = list(l3_cope_name = "group"),
    write_threshold_images = FALSE
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")
  dir.create(plan$tasks$output_directory, recursive = TRUE)
  write_fwe_test_result(plan)

  refreshed <- refresh_fwe_plan(plan)

  expect_identical(refreshed$tasks$status, "complete")
  expect_true(all(refreshed$artifacts$exists))
})

test_that("scheduler execution records submitted job identifiers", {
  root <- tempfile("fwe_execute_scheduler_")
  gpa <- make_fwe_execution_test_gpa(root)
  worker <- make_fake_ptfce_worker(root)
  spec <- fwe_spec("scheduler", targets = list(l3_cope_name = "group"),
    compute = list(cpus_per_job = 2L, memgb_per_command = 3, time_per_command = "00:20:00"))
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  testthat::local_mocked_bindings(
    cluster_submit_shell_jobs = function(job_list, log_file, ...) {
      settings <- list(...)
      expect_equal(settings$cpus_per_job, 2L)
      expect_equal(settings$memgb_per_command, 3)
      expect_identical(settings$time_per_command, "00:20:00")
      utils::write.csv(
        data.frame(
          job_number = seq_along(job_list),
          job_id = rep("12345", length(job_list)),
          cmd = job_list,
          stringsAsFactors = FALSE
        ),
        log_file,
        row.names = FALSE
      )
      "12345"
    },
    .package = "fmri.pipeline"
  )

  run <- run_fwe_plan(plan, scheduler = "slurm", worker_script = worker,
    tracking_sqlite_db = file.path(root, "tracking.sqlite"))

  expect_identical(run$execution$execution_status, "submitted")
  expect_identical(run$execution$job_id, "12345")
  expect_identical(run$job_ids, "12345")
  tracked <- fmri.pipeline:::submit_sqlite_query("SELECT job_id, job_name FROM job_tracking",
    file.path(root, "tracking.sqlite"), return_result = TRUE)
  expect_identical(tracked$job_id, "12345")
  expect_identical(tracked$job_name, "fwe_scheduler")
  expect_s3_class(run$result, "fwe_result")
  expect_false(run$result$complete)
  expect_true(file.exists(run$result$manifest_file))
  expect_true(all(run$result$manifest$execution_status == "submitted"))
  expect_true(all(run$result$manifest$job_id == "12345"))
})

test_that("local execution safely applies named environment variables", {
  expect_identical(
    fmri.pipeline:::normalize_fwe_run_environment(c(TEST_VALUE = "a b")),
    "TEST_VALUE='a b'"
  )
  expect_error(
    fmri.pipeline:::normalize_fwe_run_environment(c("BAD-NAME" = "x")),
    "invalid environment variable name"
  )

  variable <- "FMRI_PIPELINE_SCOPED_ENV_TEST"
  withr::local_envvar(stats::setNames("before", variable))
  observed <- fmri.pipeline:::with_fwe_run_environment(
    stats::setNames("during", variable),
    Sys.getenv(variable)
  )
  expect_identical(observed, "during")
  expect_identical(Sys.getenv(variable), "before")
})

test_that("result collection updates a partial persistent manifest", {
  root <- tempfile("fwe_collect_partial_")
  gpa <- make_fwe_execution_test_gpa(root)
  spec <- fwe_spec(
    "partial", targets = list(l3_cope_name = "age"),
    write_threshold_images = FALSE
  )
  plan <- plan_fwe_correction(gpa, spec, source = "setup")

  partial <- collect_fwe_results(plan)
  expect_s3_class(partial, "fwe_result")
  expect_false(partial$complete)
  expect_true(all(partial$manifest$artifact_status == "missing_required"))
  expect_error(
    collect_fwe_results(partial, require_complete = TRUE),
    "Required FWE artifacts are missing or unverified"
  )

  dir.create(plan$tasks$output_directory, recursive = TRUE)
  write_fwe_test_result(plan)
  complete <- collect_fwe_results(partial, require_complete = TRUE)

  expect_true(complete$complete)
  expect_true(all(complete$manifest$exists))
  expect_true(all(complete$manifest$artifact_status == "available"))
  expect_equal(
    read_fwe_artifact_manifest(plan$output_directory),
    complete$manifest,
    ignore_attr = TRUE
  )
})
