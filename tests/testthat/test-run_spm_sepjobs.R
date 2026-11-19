test_that("run_spm_sepjobs final tracking reuses one terminal job status", {
  spm_path <- source_tree_file("R", "spm_utils.R")
  lines <- readLines(spm_path, warn = FALSE)

  expect_true(any(grepl("if [ $job_failed -ne 0 ]; then", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"FAILED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("job_status=\\\"COMPLETED\\\"", lines, fixed = TRUE)))
  expect_true(any(grepl("--status", lines) & grepl("\\$\\{job_status\\}", lines)))
  expect_true(any(grepl("if [ \\\"${job_status}\\\" = \\\"FAILED\\\" ]; then", lines, fixed = TRUE)))
})


test_that("SPM runner stops at failures and validates before marking completion", {
  td <- tempfile("spm_runner_")
  dir.create(td)
  runner <- system.file("scripts", "spm_runner.sh", package = "fmri.pipeline")
  fake <- file.path(td, "matlab")
  # Execute script bodies in bash to exercise the generated stage sequencing.
  writeLines(c("#!/bin/bash", "expression=\"$1\"",
    "path=$(printf '%s' \"$expression\" | sed -n \"s/.*run('\\([^']*\\)').*/\\1/p\")",
    "bash \"$path\""), fake)
  Sys.chmod(fake, "0755")
  for (name in c("setup.m", "fit.m", "contrasts.m")) {
    writeLines(paste("echo", name, ">>", shQuote(file.path(td, "stages"))), file.path(td, name))
  }
  writeLines("true", file.path(td, "setup_spm_contrasts.sh"))
  writeLines("token", file.path(td, "spm_setup_token.txt"))
  writeLines("output", file.path(td, "output.nii"))
  writeLines(file.path(td, "output.nii"), file.path(td, "spm_expected_outputs.txt"))
  writeLines(paste("cp", shQuote(file.path(td, "spm_setup_token.txt")),
    shQuote(file.path(td, "spm_validated_token.txt"))), file.path(td, "validate_spm_outputs.m"))
  job <- file.path(td, "job.sh")
  write_job <- function(partial = FALSE) writeLines(c("#!/bin/bash",
    paste0("matlab_cmd=", shQuote(fake)), "matlab_args=", "run_setup=1", "run_glm=1",
    paste0("run_contrasts=", as.integer(!partial)), "spm_level=1", "cleanup_tmp=0",
    "spm_setup_script=setup.m", "spm_glm_script=fit.m", "spm_contrast_script=contrasts.m",
    "write_job_manifest() { :; }", paste("source", shQuote(runner)),
    paste("spm_runner", shQuote(td))), job)
  run <- function() suppressWarnings(system2("bash", shQuote(job), stdout = TRUE, stderr = TRUE))
  write_job()
  expect_null(attr(run(), "status"))
  expect_true(file.exists(file.path(td, ".spm_complete")))
  unlink(file.path(td, "stages"))
  writeLines("exit 7", file.path(td, "fit.m"))
  expect_equal(attr(run(), "status"), 7L)
  expect_equal(readLines(file.path(td, "stages")), "setup.m")
  expect_true(file.exists(file.path(td, ".spm_fail")))
  expect_false(file.exists(file.path(td, ".spm_complete")))
  unlink(file.path(td, "setup.m"))
  expect_equal(attr(run(), "status"), 1L)
  writeLines("true", file.path(td, "setup.m"))
  writeLines("true", file.path(td, "fit.m"))
  unlink(file.path(td, "output.nii"))
  expect_equal(attr(run(), "status"), 1L)
  expect_false(file.exists(file.path(td, ".spm_complete")))
  write_job(partial = TRUE)
  expect_null(attr(run(), "status"))
  expect_true(file.exists(file.path(td, ".spm_stage_complete")))
  expect_false(file.exists(file.path(td, ".spm_complete")))
})


test_that("SPM submission refreshes stale queue status and refuses stale setup", {
  td <- tempfile("spm_submit_")
  out <- file.path(td, "model")
  make_spm_completed_fixture(out)
  unlink(file.path(out, "con_0001.nii"))
  gpa <- create_mock_gpa(output_directory = td)
  gpa$parallel$spm <- list(l1_spm_time = "00:10:00", l1_spm_memgb = 1,
    l1_spm_cpus_per_job = 1, l1_spm_runs_per_cpu = 1)
  gpa$parallel$compute_environment <- list(global = character(), spm = "true", r = character())
  gpa$glm_settings <- list(spm = list(run_l1_setup = TRUE, run_l1_glm = TRUE, run_l1_contrasts = TRUE))
  gpa$l1_model_setup <- structure(list(spm = data.frame(spm_dir = out, spm_complete = TRUE, to_run = FALSE)), class = "l1_setup")
  submitted <- character()
  local_mocked_bindings(cluster_job_submit = function(job_script, ...) {
    submitted <<- c(submitted, job_script)
    "test-job"
  }, .package = "fmri.pipeline")
  expect_equal(run_spm_sepjobs(gpa), "test-job")
  expect_length(submitted, 1L)
  expect_equal(system2("bash", c("-n", shQuote(submitted))), 0L)
  fmri.pipeline:::spm_invalidate_setup(out)
  expect_error(run_spm_sepjobs(gpa), "setup is missing or stale")
})
