test_that("cancelled Slurm states are terminal failures", {
  calls <- 0L
  current_state <- "CANCELLED"
  local_mocked_bindings(slurm_job_status = function(...) {
    calls <<- calls + 1L
    if (calls > 1L) stop("cancelled job was polled again")
    data.frame(State = current_state)
  }, .package = "fmri.pipeline")
  for (state in c("CANCELLED", "CANCELLED by 12345", "CANCELLED+", "FAILED", "TIMEOUT")) {
    calls <- 0L
    current_state <- state
    expect_output(expect_false(wait_for_job("123", scheduler = "slurm")), "failed")
    expect_equal(calls, 1L)
  }
  calls <- 0L
  current_state <- "COMPLETED"
  expect_true(wait_for_job("123", scheduler = "slurm"))
})

test_that("waiting continues for running siblings after a cancellation", {
  calls <- 0L
  local_mocked_bindings(slurm_job_status = function(...) {
    calls <<- calls + 1L
    if (calls > 2L) stop("terminal jobs were polled again")
    data.frame(State = c("CANCELLED by 12345", if (calls == 1L) "RUNNING" else "COMPLETED"))
  }, .package = "fmri.pipeline")
  expect_output(expect_false(wait_for_job(c("123", "124"), scheduler = "slurm",
    repolling_interval = 0.1)), "failed")
  expect_equal(calls, 2L)
})

test_that("nonfatal submission errors warn and return NULL without tracking a job", {
  script <- tempfile()
  writeLines("exit 0", script)
  on.exit(unlink(script))
  local_mocked_bindings(
    system2 = function(command, args, stdout, stderr, ...) {
      writeLines("scheduler rejected submission", stderr)
      1L
    },
    insert_tracked_job = function(...) stop("must not track a failed submission"),
    add_tracked_job_parent = function(...) stop("must not track a failed submission"),
    .package = "fmri.pipeline"
  )
  for (scheduler in c("slurm", "torque")) {
    expect_warning(result <- cluster_job_submit(script, scheduler = scheduler,
      fail_on_error = FALSE, echo = FALSE, wait_jobs = "123"), "scheduler rejected submission")
    expect_null(result)
    expect_error(cluster_job_submit(script, scheduler = scheduler,
      fail_on_error = TRUE, echo = FALSE), "scheduler rejected submission")
  }
})

test_that("successful submission preserves the job id and command metadata", {
  script <- tempfile()
  writeLines("exit 0", script)
  on.exit(unlink(script))
  local_mocked_bindings(system2 = function(command, args, stdout, stderr, ...) {
    writeLines("Submitted batch job 12345", stdout)
    0L
  }, .package = "fmri.pipeline")
  result <- cluster_job_submit(script, scheduler = "slurm", echo = FALSE)
  expect_equal(as.character(result), "12345")
  expect_match(attr(result, "cmd"), "sbatch")
})

test_that("unit tests reject unmocked scheduler submission", {
  script <- tempfile()
  writeLines("exit 0", script)
  on.exit(unlink(script))
  for (scheduler in c("slurm", "torque")) {
    expect_error(cluster_job_submit(script, scheduler = scheduler, echo = FALSE),
      "Unmocked scheduler command in unit test")
  }
})
