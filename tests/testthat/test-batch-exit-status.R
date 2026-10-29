run_status_batch <- function(jobs, fork, pre = NULL, post = NULL) {
  directory <- tempfile("batch status ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  cluster_submit_shell_jobs(jobs, cpus_per_job = length(jobs), fork_jobs = fork,
    pre = pre, post = post, scheduler = "local", job_out_dir = directory,
    log_file = file.path(directory, "jobs.csv"), debug = TRUE)
  suppressWarnings(system2("sh", shQuote(file.path(directory, "job_1.bash")),
    stdout = FALSE, stderr = FALSE))
}

test_that("generated batches preserve failures and wait for all tasks", {
  skip_on_os("windows")
  for (fork in c(TRUE, FALSE)) {
    sibling <- tempfile()
    cleanup <- tempfile()
    skipped <- tempfile()
    on.exit(unlink(c(sibling, cleanup, skipped)), add = TRUE)
    # A later command in the same task, another task, and cleanup all succeed.
    jobs <- list(c("sh -c 'exit 7'", paste("touch", shQuote(skipped))),
      c("sleep 0.1", paste("touch", shQuote(sibling))))
    status <- run_status_batch(jobs, fork, post = c(paste("touch", shQuote(cleanup)), "exit 0"))
    expect_equal(status, 7L)
    expect_false(file.exists(skipped))
    expect_true(file.exists(sibling))
    expect_true(file.exists(cleanup))
    expect_equal(run_status_batch(c("sh -c 'exit 9'", "true"), fork), 9L)
    expect_equal(run_status_batch(c("true", "sh -c 'exit 8'"), fork), 8L)
    expect_equal(run_status_batch(c("true", "true"), fork), 0L)
    expect_equal(run_status_batch("true", fork, post = c("false", "true")), 1L)
    expect_equal(run_status_batch("sh -c 'exit 7'", fork, post = "false"), 7L)
  }
})

test_that("batch setup failures stop execution and setup variables reach tasks", {
  skip_on_os("windows")
  for (fork in c(TRUE, FALSE)) {
    marker <- tempfile()
    on.exit(unlink(marker), add = TRUE)
    expect_equal(run_status_batch(paste("touch", shQuote(marker)), fork,
      pre = c("sh -c 'exit 6'", "true")), 6L)
    expect_false(file.exists(marker))
    expect_equal(run_status_batch('test "$batch_test_value" = available', fork,
      pre = "batch_test_value=available"), 0L)
  }
})
