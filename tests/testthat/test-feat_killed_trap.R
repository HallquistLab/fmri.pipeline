test_that("interrupting a running FEAT child records failure and retains its reservation", {
  skip_on_os("windows")
  root <- withr::local_tempdir()
  fit <- make_valid_feat(root)
  fsf <- file.path(root, "queued.fsf")
  out <- file.path(root, "queued.feat")
  feat_reserve(fsf, "owner")
  script <- file.path(root, "interrupt.sh")
  runner <- system.file("bash", "feat_runner.sh", package = "fmri.pipeline")
  # Signal the PID recorded by the same 'feat_runner ... &' used by the scheduler.
  writeLines(c(
    paste("source", shQuote(runner)),
    "feat_level=1", "feat_binary=estimator", "feat_use_retry=0",
    "write_job_manifest() { :; }", "Rscript() { :; }",
    'estimator() { touch "$odir/started"; exec sleep 10; }',
    paste("mkdir", shQuote(out)),
    paste("feat_runner", shQuote(fsf), shQuote(out), "owner keep &"),
    "worker=$!",
    paste0("for ((i=0; i<200; i++)); do [ -e ", shQuote(file.path(out, "started")), " ] && break; sleep 0.02; done"),
    "kill -TERM \"$worker\"", "wait \"$worker\""
  ), script)
  expect_equal(system2("bash", shQuote(script)), 143L)
  expect_true(file.exists(file.path(out, ".feat_fail")))
  expect_false(file.exists(file.path(out, ".feat_complete")))
  expect_true(dir.exists(feat_lock_dir(fsf)))
  expect_true(get_feat_status(fit$out)$feat_complete)
})
