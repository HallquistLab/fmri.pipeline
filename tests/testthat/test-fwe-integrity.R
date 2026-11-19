test_that("configuration changes, changed inputs, and corrupt outputs invalidate completion", {
  root <- tempfile("fwe_integrity_")
  gpa <- make_fwe_execution_test_gpa(root)
  spec <- fwe_spec("same_name", list(l3_cope_name = "group"))
  plan <- plan_fwe_correction(gpa, spec, source = "setup")
  write_fwe_test_result(plan)
  expect_identical(refresh_fwe_plan(plan)$tasks$status, "complete")

  one_sided <- spec
  one_sided$method$two_sided <- FALSE
  changed <- plan_fwe_correction(gpa, one_sided, source = "setup")
  expect_identical(changed$tasks$status, "ready")
  expect_false(collect_fwe_results(changed, persist = FALSE)$complete)
  expect_equal(nrow(fwe_result_artifacts(changed)), 0L)

  write_fwe_test_image(plan$tasks$zstat_file, 2)
  expect_identical(refresh_fwe_plan(plan)$tasks$status, "ready")
  plan <- refresh_fwe_plan(plan)
  write_fwe_test_result(plan)
  expect_identical(refresh_fwe_plan(plan)$tasks$status, "complete")
  cat("corrupt", file = plan$artifacts$artifact_file[1L], append = TRUE)
  expect_identical(refresh_fwe_plan(plan)$tasks$status, "ready")
})

test_that("unverified and empty outputs never count as successful inference", {
  gpa <- make_fwe_execution_test_gpa(tempfile("fwe_empty_"))
  plan <- plan_fwe_correction(gpa, fwe_spec("empty", list(l3_cope_name = "age")), source = "setup")
  dir.create(plan$tasks$output_directory, recursive = TRUE)
  file.create(plan$artifacts$artifact_file)
  expect_identical(refresh_fwe_plan(plan)$tasks$status, "ready")
  result <- collect_fwe_results(plan, persist = FALSE)
  expect_false(result$complete)
  expect_true(all(result$manifest$artifact_status == "unverified"))
  expect_false(any(result$manifest$can_use_as_mask))
  expect_error(collect_fwe_results(plan, require_complete = TRUE), "unverified")
  expect_error(fmri.pipeline:::validate_ptfce_outputs(plan$tasks, plan$artifacts), "empty")
})

test_that("a failed forced rerun invalidates the previous success receipt", {
  root <- tempfile("fwe_failed_rerun_")
  gpa <- make_fwe_execution_test_gpa(root)
  plan <- plan_fwe_correction(gpa, fwe_spec("rerun", list(l3_cope_name = "age")), source = "setup")
  write_fwe_test_result(plan)
  worker <- make_fake_ptfce_worker(root, fail = TRUE)
  run <- run_fwe_plan(plan, scheduler = "local", force = TRUE, worker_script = worker, stop_on_error = FALSE)
  expect_identical(run$execution$execution_status, "failed")
  expect_false(run$result$complete)
  expect_false(file.exists(fmri.pipeline:::fwe_completion_file(plan$tasks)))
})

test_that("unused minimum cluster options are explicitly deprecated", {
  expect_warning(spec <- fwe_spec("old", list(l3_model = "group"), min_cluster_voxels = 10), "deprecated")
  expect_null(spec$method$min_cluster_voxels)
})
