test_that("L3 inference round trips through the main configuration", {
  gpa <- create_mock_gpa(n_subjects = 8)
  gpa$glm_software <- "fsl"
  config <- list(l3_models = list(group = list(model_formula = "~ age", l3_input_mode = "subject_rows",
    contrasts = list(diagonal = TRUE), inference = list(method = "ptfce", fwe_alpha = c(0.05, 0.0005),
      two_sided = FALSE, compute = list(scheduler = "local", memgb_per_command = 4)))))
  path <- tempfile(fileext = ".yaml")
  yaml::write_yaml(config, path)
  gpa <- suppressWarnings(build_l3_models(gpa, from_spec_file = path))
  spec <- gpa$l3_models$models$group$inference
  expect_s3_class(spec, "fwe_spec")
  expect_identical(spec$targets$l3_model, "group")
  expect_false(spec$method$two_sided)
  expect_equal(spec$compute$memgb_per_command, 4)
  gpa$fwe_correction <- list(additional = fwe_spec("additional", list(l3_model = "group"), two_sided = TRUE))
  out <- tempfile(fileext = ".yaml")
  suppressWarnings(export_glm_config(gpa, out))
  restored <- suppressWarnings(build_l3_models(gpa, from_spec_file = out))
  expect_equal(restored$l3_models$models$group$inference, spec)
  expect_equal(restored$fwe_correction, gpa$fwe_correction)

  # Absence of inference in an imported configuration disables it.
  config$l3_models$group$inference <- NULL
  yaml::write_yaml(config, path)
  disabled <- suppressWarnings(build_l3_models(gpa, from_spec_file = path))
  expect_null(disabled$l3_models$models$group$inference)
  expect_length(disabled$fwe_correction, 0L)
})

test_that("model inference rejects unrelated targets and unavailable methods", {
  expect_error(fmri.pipeline:::normalize_l3_inference(list(targets = list(l3_model = "other")), "group"), "own L3 model")
  expect_error(fmri.pipeline:::normalize_l3_inference(list(method = "afni_3dclustsim_permutation"), "group"), "No FWE planning adapter")
  expect_error(fmri.pipeline:::normalize_l3_inference(list(targets = list(l3_model = "group", l3_cope_name = "missing")), "group", "age"), "subset")
})

test_that("pipeline correction selects requested models and validates their backend", {
  gpa <- create_mock_gpa(include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$l3_models$models$l3_model1$inference <- list(method = "ptfce")
  specs <- fmri.pipeline:::pipeline_fwe_specs(gpa)
  expect_named(specs, "l3_model1_ptfce")
  expect_length(fmri.pipeline:::pipeline_fwe_specs(gpa, character(0)), 0L)
  expect_silent(fmri.pipeline:::validate_pipeline_fwe(gpa, specs))
  gpa$scheduler <- "slurm"
  specs[[1L]]$compute$scheduler <- "torque"
  expect_error(fmri.pipeline:::validate_pipeline_fwe(gpa, specs, automatic = TRUE), "pipeline scheduler")
  expect_silent(fmri.pipeline:::validate_pipeline_fwe(gpa, specs))
  specs[[1L]]$compute$scheduler <- "local"
  expect_silent(fmri.pipeline:::validate_pipeline_fwe(gpa, specs, automatic = TRUE))
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  expect_error(fmri.pipeline:::validate_pipeline_fwe(gpa, specs), "requires an FSL L3")
})

test_that("correction-only execution collects results without estimating GLMs", {
  root <- tempfile("fwe_only_")
  gpa <- make_fwe_execution_test_gpa(root)
  gpa$glm_software <- "fsl"
  gpa$l3_models <- list(models = list(l3_group = list(name = "l3_group",
    inference = list(method = "ptfce", targets = list(l3_model = "l3_group", l3_cope_name = "age")))))
  worker <- make_fake_ptfce_worker(root)
  run <- run_fwe_correction(gpa, scheduler = "local", worker_script = worker)
  expect_identical(run$fwe_status$l3_group_ptfce, "complete")
  expect_true(run$fwe_runs$l3_group_ptfce$result$complete)
  skipped <- run_fwe_correction(run, scheduler = "local", worker_script = worker)
  expect_identical(skipped$fwe_runs$l3_group_ptfce$execution$execution_status, "skipped_complete")
})

test_that("failed inference persists separate status while leaving GLM state intact", {
  root <- tempfile("fwe_pipeline_failure_")
  gpa <- make_fwe_execution_test_gpa(root)
  gpa$glm_software <- "fsl"
  gpa$output_locations <- list(object_cache = file.path(root, "gpa.rds"))
  gpa$l3_models <- list(models = list(l3_group = list(name = "l3_group",
    inference = list(method = "ptfce", targets = list(l3_model = "l3_group", l3_cope_name = "age")))))
  expect_error(run_fwe_correction(gpa, scheduler = "local", worker_script = make_fake_ptfce_worker(root, fail = TRUE)), "FWE execution failed")
  saved <- readRDS(gpa$output_locations$object_cache)
  expect_identical(saved$fwe_status$l3_group_ptfce, "failed")
  expect_equal(saved$l3_model_setup, gpa$l3_model_setup)
  expect_error(refresh_fwe_correction(saved, require_complete = TRUE), "Group estimation completed")
  expect_identical(readRDS(gpa$output_locations$object_cache)$fwe_status$l3_group_ptfce, "incomplete")
})
