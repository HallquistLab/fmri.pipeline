# Set RUN_FEAT_INTEGRATION=true and install the package for the Rscript children.
test_that("FSL estimates and revalidates the L1, L2 and L3 chain", {
  skip_if(Sys.getenv("RUN_FEAT_INTEGRATION") != "true")
  skip_if_not(nzchar(Sys.which("feat")))
  root <- withr::local_tempdir()
  config <- file.path(root, "fsl_sub.yml")
  writeLines("method: shell", config)
  withr::local_envvar(c(FSLSUB_CONF = config, TMPDIR = root))
  gpa <- create_mock_gpa(n_subjects = 3, n_runs = 2, include_l1_models = TRUE, output_directory = root)
  gpa$lgr_threshold <- "warn"
  gpa$multi_run <- TRUE
  gpa$use_preconvolve <- TRUE
  gpa$glm_settings <- list(fsl = list())
  gpa$output_locations$feat_l1_directory <- file.path(root, "l1", "{id}", "{l1_model}")
  gpa$output_locations$feat_l2_directory <- file.path(root, "l2", "{id}", "{l1_contrast}")
  gpa$output_locations$feat_l3_directory <- file.path(root, "l3")
  gpa$output_locations$feat_l3_fsf <- "group.fsf"
  gpa$l1_models$models$model1$regressors <- "task"
  gpa$l1_models$models$model1$contrasts <- matrix(1, 1, dimnames = list("task", "task"))
  gpa$confound_settings$l1_confound_regressors <- NULL
  gpa$subject_data$exclude_subject <- FALSE
  gpa$run_data$exclude_subject <- FALSE
  gpa$additional$feat_l1_args <- list(thresh = 0)

  set.seed(102)
  task <- sin(seq_len(60) / 6)
  timing <- file.path(root, "task.txt")
  writeLines(as.character(task), timing)
  runs <- gpa$run_data
  for (i in seq_len(nrow(runs))) {
    run <- runs$run_number[i]
    key <- paste0("run_number", run)
    nii <- file.path(root, paste0("input", i, ".nii.gz"))
    bold <- array(rnorm(9^3 * 60, 100, 2), c(9, 9, 9, 60))
    RNifti::writeNifti(sweep(bold, 4, 2 * task, "+"), nii)
    bdm <- structure(list(run_nifti = nii, run_volumes = 60L, drop_volumes = 0L,
      runs_to_output = run, tr = 1, design_convolved = setNames(list(data.frame(task)), key),
      timing_files = list(convolved = matrix(timing, 1, dimnames = list(key, "task")))), class = "bdm")
    fit <- fsl_l1_model(id = runs$id[i], session = 1L, d_obj = bdm, gpa = gpa,
      model_name = "model1", nvoxels = 729L, execute_feat = TRUE)
    expect_true(fit$feat_complete)
    runs$cope_file[i] <- file.path(fit$feat_dir, "stats", "cope1.nii.gz")
  }
  runs$l1_model <- "model1"
  runs$l1_cope_name <- "task"
  runs$l1_cope_number <- 1L
  reference <- file.path(dirname(dirname(runs$cope_file[1L])), "example_func.nii.gz")
  gpa$additional$feat_l2_args <- list(thresh = 0, featwatcher_yn = 0, regstandard = shQuote(reference, type = "cmd"))
  gpa$additional$feat_l3_args <- gpa$additional$feat_l2_args
  model <- function(data, level) {
    create_new_hi_model(data = data, level = level, spec_list = list(
      name = "mean", model_formula = "~ 1", diagonal = TRUE, num2fac = character(),
      covariate_transform = c(), reference_level = c()), lg = lgr::get_logger("test"))
  }
  gpa$l2_models$models$mean <- model(gpa$run_data, 2L)
  gpa$l3_models$models$mean <- model(gpa$subject_data, 3L)
  # Execute the generated scheduler scripts locally, including their validation and markers.
  gpa$scheduler <- "slurm"
  gpa$parallel$fsl <- list(l2_feat_time = "00:10:00", l2_feat_memgb = 1,
    l3_feat_time = "00:10:00", l3_feat_memgb = 1, l3_feat_cpusperjob = 1)
  withr::local_envvar(c(SLURM_SUBMIT_DIR = root))
  local_mocked_bindings(cluster_job_submit = function(batch_file, ...) {
    stopifnot(system2("bash", shQuote(batch_file)) == 0L)
    "1"
  }, .package = "fmri.pipeline")
  subjects <- dplyr::bind_rows(lapply(split(runs, runs$id), fsl_l2_model, l2_model = "mean", gpa = gpa))
  gpa$l2_model_setup <- structure(list(fsl = subjects), class = "l2_setup")
  run_feat_sepjobs(gpa, level = 2L)
  expect_identical(run_feat_sepjobs(gpa, level = 2L), character())
  subjects$cope_file <- file.path(subjects$feat_dir, "cope1.feat", "stats", "cope1.nii.gz")
  subjects$l2_cope_name <- "EV_Intercept"
  subjects$l3_model <- "mean"
  group <- fsl_l3_model(subjects, gpa)
  gpa$l3_model_setup <- structure(list(fsl = group), class = "l3_setup")
  run_feat_sepjobs(gpa, level = 3L)
  status <- get_feat_status(group$feat_dir, group$feat_fsf)
  expect_true(status$feat_complete)
  expect_identical(status$feat_requested_estimator, "flame12")
  expect_identical(status$feat_actual_estimator, "flame12")
  expect_identical(status$feat_flame1_fallback_count, 0L)
  expect_false(status$feat_actual_outlier_deweighting)
  expect_identical(get_feat_status(subjects$feat_dir[1L])$feat_actual_estimator, "fe")
  expect_identical(run_feat_sepjobs(gpa, level = 3L), character())
  write("0", timing, append = TRUE)
  expect_false(get_feat_status(group$feat_dir, group$feat_fsf)$feat_complete)
})
