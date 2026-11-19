afni_test_fixture <- function(path, estimation = "OLS") {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  set.seed(724)
  n <- 50L
  x <- rnorm(n * 2L)
  images <- file.path(path, c("run one.nii.gz", "run two.nii.gz"))
  for (i in 1:2) {
    rows <- seq_len(n) + (i - 1L) * n
    # Equal time series in all voxels, with known task effect and run offset.
    y <- 100 + 20 * i + 2 * x[rows] + rnorm(n, sd = .4)
    RNifti::writeNifti(array(rep(y, each = 64L), c(4, 4, 4, n)), images[i])
  }
  mask <- file.path(path, "mask.nii.gz")
  RNifti::writeNifti(array(1, c(4, 4, 4)), mask)
  gpa <- create_mock_gpa(n_subjects = 1L, include_l1_models = TRUE, output_directory = path)
  gpa$glm_software <- "afni"
  gpa$run_data$run_number <- c(2L, 5L)
  gpa$glm_settings <- list(afni = list(estimation = estimation, polort = 0L, mask = mask))
  gpa$output_locations$afni_l1_directory <- file.path(path, "afni", "sub-{id}", "ses-{session}", "{l1_model}")
  gpa$l1_models$models$model1$regressors <- "task:a"
  gpa$l1_models$models$model1$contrasts <- matrix(c(1, -1), 2, 1,
    dimnames = list(c("positive task", "negative task"), "task:a"))
  d <- structure(list(runs_to_output = c(2L, 5L), run_volumes = c(n, n), drop_volumes = c(0L, 0L),
    tr = 2, run_nifti = images, design_convolved = list(
      run_number2 = data.frame("task:a" = x[seq_len(n)], check.names = FALSE),
      run_number5 = data.frame("task:a" = x[seq_len(n) + n], check.names = FALSE))), class = "bdm")
  list(gpa = gpa, d = d)
}

test_that("AFNI design aligns nonconsecutive runs, drops and block nuisance columns", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  f$d$drop_volumes <- c(2L, 0L)
  cf <- file.path(tmp, c("conf1.txt", "conf2.txt"))
  write.table(cbind(1, seq_len(52)), cf[1], row.names = FALSE, col.names = FALSE)
  write.table(seq_len(50), cf[2], row.names = FALSE, col.names = FALSE)
  cen <- file.path(tmp, "censor.txt")
  writeLines(as.character(c(1, 1, 0, rep(1, 49))), cen)
  d <- fmri.pipeline:::afni_l1_design(f$d, "task:a", cf, c(cen, NA_character_))
  expect_equal(d$starts, c(0, 50))
  expect_equal(dim(d$design), c(100, 3))
  expect_equal(d$censor[1], 0)
  expect_true(all(d$design[51:100, 2] == 0))
  expect_true(all(d$design[1:50, 3] == 0))
  expect_equal(d$design[1:50, 2], 3:52 - mean(4:52))
  writeLines("1", cf[1])
  expect_error(fmri.pipeline:::afni_l1_design(f$d, "task:a", cf), "expected 50 retained")
  writeLines(rep("2", 50), cen)
  expect_error(fmri.pipeline:::afni_l1_design(f$d, "task:a", censor_files = c(cen, NA_character_)), "binary column")
})

test_that("AFNI-only backend resolution supplies its own L3 inputs", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "afni"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  b <- fmri.pipeline:::get_glm_backends(gpa)$afni
  expect_true(b$runs_l1)
  expect_false(b$runs_l2)
  req <- fmri.pipeline:::resolve_model_l3_requirements(gpa)
  expect_equal(req$producer_backend, "afni")
  expect_equal(req$producer_level, 1L)
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  expect_equal(fmri.pipeline:::get_effective_model_backends(gpa, 3L, type = "producer")$l3_model1, "afni")
})

test_that("AFNI status requires all outputs and invalidates changed specifications", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  out <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
  expect_false(out$afni_complete)
  manifest <- readLines(file.path(out$afni_dir, "expected_outputs.txt"))
  for (file in manifest) writeLines("dummy", file)
  file.create(file.path(out$afni_dir, ".afni_complete"))
  expect_true(get_afni_glm_status(out$afni_dir)$afni_complete)
  out2 <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
  expect_true(out2$afni_complete)
  unlink(tail(manifest, 1))
  expect_false(get_afni_glm_status(out$afni_dir)$afni_complete)
  writeLines("dummy", tail(manifest, 1))
  f$d$design_convolved$run_number2[1, 1] <- 42
  out3 <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
  expect_false(out3$afni_complete)
  expect_false(file.exists(file.path(out$afni_dir, ".afni_complete")))
})

test_that("AFNI OLS and REML scripts execute and expose correct named contrasts", {
  for (exe in c("3dDeconvolve", "3dREMLfit", "3dTcat", "3dcopy", "3dbucket")) {
    skip_if(!nzchar(Sys.which(exe)), paste(exe, "not installed"))
  }
  tmp <- withr::local_tempdir()
  for (est in c("OLS", "REML")) {
    f <- afni_test_fixture(file.path(tmp, est), est)
    out <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
    status <- system2("bash", shQuote(out$afni_script))
    expect_equal(status, 0, info = paste(readLines(file.path(out$afni_dir, "afni_glm.log")), collapse = "\n"))
    if (status != 0L) next
    expect_true(get_afni_glm_status(out$afni_dir)$afni_complete)
    con <- as.numeric(RNifti::readNifti(file.path(out$afni_dir, "con_0001.nii.gz")))[1]
    neg <- as.numeric(RNifti::readNifti(file.path(out$afni_dir, "con_0002.nii.gz")))[1]
    expect_equal(neg, -con, tolerance = 1e-5)
    expect_equal(con, 2, tolerance = .1)
    if (est == "OLS") {
      y <- unlist(lapply(f$d$run_nifti, function(p) RNifti::readNifti(p)[1, 1, 1, ]))
      x <- unlist(lapply(f$d$design_convolved, function(d) d[[1]]))
      m <- lm(y ~ x + factor(rep(1:2, each = 50)))
      expect_equal(con, unname(coef(m)[2]), tolerance = 1e-5)
      t <- as.numeric(RNifti::readNifti(file.path(out$afni_dir, "tstat_0001.nii.gz")))[1]
      expect_equal(t, unname(summary(m)$coefficients[2, "t value"]), tolerance = 1e-5)
    }
  }
})

test_that("AFNI scheduler wrapper preserves failures and vector dependencies", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  gpa <- f$gpa
  gpa$scheduler <- "slurm"
  gpa$batch_run <- list(batch_directory = tmp)
  script <- file.path(tmp, "fail.sh")
  writeLines(c("#!/bin/bash", "exit 19"), script)
  captured <- NULL
  testthat::local_mocked_bindings(
    cluster_job_submit = function(script, ...) { captured <<- script; "123" },
    get_compute_environment = function(...) character(0), .package = "fmri.pipeline")
  job <- fmri.pipeline:::submit_afni_glm_job(gpa, script, "test", "1:00:00", "2", 2L, c("10", "11"))
  expect_equal(job, "123")
  lines <- readLines(captured)
  expect_true(any(grepl("afterok:10:11", lines, fixed = TRUE)))
  expect_true(any(grepl("--cpus-per-task=2", lines, fixed = TRUE)))
  withr::local_envvar(c(SLURM_JOB_ID = "123"))
  expect_equal(system2("bash", shQuote(captured)), 19L)
})

test_that("AFNI group OLS respects covariate order and signed contrast maps", {
  for (exe in c("3dDeconvolve", "3dTcat", "3dcopy", "3dbucket", "3dmask_tool")) {
    skip_if(!nzchar(Sys.which(exe)), paste(exe, "not installed"))
  }
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  gpa <- f$gpa
  gpa$l3_models <- structure(list(models = list(group = list(l3_input_mode = "per_session"))), class = "hi_model_set")
  gpa$output_locations$afni_l3_directory <- file.path(tmp, "group", "{l3_model}", "{l1_cope_name}")
  ids <- paste0("s", 1:8)
  age <- c(-3, 2, -2, 1, 4, -1, 3, 0)
  y <- -5 + .7 * age + c(.2, -.3, .1, -.2, .3, .2, -.1, -.2)
  input <- data.frame(id = ids, session = 1L, l1_model = "model1", l3_model = "group",
    l1_cope_name = "positive task", mask_file = gpa$glm_settings$afni$mask,
    InputFile = file.path(tmp, paste0(ids, ".nii.gz")))
  for (i in seq_along(y)) RNifti::writeNifti(array(y[i], c(4, 4, 4)), input$InputFile[i])
  fit <- list(metadata = input[, c("id", "session")], model_matrix = cbind("(Intercept)" = 1, age = age),
    contrasts = matrix(c(0, 1), 1, dimnames = list("age slope", c("(Intercept)", "age"))))
  testthat::local_mocked_bindings(respecify_l3_model = function(...) fit, .package = "fmri.pipeline")
  input <- input[c(8, 3, 2, 6, 4, 1, 7, 5), ]
  gpa$glm_settings$afni$mask <- NULL # exercise source-mask intersection
  out <- fmri.pipeline:::afni_l3_model(input, gpa)
  status <- system2("bash", shQuote(out$afni_script))
  expect_equal(status, 0L, info = paste(readLines(file.path(out$afni_dir, "afni_glm.log")), collapse = "\n"))
  if (status == 0L) {
    beta <- as.numeric(RNifti::readNifti(file.path(out$afni_dir, "con_0001.nii.gz")))[1]
    expect_equal(beta, unname(coef(lm(y ~ age))[2]), tolerance = 1e-5)
    expect_true(fmri.pipeline:::get_afni_l3_status(out$output_file)$afni_complete)
  }
})

test_that("AFNI harvest uses manifest names and excludes incomplete/excluded subjects", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  out <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
  for (file in readLines(file.path(out$afni_dir, "expected_outputs.txt"))) writeLines("dummy", file)
  file.create(file.path(out$afni_dir, ".afni_complete"))
  gpa <- f$gpa
  gpa$l1_model_setup <- list(afni = out)
  gpa$l3_models <- create_mock_l3_models()
  inputs <- fmri.pipeline:::harvest_l3_inputs_afni(gpa)
  expect_length(inputs$l3_model1, 2L)
  expect_setequal(vapply(inputs$l3_model1, function(x) x$l1_cope_name, character(1)),
    c("positive task", "negative task"))
  file.create(file.path(out$afni_dir, ".afni_failed"))
  expect_null(fmri.pipeline:::harvest_l3_inputs_afni(gpa))
  unlink(file.path(out$afni_dir, ".afni_failed"))
  gpa$subject_data$exclude_subject <- TRUE
  expect_null(fmri.pipeline:::harvest_l3_inputs_afni(gpa))
})

test_that("AFNI excludes flagged runs before concatenation and validates image grids", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  f$gpa$run_data$exclude_run[1] <- TRUE
  out <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
  spec <- readRDS(file.path(out$afni_dir, "afni_setup.rds"))
  expect_equal(nrow(spec$design), 50L)
  expect_equal(spec$inputs, f$d$run_nifti[2])
  f$gpa$run_data$exclude_run[] <- TRUE
  expect_null(afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1"))
  f$gpa$run_data$exclude_run[] <- FALSE
  RNifti::writeNifti(array(1, c(3, 4, 4)), f$gpa$glm_settings$afni$mask)
  expect_error(afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1"), "spatial geometry")
})

test_that("AFNI execution applies leading volume drops and censoring to the same samples", {
  for (exe in c("3dDeconvolve", "3dTcat", "3dcopy", "3dbucket")) skip_if(!nzchar(Sys.which(exe)))
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  f$d$drop_volumes[1] <- 2L
  original <- RNifti::readNifti(f$d$run_nifti[1])
  padded <- abind::abind(array(999, c(4, 4, 4, 2)), original, along = 4L)
  RNifti::writeNifti(padded, f$d$run_nifti[1])
  censor <- file.path(tmp, "censor.1D")
  keep <- rep(1L, 50L)
  keep[c(3, 9, 12)] <- 0L
  writeLines(as.character(c(1, 1, keep)), censor)
  out <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1", censor_files = c(censor, NA_character_))
  status <- system2("bash", shQuote(out$afni_script))
  expect_equal(status, 0L)
  if (status == 0L) {
    x <- unlist(lapply(f$d$design_convolved, function(d) d[[1]]))
    y <- c(original[1, 1, 1, ], RNifti::readNifti(f$d$run_nifti[2])[1, 1, 1, ])
    m <- lm(y ~ x + factor(rep(1:2, each = 50)), subset = c(keep, rep(1L, 50L)) == 1L)
    beta <- as.numeric(RNifti::readNifti(file.path(out$afni_dir, "con_0001.nii.gz")))[1]
    expect_equal(beta, unname(coef(m)[2]), tolerance = 1e-5)
  }
})

test_that("public L1 setup routes a built pipeline design into AFNI", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  gpa <- f$gpa
  gpa$l1_models <- NULL
  gpa$drop_volumes <- 0L
  gpa$trial_data$run_number <- ifelse(gpa$trial_data$run_number == 1L, 2L, 5L)
  gpa$run_data$run_nifti <- f$d$run_nifti
  gpa$run_data$run_volumes <- 50L
  gpa$run_data$run_nifti_present <- TRUE
  gpa$run_data$tr <- 2
  gpa$run_data$nvoxels <- 64L
  gpa$log_json <- gpa$log_txt <- FALSE
  gpa$use_preconvolve <- TRUE
  gpa$additional <- list(bdm_args = list(plot = FALSE))
  spec <- file.path(tmp, "model.yaml")
  yaml::write_yaml(list(onsets = "onset", durations = "duration",
    events = list(cue = list(onset = "onset", duration = "duration")),
    signals = list(cue = list(event = "cue", value_fixed = 1, normalization = "none")),
    l1_models = list(model1 = list(signals = "cue", contrasts = list(diagonal = TRUE)))), spec)
  gpa <- build_l1_models(gpa, from_spec_file = spec)
  gpa <- setup_l1_models(gpa, "model1")
  expect_s3_class(gpa$l1_model_setup, "l1_setup")
  expect_equal(nrow(gpa$l1_model_setup$afni), 1L)
  expect_equal(nrow(gpa$l1_model_setup$fsl), 0L)
  if (nrow(gpa$l1_model_setup$afni)) {
    expect_true(file.exists(gpa$l1_model_setup$afni$afni_script))
    expect_false(gpa$l1_model_setup$afni$afni_complete)
  }
})

test_that("saved AFNI placeholders upgrade while custom backend definitions survive", {
  specs <- fmri.pipeline:::default_glm_backend_specs()
  specs$afni$l1_setup <- "__not_implemented__"
  specs$afni$l3_setup <- "afni_3dlmer_setup"
  upgraded <- fmri.pipeline:::upgrade_glm_backend_specs(specs)
  expect_equal(upgraded$afni$l1_setup, "afni_l1_model")
  specs$afni$l3_setup <- "custom_group_setup"
  expect_identical(fmri.pipeline:::upgrade_glm_backend_specs(specs), specs)
})

test_that("default AFNI masking works on the retained time-series inputs", {
  for (exe in c("3dDeconvolve", "3dTcat", "3dAutomask", "3dmask_tool", "3dbucket")) skip_if(!nzchar(Sys.which(exe)))
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  f$gpa$glm_settings$afni$mask <- NULL
  grid <- expand.grid(x = 1:16, y = 1:16, z = 1:16)
  # Spatial variation avoids a degenerate histogram in AFNI clip-level estimation.
  brain <- as.numeric((grid$x - 8)^2 + (grid$y - 8)^2 + (grid$z - 8)^2 < 36) *
    (1 + (grid$x + grid$y + grid$z) / 32)
  for (i in 1:2) {
    y <- RNifti::readNifti(f$d$run_nifti[i])[1, 1, 1, ]
    RNifti::writeNifti(array(outer(brain, y), c(16, 16, 16, 50)), f$d$run_nifti[i])
  }
  out <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
  status <- system2("bash", shQuote(out$afni_script))
  expect_equal(status, 0L, info = paste(readLines(file.path(out$afni_dir, "afni_glm.log")), collapse = "\n"))
  if (status == 0L) {
    mask <- as.numeric(RNifti::readNifti(file.path(out$afni_dir, "mask.nii.gz")))
    expect_true(any(mask > 0))
    expect_equal(mask[1L], 0)
    expect_lt(sum(mask > 0), length(mask))
  }
})

test_that("AFNI runner refreshes completion and honors model selection and rerun", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  out <- afni_l1_model("sub1", 1L, d_obj = f$d, gpa = f$gpa, model_name = "model1")
  gpa <- f$gpa
  gpa$l1_model_setup <- structure(list(afni = out), class = "l1_setup")
  submitted <- character(0)
  testthat::local_mocked_bindings(submit_afni_glm_job = function(gpa, script, ...) {
    submitted <<- c(submitted, script)
    "job1"
  }, .package = "fmri.pipeline")
  expect_equal(run_afni_sepjobs(gpa, model_names = "model1"), "job1")
  expect_length(run_afni_sepjobs(gpa, model_names = "another_model"), 0L)
  for (file in readLines(file.path(out$afni_dir, "expected_outputs.txt"))) writeLines("dummy", file)
  file.create(file.path(out$afni_dir, ".afni_complete"))
  expect_length(run_afni_sepjobs(gpa), 0L) # ignores stale FALSE in cached setup table
  expect_equal(run_afni_sepjobs(gpa, rerun = TRUE), "job1")
  unlink(file.path(out$afni_dir, "expected_outputs.txt"))
  expect_false(fmri.pipeline:::get_afni_l3_status(out$output_file)$afni_complete)
  expect_error(run_afni_sepjobs(gpa), "Missing AFNI GLM output manifest")
})

test_that("public group setup consumes AFNI L1 manifests without an L2 model", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  gpa <- f$gpa
  class(gpa$l1_models) <- c("l1_model_set", "list")
  gpa$multi_run <- TRUE
  gpa$log_json <- gpa$log_txt <- FALSE
  gpa$subject_data <- data.frame(id = paste0("sub", 1:6), session = 1L, exclude_subject = FALSE)
  gpa$run_data <- do.call(rbind, lapply(gpa$subject_data$id, function(id) {
    r <- f$gpa$run_data
    r$id <- id
    r$exclude_subject <- FALSE
    r
  }))
  gpa$output_locations$afni_l3_directory <- file.path(tmp, "group", "{l3_model}", "{l1_cope_name}")
  setup <- lapply(1:6, function(i) {
    dir <- file.path(tmp, paste0("sub", i))
    dir.create(dir)
    coef <- file.path(dir, "con_0001.nii.gz")
    mask <- file.path(dir, "mask.nii.gz")
    RNifti::writeNifti(array(i, c(4, 4, 4)), coef)
    RNifti::writeNifti(array(1, c(4, 4, 4)), mask)
    write.table(data.frame(contrast_number = 1L, contrast_name = "positive task", coef_file = coef),
      file.path(dir, "contrasts.tsv"), sep = "\t", row.names = FALSE)
    writeLines(c(coef, mask), file.path(dir, "expected_outputs.txt"))
    file.create(file.path(dir, ".afni_complete"))
    data.frame(id = paste0("sub", i), session = 1L, l1_model = "model1", afni_dir = dir)
  })
  gpa$l1_model_setup <- structure(list(afni = do.call(rbind, setup)), class = c("l1_setup", "list"))
  spec <- file.path(tmp, "group.yaml")
  yaml::write_yaml(list(l3_models = list(group_mean = list(model_formula = "~ 1",
    level = 3L, l3_input_mode = "per_session", contrasts = list(diagonal = TRUE)))), spec)
  gpa <- build_l3_models(gpa, from_spec_file = spec)
  gpa <- setup_l3_models(gpa, l3_model_names = "group_mean", backend = "afni")
  expect_true(gpa$l3_setup_status$success)
  expect_equal(nrow(gpa$l3_model_setup$afni), 1L)
  expect_equal(gpa$l3_model_setup$afni$l1_cope_name, "positive task")
  expect_true(file.exists(gpa$l3_model_setup$afni$afni_script))
  if (all(nzchar(Sys.which(c("3dDeconvolve", "3dTcat", "3dcopy", "3dbucket"))))) {
    out <- gpa$l3_model_setup$afni
    expect_equal(system2("bash", shQuote(out$afni_script)), 0L)
    gpa <- refresh_glm_status(gpa, level = 3L)
    expect_true(gpa$l3_model_setup$afni$afni_complete)
    beta <- as.numeric(RNifti::readNifti(file.path(out$afni_dir, "con_0001.nii.gz")))[1]
    expect_equal(beta, mean(1:6), tolerance = 1e-5)
  }
})

test_that("older pipeline objects get dedicated AFNI output directories", {
  tmp <- withr::local_tempdir()
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l3_models = TRUE, output_directory = tmp)
  l1 <- fmri.pipeline:::get_output_directory(id = "sub1", session = 1L, l1_model = "model1",
    gpa = gpa, glm_software = "afni")
  l3 <- fmri.pipeline:::get_output_directory(l1_model = "model1", l3_model = "l3_model1",
    l1_contrast = "test", what = "l3", gpa = gpa, glm_software = "afni")
  expect_equal(as.character(l1), file.path(tmp, "afni_l1", "sub-sub1", "ses-1", "model1"))
  expect_equal(as.character(l3), file.path(tmp, "afni_l3", "L3m-l3_model1", "L1m-model1", "l1c-test"))
})


test_that("AFNI drops spike nuisance columns for volumes already censored", {
  tmp <- withr::local_tempdir()
  f <- afni_test_fixture(tmp)
  nuisance <- file.path(tmp, "spikes.1D")
  censor <- file.path(tmp, "censor.1D")
  conf <- cbind(censored_spike = as.numeric(1:50 == 3L), retained_spike = as.numeric(1:50 == 9L))
  write.table(conf, nuisance, row.names = FALSE, col.names = FALSE)
  writeLines(as.character(as.numeric(1:50 != 3L)), censor)
  d <- fmri.pipeline:::afni_l1_design(f$d, "task:a", c(nuisance, NA_character_), c(censor, NA_character_))
  expect_equal(ncol(d$design), 2L)
  expect_equal(d$design[1:50, 2], conf[, 2] - 1/49)
  expect_equal(sum(d$design[which(d$censor == 1L), 2]), 0, tolerance = 1e-12)
  expect_true(all(d$design[51:100, 2] == 0))
})
