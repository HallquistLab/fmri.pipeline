# Unit tests for setup_l1_models, setup_l2_models, setup_l3_models
# ==================================================================
# These tests focus on input validation and helper functions that can be
# tested in isolation without requiring external data or cluster resources.

library(testthat)

# ==============================================================================
# Mock factories specific to setup_l*_models testing
# ==============================================================================

#' Create a mock gpa suitable for setup_l1_models input validation
create_mock_gpa_for_l1_setup <- function(
    n_subjects = 2,
    n_runs = 2,
    include_l1_models = TRUE
) {
  gpa <- create_mock_gpa(
    n_subjects = n_subjects,
    n_runs = n_runs,
    include_l1_models = include_l1_models
  )
  
  # Add required columns for setup_l1_models
  gpa$run_data$run_nifti_present <- TRUE
  gpa$run_data$tr <- 1.0
  gpa$run_data$nvoxels <- 100000L
  
  # Add output locations
  gpa$output_locations$setup_l1_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l1_log_txt <- tempfile(fileext = ".txt")
  
  # Add logging threshold
  gpa$lgr_threshold <- "warn"  # Suppress debug messages during testing
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  
  # Add parallel settings
  gpa$parallel$l1_setup_cores <- 1L
  
  return(gpa)
}

#' Create a mock gpa suitable for setup_l2_models input validation
create_mock_gpa_for_l2_setup <- function(
    n_subjects = 2,
    n_runs = 2
) {
  gpa <- create_mock_gpa(
    n_subjects = n_subjects,
    n_runs = n_runs,
    include_l1_models = TRUE,
    include_l2_models = TRUE
  )
  
  # Add required columns for setup_l2_models
  gpa$run_data$exclude_run <- FALSE
  gpa$run_data$exclude_subject <- FALSE
  
  # Add output locations
  gpa$output_locations$setup_l2_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l2_log_txt <- tempfile(fileext = ".txt")
  gpa$output_locations$setup_l2_l1_cope_validity_tsv <- tempfile(fileext = ".tsv")
  gpa$output_locations$feat_l2_directory <- file.path(
    gpa$output_directory, "feat_l2", "sub-{id}", "ses-{session}",
    "L1m-{l1_model}", "l1c-{l1_cope_label}", "L2m-{l2_model}"
  )
  gpa$output_locations$feat_l2_id_scope_directory <- file.path(
    gpa$output_directory, "feat_l2", "sub-{id}",
    "L1m-{l1_model}", "l1c-{l1_cope_label}", "L2m-{l2_model}"
  )
  
  # Add logging threshold
  gpa$lgr_threshold <- "warn"
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  
  # Add parallel settings
  gpa$parallel$l2_setup_cores <- 1L
  gpa$glm_settings <- list(fsl = list(force_l2_creation = TRUE))
  gpa$additional <- list(feat_l2_args = list())
  gpa$multi_run <- TRUE
  gpa$subject_data$exclude_subject <- FALSE
  
  # Add l1_model_setup (would normally be created by setup_l1_models)
  gpa$l1_model_setup <- create_mock_l1_model_setup(gpa)
  
  # Set correct class for l1_models and l2_models
  class(gpa$l1_models) <- c("l1_model_set", "list")
  class(gpa$l2_models) <- c("hi_model_set", "list")
  gpa$l1_cope_names <- list(model1 = "event_vs_baseline")
  
  return(gpa)
}

#' Create a mock gpa suitable for setup_l3_models input validation
create_mock_gpa_for_l3_setup <- function(
    n_subjects = 3,
    n_runs = 2,
    multi_run = TRUE
) {
  gpa <- create_mock_gpa(
    n_subjects = n_subjects,
    n_runs = n_runs,
    include_l1_models = TRUE,
    include_l2_models = multi_run,
    include_l3_models = TRUE
  )
  
  # Add required columns
  gpa$run_data$exclude_run <- FALSE
  gpa$run_data$exclude_subject <- FALSE
  gpa$subject_data$exclude_subject <- FALSE
  
  # Add output locations
  gpa$output_locations$setup_l3_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l3_log_txt <- tempfile(fileext = ".txt")
  
  # Add settings
  gpa$lgr_threshold <- "warn"
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  gpa$multi_run <- multi_run
  gpa$glm_software <- "fsl"
  
  # Add l1_model_setup and l2_model_setup
  gpa$l1_model_setup <- create_mock_l1_model_setup(gpa)
  if (multi_run) {
    gpa$l2_model_setup <- create_mock_l2_model_setup(gpa)
  }
  
  # Set correct classes
  class(gpa$l1_models) <- c("l1_model_set", "list")
  if (multi_run) {
    class(gpa$l2_models) <- c("hi_model_set", "list")
  }
  class(gpa$l3_models) <- c("hi_model_set", "list")
  gpa$l1_cope_names <- list(model1 = "event_vs_baseline")
  
  return(gpa)
}

#' Create a mock l1_model_setup object
create_mock_l1_model_setup <- function(gpa) {
  n_rows <- nrow(gpa$run_data)
  feat_dirs <- replicate(n_rows, tempfile(fileext = ".feat"))
  feat_fsfs <- replicate(n_rows, tempfile(fileext = ".fsf"))
  for (ii in seq_len(n_rows)) {
    dir.create(file.path(feat_dirs[ii], "stats"), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(feat_dirs[ii], ".feat_complete"))
    file.create(file.path(feat_dirs[ii], "stats", "cope1.nii.gz"))
    file.create(file.path(feat_dirs[ii], "stats", "varcope1.nii.gz"))
    writeLines("set fmri(outputdir) \"dummy\"", feat_fsfs[ii])
  }
  
  fsl_setup <- data.frame(
    id = gpa$run_data$id,
    session = gpa$run_data$session,
    run_number = gpa$run_data$run_number,
    l1_model = rep("model1", n_rows),
    feat_fsf = feat_fsfs,
    feat_dir = feat_dirs,
    feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  
  result <- list(fsl = fsl_setup, spm = list(), afni = list())
  class(result) <- c("l1_setup", "list")
  return(result)
}

#' Create a mock l2_model_setup object
create_mock_l2_model_setup <- function(gpa) {
  subjects <- unique(gpa$subject_data[, c("id", "session")])
  n_rows <- nrow(subjects)
  feat_dirs <- replicate(n_rows, tempfile(fileext = ".gfeat"))
  feat_fsfs <- replicate(n_rows, tempfile(fileext = ".fsf"))
  
  fsl_setup <- data.frame(
    id = subjects$id,
    session = subjects$session,
    l1_model = rep("model1", n_rows),
    l1_cope_number = 1L,
    l1_cope_name = "event_vs_baseline",
    l2_model = rep("l2_model1", n_rows),
    l2_scope = "id_session",
    l2_input_mode = "cope_files",
    l2_passthrough = FALSE,
    n_l2_copes = 1L,
    n_input_files = 2L,
    passthrough_cope_file = NA_character_,
    feat_fsf = feat_fsfs,
    feat_dir = feat_dirs,
    feat_complete = TRUE,
    feat_dir_exists = TRUE,
    stringsAsFactors = FALSE
  )
  fsl_setup$cope_list <- lapply(seq_len(n_rows), function(ii) {
    data.frame(
      id = subjects$id[ii],
      session = subjects$session[ii],
      l2_cope_number = 1L,
      l2_cope_name = "intercept",
      stringsAsFactors = FALSE
    )
  })
  
  result <- list(fsl = fsl_setup, spm = list(), afni = list())
  class(result) <- c("l2_setup", "list")
  return(result)
}

# ==============================================================================
# Tests for setup_l1_models input validation
# ==============================================================================

test_that("setup_l1_models rejects invalid gpa class", {
  expect_error(
    setup_l1_models(list()),
    "glm_pipeline_arguments"
  )
})

test_that("setup_l1_models rejects invalid l1_model_names type", {
  gpa <- create_mock_gpa_for_l1_setup()
  
  expect_error(
    setup_l1_models(gpa, l1_model_names = 123),
    "character"
  )
})

test_that("setup_l1_models rejects unknown model names", {
  gpa <- create_mock_gpa_for_l1_setup()
  
  expect_error(
    setup_l1_models(gpa, l1_model_names = c("nonexistent_model")),
    "subset"
  )
})

test_that("setup_l1_models requires run_data to be a data.frame", {
  gpa <- create_mock_gpa_for_l1_setup()
  gpa$run_data <- NULL
  
  expect_error(
    setup_l1_models(gpa),
    "data.frame"
  )
})

test_that("setup_l1_models requires specific columns in run_data", {
  gpa <- create_mock_gpa_for_l1_setup()
  gpa$run_data$run_nifti <- NULL  # Remove required column
  
  expect_error(
    setup_l1_models(gpa),
    "run_nifti"
  )
})

test_that("setup_l1_models defaults to all models when l1_model_names is NULL", {
  gpa <- create_mock_gpa_for_l1_setup()
  # NULL l1_model_names defaults to all models (line 33 in setup_l1_models.R)
  # The function handles missing files gracefully, so we just verify it doesn't error on validation
  # It will return a gpa with updated l1_model_setup, even if empty
  result <- setup_l1_models(gpa, l1_model_names = NULL)
  expect_s3_class(result, "glm_pipeline_arguments")
})

test_that("setup_l1_models drops runs with missing event rows before building L1 designs", {
  gpa <- create_mock_gpa_for_l1_setup(n_subjects = 1, n_runs = 2)

  gpa$subject_data$mr_dir <- file.path(gpa$output_directory, gpa$subject_data$id)
  gpa$run_data$run_nifti <- "bold.nii.gz"
  dirs <- unique(gpa$run_data$mr_dir)
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  for (nii in file.path(gpa$run_data$mr_dir, gpa$run_data$run_nifti)) {
    RNifti::writeNifti(array(0, dim = c(2, 2, 2, 100)), nii)
  }

  gpa$glm_software <- character(0)
  gpa$level_backends <- list(l1 = character(0), l2 = character(0), l3 = character(0))
  gpa$additional <- list(bdm_args = list(plot = FALSE))
  gpa$use_preconvolve <- FALSE
  gpa$glm_settings <- list()
  gpa$output_locations$feat_l1_directory <- file.path(
    gpa$output_directory, "feat_l1", "sub-{id}", "ses-{session}", "{l1_model}"
  )
  gpa$l1_models$models$model1$signals <- "parametric_signal"
  gpa$l1_models$signals$parametric_signal$value <- subset(
    gpa$l1_models$signals$parametric_signal$value,
    run_number == 1
  )

  event_data <- subset(
    gpa$trial_data,
    run_number == 1,
    select = c("id", "session", "run_number", "trial", "onset", "duration", "event")
  )
  gpa$l1_models$events <- list(
    stimulus = list(name = "stimulus", data = event_data)
  )

  result <- setup_l1_models(gpa, l1_model_names = "model1")

  expect_s3_class(result, "glm_pipeline_arguments")
  expect_equal(result$l1_model_setup$metadata$run_number, 1)
})

# ==============================================================================
# Tests for setup_l2_models input validation
# ==============================================================================

test_that("setup_l2_models rejects invalid gpa class", {
  expect_error(
    setup_l2_models(list()),
    "glm_pipeline_arguments"
  )
})

test_that("setup_l2_models rejects missing l2_models", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l2_models <- NULL
  
  expect_error(
    setup_l2_models(gpa),
    "hi_model_set"
  )
})

test_that("setup_l2_models rejects missing l1_models", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_models <- NULL
  
  expect_error(
    setup_l2_models(gpa),
    "l1_model_set"
  )
})

test_that("setup_l2_models rejects invalid l1_model_names", {
  gpa <- create_mock_gpa_for_l2_setup()
  
  expect_error(
    setup_l2_models(gpa, l1_model_names = c("nonexistent")),
    "subset"
  )
})

test_that("setup_l2_models rejects invalid l2_model_names", {
  gpa <- create_mock_gpa_for_l2_setup()
  
  expect_error(
    setup_l2_models(gpa, l2_model_names = c("nonexistent")),
    "subset"
  )
})

test_that("setup_l2_models requires l1_model_setup", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_model_setup <- NULL
  
  expect_error(
    setup_l2_models(gpa),
    "l1_model_setup|setup_l1_models"
  )
})

test_that("setup_l2_models uses l2_scope to group L2 inputs", {
  gpa <- create_mock_gpa_for_l2_setup(n_subjects = 2, n_runs = 2)

  run_data_s2 <- gpa$run_data
  run_data_s2$session <- 2L
  gpa$run_data <- dplyr::bind_rows(gpa$run_data, run_data_s2) %>%
    dplyr::arrange(id, session, run_number)

  subject_data_s2 <- gpa$subject_data
  subject_data_s2$session <- 2L
  gpa$subject_data <- dplyr::bind_rows(gpa$subject_data, subject_data_s2) %>%
    dplyr::arrange(id, session)

  l1_setup_s2 <- gpa$l1_model_setup$fsl
  l1_setup_s2$session <- 2L
  gpa$l1_model_setup$fsl <- dplyr::bind_rows(gpa$l1_model_setup$fsl, l1_setup_s2) %>%
    dplyr::arrange(id, session, run_number)

  gpa$run_data$predictor <- seq_len(nrow(gpa$run_data)) / 10

  l2_session <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_session",
      model_formula = "~ predictor",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id_session"
    ),
    lg = lgr::get_logger("test")
  )

  l2_id <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_id",
      model_formula = "~ predictor",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id"
    ),
    lg = lgr::get_logger("test")
  )

  gpa$l2_models$models <- list(
    l2_session = l2_session,
    l2_id = l2_id
  )
  class(gpa$l2_models) <- c("hi_model_set", "list")
  gpa$l1_cope_names <- list(model1 = "cope1")

  fake_status <- function(feat_dir, feat_fsf, lg = NULL) {
    data.frame(
      feat_complete = TRUE,
      feat_failed = FALSE,
      stringsAsFactors = FALSE
    )
  }

  fake_l2_setup <- function(l1_df, l2_model, gpa) {
    data.frame(
      id = l1_df$id[1L],
      session = if (dplyr::n_distinct(l1_df$session) == 1L) l1_df$session[1L] else 0L,
      l1_model = l1_df$l1_model[1L],
      l1_cope_name = l1_df$l1_cope_name[1L],
      l2_model = l2_model,
      feat_fsf = tempfile(fileext = ".fsf"),
      feat_dir = tempfile(pattern = "feat_l2_"),
      feat_complete = TRUE,
      n_input_rows = nrow(l1_df),
      n_sessions_in_call = dplyr::n_distinct(l1_df$session),
      stringsAsFactors = FALSE
    )
  }

  gpa$glm_software <- "fsl"
  gpa$glm_backend_specs <- list(
    fsl = list(
      name = "fsl",
      runs_l1 = TRUE, runs_l2 = TRUE, runs_l3 = TRUE,
      multi_run_strategy = "explicit_l2",
      produced_artifacts = c("run_level_contrasts", "subject_session_contrasts", "group_level_stats"),
      l1_setup = NULL,
      l2_setup = fake_l2_setup,
      l3_setup = NULL,
      l1_status = fake_status,
      l2_status = fake_status,
      l3_status = NULL,
      l1_status_inputs = c("feat_dir", "feat_fsf"),
      l2_status_inputs = c("feat_dir", "feat_fsf"),
      l3_status_inputs = character(0),
      output_dir = NULL,
      l1_run = NULL,
      l2_run = NULL,
      l3_run = NULL
    )
  )

  result <- setup_l2_models(
    gpa,
    l1_model_names = "model1",
    l2_model_names = c("l2_session", "l2_id"),
    backend = "fsl"
  )

  expect_s3_class(result$l2_model_setup, "l2_setup")
  expect_true(is.data.frame(result$l2_model_setup$fsl))

  by_model <- split(result$l2_model_setup$fsl, result$l2_model_setup$fsl$l2_model)

  expect_equal(nrow(by_model$l2_session), 4L)
  expect_true(all(by_model$l2_session$n_sessions_in_call == 1L))

  expect_equal(nrow(by_model$l2_id), 2L)
  expect_true(all(by_model$l2_id$n_sessions_in_call == 2L))
  expect_true(all(by_model$l2_id$session == 0L))
})

test_that("build_l1_cope_validity classifies valid and invalid L1 copes", {
  gpa <- create_mock_gpa_for_l2_setup(n_subjects = 1, n_runs = 1)

  contrasts <- rbind(
    copeA = c(regA = 1, regB = 0),
    copeB = c(regA = 0, regB = 1),
    copeAB = c(regA = 1, regB = 1)
  )
  gpa$l1_models$models$model1$contrasts <- contrasts
  gpa$l1_cope_names <- list(model1 = rownames(contrasts))

  timing_dir <- tempfile("l1_validity_")
  dir.create(timing_dir, recursive = TRUE)
  zero_file <- file.path(timing_dir, "run_number1_regA.1D")
  nonzero_file <- file.path(timing_dir, "run_number1_regB.1D")
  writeLines(rep("0", 5), zero_file)
  writeLines(c("0", "1", "0.5", "0"), nonzero_file)
  writeLines(c("0\t1\tNA", "2\t1\tNA"), file.path(timing_dir, "run1_regA_FSL3col.txt"))

  fsf_file <- tempfile(fileext = ".fsf")
  writeLines(c(
    "set fmri(evtitle1) \"regA\"",
    "set fmri(shape1) 2",
    sprintf("set fmri(custom1) \"%s\"", zero_file),
    "set fmri(evtitle2) \"regB\"",
    "set fmri(shape2) 2",
    sprintf("set fmri(custom2) \"%s\"", nonzero_file)
  ), fsf_file)

  row <- data.frame(
    id = "sub1",
    session = 1L,
    run_number = 1L,
    l1_model = "model1",
    feat_fsf = fsf_file,
    feat_dir = file.path(timing_dir, "run1.feat"),
    feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  dir.create(file.path(row$feat_dir, "stats"), recursive = TRUE)
  for (cc in seq_len(nrow(contrasts))) {
    file.create(file.path(row$feat_dir, "stats", paste0("cope", cc, ".nii.gz")))
    file.create(file.path(row$feat_dir, "stats", paste0("varcope", cc, ".nii.gz")))
  }

  validity <- fmri.pipeline:::build_l1_cope_validity(row, gpa, lg = lgr::get_logger("test"))
  validity_by_cope <- split(validity, validity$l1_cope_name)

  expect_equal(validity_by_cope$copeA$validity_state, "invalid_unexpected_empty")
  expect_false(validity_by_cope$copeA$valid_for_higher_level)
  expect_equal(validity_by_cope$copeA$empty_evs, "regA")
  expect_identical(validity_by_cope$copeA$validity_reason, "all_missing_parametric_modulator")
  expect_equal(validity_by_cope$copeB$validity_state, "valid")
  expect_true(validity_by_cope$copeB$valid_for_higher_level)
  expect_equal(validity_by_cope$copeAB$validity_state, "invalid_unexpected_empty")
  expect_false(validity_by_cope$copeAB$valid_for_higher_level)
  expect_equal(validity_by_cope$copeAB$empty_evs, "regA")
  expect_identical(validity_by_cope$copeAB$validity_reason, "all_missing_parametric_modulator")

  expected_fsf <- tempfile(fileext = ".fsf")
  writeLines(c(
    "set fmri(evtitle1) \"regA\"",
    "set fmri(shape1) 10",
    "set fmri(evtitle2) \"regB\"",
    "set fmri(shape2) 2",
    sprintf("set fmri(custom2) \"%s\"", nonzero_file)
  ), expected_fsf)
  expected_row <- row
  expected_row$feat_fsf <- expected_fsf
  expected_validity <- fmri.pipeline:::build_l1_cope_validity(expected_row, gpa, lg = lgr::get_logger("test"))
  expect_equal(
    expected_validity$validity_state[expected_validity$l1_cope_name == "copeA"],
    "invalid_expected_empty"
  )

  missing_row <- row
  unlink(file.path(row$feat_dir, "stats", "varcope2.nii.gz"))
  missing_validity <- fmri.pipeline:::build_l1_cope_validity(missing_row, gpa, lg = lgr::get_logger("test"))
  expect_equal(
    missing_validity$validity_state[missing_validity$l1_cope_name == "copeB"],
    "invalid_missing_file"
  )
})

mark_mock_regA_empty <- function(gpa, empty_run_numbers) {
  for (ii in seq_len(nrow(gpa$l1_model_setup$fsl))) {
    rr <- gpa$l1_model_setup$fsl$run_number[ii]
    timing_dir <- tempfile(sprintf("l1_run%s_", rr))
    dir.create(timing_dir, recursive = TRUE)

    regA_file <- file.path(timing_dir, sprintf("run_number%s_regA.1D", rr))
    regB_file <- file.path(timing_dir, sprintf("run_number%s_regB.1D", rr))
    if (rr %in% empty_run_numbers) {
      writeLines(rep("0", 5), regA_file)
      writeLines(c("0\t1\tNA", "2\t1\tNA"), file.path(timing_dir, sprintf("run%s_regA_FSL3col.txt", rr)))
    } else {
      writeLines(c("0", "1", "0.5", "0"), regA_file)
    }
    writeLines(c("0", "1", "0.5", "0"), regB_file)

    fsf_file <- tempfile(fileext = ".fsf")
    writeLines(c(
      "set fmri(evtitle1) \"regA\"",
      "set fmri(shape1) 2",
      sprintf("set fmri(custom1) \"%s\"", regA_file),
      "set fmri(evtitle2) \"regB\"",
      "set fmri(shape2) 2",
      sprintf("set fmri(custom2) \"%s\"", regB_file)
    ), fsf_file)
    gpa$l1_model_setup$fsl$feat_fsf[ii] <- fsf_file
  }
  gpa
}

ensure_mock_l1_stats <- function(gpa, n_copes) {
  for (feat_dir in gpa$l1_model_setup$fsl$feat_dir) {
    dir.create(file.path(feat_dir, "stats"), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(feat_dir, ".feat_complete"))
    for (cc in seq_len(n_copes)) {
      file.create(file.path(feat_dir, "stats", paste0("cope", cc, ".nii.gz")))
      file.create(file.path(feat_dir, "stats", paste0("varcope", cc, ".nii.gz")))
    }
  }
  gpa
}

test_that("setup_l2_models creates per-cope L2 analyses and refits design after dropping invalid runs", {
  gpa <- create_mock_gpa_for_l2_setup(n_subjects = 1, n_runs = 3)
  gpa$run_data$predictor <- c(10.1, 20.2, 30.3)

  contrasts <- rbind(
    copeA = c(regA = 1, regB = 0),
    copeB = c(regA = 0, regB = 1)
  )
  gpa$l1_models$models$model1$contrasts <- contrasts
  gpa$l1_models$n_contrasts <- c(model1 = 2L)
  gpa$l1_cope_names <- list(model1 = rownames(contrasts))
  gpa <- ensure_mock_l1_stats(gpa, n_copes = 2L)

  l2_model <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_model1",
      model_formula = "~ predictor",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id_session"
    ),
    lg = lgr::get_logger("test")
  )
  gpa$l2_models$models <- list(l2_model1 = l2_model)
  class(gpa$l2_models) <- c("hi_model_set", "list")

  gpa <- mark_mock_regA_empty(gpa, empty_run_numbers = 2L)

  result <- setup_l2_models(
    gpa,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    backend = "fsl"
  )

  expect_s3_class(result$l2_model_setup, "l2_setup")
  expect_s3_class(result$l1_cope_validity, "l1_cope_validity")
  expect_true("l2_model" %in% names(result$l1_cope_validity$fsl))
  expect_true(all(result$l1_cope_validity$fsl$l2_model == "l2_model1"))

  l2_setup <- result$l2_model_setup$fsl
  expect_equal(nrow(l2_setup), 2L)
  expect_true(all(l2_setup$l2_input_mode == "cope_files"))

  row_copeA <- l2_setup[l2_setup$l1_cope_name == "copeA", , drop = FALSE]
  row_copeB <- l2_setup[l2_setup$l1_cope_name == "copeB", , drop = FALSE]

  expect_equal(row_copeA$n_input_files, 2L)
  expect_equal(row_copeB$n_input_files, 3L)
  expect_match(row_copeA$feat_dir, file.path("L1m-model1", "l1c-01_copeA", "L2m-l2_model1.gfeat"), fixed = TRUE)
  expect_match(row_copeB$feat_dir, file.path("L1m-model1", "l1c-02_copeB", "L2m-l2_model1.gfeat"), fixed = TRUE)
  expect_match(row_copeA$feat_fsf, file.path("L1m-model1", "l1c-01_copeA", "L2m-l2_model1.fsf"), fixed = TRUE)

  fsf_A <- readLines(row_copeA$feat_fsf, warn = FALSE)
  fsf_B <- readLines(row_copeB$feat_fsf, warn = FALSE)

  expect_true(any(grepl("set fmri\\(inputtype\\) 2$", fsf_A)))
  expect_false(any(grepl("set fmri\\(ncopeinputs\\)", fsf_A)))
  expect_equal(sum(grepl("^set feat_files\\(", fsf_A)), 2L)
  expect_equal(sum(grepl("^set feat_files\\(", fsf_B)), 3L)

  predictor_lines_A <- fsf_A[grepl("^set fmri\\(evg[0-9]+\\.2\\)", fsf_A)]
  predictor_lines_B <- fsf_B[grepl("^set fmri\\(evg[0-9]+\\.2\\)", fsf_B)]
  expect_match(paste(predictor_lines_A, collapse = "\n"), "10\\.1")
  expect_match(paste(predictor_lines_A, collapse = "\n"), "30\\.3")
  expect_false(any(grepl("20\\.2", predictor_lines_A)))
  expect_match(paste(predictor_lines_B, collapse = "\n"), "20\\.2")
})

test_that("setup_l2_models skips subject/cope when no runs survive validity filtering", {
  gpa <- create_mock_gpa_for_l2_setup(n_subjects = 1, n_runs = 3)
  gpa$run_data$predictor <- c(10.1, 20.2, 30.3)

  contrasts <- rbind(
    copeA = c(regA = 1, regB = 0),
    copeB = c(regA = 0, regB = 1)
  )
  gpa$l1_models$models$model1$contrasts <- contrasts
  gpa$l1_models$n_contrasts <- c(model1 = 2L)
  gpa$l1_cope_names <- list(model1 = rownames(contrasts))
  gpa <- ensure_mock_l1_stats(gpa, n_copes = 2L)

  l2_model <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_model1",
      model_formula = "~ predictor",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id_session"
    ),
    lg = lgr::get_logger("test")
  )
  gpa$l2_models$models <- list(l2_model1 = l2_model)
  class(gpa$l2_models) <- c("hi_model_set", "list")

  gpa <- mark_mock_regA_empty(gpa, empty_run_numbers = 1:3)

  result <- setup_l2_models(
    gpa,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    backend = "fsl"
  )

  l2_setup <- result$l2_model_setup$fsl
  expect_false("copeA" %in% l2_setup$l1_cope_name)
  expect_true("copeB" %in% l2_setup$l1_cope_name)
})

test_that("setup_l2_models passes a single surviving run directly through intercept-only L2", {
  gpa <- create_mock_gpa_for_l2_setup(n_subjects = 1, n_runs = 2)

  contrasts <- matrix(1, nrow = 1, ncol = 1, dimnames = list("copeA", "regA"))
  gpa$l1_models$models$model1$contrasts <- contrasts
  gpa$l1_models$n_contrasts <- c(model1 = 1L)
  gpa$l1_cope_names <- list(model1 = rownames(contrasts))

  l2_model <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_model1",
      model_formula = "~1",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id_session"
    ),
    lg = lgr::get_logger("test")
  )
  gpa$l2_models$models <- list(l2_model1 = l2_model)
  class(gpa$l2_models) <- c("hi_model_set", "list")

  root <- tempfile("l2_passthrough_")
  dir.create(root, recursive = TRUE)
  feat_dirs <- file.path(root, paste0("run", seq_len(2), ".feat"))
  feat_fsfs <- file.path(root, paste0("run", seq_len(2), ".fsf"))
  for (ii in seq_along(feat_dirs)) {
    dir.create(file.path(feat_dirs[ii], "stats"), recursive = TRUE)
    file.create(file.path(feat_dirs[ii], ".feat_complete"))
    file.create(file.path(feat_dirs[ii], "stats", "cope1.nii.gz"))
    file.create(file.path(feat_dirs[ii], "stats", "varcope1.nii.gz"))
    writeLines("set fmri(outputdir) \"dummy\"", feat_fsfs[ii])
  }
  gpa$l1_model_setup$fsl$feat_dir <- feat_dirs
  gpa$l1_model_setup$fsl$feat_fsf <- feat_fsfs
  gpa <- mark_mock_regA_empty(gpa, empty_run_numbers = 2L)
  gpa$l1_model_setup$fsl$feat_dir <- feat_dirs

  result <- setup_l2_models(
    gpa,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    backend = "fsl"
  )

  l2_setup <- result$l2_model_setup$fsl
  expect_equal(nrow(l2_setup), 1L)
  expect_true(l2_setup$l2_passthrough[1])
  expect_equal(l2_setup$l2_input_mode[1], "l1_cope_file_passthrough")
  expect_equal(l2_setup$n_input_files[1], 1L)
  expect_false(l2_setup$to_run[1])
  expect_equal(l2_setup$passthrough_cope_file[1], file.path(feat_dirs[1], "stats", "cope1.nii.gz"))

  cfg <- data.frame(
    id = l2_setup$id[1],
    session = l2_setup$session[1],
    l1_model = "model1",
    l2_model = "l2_model1",
    l3_model = "l3_model1",
    l3_input_mode = "per_session",
    l1_cope_number = 1L,
    l2_cope_number = 1L,
    l1_cope_name = "copeA",
    l2_cope_name = l2_setup$cope_list[[1]]$l2_cope_name[1],
    stringsAsFactors = FALSE
  )
  inputs <- fmri.pipeline:::get_feat_l3_inputs(result, cfg, lg = lgr::get_logger("test"))
  expect_equal(inputs[[1]]$cope_file[1], l2_setup$passthrough_cope_file[1])
})

test_that("setup_l2_models skips single surviving run for non-pass-through L2 design", {
  gpa <- create_mock_gpa_for_l2_setup(n_subjects = 1, n_runs = 2)
  gpa$run_data$predictor <- c(10.1, 20.2)

  contrasts <- matrix(1, nrow = 1, ncol = 1, dimnames = list("copeA", "regA"))
  gpa$l1_models$models$model1$contrasts <- contrasts
  gpa$l1_models$n_contrasts <- c(model1 = 1L)
  gpa$l1_cope_names <- list(model1 = rownames(contrasts))

  l2_model <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_model1",
      model_formula = "~ predictor",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id_session"
    ),
    lg = lgr::get_logger("test")
  )
  gpa$l2_models$models <- list(l2_model1 = l2_model)
  class(gpa$l2_models) <- c("hi_model_set", "list")

  gpa <- mark_mock_regA_empty(gpa, empty_run_numbers = 2L)

  result <- setup_l2_models(
    gpa,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    backend = "fsl"
  )

  l2_setup <- result$l2_model_setup$fsl
  expect_true(is.null(l2_setup) || nrow(l2_setup) == 0L)
  expect_s3_class(result$l1_cope_validity, "l1_cope_validity")
  expect_equal(sum(result$l1_cope_validity$fsl$valid_for_higher_level), 1L)
})

test_that("setup_l2_models output feeds setup_l3_models with per-cope L2 inputs", {
  gpa <- create_mock_gpa_for_l2_setup(n_subjects = 4, n_runs = 3)
  gpa$run_data$predictor <- gpa$run_data$run_number * 10 + seq_len(nrow(gpa$run_data)) / 100

  contrasts <- rbind(
    copeA = c(regA = 1, regB = 0),
    copeB = c(regA = 0, regB = 1)
  )
  gpa$l1_models$models$model1$contrasts <- contrasts
  gpa$l1_models$n_contrasts <- c(model1 = 2L)
  gpa$l1_cope_names <- list(model1 = rownames(contrasts))
  gpa <- ensure_mock_l1_stats(gpa, n_copes = 2L)

  l2_model <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_model1",
      model_formula = "~ predictor",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id_session"
    ),
    lg = lgr::get_logger("test")
  )
  gpa$l2_models$models <- list(l2_model1 = l2_model)
  class(gpa$l2_models) <- c("hi_model_set", "list")

  l3_model <- fmri.pipeline:::create_new_hi_model(
    data = gpa$subject_data,
    level = 3L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l3_model1",
      model_formula = "~1",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l3_input_mode = "per_session",
      fsl_outlier_deweighting = FALSE
    ),
    lg = lgr::get_logger("test")
  )
  gpa$l3_models <- structure(list(models = list(l3_model1 = l3_model)), class = c("hi_model_set", "list"))
  gpa$output_locations$setup_l3_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l3_log_txt <- tempfile(fileext = ".txt")
  gpa$output_locations$feat_l3_directory <- file.path(
    gpa$output_directory, "feat_l3", "L3m-{l3_model}",
    "L1m-{l1_model}", "L2m-{l2_model}", "l2c-{l2_contrast}"
  )
  gpa$output_locations$feat_l3_fsf <- "FEAT_l1c-{l1_contrast}.fsf"
  gpa$additional$feat_l3_args <- list()
  gpa$glm_settings$fsl$force_l3_creation <- TRUE

  gpa <- mark_mock_regA_empty(gpa, empty_run_numbers = 2L)

  l2_result <- setup_l2_models(
    gpa,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    backend = "fsl"
  )

  l2_setup <- l2_result$l2_model_setup$fsl
  expect_equal(nrow(l2_setup), 8L)
  expect_true(all(l2_setup$n_input_files[l2_setup$l1_cope_name == "copeA"] == 2L))
  expect_true(all(l2_setup$n_input_files[l2_setup$l1_cope_name == "copeB"] == 3L))

  for (ii in seq_len(nrow(l2_setup))) {
    dir.create(file.path(l2_setup$feat_dir[ii], "cope1.feat", "stats"), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(l2_setup$feat_dir[ii], ".feat_complete"))
    for (cc in seq_len(l2_setup$n_l2_copes[ii])) {
      file.create(file.path(l2_setup$feat_dir[ii], "cope1.feat", "stats", paste0("cope", cc, ".nii.gz")))
    }
  }

  l3_result <- setup_l3_models(
    l2_result,
    l1_model_names = "model1",
    l2_model_names = "l2_model1",
    l3_model_names = "l3_model1",
    backend = "fsl"
  )

  expect_s3_class(l3_result$l3_model_setup, "l3_setup")
  l3_setup <- l3_result$l3_model_setup$fsl
  expect_equal(nrow(l3_setup), 4L)
  expect_setequal(l3_setup$l1_cope_name, c("copeA", "copeB"))
  expect_true(all(file.exists(l3_setup$feat_fsf)))

  fsfs_by_cope <- split(l3_setup$feat_fsf, l3_setup$l1_cope_name)
  feat_inputs <- lapply(fsfs_by_cope, function(fsfs) {
    unlist(lapply(fsfs, function(fsf) {
      lines <- readLines(fsf, warn = FALSE)
      sub("^set feat_files\\([0-9]+\\) \"([^\"]+)\"$", "\\1", grep("^set feat_files\\(", lines, value = TRUE))
    }), use.names = FALSE)
  })

  expect_equal(length(feat_inputs$copeA), 8L)
  expect_equal(length(feat_inputs$copeB), 8L)
  expect_true(all(grepl(file.path("l1c-01_copeA", "L2m-l2_model1.gfeat", "cope1.feat", "stats"), feat_inputs$copeA, fixed = TRUE)))
  expect_true(all(grepl(file.path("l1c-02_copeB", "L2m-l2_model1.gfeat", "cope1.feat", "stats"), feat_inputs$copeB, fixed = TRUE)))
  expect_true(any(grepl("cope2.nii.gz", feat_inputs$copeA, fixed = TRUE)))
  expect_true(any(grepl("cope2.nii.gz", feat_inputs$copeB, fixed = TRUE)))
  expect_false(any(grepl("l1c-02_copeB", feat_inputs$copeA, fixed = TRUE)))
  expect_false(any(grepl("l1c-01_copeA", feat_inputs$copeB, fixed = TRUE)))
})

# ==============================================================================
# Tests for setup_l3_models input validation
# ==============================================================================

test_that("setup_l3_models rejects invalid gpa class", {
  expect_error(
    setup_l3_models(list()),
    "glm_pipeline_arguments"
  )
})

test_that("setup_l3_models rejects missing l3_models", {
  gpa <- create_mock_gpa_for_l3_setup()
  gpa$l3_models <- NULL
  
  expect_error(
    setup_l3_models(gpa),
    "hi_model_set"
  )
})

test_that("setup_l3_models rejects missing l1_models", {
  gpa <- create_mock_gpa_for_l3_setup()
  gpa$l1_models <- NULL
  
  expect_error(
    setup_l3_models(gpa),
    "l1_model_set"
  )
})

test_that("setup_l3_models handles invalid l3_model_names in downstream processing", {
  # Note: setup_l3_models does not directly validate l3_model_names as a subset
  # It errors when trying to process unknown models in the actual setup loop
  # The validation happens via choose_glm_models in run_glm_pipeline
  gpa <- create_mock_gpa_for_l3_setup()
  # Ensure l2 models are marked complete to get past that check
  gpa$l2_model_setup$fsl$feat_complete <- TRUE
  gpa$l1_model_setup$fsl$feat_complete <- TRUE
  
  # With a nonexistent model name, it should error somewhere in processing
  # (either on subset check or when the model can't be found)
  expect_error(
    setup_l3_models(gpa, l3_model_names = c("nonexistent"))
  )
})

test_that("setup_l3_models requires l2_models for multi_run", {
  gpa <- create_mock_gpa_for_l3_setup(multi_run = TRUE)
  gpa$l2_models <- NULL
  
  expect_error(
    setup_l3_models(gpa),
    "hi_model_set"
  )
})

test_that("setup_l3_models reports signature mismatch details when no pairs remain", {
  gpa <- create_mock_gpa_for_l3_setup(multi_run = TRUE)
  gpa$l2_models$models$l2_model1$l2_scope <- "id"
  gpa$l3_models$models$l3_model1$l3_input_mode <- "per_session"

  expect_warning(
    result <- setup_l3_models(gpa),
    "No L2\\+L3 model combinations remain"
  )

  expect_false(result$l3_setup_status$success)
  expect_match(result$l3_setup_status$reason, "requires l2_scope='id_session'")
})

test_that("get_feat_l3_inputs splits by session only for per_session mode", {
  root <- tempfile("l3_split_mode_")
  dir.create(root, recursive = TRUE)

  l2_setup <- data.frame(
    id = rep(c("sub1", "sub2"), each = 2),
    session = rep(c(1L, 2L), times = 2),
    l1_model = "model1",
    l1_cope_number = 1L,
    l1_cope_name = "cope1",
    l2_model = "l2_model1",
    l2_scope = "id_session",
    l2_input_mode = "cope_files",
    l2_passthrough = FALSE,
    n_l2_copes = 1L,
    n_input_files = 2L,
    passthrough_cope_file = NA_character_,
    feat_dir = file.path(root, paste0("feat_l2_", seq_len(4))),
    feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  l2_setup$cope_list <- lapply(seq_len(nrow(l2_setup)), function(ii) {
    data.frame(
      id = l2_setup$id[ii],
      session = l2_setup$session[ii],
      l2_cope_number = 1L,
      l2_cope_name = "cope1",
      stringsAsFactors = FALSE
    )
  })

  for (ii in seq_len(nrow(l2_setup))) {
    cope_stats <- file.path(l2_setup$feat_dir[ii], "cope1.feat", "stats")
    dir.create(cope_stats, recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(cope_stats, "cope1.nii.gz"))
  }

  gpa <- list(
    multi_run = TRUE,
    l2_model_setup = structure(list(fsl = l2_setup), class = c("l2_setup", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  cfg_base <- data.frame(
    id = rep(c("sub1", "sub2"), each = 2),
    session = rep(c(1L, 2L), times = 2),
    l1_model = "model1",
    l2_model = "l2_model1",
    l3_model = "l3_model1",
    l1_cope_number = 1L,
    l2_cope_number = 1L,
    l1_cope_name = "cope1",
    l2_cope_name = "cope1",
    stringsAsFactors = FALSE
  )

  cfg_sep <- cfg_base
  cfg_sep$l3_input_mode <- "per_session"
  sep_inputs <- fmri.pipeline:::get_feat_l3_inputs(gpa, cfg_sep, lg = lgr::get_logger("test"))
  expect_equal(length(sep_inputs), 2L)
  expect_true(all(vapply(sep_inputs, function(df) length(unique(df$session)), integer(1)) == 1L))

  cfg_pool <- cfg_base
  cfg_pool$l3_input_mode <- "pooled_sessions_subject_ev"
  pooled_inputs <- fmri.pipeline:::get_feat_l3_inputs(gpa, cfg_pool, lg = lgr::get_logger("test"))
  expect_equal(length(pooled_inputs), 1L)
  expect_equal(nrow(pooled_inputs[[1L]]), nrow(cfg_pool))
})

test_that("get_feat_l3_inputs resolves cope-files L2 outputs from cope1.feat", {
  root <- tempfile("l3_cope_files_")
  dir.create(root, recursive = TRUE)

  l2_setup <- data.frame(
    id = "sub1",
    session = 1L,
    l1_model = "model1",
    l1_cope_number = 2L,
    l1_cope_name = "copeB",
    l2_model = "l2_model1",
    l2_scope = "id_session",
    l2_input_mode = "cope_files",
    l2_passthrough = FALSE,
    n_l2_copes = 1L,
    n_input_files = 2L,
    passthrough_cope_file = NA_character_,
    feat_dir = file.path(root, "feat_l2_single_cope"),
    feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  l2_setup$cope_list <- list(data.frame(
    id = "sub1",
    session = 1L,
    l2_cope_number = 1L,
    l2_cope_name = "cope1",
    stringsAsFactors = FALSE
  ))

  cope_stats <- file.path(l2_setup$feat_dir, "cope1.feat", "stats")
  dir.create(cope_stats, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(cope_stats, "cope1.nii.gz"))

  gpa <- list(
    multi_run = TRUE,
    l2_model_setup = structure(list(fsl = l2_setup), class = c("l2_setup", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  cfg <- data.frame(
    id = "sub1",
    session = 1L,
    l1_model = "model1",
    l2_model = "l2_model1",
    l3_model = "l3_model1",
    l3_input_mode = "per_session",
    l1_cope_number = 2L,
    l2_cope_number = 1L,
    l1_cope_name = "copeB",
    l2_cope_name = "cope1",
    stringsAsFactors = FALSE
  )

  inputs <- fmri.pipeline:::get_feat_l3_inputs(gpa, cfg, lg = lgr::get_logger("test"))
  expect_equal(length(inputs), 1L)
  expect_match(
    inputs[[1L]]$cope_file[1],
    file.path("cope1.feat", "stats", "cope1.nii.gz"),
    fixed = TRUE
  )
})

test_that("get_spm_l3_inputs splits by session only for per_session mode", {
  root <- tempfile("spm_l3_split_mode_")
  dir.create(root, recursive = TRUE)

  spm_setup <- data.frame(
    id = rep(c("sub1", "sub2"), each = 2),
    session = rep(c(1L, 2L), times = 2),
    l1_model = "model1",
    spm_dir = file.path(root, paste0("spm_l1_", seq_len(4))),
    spm_complete = TRUE,
    spm_contrast_exists = TRUE,
    stringsAsFactors = FALSE
  )
  for (ii in seq_len(nrow(spm_setup))) {
    dir.create(spm_setup$spm_dir[ii], recursive = TRUE, showWarnings = FALSE)
    make_spm_completed_fixture(spm_setup$spm_dir[ii], numbers = 3L)
  }

  gpa <- list(
    l1_model_setup = structure(list(spm = spm_setup), class = c("l1_setup", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  cfg_base <- data.frame(
    id = rep(c("sub1", "sub2"), each = 2),
    session = rep(c(1L, 2L), times = 2),
    l1_model = "model1",
    l3_model = "l3_model1",
    l1_cope_number = 1L,
    l1_cope_name = "cope1",
    stringsAsFactors = FALSE
  )

  cfg_sep <- cfg_base
  cfg_sep$l3_input_mode <- "per_session"
  sep_inputs <- fmri.pipeline:::get_spm_l3_inputs(gpa, cfg_sep, lg = lgr::get_logger("test"))
  expect_equal(length(sep_inputs), 2L)
  expect_true(all(vapply(sep_inputs, function(df) length(unique(df$session)), integer(1)) == 1L))

  cfg_pool <- cfg_base
  cfg_pool$l3_input_mode <- "pooled_sessions_subject_ev"
  pooled_inputs <- fmri.pipeline:::get_spm_l3_inputs(gpa, cfg_pool, lg = lgr::get_logger("test"))
  expect_equal(length(pooled_inputs), 1L)
  expect_equal(nrow(pooled_inputs[[1L]]), nrow(cfg_pool))
  expect_true(all(basename(pooled_inputs[[1L]]$con_file) == "con_0003.nii"))
})

# ==============================================================================
# Tests for enforce_glms_complete helper function
# ==============================================================================

test_that("enforce_glms_complete rejects invalid gpa", {
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(list(), level = 1L, lg = lg),
    "glm_pipeline_arguments"
  )
})

test_that("enforce_glms_complete rejects invalid level", {
  gpa <- create_mock_gpa()
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 0L, lg = lg),
    "level"
  )
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 4L, lg = lg),
    "level"
  )
})

test_that("enforce_glms_complete errors when l1_model_setup is missing", {
  gpa <- create_mock_gpa()
  gpa$l1_model_setup <- NULL
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg),
    "l1_model_setup"
  )
})

test_that("enforce_glms_complete errors when all feat runs incomplete", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_model_setup$fsl$feat_complete <- FALSE
  gpa$glm_software <- "fsl"
  lg <- lgr::get_logger("test")
  
  expect_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg),
    "incomplete"
  )
})

test_that("enforce_glms_complete warns when some feat runs incomplete", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$l1_model_setup$fsl$feat_complete[1] <- FALSE
  gpa$glm_software <- "fsl"
  lg <- lgr::get_logger("test")
  lg$set_threshold("warn")
  
  # Should not error, just warn
  expect_no_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg)
  )
})

test_that("enforce_glms_complete passes with all complete runs", {
  gpa <- create_mock_gpa_for_l2_setup()
  gpa$glm_software <- "fsl"
  lg <- lgr::get_logger("test")
  
  expect_no_error(
    fmri.pipeline:::enforce_glms_complete(gpa, level = 1L, lg = lg)
  )
})

# ==============================================================================
# Tests for respecify_l2_models_by_subject helper function
# ==============================================================================

test_that("respecify_l2_models_by_subject rejects invalid model object", {
  data <- data.frame(id = c("s1", "s2"), session = c(1, 1), value = c(1, 2))
  
  expect_error(
    fmri.pipeline:::respecify_l2_models_by_subject(list(), data),
    "l1_model_spec|hi_model_spec"
  )
})

test_that("respecify_l2_models_by_subject supports id and id_session grouping", {
  data <- data.frame(
    id = rep(c("s1", "s2"), each = 4),
    session = rep(rep(1:2, each = 2), times = 2),
    run_number = rep(1:2, times = 4),
    predictor = seq(0.1, 0.8, by = 0.1)
  )

  mobj <- fmri.pipeline:::create_new_hi_model(
    data = data,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "respecify_scope_test",
      model_formula = "~ predictor",
      diagonal = TRUE,
      num2fac = character(0),
      covariate_transform = c(),
      reference_level = c(),
      l2_scope = "id_session"
    ),
    lg = lgr::get_logger("test")
  )

  by_id_session <- fmri.pipeline:::respecify_l2_models_by_subject(
    mobj,
    data,
    split_on = c("id", "session")
  )
  expect_equal(nrow(by_id_session$by_subject), 4L)
  expect_true(all(by_id_session$by_subject$grouping_scope == "id_session"))

  by_id <- fmri.pipeline:::respecify_l2_models_by_subject(
    mobj,
    data,
    split_on = "id",
    aggregated_session = 0L
  )
  expect_equal(nrow(by_id$by_subject), 2L)
  expect_true(all(by_id$by_subject$session == 0L))
  expect_true(all(by_id$by_subject$grouping_scope == "id"))
})

test_that("respecify_l2_models_by_subject requires id and session columns", {
  # Create minimal model object
  mobj <- list(
    lmfit = lm(value ~ 1, data = data.frame(value = 1:3)),
    contrasts = list(intercept = c(intercept = 1))
  )
  class(mobj) <- c("hi_model_spec", "list")
  
  data <- data.frame(x = 1:3)
  
  expect_error(
    fmri.pipeline:::respecify_l2_models_by_subject(mobj, data),
    "id|session"
  )
})

# ==============================================================================
# Tests for choose_glm_models helper function
# ==============================================================================

test_that("choose_glm_models returns all models when 'all' specified", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  result <- fmri.pipeline:::choose_glm_models(gpa, "all", level = 1)
  
  expect_equal(result, names(gpa$l1_models$models))
})

test_that("choose_glm_models returns NULL when 'none' specified", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  result <- fmri.pipeline:::choose_glm_models(gpa, "none", level = 1)
  
  expect_null(result)
})

test_that("choose_glm_models returns specified model names", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  model_name <- names(gpa$l1_models$models)[1]
  
  result <- fmri.pipeline:::choose_glm_models(gpa, model_name, level = 1)
  
  expect_equal(result, model_name)
})

test_that("choose_glm_models rejects invalid model names", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  expect_error(
    fmri.pipeline:::choose_glm_models(gpa, "nonexistent_model", level = 1),
    "subset"
  )
})

test_that("choose_glm_models validates level argument", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  # Level 0 is invalid - should error on checkmate validation
  expect_error(
    fmri.pipeline:::choose_glm_models(gpa, "all", level = 0),
    "Element 1 is not >= 1"
  )
  
  # Level 4 is invalid - should error on checkmate validation
  expect_error(
    fmri.pipeline:::choose_glm_models(gpa, "all", level = 4),
    "Element 1 is not <= 3"
  )
})

test_that("choose_glm_models handles missing model object gracefully", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  # Remove l2_models
  gpa$l2_models <- NULL
  
  # When the model object doesn't exist, the function should handle it
  # (either error or return NULL)
  result <- tryCatch(
    fmri.pipeline:::choose_glm_models(gpa, "all", level = 2),
    error = function(e) "error"
  )
  
  # Either returns NULL (no models) or errors
  expect_true(is.null(result) || identical(result, "error"))
})

# ==============================================================================
# Tests for model setup data structures
# ==============================================================================

test_that("l1_setup object has correct structure", {
  gpa <- create_mock_gpa_for_l2_setup()
  setup <- gpa$l1_model_setup
  
  expect_s3_class(setup, "l1_setup")
  expect_true("fsl" %in% names(setup))
  expect_s3_class(setup$fsl, "data.frame")
})

test_that("l1_setup$fsl has required columns", {
  gpa <- create_mock_gpa_for_l2_setup()
  setup <- gpa$l1_model_setup
  
  required_cols <- c("id", "session", "run_number", "l1_model", "feat_fsf", "feat_dir", "feat_complete")
  expect_true(all(required_cols %in% names(setup$fsl)))
})

test_that("l2_setup object has correct structure", {
  gpa <- create_mock_gpa_for_l3_setup()
  setup <- gpa$l2_model_setup
  
  expect_s3_class(setup, "l2_setup")
  expect_true("fsl" %in% names(setup))
  expect_s3_class(setup$fsl, "data.frame")
})

test_that("l2_setup$fsl has required columns", {
  gpa <- create_mock_gpa_for_l3_setup()
  setup <- gpa$l2_model_setup
  
  required_cols <- c(
    "id", "session", "l1_model", "l1_cope_number", "l1_cope_name",
    "l2_model", "l2_scope", "l2_input_mode", "l2_passthrough",
    "n_l2_copes", "n_input_files", "passthrough_cope_file",
    "feat_fsf", "feat_dir", "feat_complete", "cope_list"
  )
  expect_true(all(required_cols %in% names(setup$fsl)))
})

# ==============================================================================
# Tests for run/subject exclusion logic
# ==============================================================================

test_that("run_data exclusion flags are properly structured", {
  gpa <- create_mock_gpa_for_l2_setup()
  
  expect_true("exclude_run" %in% names(gpa$run_data))
  expect_true("exclude_subject" %in% names(gpa$run_data))
  expect_type(gpa$run_data$exclude_run, "logical")
  expect_type(gpa$run_data$exclude_subject, "logical")
})

test_that("subject_data exclusion flags are properly structured", {
  gpa <- create_mock_gpa_for_l3_setup()
  
  expect_true("exclude_subject" %in% names(gpa$subject_data))
  expect_type(gpa$subject_data$exclude_subject, "logical")
})

# ==============================================================================
# Tests for multi_run vs single_run configuration
# ==============================================================================

test_that("multi_run gpa requires l2_models", {
  gpa <- create_mock_gpa_for_l3_setup(multi_run = TRUE)
  
  expect_true(gpa$multi_run)
  expect_s3_class(gpa$l2_models, "hi_model_set")
})

test_that("single_run gpa can work without l2_models", {
  gpa <- create_mock_gpa_for_l3_setup(multi_run = FALSE)
  gpa$l2_models <- NULL  # Remove l2_models
  
  expect_false(gpa$multi_run)
  expect_null(gpa$l2_models)
})

# ==============================================================================
# Tests for model name defaulting behavior
# ==============================================================================

test_that("NULL l1_model_names defaults to all models", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  all_models <- names(gpa$l1_models$models)
  
  # The setup functions should use all models when NULL is passed
  # We test this indirectly through the choose_glm_models helper
  result <- fmri.pipeline:::choose_glm_models(gpa, "all", level = 1)
  expect_equal(result, all_models)
})

test_that("empty model names vector is handled", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  
  result <- fmri.pipeline:::choose_glm_models(gpa, NULL, level = 1)
  expect_null(result)
})
