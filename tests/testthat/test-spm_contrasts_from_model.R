run_spm_contrast_script <- function(script, arguments) {
  system2(
    rscript_executable(),
    c("--no-save", "--no-restore", script, arguments),
    stdout = TRUE,
    stderr = TRUE
  )
}

test_that("generate_spm_contrasts_from_model maps L1 contrasts to SPM columns", {
  skip_if_not_installed("R.matlab")

  tmp_dir <- tempfile("spm_contrast_model_")
  dir.create(tmp_dir, recursive = TRUE)

  # Fake SPM design columns (two runs, two conditions)
  mnames <- c(
    "Sn(1) feedback.A*bf(1)",
    "Sn(1) feedback.B*bf(1)",
    "Sn(2) feedback.A*bf(1)",
    "Sn(2) feedback.B*bf(1)",
    "Sn(1) constant",
    "Sn(2) constant"
  )
  cpos <- 1:4
  bpos <- 5:6
  npos <- integer(0)

  matfile <- file.path(tmp_dir, "design_columns.mat")
  R.matlab::writeMat(matfile, mnames = mnames, cpos = cpos, bpos = bpos, npos = npos)

  # Minimal l1_model_spec with contrasts
  cmat <- matrix(c(1, -1), nrow = 1, dimnames = list("A_gt_B", c("feedback.A", "feedback.B")))
  mobj <- list(contrasts = cmat)
  class(mobj) <- c("list", "l1_model_spec")

  # Generate spec + commands
  cmds <- fmri.pipeline:::generate_spm_contrasts_from_model(
    output_dir = tmp_dir,
    mobj = mobj,
    spm_path = tmp_dir,
    execute = FALSE,
    average_across_runs = TRUE
  )

  spec_path <- file.path(tmp_dir, "spm_contrast_spec.rds")
  expect_true(file.exists(spec_path))
  expect_true(file.exists(file.path(tmp_dir, "extract_design_columns.m")))

  # Run contrast script directly with the prewritten .mat file
  setup_script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package = "fmri.pipeline")
  expect_true(file.exists(setup_script))
  status <- run_spm_contrast_script(
    setup_script,
    c(
      "-mat_file", matfile,
      "-contrast_rds", spec_path,
      "-average_across_runs", "TRUE",
      "-spm_path", tmp_dir
    )
  )
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))

  mfile <- file.path(tmp_dir, "estimate_glm_contrasts.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  convec_line <- lines[grepl("tcon.convec", lines)][1]
  expect_true(nchar(convec_line) > 0)

  # Extract numeric weights from convec line
  weights <- as.numeric(strsplit(gsub(".*\\[|\\].*", "", convec_line), ",")[[1]])
  expect_equal(length(weights), length(mnames))
  # Average across runs => 0.5 for feedback.A in each run, -0.5 for feedback.B in each run
  expect_equal(weights[1:4], c(0.5, -0.5, 0.5, -0.5), tolerance = 1e-8)
  expect_equal(weights[5:6], c(0, 0))
})

test_that("generate_spm_contrasts_from_model errors on missing contrasts", {
  mobj <- list(contrasts = NULL)
  class(mobj) <- c("list", "l1_model_spec")
  expect_error(
    fmri.pipeline:::generate_spm_contrasts_from_model(output_dir = tempdir(), mobj = mobj),
    "missing or empty"
  )
})

test_that("generate_spm_contrasts_from_model matches compact SPM pmod labels", {
  skip_if_not_installed("R.matlab")

  tmp_dir <- tempfile("spm_contrast_pmod_")
  dir.create(tmp_dir, recursive = TRUE)

  # Compact pmod form from SPM (no spaces around x)
  mnames <- c(
    "Sn(1) choice*bf(1)",
    "Sn(1) choicexchoice_trial-pmod^1*bf(1)",
    "Sn(1) feedback*bf(1)",
    "Sn(1) constant"
  )
  cpos <- 1:3
  bpos <- 4
  npos <- integer(0)

  matfile <- file.path(tmp_dir, "design_columns.mat")
  R.matlab::writeMat(matfile, mnames = mnames, cpos = cpos, bpos = bpos, npos = npos)

  cmat <- matrix(
    c(0, 0, 1),
    nrow = 1,
    dimnames = list("EV_choice_trial", c("choice", "feedback", "choice_trial"))
  )
  mobj <- list(contrasts = cmat)
  class(mobj) <- c("list", "l1_model_spec")

  fmri.pipeline:::generate_spm_contrasts_from_model(
    output_dir = tmp_dir,
    mobj = mobj,
    spm_path = tmp_dir,
    execute = FALSE,
    average_across_runs = TRUE
  )

  spec_path <- file.path(tmp_dir, "spm_contrast_spec.rds")
  setup_script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package = "fmri.pipeline")
  status <- run_spm_contrast_script(
    setup_script,
    c(
      "-mat_file", matfile,
      "-contrast_rds", spec_path,
      "-average_across_runs", "TRUE",
      "-spm_path", tmp_dir
    )
  )
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))

  mfile <- file.path(tmp_dir, "estimate_glm_contrasts.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_true(any(grepl("EV_choice_trial", lines, fixed = TRUE)))

  convec_line <- lines[grepl("tcon.convec", lines)][1]
  weights <- as.numeric(strsplit(gsub(".*\\[|\\].*", "", convec_line), ",")[[1]])
  expect_equal(weights, c(0, 1, 0, 0), tolerance = 1e-8)
})

test_that("generate_spm_contrasts_from_model adds session-specific and difference interaction contrasts", {
  skip_if_not_installed("R.matlab")

  tmp_dir <- tempfile("spm_contrast_proj_modes_")
  dir.create(tmp_dir, recursive = TRUE)

  mnames <- c(
    "Sn(1) face_onset*bf(1)",
    "Sn(1) face_onset x face_onset_x_l2_face_emotion-pmod^1*bf(1)",
    "Sn(2) face_onset*bf(1)",
    "Sn(2) face_onset x face_onset_x_l2_face_emotion-pmod^1*bf(1)",
    "Sn(1) constant",
    "Sn(2) constant"
  )
  cpos <- 1:4
  bpos <- 5:6
  npos <- integer(0)

  matfile <- file.path(tmp_dir, "design_columns.mat")
  R.matlab::writeMat(matfile, mnames = mnames, cpos = cpos, bpos = bpos, npos = npos)

  cmat <- matrix(
    c(1, 0),
    nrow = 1,
    dimnames = list("face_main", c("face_onset", "face_onset_x_l2_face_emotion"))
  )
  mobj <- list(contrasts = cmat)
  class(mobj) <- c("list", "l1_model_spec")

  fmri.pipeline:::generate_spm_contrasts_from_model(
    output_dir = tmp_dir,
    mobj = mobj,
    spm_path = tmp_dir,
    execute = FALSE,
    average_across_runs = TRUE,
    projection_interaction_terms = "face_onset_x_l2_face_emotion",
    projection_interaction_contrast_modes = c("session_specific", "session_differences"),
    projection_interaction_run_labels = c(1L, 2L)
  )

  spec_path <- file.path(tmp_dir, "spm_contrast_spec.rds")
  setup_script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package = "fmri.pipeline")
  status <- run_spm_contrast_script(
    setup_script,
    c(
      "-mat_file", matfile,
      "-contrast_rds", spec_path,
      "-average_across_runs", "TRUE",
      "-spm_path", tmp_dir
    )
  )
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))

  mfile <- file.path(tmp_dir, "estimate_glm_contrasts.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  expect_true(any(grepl("proj_int_face_onset_x_l2_face_emotion_run1", lines, fixed = TRUE)))
  expect_true(any(grepl("proj_int_face_onset_x_l2_face_emotion_run2", lines, fixed = TRUE)))
  expect_true(any(grepl("proj_int_face_onset_x_l2_face_emotion_run2_minus_run1", lines, fixed = TRUE)))

  convec_lines <- lines[grepl("tcon.convec", lines)]
  vectors <- lapply(convec_lines, function(line) {
    as.numeric(strsplit(gsub(".*\\[|\\].*", "", line), ",")[[1]])
  })

  expect_true(any(vapply(vectors, function(v) isTRUE(all.equal(v, c(0, 1, 0, 0, 0, 0), tolerance = 1e-8)), logical(1))))
  expect_true(any(vapply(vectors, function(v) isTRUE(all.equal(v, c(0, 0, 0, 1, 0, 0), tolerance = 1e-8)), logical(1))))
  expect_true(any(vapply(vectors, function(v) isTRUE(all.equal(v, c(0, -1, 0, 1, 0, 0), tolerance = 1e-8)), logical(1))))
})

test_that("generate_spm_contrasts_from_model adds centered projection main-effect contrasts", {
  skip_if_not_installed("R.matlab")

  tmp_dir <- tempfile("spm_contrast_proj_main_")
  dir.create(tmp_dir, recursive = TRUE)

  mnames <- c(
    "Sn(1) face_onset*bf(1)",
    "Sn(2) face_onset*bf(1)",
    "Sn(3) face_onset*bf(1)",
    "Sn(4) face_onset*bf(1)",
    "Sn(1) constant",
    "Sn(2) constant",
    "Sn(3) constant",
    "Sn(4) constant"
  )
  cpos <- 1:4
  bpos <- 5:8
  npos <- integer(0)

  matfile <- file.path(tmp_dir, "design_columns.mat")
  R.matlab::writeMat(matfile, mnames = mnames, cpos = cpos, bpos = bpos, npos = npos)

  cmat <- matrix(
    1,
    nrow = 1,
    dimnames = list("face_main", "face_onset")
  )
  mobj <- list(contrasts = cmat)
  class(mobj) <- c("list", "l1_model_spec")

  fmri.pipeline:::generate_spm_contrasts_from_model(
    output_dir = tmp_dir,
    mobj = mobj,
    spm_path = tmp_dir,
    execute = FALSE,
    average_across_runs = TRUE,
    projection_main_effect_terms = "drug",
    projection_main_effect_weights = data.frame(
      drug = c(1, 1, 0, 0),
      stringsAsFactors = FALSE
    )
  )

  spec_path <- file.path(tmp_dir, "spm_contrast_spec.rds")
  setup_script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package = "fmri.pipeline")
  status <- run_spm_contrast_script(
    setup_script,
    c(
      "-mat_file", matfile,
      "-contrast_rds", spec_path,
      "-average_across_runs", "TRUE",
      "-spm_path", tmp_dir
    )
  )
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))

  mfile <- file.path(tmp_dir, "estimate_glm_contrasts.m")
  expect_true(file.exists(mfile))
  lines <- readLines(mfile)
  row_idx <- which(grepl("proj_drug_x_face_main", lines, fixed = TRUE))
  expect_equal(length(row_idx), 1L)
  convec_line <- lines[row_idx + 1L]
  weights <- as.numeric(strsplit(gsub(".*\\[|\\].*", "", convec_line), ",")[[1]])

  expect_equal(weights[1:4], c(0.5, 0.5, -0.5, -0.5), tolerance = 1e-8)
  expect_equal(weights[5:8], rep(0, 4), tolerance = 1e-8)
})


test_that("SPM contrast mapping rejects missing and non-estimable terms", {
  td <- tempfile("spm_invalid_contrast_")
  dir.create(td)
  mat <- file.path(td, "columns.mat")
  spec <- file.path(td, "contrasts.rds")
  script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package = "fmri.pipeline")
  run <- function() suppressWarnings(run_spm_contrast_script(script,
    c("-mat_file", mat, "-contrast_rds", spec)))
  cmat <- matrix(c(1, -1), 1, dimnames = list("A_minus_B", c("A", "B")))
  saveRDS(list(contrast_matrix = cmat), spec)
  R.matlab::writeMat(mat, mnames = c("Sn(1) A*bf(1)", "Sn(1) constant"),
    cpos = 1, bpos = 2, npos = integer(0), design_matrix = diag(2))
  result <- run()
  expect_equal(attr(result, "status"), 1L)
  expect_true(any(grepl("no SPM columns matched.*B", result)))
  expect_false(file.exists(file.path(td, "estimate_glm_contrasts.m")))

  R.matlab::writeMat(mat, mnames = c("Sn(1) A*bf(1)", "Sn(1) B*bf(1)", "Sn(1) constant"),
    cpos = 1:2, bpos = 3, npos = integer(0), design_matrix = cbind(1:5, 1:5, 1))
  result <- run()
  expect_equal(attr(result, "status"), 1L)
  expect_true(any(grepl("Non-estimable SPM contrasts: A_minus_B", result)))

  cmat[1, ] <- c(1, 1)
  rownames(cmat) <- "Participant's sum"
  saveRDS(list(contrast_matrix = cmat), spec)
  result <- run()
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
  expect_equal(readLines(file.path(td, "spm_expected_contrasts.txt")), "Participant's sum")
  lines <- readLines(file.path(td, "estimate_glm_contrasts.m"))
  expect_true(any(grepl("Participant''s sum", lines, fixed = TRUE)))
  expect_true(any(grepl("con.delete = 1", lines, fixed = TRUE)))

  # A small regressor scale alone must not make a full-rank contrast invalid.
  R.matlab::writeMat(mat, mnames = c("Sn(1) A*bf(1)", "Sn(1) B*bf(1)", "Sn(1) constant"),
    cpos = 1:2, bpos = 3, npos = integer(0), design_matrix = cbind(1e-9 * (1:5), c(1, 1, 2, 2, 4), 1))
  result <- run()
  expect_null(attr(result, "status"), info = paste(result, collapse = "\n"))
})
