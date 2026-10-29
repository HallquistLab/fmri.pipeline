test_that("concat_l1_confounds keeps spikes independent and aligns shared confounds", {
  skip_if_not_installed("data.table")

  tmp_dir <- tempfile("confound_concat_")
  dir.create(tmp_dir, recursive = TRUE)
  sqlite_db <- file.path(tmp_dir, "test.sqlite")

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = FALSE, output_directory = tmp_dir)
  gpa$output_locations$sqlite_db <- sqlite_db
  gpa$glm_software <- "spm"
  gpa$parallel$compute_environment <- list(global = character(0), fsl = character(0), afni = character(0), spm = character(0), r = character(0))

  run1 <- data.frame(
    csf = c(1, 2, 3),
    wm = c(4, 5, 6),
    expr1_spike_1 = c(0, 1, 0),
    stringsAsFactors = FALSE
  )
  run2 <- data.frame(
    csf = c(7, 8, 9),
    wm = c(10, 11, 12),
    expr1_spike_1 = c(0, 0, 1),
    expr1_spike_2 = c(1, 0, 0),
    stringsAsFactors = FALSE
  )

  f1 <- file.path(tmp_dir, "run1_l1_confounds.txt")
  f2 <- file.path(tmp_dir, "run2_l1_confounds.txt")
  write.table(run1, f1, row.names = FALSE, col.names = FALSE)
  write.table(run2, f2, row.names = FALSE, col.names = FALSE)

  cols1 <- data.frame(
    col_index = 1:3,
    col_name = c("csf", "wm", "expr1_spike_1"),
    is_spike = c(FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  cols2 <- data.frame(
    col_index = 1:4,
    col_name = c("csf", "wm", "expr1_spike_1", "expr1_spike_2"),
    is_spike = c(FALSE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  fmri.pipeline:::insert_df_sqlite(gpa, id = "sub1", session = 1L, run_number = 1L, data = cols1, table = "l1_confound_columns", immediate = TRUE)
  fmri.pipeline:::insert_df_sqlite(gpa, id = "sub1", session = 1L, run_number = 2L, data = cols2, table = "l1_confound_columns", immediate = TRUE)

  concat_file <- fmri.pipeline:::concat_l1_confounds(
    gpa = gpa, id = "sub1", session = 1L,
    run_numbers = c(1L, 2L),
    confound_files = c(f1, f2),
    output_dir = tmp_dir
  )

  expect_true(file.exists(concat_file))
  concat_df <- data.table::fread(concat_file, header = FALSE, data.table = FALSE)
  expect_equal(nrow(concat_df), 6)
  expect_equal(ncol(concat_df), 5)
  expected <- cbind(
    csf = c(1, 2, 3, 7, 8, 9),
    wm = c(4, 5, 6, 10, 11, 12),
    run1_spike = c(0, 1, 0, 0, 0, 0),
    run2_spike1 = c(0, 0, 0, 0, 0, 1),
    run2_spike2 = c(0, 0, 0, 1, 0, 0)
  )
  expect_equal(unname(as.matrix(concat_df)), unname(expected))
})

test_that("pooled sessions keep single spikes with repeated run numbers independent", {
  root <- tempfile("spike_sessions_")
  dir.create(root)
  gpa <- structure(list(output_locations = list(sqlite_db = file.path(root, "cache.sqlite"))),
    class = c("glm_pipeline_arguments", "list"))
  files <- file.path(root, paste0("session", 1:2, ".txt"))
  lg <- lgr::get_logger("test/spike_sessions")
  for (session in 1:2) {
    spikes <- fmri.pipeline:::compute_spike_regressors(
      data.frame(fd = c(0, 1, 0)), "fd > 0.5", lg = lg
    )
    expect_identical(dim(spikes), c(3L, 1L))
    expect_identical(colnames(spikes), "expr1_spike_1")
    confounds <- cbind(data.frame(csf = 1:3), spikes)
    metadata <- data.frame(col_index = seq_along(confounds), col_name = names(confounds),
      is_spike = grepl("_spike_", names(confounds)))
    fmri.pipeline:::insert_df_sqlite(gpa, id = "s1", session = session, run_number = 1L,
      data = metadata, table = "l1_confound_columns", immediate = TRUE)
    write.table(confounds, files[session], row.names = FALSE, col.names = FALSE)
  }
  out <- fmri.pipeline:::concat_l1_confounds(gpa, id = "s1", session = 0L,
    run_numbers = c(1L, 1L), run_sessions = 1:2, confound_files = files, output_dir = root)
  actual <- as.matrix(read.table(out))
  expected <- cbind(c(1:3, 1:3), c(0, 1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, 0))
  expect_equal(unname(actual), expected)
  # Independent spikes can remove different amplitudes at the two contaminated volumes.
  y <- c(0, 10, 0, 0, 20, 0)
  expect_equal(unname(lm.fit(actual[, -1, drop = FALSE], y)$residuals), rep(0, 6), tolerance = 1e-10)
})

test_that("spike deduplication preserves single-column matrices and no-spike results", {
  lg <- lgr::get_logger("test/spike_dimensions")
  spikes <- fmri.pipeline:::compute_spike_regressors(
    data.frame(fd = c(0, 1, 0)), c("fd > 0.5", "fd > 0.8"), lg = lg
  )
  expect_identical(dim(spikes), c(3L, 1L))
  expect_identical(colnames(spikes), "expr1_spike_1")
  expect_null(fmri.pipeline:::compute_spike_regressors(
    data.frame(fd = c(0, 0, 0)), "fd > 0.5", lg = lg
  ))
})
