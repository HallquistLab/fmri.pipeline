test_that("get_spm_status does not treat marker-only runs as complete", {
  td <- tempfile("spm_status_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines(as.character(c(Sys.time(), Sys.time())), file.path(td, ".spm_complete"))

  status <- get_spm_status(td)

  expect_false(isTRUE(status$spm_complete[[1]]))
  expect_true(isTRUE(status$spm_failed[[1]]))
  expect_false(isTRUE(status$spm_mat_exists[[1]]))
  expect_false(isTRUE(status$spm_beta_exists[[1]]))
  expect_false(isTRUE(status$spm_minimal_outputs_exist[[1]]))
})

test_that("get_spm_status rejects empty SPM.mat and beta files", {
  td <- tempfile("spm_status_")
  dir.create(td, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)

  writeLines(as.character(c(Sys.time(), Sys.time())), file.path(td, ".spm_complete"))
  file.create(file.path(td, "SPM.mat"))
  file.create(file.path(td, "beta_0001.nii"))

  status <- get_spm_status(td)

  expect_false(isTRUE(status$spm_complete[[1]]))
  expect_true(isTRUE(status$spm_failed[[1]]))
  expect_true(isTRUE(status$spm_mat_exists[[1]]))
  expect_true(isTRUE(status$spm_beta_exists[[1]]))
  expect_false(isTRUE(status$spm_minimal_outputs_exist[[1]]))
})
