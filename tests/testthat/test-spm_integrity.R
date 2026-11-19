test_that("SPM completion requires every validated output and a matching setup", {
  td <- tempfile("spm_valid_")
  manifest <- make_spm_completed_fixture(td, c("B", "A"), c(1L, 3L))
  expect_true(get_spm_status(td)$spm_complete)
  expect_equal(fmri.pipeline:::get_spm_contrast_file(td, "A"), manifest$coef_file[2])
  expect_true(is.na(fmri.pipeline:::get_spm_contrast_file(td, "unknown")))
  for (name in c(".spm_complete", "SPM.mat", "beta_0002.nii", "spmT_0003.nii", "con_0003.nii", "spm_validated_token.txt")) {
    path <- file.path(td, name)
    old <- readBin(path, "raw", n = file.info(path)$size)
    unlink(path)
    expect_false(get_spm_status(td)$spm_complete, info = name)
    writeBin(old, path)
  }
  writeLines(rep(as.character(Sys.time()), 2L), file.path(td, ".spm_fail"))
  expect_false(get_spm_status(td)$spm_complete)
  unlink(file.path(td, ".spm_fail"))
  writeLines("changed input", rownames(readRDS(file.path(td, "spm_setup.rds"))$signature$files))
  expect_false(get_spm_status(td)$spm_complete)
})

test_that("SPM setup invalidation removes stale completion and checks design specifications", {
  td <- tempfile("spm_cache_")
  make_spm_completed_fixture(td)
  saved <- readRDS(file.path(td, "spm_setup.rds"))$signature
  expect_true(fmri.pipeline:::spm_setup_current(td, saved))
  changed <- saved
  changed$spec$settings$hpf <- 40
  expect_false(fmri.pipeline:::spm_setup_current(td, changed))
  fmri.pipeline:::spm_invalidate_setup(td)
  expect_false(file.exists(file.path(td, ".spm_complete")))
  expect_false(fmri.pipeline:::spm_setup_current(td, saved))
  fmri.pipeline:::spm_save_setup(td, changed)
  expect_true(fmri.pipeline:::spm_setup_current(td, changed))
  expect_false(get_spm_status(td)$spm_complete)
})

test_that("SPM aligns original and retained image and nuisance rows", {
  td <- tempfile("spm_alignment_")
  bdm <- make_spm_input_fixture(td)
  conf <- file.path(td, "conf.txt")
  write.table(matrix(1:20, nrow = 10), conf, row.names = FALSE, col.names = FALSE)
  aligned <- fmri.pipeline:::spm_align_inputs(bdm, bdm$run_niftis, conf, TRUE, td)
  expect_equal(aligned$lengths, c(3L, 5L))
  expect_equal(aligned$offsets, c(2L, 0L))
  expect_equal(read.table(aligned$ts_files)$V1, 3:10)
  retained <- make_spm_input_fixture(file.path(td, "cropped"), c(3L, 5L), c(0L, 0L))
  bdm$run_niftis <- retained$run_niftis
  aligned <- fmri.pipeline:::spm_align_inputs(bdm, bdm$run_niftis, aligned$ts_files, TRUE, td)
  expect_equal(aligned$offsets, c(0L, 0L))
  expect_equal(nrow(read.table(aligned$ts_files)), 8L)
  writeLines(rep("0 1", 9), conf)
  expect_error(fmri.pipeline:::spm_align_inputs(bdm, bdm$run_niftis, conf, TRUE, td), "row count")
  writeLines(c("0 1", "NA 1"), conf)
  expect_error(fmri.pipeline:::spm_align_inputs(bdm, bdm$run_niftis, conf, TRUE, td), "finite numeric")
  expect_error(fmri.pipeline:::spm_align_inputs(bdm, bdm$run_niftis, paste0(conf, "missing"), TRUE, td))
  bdm$run_volumes <- c(2L, 5L)
  expect_error(fmri.pipeline:::spm_align_inputs(bdm, bdm$run_niftis, NULL, TRUE, td), "image lengths")
})

test_that("SPM keeps identically named compressed runs separate and honors leading drops", {
  td <- tempfile("spm_gunzip_")
  bdm <- make_spm_input_fixture(td, compressed = TRUE)
  out <- file.path(td, "model")
  generate_spm_mat(bdm, output_dir = out, spm_path = td, concatenate_runs = TRUE,
    nifti_tmpdir = file.path(td, "shared_tmp"))
  lines <- readLines(file.path(out, "glm_design_batch.m"))
  expect_true(any(grepl("Nt = 3;", lines, fixed = TRUE)))
  expect_true(any(grepl("num2str(i + 2)", lines, fixed = TRUE)))
  tmpdir <- readLines(file.path(out, ".nifti_tmpdir"))
  for (attempt in 1:2) {
    expect_equal(system2("bash", shQuote(file.path(out, "gunzip_commands.sh"))), 0L)
    images <- list.files(tmpdir, full.names = TRUE, pattern = "nii$")
    expect_length(images, 2L)
    expect_equal(vapply(images, function(f) as.numeric(RNifti::readNifti(f)[1]), numeric(1)), c(1, 2), ignore_attr = TRUE)
    unlink(tmpdir, recursive = TRUE)
  }
  out2 <- file.path(td, "another_model")
  generate_spm_mat(bdm, output_dir = out2, spm_path = td, nifti_tmpdir = file.path(td, "shared_tmp"))
  expect_false(identical(tmpdir, readLines(file.path(out2, ".nifti_tmpdir"))))
})

test_that("SPM BDM cache compares timing and model inputs", {
  path <- tempfile(fileext = ".RData")
  d_obj <- structure(list(design = "old"), class = "bdm")
  bdm_args <- list(events = c(1, 3), signals = "A", lg = NULL)
  save(d_obj, bdm_args, file = path)
  lg <- lgr::get_logger("test")
  expect_false(fmri.pipeline:::load_cached_l1_bdm(path, lg, bdm_args)$run_bdm)
  bdm_args$events <- c(2, 4)
  expect_true(fmri.pipeline:::load_cached_l1_bdm(path, lg, bdm_args)$run_bdm)
})


test_that("SPM L1 setup refreshes changed settings, events and confounds", {
  td <- tempfile("spm_l1_cache_")
  bdm <- make_spm_input_fixture(td, c(10L, 10L), c(0L, 0L))
  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = TRUE, output_directory = td)
  gpa$glm_software <- "spm"
  gpa$glm_settings <- list(spm = list(spm_path = td))
  gpa$parallel$compute_environment <- list(global = character(), spm = character(), r = character())
  gpa$output_locations$spm_l1_directory <- file.path(td, "{id}", "ses-{session}", "{l1_model}")
  gpa$l1_models$models <- list(model1 = structure(list(name = "model1",
    contrasts = matrix(1, 1, dimnames = list("A", "A"))), class = c("l1_model_spec", "list")))
  setup <- function(conf = NULL) spm_l1_model(id = "sub1", session = 1L, model_name = "model1",
    d_obj = bdm, run_nifti = bdm$run_niftis, l1_confound_files = conf, gpa = gpa)
  out <- setup()$spm_dir
  token <- function() readLines(file.path(out, "spm_setup_token.txt"))
  old <- token()
  setup()
  expect_identical(token(), old)
  bdm$design[[1]][1, "onset"] <- 1
  setup()
  expect_false(identical(token(), old))
  old <- token()
  gpa$glm_settings$spm$hpf <- 50
  setup()
  expect_false(identical(token(), old))
  expect_error(setup(c("missing1", "missing2")), "refusing to omit nuisance")
})


test_that("SPM confound resolution never drops requested nuisance adjustment", {
  td <- tempfile("spm_confounds_")
  dir.create(td)
  file <- file.path(td, "nuisance.txt")
  writeLines("1 2", file)
  gpa <- structure(list(glm_settings = list(spm = list(concatenate_runs = TRUE))), class = "glm_pipeline_arguments")
  resolve <- function(files) fmri.pipeline:::resolve_spm_confound_inputs(gpa, "sub1", 1L,
    1:2, files, td, lgr::get_logger("test"))
  expect_null(resolve(c("", "")))
  expect_error(resolve(c(file, "missing.txt")), "refusing to omit nuisance")
  local_mocked_bindings(concat_l1_confounds = function(...) NULL, .package = "fmri.pipeline")
  expect_error(resolve(rep(file, 2)), "Unable to concatenate SPM confounds")
})


test_that("SPM L3 setup refreshes changed contrasts and lower-level maps", {
  td <- tempfile("spm_l3_cache_")
  dir.create(td)
  maps <- file.path(td, paste0("subject", 1:4, ".nii"))
  for (map in maps) writeLines("map", map)
  input <- data.frame(id = paste0("sub", 1:4), session = 1L,
    l1_model = "model1", l1_cope_name = "A", l3_model = "group", con_file = maps)
  gpa <- create_mock_gpa(n_subjects = 4, include_l1_models = TRUE, output_directory = td)
  gpa$glm_settings <- list(spm = list(spm_path = td))
  gpa$parallel$compute_environment <- list(global = character(), spm = character(), r = character())
  gpa$output_locations$spm_l3_directory <- file.path(td, "group")
  gpa$l3_models <- list(models = list(group = list(contrasts = matrix(1, 1, dimnames = list("mean", "(Intercept)")))))
  local_mocked_bindings(respecify_l3_model = function(mobj, new_data) {
    mobj$metadata <- new_data[, c("id", "session")]
    mobj$model_matrix <- matrix(1, nrow(new_data), 1, dimnames = list(NULL, "(Intercept)"))
    mobj
  }, .package = "fmri.pipeline")
  out <- spm_l3_model(input, gpa)$spm_dir
  token <- function() readLines(file.path(out, "spm_setup_token.txt"))
  old <- token()
  spm_l3_model(input, gpa)
  expect_identical(token(), old)
  rownames(gpa$l3_models$models$group$contrasts) <- "renamed"
  spm_l3_model(input, gpa)
  expect_false(identical(token(), old))
  expect_equal(readLines(file.path(out, "spm_expected_contrasts.txt")), "renamed")
  old <- token()
  writeLines("modified map", maps[1])
  spm_l3_model(input, gpa)
  expect_false(identical(token(), old))
})
