confound_cache_fixture <- function() {
  outdir <- tempfile("confound_cache_")
  dir.create(outdir)
  confounds <- file.path(outdir, "confounds.tsv")
  write.table(data.frame(csf = seq_len(20), wm = seq_len(20)^2), confounds,
    sep = "\t", row.names = FALSE, quote = FALSE)
  motion <- file.path(outdir, "motion.par")
  write.table(cbind(0, 0, 0, seq(0, 1.9, by = 0.1), 0, 0), motion,
    row.names = FALSE, col.names = FALSE)
  image <- file.path(outdir, "bold.nii.gz")
  RNifti::writeNifti(array(seq_len(160), c(2, 2, 2, 20)), image)
  structure(list(
    output_directory = outdir,
    output_locations = list(sqlite_db = file.path(outdir, "cache.sqlite"),
      feat_sub_directory = "{gpa$output_directory}/feat_l1/sub-{id}"),
    run_data = data.frame(id = "sub1", session = 1L, run_number = 1L,
      run_nifti = image, run_volumes = 20L, drop_volumes = 0L, tr = 1,
      confound_input_file = confounds, confound_input_file_present = TRUE,
      motion_params_file = motion, motion_params_present = TRUE),
    drop_volumes = 0L,
    confound_settings = list(motion_params_colnames = c("rx", "ry", "rz", "tx", "ty", "tz"),
      l1_confound_regressors = c("csf", "wm"), exclude_run = NULL,
      run_exclusion_columns = NULL, truncate_run = NULL, spike_volumes = NULL, na_strings = "NA"),
    lgr_threshold = "fatal"
  ), class = c("glm_pipeline_arguments", "list"))
}

cached_run <- function(gpa, run_df = gpa$run_data, demean = FALSE) {
  fmri.pipeline:::get_l1_confounds(run_df = run_df, gpa = gpa, demean = demean)
}

cached_values <- function(run) as.matrix(read.table(run$l1_confound_file, header = FALSE))

test_that("confound settings, demeaning, and input contents invalidate all cache layers", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  first <- cached_run(gpa)
  expect_equal(unname(cached_values(first)), cbind(seq_len(20), seq_len(20)^2))

  # An unchanged run really reuses its file, including when passed processed data.
  stamp <- as.POSIXct("2001-01-01", tz = "UTC")
  Sys.setFileTime(first$l1_confound_file, stamp)
  same <- cached_run(gpa, first)
  expect_equal(as.numeric(file.info(same$l1_confound_file)$mtime), as.numeric(stamp))

  gpa$confound_settings$l1_confound_regressors <- c("wm", "csf")
  reordered <- cached_run(gpa, same)
  expect_equal(unname(cached_values(reordered)), cbind(seq_len(20)^2, seq_len(20)))
  centered <- cached_run(gpa, reordered, demean = TRUE)
  expect_equal(unname(colMeans(cached_values(centered))), c(0, 0))

  write.table(data.frame(csf = seq_len(20) + 100, wm = seq_len(20)^2 + 200),
    gpa$run_data$confound_input_file, sep = "\t", row.names = FALSE, quote = FALSE)
  changed <- cached_run(gpa, centered)
  expect_equal(unname(cached_values(changed)), cbind(seq_len(20)^2 + 200, seq_len(20) + 100))

  gpa$confound_settings$l1_confound_regressors <- NULL
  removed <- cached_run(gpa, changed)
  expect_true(is.na(removed$l1_confound_file))
  expect_false(file.exists(first$l1_confound_file))
  expect_null(fmri.pipeline:::read_df_sqlite(gpa, id = "sub1", session = 1L,
    run_number = 1L, table = "l1_confound_columns"))
})

test_that("exclusion and spike changes replace previous decisions and motion data", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  gpa$confound_settings$exclude_run <- "max(framewise_displacement) > 0.05"
  gpa$confound_settings$run_exclusion_columns <- "framewise_displacement"
  first <- cached_run(gpa)
  expect_true(first$exclude_run)
  gpa$confound_settings$exclude_run <- "max(framewise_displacement) > 0.5"
  second <- cached_run(gpa, first)
  expect_false(second$exclude_run)
  gpa$confound_settings$spike_volumes <- "framewise_displacement > 0.05"
  spiked <- cached_run(gpa, second)
  expect_equal(ncol(cached_values(spiked)), 21L)
  gpa$confound_settings$spike_volumes <- "framewise_displacement > 0.5"
  unspiked <- cached_run(gpa, spiked)
  expect_equal(ncol(cached_values(unspiked)), 2L)
  write.table(matrix(0, 20, 6), gpa$run_data$motion_params_file,
    row.names = FALSE, col.names = FALSE)
  gpa$confound_settings$exclude_run <- "max(framewise_displacement) > 0.05"
  expect_false(cached_run(gpa, unspiked)$exclude_run)
  gpa$confound_settings$exclude_run <- NULL
  expect_false(cached_run(gpa)$exclude_run)
})

test_that("truncation changes use the original image without dropping volumes twice", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  gpa$run_data$drop_volumes <- 2L
  gpa$confound_settings$truncate_run <- "csf > 15"
  gpa$confound_settings$run_truncation_columns <- "csf"
  first <- cached_run(gpa)
  expect_equal(first$run_volumes, 13L)
  expect_equal(cached_run(gpa, first)$run_volumes, 13L)
  stamp <- as.POSIXct("2001-01-01", tz = "UTC")
  Sys.setFileTime(first$l1_confound_file, stamp)
  cached_run(gpa, first)
  expect_equal(as.numeric(file.info(first$l1_confound_file)$mtime), as.numeric(stamp))
  gpa$confound_settings$truncate_run <- "csf > 18"
  second <- cached_run(gpa, first)
  expect_equal(second$run_volumes, 16L)
  expect_equal(unname(cached_values(second)[, 1]), 3:18)
  expect_equal(as.array(RNifti::readNifti(second$run_nifti)),
    array(seq_len(160), c(2, 2, 2, 20))[, , , 3:18], ignore_attr = TRUE)
  gpa$confound_settings$truncate_run <- NULL
  second$drop_volumes <- 0L
  restored <- cached_run(gpa, second)
  expect_equal(restored$run_nifti, gpa$run_data$run_nifti)
  expect_equal(restored$run_volumes, 20L)
})

test_that("changed images replace an existing truncated output", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  gpa$run_data$drop_volumes <- 2L
  first <- cached_run(gpa)
  RNifti::writeNifti(array(42, c(2, 2, 2, 20)), gpa$run_data$run_nifti)
  Sys.setFileTime(gpa$run_data$run_nifti, Sys.time() + 5)
  second <- cached_run(gpa, first)
  expect_equal(second$run_nifti, first$run_nifti)
  expect_true(all(RNifti::readNifti(second$run_nifti) == 42))
})

test_that("legacy caches rebuild and failed rebuilds remain invalid", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  first <- cached_run(gpa)
  manifest <- file.path(dirname(first$l1_confound_file), "run1_l1_confounds_cache.rds")
  unlink(manifest)
  writeLines("999 999", first$l1_confound_file)
  rebuilt <- cached_run(gpa)
  expect_equal(nrow(cached_values(rebuilt)), 20L)
  gpa$confound_settings$l1_confound_regressors <- c("wm", "csf")
  local_mocked_bindings(confound_manipulations = function(...) stop("interrupted rebuild"),
    .package = "fmri.pipeline")
  expect_error(cached_run(gpa), "interrupted rebuild")
  expect_null(readRDS(manifest)$signature)
})

test_that("cache invalidation is scoped to one run and repairs missing artifacts", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  first <- cached_run(gpa)
  other_run <- gpa$run_data
  other_run$run_number <- 2L
  other <- cached_run(gpa, other_run)
  stamp <- as.POSIXct("2001-01-01", tz = "UTC")
  Sys.setFileTime(other$l1_confound_file, stamp)
  gpa$confound_settings$l1_confound_regressors <- c("wm", "csf")
  cached_run(gpa, first)
  expect_equal(as.numeric(file.info(other$l1_confound_file)$mtime), as.numeric(stamp))
  expect_equal(unname(cached_values(other)), cbind(seq_len(20), seq_len(20)^2))
  other_cached <- fmri.pipeline:::read_df_sqlite(gpa, id = "sub1", session = 1L,
    run_number = 2L, table = "l1_run_calculations")
  expect_equal(nrow(other_cached), 1L)
  unlink(first$l1_confound_file)
  repaired <- cached_run(gpa)
  expect_equal(unname(cached_values(repaired)), cbind(seq_len(20)^2, seq_len(20)))
})

test_that("legacy truncated runs require source provenance", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  gpa$run_data$drop_volumes <- 2L
  first <- cached_run(gpa)
  unlink(file.path(dirname(first$l1_confound_file), "run1_l1_confounds_cache.rds"))
  expect_error(cached_run(gpa, first), "Supply the original run_nifti and run_volumes")
  expect_equal(cached_run(gpa)$run_volumes, 18L)
})

test_that("no-spike results can be cached and later gain spike regressors", {
  gpa <- confound_cache_fixture()
  on.exit(unlink(gpa$output_directory, recursive = TRUE))
  gpa$confound_settings$l1_confound_regressors <- NULL
  gpa$confound_settings$spike_volumes <- "framewise_displacement > 0.5"
  first <- cached_run(gpa)
  expect_true(is.na(first$l1_confound_file))
  manifest <- file.path(gpa$output_directory, "feat_l1", "sub-sub1", "run1_l1_confounds_cache.rds")
  stamp <- as.POSIXct("2001-01-01", tz = "UTC")
  Sys.setFileTime(manifest, stamp)
  expect_true(is.na(cached_run(gpa, first)$l1_confound_file))
  expect_equal(as.numeric(file.info(manifest)$mtime), as.numeric(stamp))
  gpa$confound_settings$spike_volumes <- "framewise_displacement > 0.05"
  expect_equal(ncol(cached_values(cached_run(gpa, first))), 19L)
})
