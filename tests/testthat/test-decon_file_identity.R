make_alignment_file <- function(name = "subA_run1_deconvolved.csv.gz") {
  root <- tempfile("alignment_")
  dir.create(root)
  path <- file.path(root, name)
  readr::write_csv(data.frame(time = 0:20, decon = 0:20, vnum = 1L, atlas_value = 1L), path)
  list(root = root, file = path,
    trials = data.frame(id = "A", run = 1L, trial = 1L, onset = 10),
    alignment = list(evt_col = "onset", time_before = -1, time_after = 1,
      pad_before = -1, pad_after = 1, output_resolution = 1))
}

test_that("legacy and BIDS identifiers preserve strings and sessions", {
  expect_equal(fmri.pipeline:::parse_decon_identity(".tmp.data.sub001_run02_bold_deconvolved.csv.gz"),
    list(id = "001", run = "02"))
  expect_equal(fmri.pipeline:::parse_decon_identity("/data/sub-ABC/ses-01/sub-ABC_ses-01_run-02_bold.nii.gz"),
    list(id = "ABC", run = "02", session = "01"))
  expect_error(fmri.pipeline:::parse_decon_identity("sub-A_sub-B_run1.csv"), "Ambiguous id")
})

test_that("legacy files align and expose a complete status manifest", {
  f <- make_alignment_file()
  on.exit(unlink(f$root, recursive = TRUE))
  result <- evt_align_decon_files(f$file, f$trials, f$alignment, tr = 1, ncpus = 1)
  expect_equal(result$decon_mean, c(9, 10, 11))
  expect_equal(unique(result$id), "A")
  expect_equal(unique(result$run), 1L)
  expect_equal(attr(result, "file_status")$status, "ok")
  expect_equal(attr(result, "file_status")$file, f$file)
})

test_that("sidecars support arbitrary filenames and leading-zero session identities", {
  f <- make_alignment_file("arbitrary_deconvolved.csv.gz")
  on.exit(unlink(f$root, recursive = TRUE))
  readr::write_csv(data.frame(id = "001", run = "01", session = "dayB"), paste0(f$file, ".metadata.csv"))
  trials <- data.frame(id = c("001", "1", "001"), run = 1L,
    session = c("dayB", "dayB", "dayA"), trial = 1L, onset = c(10, 15, 5))
  result <- evt_align_decon_files(f$file, trials, f$alignment, tr = 1, ncpus = 1)
  expect_equal(unique(result$id), "001")
  expect_equal(unique(result$session), "dayB")
  expect_equal(result$decon_mean, c(9, 10, 11))
})

test_that("explicit manifests override filenames and session ambiguity is rejected", {
  f <- make_alignment_file("subWrong_run9_deconvolved.csv.gz")
  on.exit(unlink(f$root, recursive = TRUE))
  trials <- data.frame(id = "participantA", run = 1L, session = c("a", "b"), trial = 1L, onset = c(5, 10))
  manifest <- data.frame(file = f$file, id = "participantA", run = 1L, session = "b")
  result <- evt_align_decon_files(f$file, trials, f$alignment, tr = 1, ncpus = 1, file_manifest = manifest)
  expect_equal(unique(result$id), "participantA")
  expect_equal(result$decon_mean, c(9, 10, 11))
  expect_error(evt_align_decon_files(f$file, trials, f$alignment, tr = 1, ncpus = 1,
    file_manifest = manifest[c("file", "id", "run")]), "Ambiguous trial identity")
  expect_error(evt_align_decon_files(f$file, trials, f$alignment, tr = 1, ncpus = 1,
    file_manifest = manifest[c(1L, 1L), ]), "unique")
})

test_that("unknown file identities fail with actionable manifest guidance", {
  f <- make_alignment_file("custom_deconvolved.csv.gz")
  on.exit(unlink(f$root, recursive = TRUE))
  expect_error(evt_align_decon_files(f$file, f$trials, f$alignment, tr = 1, ncpus = 1),
    "supply file_manifest")
})

test_that("worker failures are reported and partial output requires an explicit option", {
  f <- make_alignment_file()
  on.exit(unlink(f$root, recursive = TRUE))
  bad <- file.path(f$root, "subB_run2_deconvolved.csv.gz")
  readr::write_csv(data.frame(time = 0:20, decon = 0:20, atlas_value = 1), bad)
  output <- file.path(f$root, "aligned.csv.gz")
  err <- tryCatch(evt_align_decon_files(c(f$file, bad), f$trials, f$alignment,
    tr = 1, ncpus = 1, out_file = output), error = identity)
  expect_s3_class(err, "decon_alignment_error")
  expect_equal(err$file_status$status, c("ok", "error"))
  expect_match(err$file_status$message[2L], "vnum")
  expect_false(file.exists(output))
  expect_true(file.exists(paste0(output, ".files.csv")))
  expect_warning(result <- evt_align_decon_files(c(f$file, bad), f$trials, f$alignment,
    tr = 1, ncpus = 1, out_file = output, fail_on_error = FALSE), "1 deconvolution file")
  expect_equal(nrow(result), 3L)
  expect_equal(attr(result, "file_status")$status, c("ok", "error"))
  expect_true(file.exists(output))
})

test_that("all-failed batches retain their status instead of disappearing", {
  f <- make_alignment_file()
  on.exit(unlink(f$root, recursive = TRUE))
  readr::write_csv(data.frame(time = 0:20, decon = NA_real_, vnum = 1, atlas_value = 1), f$file)
  expect_warning(result <- evt_align_decon_files(f$file, f$trials, f$alignment,
    tr = 1, ncpus = 1, fail_on_error = FALSE), "No finite deconvolution")
  expect_equal(nrow(result), 0L)
  expect_equal(attr(result, "file_status")$status, "error")
})

test_that("alignment discovers old hidden files and excludes metadata sidecars", {
  f <- make_alignment_file()
  on.exit(unlink(f$root, recursive = TRUE))
  atlas <- file.path(f$root, "atlas.nii.gz")
  RNifti::writeNifti(array(1, c(2L, 2L, 2L)), atlas)
  adir <- file.path(f$root, "atlas", "deconvolved")
  dir.create(adir, recursive = TRUE)
  names <- c(".old_subA_run1_deconvolved.csv.gz", "decon_subA_run2_deconvolved.csv.gz")
  file.create(file.path(adir, names))
  file.create(file.path(adir, paste0(names, ".metadata.csv")))
  captured <- NULL
  local_mocked_bindings(R_batch_job = list(new = function(...) {
    captured <<- list(...)
    list(submit = function() invisible(NULL))
  }))
  run_decon_alignment(atlas, f$root, f$trials, list(onset = f$alignment),
    tr = 1, scheduler = "local", ncpus = 1)
  expect_setequal(basename(captured$input_objects$d_files), names)
  expect_true(captured$input_objects$fail_on_error)
})

test_that("file alignment passes validity masks and overlap policy through to extraction", {
  f <- make_alignment_file()
  on.exit(unlink(f$root, recursive = TRUE))
  readr::write_csv(data.frame(time = 0:20, decon = 0:20, vnum = 1, atlas_value = 1,
    usable = (0:20) != 10), f$file)
  f$alignment$valid_col <- "usable"
  f$alignment$overlap_policy <- "allow"
  f$alignment$neighbor_events <- "onset"
  out <- file.path(f$root, "aligned.csv")
  result <- evt_align_decon_files(f$file, f$trials, f$alignment, tr = 1, ncpus = 1, out_file = out)
  expect_true(file.exists(paste0(out, ".sources.csv")))
  clock <- attr(result, "source_clock")
  expect_equal(clock$time, as.numeric(0:20))
  expect_equal(clock$source_id, 1:21)
  expect_equal(unique(result$alignment_file), normalizePath(f$file))
  expect_equal(result$decon_source_max_epochs[result$evt_time == 0], 0L)
  expect_true(is.na(result$decon_mean[result$evt_time == 0]))
  expect_equal(result$decon_n_censored[result$evt_time == 0], 1)
  expect_true(all(result$align_overlap_policy == "allow"))
  expect_equal(attr(result, "file_status")$status, "ok")
})

test_that("an empty acquisition clock retains no-data epochs and an empty source manifest", {
  f <- make_alignment_file()
  on.exit(unlink(f$root, recursive = TRUE))
  readr::write_csv(data.frame(time = NA_real_, decon = 0:20, vnum = 1, atlas_value = 1), f$file)
  result <- evt_align_decon_files(f$file, f$trials, f$alignment, tr = 1, ncpus = 1)
  expect_equal(nrow(result), 3L)
  expect_true(all(result$epoch_status == "no_data"))
  expect_true(all(is.na(result$decon_mean)))
  expect_equal(nrow(attr(result, "source_clock")), 0L)
  expect_equal(attr(result, "file_status")$status, "ok")
})
