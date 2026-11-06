test_that("the retired external Bush 2011 backend is rejected", {
  expect_error(
    voxelwise_deconvolution(
      niftis = character(0),
      TR = 1,
      algorithm = "bush2011_external"
    ),
    "one of"
  )
})

make_decon_fixture <- function() {
  root <- tempfile("decon_fixture_")
  dir.create(root)
  atlas <- array(seq_len(8L), c(2L, 2L, 2L))
  mask <- array(0, dim(atlas))
  mask[c(2L, 5L, 7L)] <- 1
  bold <- array(0, c(dim(atlas), 30L))
  for (i in seq_len(8L)) bold[i + 8L * (0:29)] <- 100 + sin((0:29) * i / 8)
  bold[5L + 8L * (0:29)] <- 100 # constant voxel inside the mask
  bold[7L + 8L * (0:29)] <- NA_real_
  paths <- file.path(root, c("atlas.nii.gz", "mask.nii.gz", "bold.nii.gz"))
  RNifti::writeNifti(atlas, paths[1L])
  RNifti::writeNifti(mask, paths[2L])
  RNifti::writeNifti(bold, paths[3L])
  list(root = root, atlas = paths[1L], mask = paths[2L], bold = paths[3L], data = bold)
}

test_that("internal extraction uses the same masked indices as voxel metadata", {
  f <- make_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  local_mocked_bindings(run_afni_command = function(...) stop("AFNI must not be called"))

  expect_warning(
    voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas, mask = f$mask,
      nprocs = 1, algorithm = "reglin", out_dir = f$root,
      out_file_expression = expression("sub1_run1"),
      decon_settings = list(kernel = c(0, 0.5, 0.5), lambda_grid = c(0.01, 0.1))),
    "2 constant or nonfinite voxel signals"
  )
  d <- data.table::fread(file.path(f$root, "atlas", "deconvolved", "sub1_run1_deconvolved.csv.gz"))
  original <- data.table::fread(file.path(f$root, "atlas", "original", "sub1_run1_original.csv.gz"))
  expect_equal(nrow(d), 90L)
  expect_equal(unique(d$atlas_value), c(2L, 5L, 7L))
  expect_equal(unique(d$decon_status), c("ok", "constant", "nonfinite"))
  expect_true(all(is.finite(d$decon[d$vnum == 1L])))
  expect_true(all(is.na(d$decon[d$vnum != 1L])))
  y <- as.numeric(scale(matrix(f$data, nrow = 8L)[2L, ]))
  expect_equal(original$BOLD_z[original$vnum == 1L], y, tolerance = 1e-7)
  reference <- deconvolve_reglin(c(y, rep(0, 3L)), kernel = c(0, 0.5, 0.5),
    lambda_grid = c(0.01, 0.1), normalize = FALSE)[seq_len(30L), 1L]
  expect_equal(d$decon[d$vnum == 1L], reference, tolerance = 1e-7)
  expect_false(file.exists(file.path(f$root, "deconvolve_errors")))
})

test_that("all invalid voxel signals retain their output rows", {
  f <- make_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  RNifti::writeNifti(array(100, c(2L, 2L, 2L, 30L)), f$bold)
  local_mocked_bindings(
    run_afni_command = function(...) stop("AFNI must not be called"),
    deconvolve_reglin = function(...) stop("Invalid signals must not reach the solver")
  )
  expect_warning(
    voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas, mask = f$mask,
      nprocs = 1, algorithm = "reglin", out_dir = f$root,
      out_file_expression = expression("sub1_run1"),
      decon_settings = list(kernel = c(0, 0.5, 0.5))),
    "3 constant or nonfinite voxel signals"
  )
  d <- data.table::fread(file.path(f$root, "atlas", "deconvolved", "sub1_run1_deconvolved.csv.gz"))
  expect_equal(nrow(d), 90L)
  expect_true(all(is.na(d$decon)))
  expect_true(all(d$decon_status == "constant"))
})

test_that("a single selected voxel keeps its signal and checked world coordinates", {
  f <- make_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  mask <- array(0, c(2L, 2L, 2L))
  mask[2L] <- 1
  RNifti::writeNifti(mask, f$mask)
  affine <- diag(4)
  affine[1:3, 4L] <- c(10, 20, 30)
  attr(affine, "code") <- 1L
  for (path in c(f$atlas, f$mask, f$bold)) {
    image <- RNifti::readNifti(path)
    RNifti::sform(image) <- affine
    RNifti::writeNifti(image, path)
  }
  expect_no_warning(voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas,
    mask = f$mask, nprocs = 1, algorithm = "reglin", out_dir = f$root,
    out_file_expression = expression("sub1_run1"), afni_dir = "/unused/afni",
    decon_settings = list(kernel = c(0, 0.5, 0.5), lambda = 0.1)))
  d <- data.table::fread(file.path(f$root, "atlas", "deconvolved", "sub1_run1_deconvolved.csv.gz"))
  expect_equal(nrow(d), 30L)
  expect_equal(unique(d$atlas_value), 2L)
  expect_equal(unique(d[, c("x", "y", "z")]), data.table::data.table(x = 11, y = 20, z = 30))
  expect_true(all(is.finite(d$decon)))
})

test_that("spatially incompatible masks and inputs fail before extraction", {
  f <- make_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  local_mocked_bindings(run_afni_command = function(...) stop("Unexpected extraction"))
  run <- function() voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas,
    mask = f$mask, nprocs = 1, algorithm = "reglin", out_dir = f$root,
    out_file_expression = expression("sub1_run1"), decon_settings = list(kernel = c(0, 0.5, 0.5)))
  RNifti::writeNifti(array(1, c(3L, 2L, 2L)), f$mask)
  expect_error(run(), "incompatible spatial geometry")
  RNifti::writeNifti(array(1, c(2L, 2L, 2L)), f$mask)
  RNifti::writeNifti(array(1, c(3L, 2L, 2L, 30L)), f$bold)
  expect_error(run(), "incompatible spatial geometry")
})

test_that("default output names and identity sidecars align without filename parsing", {
  f <- make_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  expect_warning(voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas,
    mask = f$mask, nprocs = 1, algorithm = "reglin", out_dir = f$root,
    add_metadata = data.frame(id = "001", run = 1L, session = "B"),
    decon_settings = list(kernel = c(0, 0.5, 0.5), lambda = 0.1)),
    "2 constant or nonfinite voxel signals")
  files <- list.files(file.path(f$root, "atlas", "deconvolved"),
    pattern = "_deconvolved\\.csv\\.gz$", full.names = TRUE)
  expect_length(files, 1L)
  expect_false(startsWith(basename(files), "."))
  expect_true(file.exists(paste0(files, ".metadata.csv")))
  result <- evt_align_decon_files(files,
    data.frame(id = "001", run = 1L, session = "B", trial = 1L, onset = 10),
    list(evt_col = "onset", time_before = -1, time_after = 1), tr = 1, ncpus = 1)
  expect_equal(unique(result$id), "001")
  expect_equal(unique(result$session), "B")
  expect_true(all(is.finite(result$decon_mean[result$atlas_value == 2])))
  expect_equal(attr(result, "file_status")$status, "ok")
})
