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
  reference <- deconvolve_reglin(y, kernel = c(0, 0.5, 0.5),
    lambda_grid = c(0.01, 0.1), normalize = FALSE, demean = FALSE)[seq_len(30L), 1L]
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


test_that("voxelwise linear scaling is explicit, recorded and preserves amplitude", {
  f <- make_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  mask <- array(0, c(2, 2, 2)); mask[2] <- 1
  RNifti::writeNifti(mask, f$mask)
  raw <- matrix(f$data, nrow = 8)[2, ]
  for (mode in c("none", "demean", "psc", "zscore")) {
    voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas, mask = f$mask,
      nprocs = 1, algorithm = "reglin", out_dir = f$root,
      out_file_expression = substitute(NAME, list(NAME = mode)),
      decon_settings = list(kernel = c(0, 0.5, 0.5), lambda = 0.1, scaling = mode))
    file <- file.path(f$root, "atlas", "deconvolved", paste0(mode, "_deconvolved.csv.gz"))
    d <- data.table::fread(file)
    transform <- data.table::fread(paste0(file, ".scaling.csv"))
    expect_equal(transform$scaling, mode)
    y <- (raw - transform$center) / transform$scale
    reference <- deconvolve_reglin(y, c(0, 0.5, 0.5), lambda = 0.1,
      scaling = "none")[1:30, 1]
    expect_equal(d$decon, reference, tolerance = 1e-7)
    original <- data.table::fread(file.path(f$root, "atlas", "original", paste0(mode, "_original.csv.gz")))
    expect_equal(original$BOLD_input, y, tolerance = 1e-7)
    expect_equal("BOLD_z" %in% names(original), mode == "zscore")
    expect_equal(data.table::fread(paste0(file, ".metadata.csv"))$scaling, mode)
  }
  # Explicit tail_padding is respected when requested for reglin.
  voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas, mask = f$mask,
    nprocs = 1, algorithm = "reglin", out_dir = f$root,
    out_file_expression = expression("padded"),
    decon_settings = list(kernel = c(0, 0.5, 0.5), lambda = 0.1, scaling = "none", tail_padding = 3L))
  padded_fit <- data.table::fread(file.path(f$root, "atlas", "deconvolved", "padded_deconvolved.csv.gz"))
  padded_ref <- deconvolve_reglin(c(raw, 0, 0, 0), c(0, 0.5, 0.5), lambda = 0.1,
    scaling = "none")[1:30, 1]
  expect_equal(padded_fit$decon, padded_ref, tolerance = 1e-7)
  # The legacy opt-out must reach the solver without an earlier hidden z-score.
  voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas, mask = f$mask,
    nprocs = 1, algorithm = "reglin", out_dir = f$root,
    out_file_expression = expression("legacy"),
    decon_settings = list(kernel = c(0, 0.5, 0.5), lambda = 0.1, normalize = FALSE, demean = FALSE))
  legacy <- data.table::fread(file.path(f$root, "atlas", "deconvolved", "legacy_deconvolved.csv.gz"))
  raw_fit <- data.table::fread(file.path(f$root, "atlas", "deconvolved", "none_deconvolved.csv.gz"))
  expect_equal(legacy$decon, raw_fit$decon)
})

test_that("voxelwise Bush inputs are standardized exactly once before padding", {
  f <- make_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE))
  mask <- array(0, c(2, 2, 2)); mask[2] <- 1
  RNifti::writeNifti(mask, f$mask)
  raw <- matrix(f$data, nrow = 8)[2, ]
  local_mocked_bindings(deconvolve_nlreg = function(BOLDobs, normalize, nev_lr, epsilon, beta, ...) {
    expect_equal(c(nev_lr, epsilon, beta), c(0.01, 0.005, 60))
    expect_false(normalize)
    expect_equal(as.numeric(BOLDobs), c(as.numeric(scale(raw)), 0, 0, 0))
    BOLDobs
  })
  voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas, mask = f$mask,
    nprocs = 1, algorithm = "bush2011", out_dir = f$root,
    decon_settings = list(kernel = c(0, 0.5, 0.5)))
  local_mocked_bindings(deconvolve_nlreg_resample = function(bold_obs, beta, normalize, ...) {
    expect_equal(beta, 17)
    expect_false(normalize)
    expect_equal(bold_obs, c(as.numeric(scale(raw)), 0, 0, 0))
    list(NEVmean = bold_obs)
  })
  voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas, mask = f$mask,
    nprocs = 1, algorithm = "bush2015", out_dir = f$root,
    out_file_expression = expression("resampled"),
    decon_settings = list(kernel = c(0, 0.5, 0.5), beta = 17))
  for (algorithm in c("bush2011", "bush2015")) {
    expect_error(voxelwise_deconvolution(f$bold, TR = 1, atlas_files = f$atlas,
      nprocs = 1, algorithm = algorithm,
      decon_settings = list(kernel = c(0, 0.5, 0.5), scaling = "none")),
      "requires scaling = 'zscore'")
  }
})

test_that("Bush resampling uses the same standardized observations in its base and residual fits", {
  set.seed(2026)
  calls <- list()
  local_mocked_bindings(deconvolve_nlreg = function(BOLDobs, kernel, ..., normalize, trim_kernel) {
    calls[[length(calls) + 1L]] <<- BOLDobs
    expect_false(normalize)
    matrix(rnorm(nrow(BOLDobs) + length(kernel) - 1L), ncol = 1)
  })
  y <- 100 + sin(1:30)
  deconvolve_nlreg_resample(y, c(0.1, 0.6, 0.3), n_resample = 3)
  expect_equal(as.numeric(calls[[1]]), as.numeric(scale(y)))
  calls <- list()
  deconvolve_nlreg_resample(y, c(0.1, 0.6, 0.3), n_resample = 3, normalize = FALSE)
  expect_equal(as.numeric(calls[[1]]), y)
})

test_that("invalid PSD tuning settings fail before voxelwise workers run", {
  expect_error(voxelwise_deconvolution(character(), TR = 1, algorithm = "reglin",
    decon_settings = list(kernel = c(0.1, 0.6, 0.3), backend = "wiener", wiener_shrinkage = "psd")),
    "Automatic GCV tuning is not supported")
})
