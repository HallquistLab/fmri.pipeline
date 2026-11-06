test_that("native extraction preserves requested voxel order and all volumes", {
  image <- array(seq_len(2L * 3L * 2L * 4L) / 7, c(2L, 3L, 2L, 4L))
  indices <- c(12L, 2L, 7L, 1L)
  for (extension in c(".nii", ".nii.gz")) {
    path <- tempfile(fileext = extension)
    RNifti::writeNifti(image, path)
    expected <- matrix(as.numeric(RNifti::readNifti(path)), nrow = 12L)[indices, , drop = FALSE]
    expect_equal(fmri.pipeline:::extract_nifti_voxel_ts(path, indices), expected)
    expect_equal(dim(fmri.pipeline:::extract_nifti_voxel_ts(path, 7L)), c(1L, 4L))
    unlink(path)
  }
})

test_that("native extraction applies on-disk integer scaling and retains missing data", {
  path <- tempfile(fileext = ".nii")
  on.exit(unlink(path))
  image <- array(seq_len(48L), c(2L, 2L, 2L, 6L))
  RNifti::writeNifti(image, path, datatype = "int16")
  # Set scaling on the stored integer data without rescaling the payload.
  # NIfTI-1 stores float32 scl_slope and scl_inter at byte offsets 112 and 116.
  con <- file(path, "r+b")
  seek(con, 112L, rw = "write")
  writeBin(c(2, -10), con, size = 4L, endian = .Platform$endian)
  close(con)
  expect_equal(RNifti::niftiHeader(path)$scl_slope, 2)
  expected <- matrix(2 * as.numeric(image) - 10, nrow = 8L)[c(8L, 2L), , drop = FALSE]
  expect_equal(fmri.pipeline:::extract_nifti_voxel_ts(path, c(8L, 2L)), expected)

  image <- array(as.numeric(image), dim(image))
  image[c(2L, 8L, 16L)] <- c(NA_real_, NaN, Inf)
  RNifti::writeNifti(image, path)
  expected <- matrix(as.numeric(RNifti::readNifti(path)), nrow = 8L)[c(8L, 2L), , drop = FALSE]
  expect_equal(fmri.pipeline:::extract_nifti_voxel_ts(path, c(8L, 2L)), expected)
})

test_that("native extraction validates indices and handles a single volume", {
  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path))
  RNifti::writeNifti(array(seq_len(8L), c(2L, 2L, 2L)), path)
  expect_equal(fmri.pipeline:::extract_nifti_voxel_ts(path, c(8L, 1L)), matrix(c(8, 1), ncol = 1L))
  for (bad in list(0, 9, -1, NA_real_, Inf, 1.5)) {
    expect_error(fmri.pipeline:::extract_nifti_voxel_ts(path, bad), "one-based indices")
  }
})
