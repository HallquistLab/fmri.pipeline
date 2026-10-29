# Header-only geometry lookup: prefer the sform, then the qform, then the
# voxel-size fallback used by RNifti for images without either transform.
roi_image_geometry <- function(image) {
  header <- RNifti::niftiHeader(image)
  list(
    dim = as.integer(header$dim[2:4]),
    volumes = if (header$dim[1L] > 3L) prod(header$dim[5L:(header$dim[1L] + 1L)]) else 1,
    spacing = as.numeric(header$pixdim[2:4]),
    units = bitwAnd(as.integer(header$xyzt_units), 7L),
    affine = matrix(as.numeric(RNifti::xform(header, useQuaternionFirst = FALSE)), 4L, 4L)
  )
}

validate_roi_geometry <- function(mask, image, mask_file, image_file, tolerance = 1e-4) {
  close <- function(x, y) {
    length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
      all(abs(x - y) <= tolerance)
  }
  compatible <- identical(mask$dim, image$dim) && image$volumes == 1 &&
    identical(mask$units, image$units) && close(mask$spacing, image$spacing) &&
    close(mask$affine, image$affine)
  if (!isTRUE(compatible)) {
    stop("ROI mask and statistic image have incompatible spatial geometry (dimensions, voxel sizes, ",
      "units, or voxel-to-world transform). Resample/register the mask to the statistic image grid ",
      "before extraction. Mask: ", mask_file, "; image: ", image_file, call. = FALSE)
  }
  invisible(TRUE)
}
