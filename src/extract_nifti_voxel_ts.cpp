#define RNIFTI_NIFTILIB_VERSION 2

#include <Rcpp.h>
#include "RNifti.h"
#include <cmath>
#include <limits>

//' Extract selected voxel time series in a specified order
//'
//' @param infile Path to a scalar NIfTI image, optionally gzipped.
//' @param voxel_indices One-based linear indices within the spatial volume.
//' @return A voxels by volumes numeric matrix, in the supplied index order.
//' @details Uses the RNifti C++ API, including NIfTI slope/intercept scaling,
//'   without creating a full double-precision R array or intermediate text file.
//'   Spatial alignment with the atlas must be checked by the caller.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix extract_nifti_voxel_ts(std::string infile, Rcpp::NumericVector voxel_indices) {
  RNifti::NiftiImage image(infile, true);
  if (image.isNull()) Rcpp::stop("Could not read NIfTI image: %s", infile.c_str());
  const auto dims = image.dim();
  if (dims.size() < 2) Rcpp::stop("Input image must have at least two spatial dimensions");
  for (size_t d = 4; d < dims.size(); ++d) {
    if (dims[d] != 1) Rcpp::stop("Input must be a scalar 3D or 4D NIfTI image");
  }
  size_t n_spatial = static_cast<size_t>(dims[0]) * dims[1];
  if (dims.size() >= 3) n_spatial *= dims[2];
  const size_t n_time = image.nVoxels() / n_spatial;
  if (n_time > static_cast<size_t>(std::numeric_limits<int>::max()) ||
      voxel_indices.size() > std::numeric_limits<int>::max()) {
    Rcpp::stop("Selected time series exceed R matrix dimensions");
  }

  std::vector<size_t> indices(voxel_indices.size());
  for (R_xlen_t v = 0; v < voxel_indices.size(); ++v) {
    const double index = voxel_indices[v];
    if (!std::isfinite(index) || index != std::floor(index) || index < 1 || index > n_spatial) {
      Rcpp::stop("voxel_indices must be finite, whole, one-based indices within the spatial volume");
    }
    indices[v] = static_cast<size_t>(index - 1);
  }

  const RNifti::NiftiImageData data = image.data();
  if (data.isComplex() || data.isRgb()) Rcpp::stop("Input NIfTI data must be real-valued scalar data");
  if (data.length() != image.nVoxels()) Rcpp::stop("NIfTI image data are unavailable");
  Rcpp::NumericMatrix out(voxel_indices.size(), n_time);
  for (size_t t = 0; t < n_time; ++t) {
    Rcpp::checkUserInterrupt();
    const size_t offset = t * n_spatial;
    for (size_t v = 0; v < indices.size(); ++v) {
      out(v, t) = static_cast<double>(data[offset + indices[v]]);
    }
  }
  return out;
}
