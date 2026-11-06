# Extract selected voxel time series in a specified order

Extract selected voxel time series in a specified order

## Usage

``` r
extract_nifti_voxel_ts(infile, voxel_indices)
```

## Arguments

- infile:

  Path to a scalar NIfTI image, optionally gzipped.

- voxel_indices:

  One-based linear indices within the spatial volume.

## Value

A voxels by volumes numeric matrix, in the supplied index order.

## Details

Uses the RNifti C++ API, including NIfTI slope/intercept scaling,
without creating a full double-precision R array or intermediate text
file. Spatial alignment with the atlas must be checked by the caller.
