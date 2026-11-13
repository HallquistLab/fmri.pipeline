# Internal helpers to check and resolve symbols from the standalone `medusa` package
medusa_is_installed <- function() {
  requireNamespace("medusa", quietly = TRUE)
}

medusa_get_export <- function(name) {
  getExportedValue("medusa", name)
}

# Internal helper to warn and forward deprecated MEDuSA / mixed_by calls to `medusa`
forward_to_medusa <- function(fn_name, ...) {
  .Deprecated(
    new = paste0("medusa::", fn_name),
    package = "fmri.pipeline",
    msg = sprintf(
      "'%s()' has moved to the 'medusa' package and is deprecated in 'fmri.pipeline'. Please use 'medusa::%s()' instead.",
      fn_name, fn_name
    )
  )
  if (!medusa_is_installed()) {
    stop(
      sprintf(
        "'%s()' has moved to the 'medusa' package, which is not installed. Install 'medusa' (e.g., remotes::install_github('HallquistLab/medusa')) and call 'medusa::%s()'.",
        fn_name, fn_name
      ),
      call. = FALSE
    )
  }
  medusa_get_export(fn_name)(...)
}

forward_to_medusa_r6 <- function(...) {
  .Deprecated(
    new = "medusa::fmri_ts",
    package = "fmri.pipeline",
    msg = "'fmri_ts' has moved to the 'medusa' package and is deprecated in 'fmri.pipeline'. Please use 'medusa::fmri_ts$new()' instead."
  )
  if (!medusa_is_installed()) {
    stop(
      "'fmri_ts' has moved to the 'medusa' package, which is not installed. Install 'medusa' (e.g., remotes::install_github('HallquistLab/medusa')) and call 'medusa::fmri_ts$new()'.",
      call. = FALSE
    )
  }
  medusa_get_export("fmri_ts")$new(...)
}

#' Deprecated MEDuSA deconvolution, event-alignment, and multilevel functions
#'
#' @description
#' All MEDuSA deconvolution, event-locked interpolation, and split-wise
#' multilevel modeling functions have moved to the standalone `medusa` package.
#' Calling these legacy entry points from `fmri.pipeline` emits a deprecation
#' warning and forwards to `medusa` when installed.
#'
#' @param ... Arguments forwarded directly to the corresponding function in
#'   `medusa`.
#'
#' @return The return value of the corresponding function in `medusa`.
#' @name medusa-deprecated
#' @aliases fmri.pipeline-deprecated
#' @keywords internal
NULL

#' @rdname medusa-deprecated
#' @export
voxelwise_deconvolution <- function(...) {
  forward_to_medusa("voxelwise_deconvolution", ...)
}

#' @rdname medusa-deprecated
#' @export
deconvolve_nlreg <- function(...) {
  forward_to_medusa("deconvolve_nlreg", ...)
}

#' @rdname medusa-deprecated
#' @export
deconvolve_nlreg_resample <- function(...) {
  forward_to_medusa("deconvolve_nlreg_resample", ...)
}

#' @rdname medusa-deprecated
#' @export
deconvolve_reglin <- function(...) {
  forward_to_medusa("deconvolve_reglin", ...)
}

#' @rdname medusa-deprecated
#' @export
event_lock_ts <- function(...) {
  forward_to_medusa("event_lock_ts", ...)
}

#' @rdname medusa-deprecated
#' @export
get_medusa_interpolated_ts <- function(...) {
  forward_to_medusa("get_medusa_interpolated_ts", ...)
}

#' @rdname medusa-deprecated
#' @export
get_medusa_compression_score <- function(...) {
  forward_to_medusa("get_medusa_compression_score", ...)
}

#' @rdname medusa-deprecated
#' @export
compress_mts_pca <- function(...) {
  forward_to_medusa("compress_mts_pca", ...)
}

#' @rdname medusa-deprecated
#' @export
evt_align_decon_files <- function(...) {
  forward_to_medusa("evt_align_decon_files", ...)
}

#' @rdname medusa-deprecated
#' @export
run_decon_alignment <- function(...) {
  forward_to_medusa("run_decon_alignment", ...)
}

#' @rdname medusa-deprecated
#' @export
mixed_by <- function(...) {
  forward_to_medusa("mixed_by", ...)
}

#' @rdname medusa-deprecated
#' @export
meta_mixed_by <- function(...) {
  forward_to_medusa("meta_mixed_by", ...)
}

#' @rdname medusa-deprecated
#' @export
jlmer <- function(...) {
  forward_to_medusa("jlmer", ...)
}

#' @rdname medusa-deprecated
#' @format An object providing `$new(...)` that emits a deprecation warning and
#'   delegates to `medusa::fmri_ts$new(...)`.
#' @export
fmri_ts <- structure(
  list(
    classname = "fmri_ts",
    new = function(...) {
      forward_to_medusa_r6(...)
    }
  ),
  class = "medusa_deprecated_r6"
)
