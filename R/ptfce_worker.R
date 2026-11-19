ptfce_require_dependencies <- function() {
  if (!requireNamespace("pTFCE", quietly = TRUE) || utils::packageVersion("pTFCE") < "0.4.0") {
    stop("Install HallquistLab/pTFCE (>= 0.4.0) on the execution host for native sided enhancement.", call. = FALSE)
  }
}

ptfce_fwe_threshold <- function(resels, alpha, two_sided) {
  ptfce_require_dependencies()
  checkmate::assert_number(resels, lower = .Machine$double.eps, finite = TRUE)
  checkmate::assert_number(alpha, lower = 1e-10, upper = 1 - 1e-10)
  checkmate::assert_flag(two_sided)
  pTFCE::fwe.p2z(resels, FWEP = if (two_sided) alpha / 2 else alpha)
}

ptfce_threshold_from_table <- function(table, alpha) {
  checkmate::assert_names(names(table), must.include = c("p_value", "two_sided", "z_ptfce", "number_of_resels"))
  sides <- unique(as.logical(table$two_sided))
  checkmate::assert_flag(sides)
  hit <- which(abs(table$p_value - alpha) <= abs(alpha) * 1e-12)
  if (length(hit) > 1L) stop("Duplicate pTFCE alpha values in threshold table.")
  threshold <- if (length(hit)) table$z_ptfce[hit] else {
    ptfce_fwe_threshold(table$number_of_resels[1L], alpha, sides)
  }
  checkmate::assert_number(threshold, lower = 0, finite = TRUE)
  list(z = threshold, two_sided = sides)
}

ptfce_read_mask <- function(file, image_file) {
  geometry <- roi_image_geometry(file)
  if (geometry$volumes != 1L) stop("FWE mask must be a 3D image: ", file)
  validate_roi_geometry(geometry, roi_image_geometry(image_file), file, image_file)
  values <- as.array(RNifti::readNifti(file))
  if (any(!is.finite(values)) || any(!values %in% c(0, 1)) || !any(values == 1)) {
    stop("FWE mask must be nonempty, finite, and binary: ", file)
  }
  values == 1
}

ptfce_smoothness <- function(file, mask) {
  values <- utils::read.table(file, col.names = c("key", "value"), nrows = 3L)
  if (!setequal(values$key, c("DLH", "VOLUME", "RESELS"))) stop("Invalid FSL smoothness file: ", file)
  values <- stats::setNames(values$value, values$key)
  if (any(!is.finite(values)) || any(values <= 0)) stop("Invalid FSL smoothness values: ", file)
  # FSL RESELS is voxels per resel, not the number of resels. Reuse the
  # model's stationary smoothness estimate, but use the actual search volume.
  volume <- sum(mask)
  list(V = volume, Rd = unname(values["DLH"]) * volume, resels = unname(values["RESELS"]))
}

run_ptfce_image <- function(zstat, mask, out_dir = dirname(zstat), fwe_alpha = 0.05,
                            two_sided = TRUE, write_threshold_images = TRUE,
                            fsl_smoothest = NULL, model_mask = NULL,
                            residuals = NULL, dof = NULL, verbose = FALSE) {
  ptfce_require_dependencies()
  checkmate::assert_file_exists(zstat)
  checkmate::assert_file_exists(mask)
  checkmate::assert_numeric(fwe_alpha, lower = 1e-10, upper = 1 - 1e-10,
    any.missing = FALSE, unique = TRUE, min.len = 1L)
  checkmate::assert_flag(two_sided)
  checkmate::assert_flag(write_threshold_images)
  inside <- ptfce_read_mask(mask, zstat)
  if (!is.null(model_mask)) {
    model_inside <- ptfce_read_mask(model_mask, zstat)
    if (any(inside & !model_inside)) stop("Correction mask extends outside the model mask.")
  }
  Z <- RNifti::readNifti(zstat)
  if (any(!is.finite(Z[inside]))) stop("Non-finite Z statistics inside the correction mask.")
  smoothness <- list()
  if (!is.null(residuals)) {
    checkmate::assert_file_exists(residuals)
    checkmate::assert_number(dof, lower = 1, finite = TRUE)
    residual_geometry <- roi_image_geometry(residuals)
    residual_geometry$volumes <- 1L
    validate_roi_geometry(roi_image_geometry(mask), residual_geometry, mask, residuals)
    smoothness <- list(residual = residuals, dof = dof)
  } else if (!is.null(fsl_smoothest)) {
    smoothness <- ptfce_smoothness(fsl_smoothest, inside)
  }
  result <- do.call(pTFCE::ptfce, c(list(img = Z, mask = inside, verbose = verbose,
    sided = if (two_sided) "two.sided" else "positive", fwe_p = fwe_alpha), smoothness))
  output <- result$Z
  if (any(!is.finite(output[inside]))) stop("pTFCE produced non-finite values inside the correction mask.")
  checkmate::assert_numeric(result$fwe.Z, len = length(fwe_alpha), lower = 0, finite = TRUE, any.missing = FALSE)
  if (!two_sided) output[output < 0] <- 0
  output[!inside] <- 0
  parts <- fwe_nifti_parts(zstat)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(out_dir, paste0(parts$stem, "_ptfce"))
  labels <- vapply(fwe_alpha, format_fwe_probability, character(1L))
  if (anyDuplicated(labels)) stop("FWE alpha values are too close to give distinct output names.")
  table_file <- paste0(base, "_zthresh.csv")
  if (file.exists(table_file) && !file.remove(table_file)) stop("Cannot replace ", table_file)
  thresholds <- lapply(seq_along(fwe_alpha), function(i) {
    z <- unname(result$fwe.Z[i])
    significant <- inside & abs(output) >= z & output != 0
    if (write_threshold_images) {
      thresholded <- output
      thresholded[!significant] <- 0
      RNifti::writeNifti(thresholded, paste0(base, "_fwep_", labels[i], parts$extension), datatype = "float")
    }
    data.frame(z_ptfce = z, p_value = fwe_alpha[i], two_sided = two_sided,
      number_of_resels = result$number_of_resels, search_voxels = sum(inside),
      significant_voxels = sum(significant), ptfce_version = as.character(utils::packageVersion("pTFCE")),
      write_threshold_images = write_threshold_images)
  })
  RNifti::writeNifti(output, paste0(base, parts$extension), datatype = "float")
  utils::write.csv(do.call(rbind, thresholds), table_file, row.names = FALSE)
  invisible(do.call(rbind, thresholds))
}

ptfce_cli <- function(args) {
  if (!length(args) || "--help" %in% args) {
    cat("Usage: ptfce_zstat.R --zstat image --mask mask [--model_mask mask] [--fsl_smoothest file]\n",
      "  [--residuals image --dof number_or_file] [--out_dir directory]\n",
      "  [--fwep alpha ...] [--twosided | --onesided] [--write_thresh_imgs] [--verbose]\n")
    cat("   or: ptfce_zstat.R --gfeat directory [--fwep alpha ...] [--twosided | --onesided]",
      "[--write_thresh_imgs] [--force]\n")
    return(invisible(NULL))
  }
  options <- list(two_sided = TRUE, write_threshold_images = FALSE)
  flags <- c("--twosided" = "two_sided", "--onesided" = "two_sided",
    "--write_thresh_imgs" = "write_threshold_images", "--verbose" = "verbose")
  values <- c("--zstat" = "zstat", "--mask" = "mask", "--model_mask" = "model_mask",
    "--out_dir" = "out_dir", "--fsl_smoothest" = "fsl_smoothest", "--residuals" = "residuals", "--dof" = "dof",
    "--gfeat" = "gfeat_dir")
  force <- "--force" %in% args
  args <- args[args != "--force"]
  i <- 1L
  while (i <= length(args)) {
    flag <- args[i]
    if (flag %in% names(flags)) {
      options[[flags[[flag]]]] <- flag != "--onesided"
      i <- i + 1L
    } else if (flag == "--fwep") {
      end <- i + 1L
      while (end <= length(args) && !startsWith(args[end], "--")) end <- end + 1L
      if (end == i + 1L) stop("Missing values for --fwep")
      options$fwe_alpha <- as.numeric(args[seq.int(i + 1L, end - 1L)])
      i <- end
    } else if (flag %in% names(values)) {
      if (i == length(args) || startsWith(args[i + 1L], "--")) stop("Missing value for ", flag)
      options[[values[[flag]]]] <- args[i + 1L]
      i <- i + 2L
    } else stop("Unknown argument: ", flag)
  }
  if (!is.null(options$gfeat_dir)) {
    allowed <- c("gfeat_dir", "fwe_alpha", "two_sided", "write_threshold_images")
    if (!all(names(options) %in% allowed)) stop("Use --gfeat with directory-mode options only (see --help).")
    spec <- ptfce_spec$new(gfeat_dir = options$gfeat_dir, fwe_p = options$fwe_alpha %||% 0.05,
      two_sided = options$two_sided, write_thresh_imgs = options$write_threshold_images)
    spec$run(force = force)
    return(invisible(NULL))
  }
  if (force) stop("--force applies only with --gfeat.")
  if (!all(c("zstat", "mask") %in% names(options))) stop("--zstat and --mask are required.")
  if (!is.null(options$dof)) {
    if (file.exists(options$dof)) options$dof <- readLines(options$dof, warn = FALSE)
    options$dof <- as.numeric(options$dof)
  }
  do.call(run_ptfce_image, options)
}
