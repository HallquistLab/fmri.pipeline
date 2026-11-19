# AFNI time-series estimation uses the pipeline's convolved regressors. AFNI
# supplies run-specific polynomial baselines and REML temporal prewhitening.
afni_glm_defaults <- function() {
  list(estimation = "REML", polort = "A", mask = NULL, censor_column = "afni_censor_file",
       write_residuals = TRUE, force_l1_creation = FALSE, force_l3_creation = FALSE)
}

afni_glm_settings <- function(gpa) {
  x <- populate_defaults(gpa$glm_settings$afni, afni_glm_defaults())
  x$estimation <- toupper(x$estimation)
  checkmate::assert_choice(x$estimation, c("REML", "OLS"))
  if (!identical(x$polort, "A")) checkmate::assert_integerish(x$polort, lower = 0, len = 1L)
  checkmate::assert_string(x$mask, null.ok = TRUE)
  checkmate::assert_string(x$censor_column)
  checkmate::assert_flag(x$write_residuals)
  checkmate::assert_flag(x$force_l1_creation)
  checkmate::assert_flag(x$force_l3_creation)
  x
}

# Align an auxiliary time series with a BDM run. Accept either the retained
# samples or the original samples with exactly drop_volumes leading rows.
afni_align_rows <- function(x, n, drop, what) {
  x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x))) stop(what, " must contain finite numeric values.")
  if (nrow(x) == n) return(x)
  if (drop > 0L && nrow(x) == n + drop) return(x[seq_len(n) + drop, , drop = FALSE])
  stop(what, " has ", nrow(x), " rows; expected ", n, " retained rows or ", n + drop, " original rows.")
}

afni_write_1d <- function(x, path) {
  utils::write.table(x, path, row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Build task and block-diagonal nuisance matrices in precisely the BDM run order.
afni_l1_design <- function(d_obj, regressors, confound_files = NULL, censor_files = NULL) {
  nr <- length(d_obj$runs_to_output)
  n <- as.integer(d_obj$run_volumes)
  drops <- rep_len(d_obj$drop_volumes %||% 0L, nr)
  checkmate::assert_integerish(n, lower = 1L, len = nr, any.missing = FALSE)
  checkmate::assert_integerish(drops, lower = 0L, len = nr, any.missing = FALSE)
  checkmate::assert_character(regressors, min.len = 1L, any.missing = FALSE, unique = TRUE)
  if (is.null(confound_files)) confound_files <- rep(NA_character_, nr)
  if (is.null(censor_files)) censor_files <- rep(NA_character_, nr)
  stopifnot(length(confound_files) == nr, length(censor_files) == nr)
  starts <- c(0L, head(cumsum(n), -1L))
  task <- matrix(0, sum(n), length(regressors), dimnames = list(NULL, regressors))
  nuisance <- matrix(numeric(0), sum(n), 0L)
  censor <- rep(1L, sum(n))
  for (i in seq_len(nr)) {
    key <- paste0("run_number", d_obj$runs_to_output[i])
    r <- d_obj$design_convolved[[key]]
    if (is.null(r)) stop("Missing convolved design for ", key)
    if (!all(regressors %in% colnames(r))) stop("Missing task regressors in ", key)
    rows <- starts[i] + seq_len(n[i])
    task[rows, ] <- afni_align_rows(r[, regressors, drop = FALSE], n[i], 0L, key)
    cf <- censor_files[i]
    if (!is.na(cf) && nzchar(cf)) {
      checkmate::assert_file_exists(cf)
      cen <- afni_align_rows(data.table::fread(cf, header = FALSE, data.table = FALSE),
                            n[i], drops[i], paste("Censor for", key))
      if (ncol(cen) != 1L || any(!cen %in% c(0, 1))) stop("AFNI censor files must have one binary column (1=keep).")
      censor[rows] <- cen[, 1L]
    }
    if (!any(censor[rows] == 1L)) stop("All volumes are censored in ", key)
    cf <- confound_files[i]
    if (!is.na(cf) && nzchar(cf)) {
      checkmate::assert_file_exists(cf)
      conf <- afni_align_rows(data.table::fread(cf, header = FALSE, data.table = FALSE),
                             n[i], drops[i], paste("Confounds for", key))
      # Constants are already represented by AFNI's run baseline.
      # Center on observations retained by censoring. A spike for an already
      # censored volume otherwise becomes a duplicate of the run baseline.
      observed <- censor[rows] == 1L
      conf <- sweep(conf, 2L, colMeans(conf[observed, , drop = FALSE]), "-")
      conf <- conf[, colSums(abs(conf[observed, , drop = FALSE])) > 1e-10, drop = FALSE]
      block <- matrix(0, sum(n), ncol(conf))
      block[rows, ] <- conf
      if (ncol(conf)) colnames(block) <- paste0("run", d_obj$runs_to_output[i], "_nuis", seq_len(ncol(conf)))
      nuisance <- cbind(nuisance, block)
    }
  }
  list(design = cbind(task, nuisance), n_task = ncol(task), censor = censor, starts = starts)
}

#' Set up an AFNI first-level GLM
#'
#' Runs within a subject/session share task coefficients and have independent
#' polynomial baselines and nuisance coefficients. The pipeline's convolved
#' design is passed to 3dDeconvolve; 3dREMLfit is the default estimator.
#' @param id,session Subject and session identifiers.
#' @param l1_confound_files Optional nuisance files, one per BDM run.
#' @param d_obj A build_design_matrix (bdm) object.
#' @param gpa A glm_pipeline_arguments object.
#' @param model_name Name of the L1 model.
#' @param run_nifti Optional input images when absent from d_obj.
#' @param censor_files Optional binary AFNI censor files, one per BDM run.
#' @return A one-row AFNI setup/status data.frame.
#' @export
afni_l1_model <- function(id, session, l1_confound_files = NULL, d_obj, gpa,
                          model_name, run_nifti = NULL, censor_files = NULL) {
  checkmate::assert_class(d_obj, "bdm")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_scalar(id)
  checkmate::assert_scalar(session)
  checkmate::assert_string(model_name)
  settings <- afni_glm_settings(gpa)
  model <- gpa$l1_models$models[[model_name]]
  if (is.null(model)) stop("Unknown AFNI L1 model: ", model_name)
  files <- d_obj$run_nifti %||% run_nifti
  checkmate::assert_character(files, min.len = 1L, any.missing = FALSE)
  checkmate::assert_character(l1_confound_files, len = length(files), null.ok = TRUE)
  checkmate::assert_character(censor_files, len = length(files), null.ok = TRUE)
  metadata <- gpa$run_data[gpa$run_data$id == id & gpa$run_data$session == session, , drop = FALSE]
  matched <- match(d_obj$runs_to_output, metadata$run_number)
  if (anyNA(matched)) stop("AFNI BDM runs do not match gpa$run_data.")
  metadata <- metadata[matched, , drop = FALSE]
  keep <- !metadata$exclude_run
  if ("exclude_subject" %in% names(metadata)) keep <- keep & !metadata$exclude_subject
  subjects <- gpa$subject_data[gpa$subject_data$id == id & gpa$subject_data$session == session, , drop = FALSE]
  if ("exclude_subject" %in% names(subjects) && any(subjects$exclude_subject %in% TRUE)) return(NULL)
  keep <- which(keep %in% TRUE)
  if (!length(keep)) return(NULL)
  d_obj$drop_volumes <- rep_len(d_obj$drop_volumes %||% 0L, length(files))[keep]
  d_obj$runs_to_output <- d_obj$runs_to_output[keep]
  d_obj$run_volumes <- d_obj$run_volumes[keep]
  files <- files[keep]
  if (!is.null(l1_confound_files)) l1_confound_files <- l1_confound_files[keep]
  if (!is.null(censor_files)) censor_files <- censor_files[keep]
  checkmate::assert_file_exists(files)
  if (length(files) != length(d_obj$run_volumes)) stop("AFNI inputs and BDM runs must align.")
  files <- normalizePath(files, mustWork = TRUE)
  inputs <- files
  drops <- rep_len(d_obj$drop_volumes %||% 0L, length(files))
  for (i in seq_along(files)) {
    dims <- RNifti::niftiHeader(files[i])$dim
    actual <- if (dims[1L] >= 4L) dims[5L] else 1L
    n <- d_obj$run_volumes[i]
    if (actual != n) {
      if (drops[i] > 0L && actual == n + drops[i]) {
        inputs[i] <- paste0(files[i], "[", drops[i], "..", actual - 1L, "]")
      } else stop("AFNI image length disagrees with BDM for ", files[i])
    }
  }
  design <- afni_l1_design(d_obj, model$regressors, l1_confound_files, censor_files)
  outdir <- get_output_directory(id = id, session = session, l1_model = model_name,
                                 gpa = gpa, glm_software = "afni", create_if_missing = TRUE)
  outdir <- normalizePath(outdir, mustWork = TRUE)
  generated <- afni_write_glm(outdir, inputs, design$design, model$contrasts,
    n_task = design$n_task, censor = design$censor, starts = design$starts,
    polort = settings$polort, tr = d_obj$tr, estimation = settings$estimation,
    mask = settings$mask, write_residuals = settings$write_residuals,
    cpus = gpa$parallel$afni$l1_glm_cpus %||% 1L,
    force = settings$force_l1_creation, input_files = c(files, l1_confound_files, censor_files))
  data.frame(id = id, session = session, l1_model = model_name,
    afni_dir = outdir, afni_script = generated$script, output_file = generated$output_file,
    get_afni_glm_status(outdir), stringsAsFactors = FALSE)
}

# Run these checks inside the configured AFNI environment (including modules),
# rather than checking the login shell's PATH.
test_afni_glm_compute_environment <- function(gpa) {
  required <- c("3dDeconvolve", "3dTcat", "3dcopy", "3dbucket", "3dAutomask", "3dmask_tool")
  if (afni_glm_settings(gpa)$estimation == "REML") required <- c(required, "3dREMLfit")
  script <- tempfile(fileext = ".sh")
  on.exit(unlink(script), add = TRUE)
  writeLines(c("#!/bin/bash", "set -e", get_compute_environment(gpa, "afni"),
    sprintf("command -v %s >/dev/null || { echo 'Required AFNI executable missing: %s'; exit 1; }", required, required)), script)
  output <- system2("bash", shQuote(script), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(output, "status")) && attr(output, "status") != 0L) {
    stop("AFNI GLM compute environment failed: ", paste(output, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}
