#' Refresh GLM backend status for a given analysis level
#'
#' @param gpa a \code{glm_pipeline_arguments} object
#' @param level the level of analysis to be refreshed (1, 2, or 3)
#' @param lg an optional lgr logger object used for logging
#' @param glm_software optional backend(s) to refresh. If NULL, refreshes all
#'   backends in the gpa object.
#' @return a modified copy of \code{gpa} with backend status columns refreshed
#'
#' @export
#' @importFrom dplyr select all_of
#' @importFrom purrr pmap_dfr
refresh_glm_status <- function(gpa, level = 1L, lg = NULL, glm_software = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)
  checkmate::assert_character(glm_software, null.ok = TRUE)

  setup_name <- paste0("l", level, "_model_setup")
  if (!setup_name %in% names(gpa)) {
    return(gpa)
  }

  if (is.null(lg)) lg <- lgr::get_logger()

  gpa_for_backends <- gpa
  if (!is.null(glm_software)) {
    gpa_for_backends$glm_software <- unique(tolower(glm_software))
    level_backends <- normalize_level_backend_config(gpa_for_backends)
    level_backends[[paste0("l", level)]] <- unique(tolower(glm_software))
    gpa_for_backends$level_backends <- level_backends
  }

  glm_backends <- get_glm_backends(gpa_for_backends, must_exist = FALSE, level = level)
  if (length(glm_backends) == 0L) {
    lg$warn("No GLM backends registered for status refresh.")
    return(gpa)
  }

  status_fn_name <- paste0("l", level, "_status")
  status_inputs_name <- paste0("l", level, "_status_inputs")

  for (backend_name in names(glm_backends)) {
    backend <- glm_backends[[backend_name]]
    if (is.null(backend)) next

    backend_tables <- gpa[[setup_name]]
    if (is.null(backend_tables) || !backend_name %in% names(backend_tables)) {
      lg$debug("No %s table for backend '%s' in gpa$%s. Skipping.", backend_name, setup_name)
      next
    }

    orig <- backend_tables[[backend_name]]
    if (is.null(orig) || (is.data.frame(orig) && nrow(orig) == 0L)) {
      lg$warn("Could not find populated $%s object in gpa$%s$%s", backend_name, setup_name, backend_name)
      next
    }

    status_fn <- backend[[status_fn_name]]
    if (is.null(status_fn)) {
      lg$debug("No %s defined for backend '%s'. Skipping refresh.", status_fn_name, backend_name)
      next
    }

    if (isTRUE(attr(status_fn, "glm_backend_not_implemented"))) {
      lg$debug("Backend '%s' does not implement %s. Skipping refresh.", backend_name, status_fn_name)
      next
    }

    status_inputs <- backend[[status_inputs_name]]
    if (is.null(status_inputs) || length(status_inputs) == 0L) {
      lg$warn("No %s defined for backend '%s'. Skipping refresh.", status_inputs_name, backend_name)
      next
    }

    status_inputs <- unique(status_inputs)
    missing_cols <- setdiff(status_inputs, names(orig))
    if (length(missing_cols) > 0L) {
      lg$warn(
        "Missing required status inputs for backend '%s': %s",
        backend_name, paste(missing_cols, collapse = ", ")
      )
      next
    }

    lg$info(
      "Found existing %s field. Refreshing status of L%d %s execution and outputs.",
      setup_name, level, toupper(backend_name)
    )

    refresh <- orig %>%
      dplyr::select(dplyr::all_of(status_inputs)) %>%
      purrr::pmap_dfr(status_fn, lg = lg)

    gpa[[setup_name]][[backend_name]][, names(refresh)] <- refresh
  }

  return(gpa)
}

#' Helper to check whether expected SPM outputs exist
#'
#' @param spm_dir Directory containing SPM outputs
#' @param lg Optional logger
#' @param prefix Optional prefix for column names
#' @return A data.frame with SPM status columns including completion status,
#'   timing information, and output file existence checks.
#' @export
#' @importFrom anytime anytime
get_spm_status <- function(spm_dir, lg = NULL, prefix = NULL) {
  checkmate::assert_string(spm_dir)
  if (is.null(lg)) lg <- lgr::get_logger()

  spm_checks <- list()
  spm_checks$spm_dir_exists <- dir.exists(spm_dir)
  spm_checks$spm_execution_start <- as.POSIXct(NA)
  spm_checks$spm_execution_end <- as.POSIXct(NA)
  spm_checks$spm_execution_min <- NA_real_
  spm_checks$spm_complete <- FALSE
  spm_checks$spm_failed <- NA

  spm_complete_file <- file.path(spm_dir, ".spm_complete")
  spm_fail_file <- file.path(spm_dir, ".spm_fail")
  spm_checks$spm_complete_file <- spm_complete_file
  spm_checks$spm_fail_file <- spm_fail_file

  spm_mat <- file.path(spm_dir, "SPM.mat")
  spm_mat_exists <- file.exists(spm_mat)
  spm_checks$spm_mat <- spm_mat
  spm_checks$spm_mat_exists <- spm_mat_exists
  spm_checks$spm_mat_modified_date <- if (spm_mat_exists) file.info(spm_mat)$mtime else as.POSIXct(NA)

  beta_files <- list.files(spm_dir, pattern = "^beta_\\d+\\.nii(\\.gz)?$", full.names = TRUE)
  con_files <- list.files(spm_dir, pattern = "^con_\\d+\\.nii(\\.gz)?$", full.names = TRUE)
  spm_checks$spm_beta_exists <- length(beta_files) > 0L
  spm_checks$spm_beta_count <- length(beta_files)
  spm_checks$spm_contrast_exists <- length(con_files) > 0L
  spm_checks$spm_contrast_count <- length(con_files)
  spm_checks$spm_minimal_outputs_exist <- spm_nonempty_files(c(spm_mat, beta_files)) && length(beta_files) > 0L
  spm_checks$spm_outputs_valid <- spm_outputs_valid(spm_dir)
  spm_checks$spm_complete <- file.exists(spm_complete_file) && !file.exists(spm_fail_file) &&
    isTRUE(spm_checks$spm_outputs_valid)
  spm_checks$spm_failed <- file.exists(spm_fail_file) ||
    (file.exists(spm_complete_file) && !isTRUE(spm_checks$spm_outputs_valid))

  if (dir.exists(spm_dir)) {
    timing_file <- if (file.exists(spm_fail_file)) spm_fail_file else if (file.exists(spm_complete_file)) spm_complete_file else NULL

    if (!is.null(timing_file)) {
      timing <- readLines(timing_file)
      if (length(timing) > 0L) {
        timing <- anytime::anytime(timing)
        spm_checks$spm_execution_start <- timing[1L]
        if (length(timing) == 2L) {
          spm_checks$spm_execution_end <- timing[2L]
          spm_checks$spm_execution_min <- as.numeric(difftime(timing[2L], timing[1L], units = "mins"))
        } else {
          lg$warn("Did not find two timing entries in %s.", timing_file)
          lg$warn("File contents: %s", timing)
        }
      }
    }
  }

  df <- as.data.frame(spm_checks)
  if (!is.null(prefix)) {
    names(df) <- paste0(prefix, names(df))
  }
  return(df)
}
