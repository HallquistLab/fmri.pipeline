# Internal helpers for projecting L2 covariates into SPM L1 design matrices.
resolve_spm_l1_session_mode <- function(mode, lg = NULL) {
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (is.null(mode) || (is.character(mode) && length(mode) == 1L && !nzchar(mode))) {
    return("separate")
  }

  if (!is.character(mode) || length(mode) != 1L) {
    stop("spm$l1_session_mode must be NULL or a single character value.", call. = FALSE)
  }

  normalized <- tolower(trimws(mode))
  if (normalized %in% c("separate", "per_session", "id_session", "session")) {
    return("separate")
  }
  if (normalized %in% c("pooled", "id", "across_sessions", "subject")) {
    return("pooled")
  }

  stop("Unknown spm$l1_session_mode: ", mode, ". Allowed: separate, pooled", call. = FALSE)
}

resolve_spm_l2_projection_model <- function(gpa, l1_model, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_string(l1_model)
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  available <- character(0)
  if (checkmate::test_class(gpa$l2_models, "hi_model_set") &&
      !is.null(gpa$l2_models$models) &&
      length(gpa$l2_models$models) > 0L) {
    available <- names(gpa$l2_models$models)
  }

  configured <- gpa$glm_settings$spm$l2_projection_model

  if (length(available) == 0L) {
    if (!is.null(configured) &&
        is.character(configured) &&
        length(configured) == 1L &&
        nzchar(configured) &&
        !tolower(configured) %in% c("none", "off", "false")) {
      lg$warn(
        "SPM l2_projection_model is set to '%s' but no L2 models are defined. Skipping L2 projection.",
        configured
      )
    }
    return(NULL)
  }

  if (is.null(configured) ||
      (is.character(configured) && length(configured) == 1L && !nzchar(configured))) {
    if (l1_model %in% available) {
      lg$info(
        "Using L2 model '%s' for SPM L1 projection (matched by L1 model name).",
        l1_model
      )
      return(l1_model)
    }
    if (length(available) == 1L) return(available[1L])

    chosen <- available[1L]
    lg$warn(
      "Multiple L2 models available for SPM L1 projection (%s). Using '%s'. Set glm_settings$spm$l2_projection_model to choose explicitly.",
      paste(available, collapse = ", "),
      chosen
    )
    return(chosen)
  }

  if (!is.character(configured) || length(configured) != 1L) {
    stop(
      "glm_settings$spm$l2_projection_model must be NULL or a single model name.",
      call. = FALSE
    )
  }

  if (tolower(configured) %in% c("none", "off", "false")) {
    return(NULL)
  }

  if (!configured %in% available) {
    stop(
      "glm_settings$spm$l2_projection_model='", configured,
      "' is not available. Choose one of: ",
      paste(available, collapse = ", "),
      call. = FALSE
    )
  }

  configured
}

resolve_spm_l2_projection_interactions <- function(interactions, available_terms, lg = NULL) {
  checkmate::assert_character(available_terms, any.missing = FALSE, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (length(available_terms) == 0L) return(character(0))

  if (is.null(interactions)) return(character(0))

  if (is.logical(interactions)) {
    if (length(interactions) != 1L || is.na(interactions)) {
      stop("spm$l2_projection_interactions must be a scalar logical, character vector, or NULL.", call. = FALSE)
    }
    return(if (isTRUE(interactions)) available_terms else character(0))
  }

  if (!is.character(interactions)) {
    stop("spm$l2_projection_interactions must be a scalar logical, character vector, or NULL.", call. = FALSE)
  }

  interactions <- interactions[!is.na(interactions) & nzchar(interactions)]
  if (length(interactions) == 0L) return(character(0))

  lower <- tolower(interactions)
  if (length(lower) == 1L && lower %in% c("none", "off", "false", "no", "0")) {
    return(character(0))
  }
  if (any(lower %in% c("all", "true", "yes", "1"))) {
    return(available_terms)
  }

  bad <- setdiff(interactions, available_terms)
  if (length(bad) > 0L) {
    stop(
      "Unknown spm$l2_projection_interactions term(s): ",
      paste(bad, collapse = ", "),
      ". Available projected L2 terms: ",
      paste(available_terms, collapse = ", "),
      call. = FALSE
    )
  }

  unique(interactions)
}

resolve_spm_l2_projection_interaction_contrast_modes <- function(modes, legacy_unit_contrasts = TRUE, lg = NULL) {
  checkmate::assert_flag(legacy_unit_contrasts)
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  allowed <- c("pooled", "session_specific", "session_differences")

  if (is.null(modes)) {
    return(if (isTRUE(legacy_unit_contrasts)) "pooled" else character(0))
  }

  if (is.logical(modes)) {
    if (length(modes) != 1L || is.na(modes)) {
      stop(
        "spm$l2_projection_interaction_contrast_modes must be scalar logical, character vector, or NULL.",
        call. = FALSE
      )
    }
    return(if (isTRUE(modes)) "pooled" else character(0))
  }

  if (!is.character(modes)) {
    stop(
      "spm$l2_projection_interaction_contrast_modes must be scalar logical, character vector, or NULL.",
      call. = FALSE
    )
  }

  modes <- modes[!is.na(modes) & nzchar(modes)]
  if (length(modes) == 0L) return(character(0))

  lower <- tolower(modes)
  if (length(lower) == 1L && lower %in% c("none", "off", "false", "no", "0")) {
    return(character(0))
  }
  if (any(lower %in% c("all", "true", "yes", "1"))) {
    return(allowed)
  }

  normalized <- vapply(lower, function(x) {
    if (x %in% c("pooled", "condition_means", "condition-means", "main_effect")) return("pooled")
    if (x %in% c("session_specific", "session-specific", "cell_means", "cell-means", "by_session")) return("session_specific")
    if (x %in% c("session_difference", "session-difference", "session_differences", "session-differences", "pairwise_diff", "pairwise_differences")) return("session_differences")
    x
  }, character(1))

  bad <- setdiff(normalized, allowed)
  if (length(bad) > 0L) {
    stop(
      "Unknown spm$l2_projection_interaction_contrast_modes value(s): ",
      paste(bad, collapse = ", "),
      ". Allowed: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }

  unique(normalized)
}

augment_spm_l1_contrasts_for_projection <- function(mobj, interaction_terms, contrast_modes = "pooled", lg = NULL) {
  checkmate::assert_list(mobj)
  checkmate::assert_character(interaction_terms, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_character(contrast_modes, any.missing = FALSE, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (length(interaction_terms) == 0L) return(mobj)

  if (is.null(mobj$contrasts) || !inherits(mobj$contrasts, "matrix") || nrow(mobj$contrasts) == 0L) {
    return(mobj)
  }

  cmat <- mobj$contrasts
  existing_cols <- colnames(cmat)
  if (is.null(existing_cols)) {
    stop("L1 contrast matrix must have column names before adding projected interaction terms.", call. = FALSE)
  }

  missing_cols <- setdiff(interaction_terms, existing_cols)
  if (length(missing_cols) > 0L) {
    zeros <- matrix(0, nrow = nrow(cmat), ncol = length(missing_cols))
    colnames(zeros) <- missing_cols
    cmat <- cbind(cmat, zeros)
    lg$info(
      "Added %d projected interaction regressor column(s) to SPM contrast matrix.",
      length(missing_cols)
    )
  }

  if ("pooled" %in% contrast_modes) {
    for (tt in interaction_terms) {
      row_name <- paste0("proj_int_", tt)
      if (!is.null(rownames(cmat)) && row_name %in% rownames(cmat)) next
      add_row <- rep(0, ncol(cmat))
      names(add_row) <- colnames(cmat)
      add_row[[tt]] <- 1
      cmat <- rbind(cmat, add_row)
      rownames(cmat)[nrow(cmat)] <- row_name
    }
  }

  mobj$contrasts <- cmat
  mobj
}

extract_spm_l2_projection <- function(gpa, l2_model_name, id, session, run_numbers,
                                      run_sessions = NULL, source_run_numbers = NULL,
                                      l1_session_mode = "separate", lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_string(l2_model_name)
  checkmate::assert_scalar(id)
  checkmate::assert_scalar(session, null.ok = TRUE)
  checkmate::assert_integerish(run_numbers, lower = 1L, min.len = 1L, any.missing = FALSE)
  checkmate::assert_integerish(run_sessions, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_integerish(source_run_numbers, null.ok = TRUE, any.missing = FALSE, lower = 1L)
  checkmate::assert_string(l1_session_mode)
  if (is.null(lg)) lg <- lgr::get_logger()
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  l1_session_mode <- resolve_spm_l1_session_mode(l1_session_mode, lg = lg)
  run_numbers <- as.integer(run_numbers)
  if (is.null(source_run_numbers)) source_run_numbers <- run_numbers
  source_run_numbers <- as.integer(source_run_numbers)
  if (is.null(run_sessions)) {
    if (is.null(session)) {
      stop("session must be provided when run_sessions is NULL.", call. = FALSE)
    }
    run_sessions <- rep.int(as.integer(session), length(run_numbers))
  }
  run_sessions <- as.integer(run_sessions)

  if (length(source_run_numbers) != length(run_numbers) || length(run_sessions) != length(run_numbers)) {
    stop("run_numbers, source_run_numbers, and run_sessions must have the same length.", call. = FALSE)
  }

  mobj <- gpa$l2_models$models[[l2_model_name]]
  if (is.null(mobj)) {
    stop("Could not find requested L2 model for SPM projection: ", l2_model_name, call. = FALSE)
  }

  if (is.null(mobj$model_matrix) || is.null(mobj$metadata)) {
    lg$warn(
      "L2 model '%s' does not include model_matrix/metadata. Skipping SPM L2 projection.",
      l2_model_name
    )
    return(NULL)
  }

  if (!is.data.frame(mobj$metadata)) {
    stop("L2 model metadata must be a data.frame for SPM projection.", call. = FALSE)
  }

  if (!all(c("id", "session", "run_number") %in% names(mobj$metadata))) {
    stop(
      "L2 model metadata must include id, session, and run_number columns for SPM projection.",
      call. = FALSE
    )
  }

  mm <- as.data.frame(mobj$model_matrix, stringsAsFactors = FALSE)
  if (nrow(mm) != nrow(mobj$metadata)) {
    stop(
      "L2 model metadata rows do not match model_matrix rows for SPM projection model '",
      l2_model_name, "'.",
      call. = FALSE
    )
  }

  regressor_cols <- setdiff(colnames(mm), "(Intercept)")
  if (length(regressor_cols) == 0L) {
    lg$info(
      "L2 model '%s' has only an intercept. No L2 regressors will be projected into SPM L1.",
      l2_model_name
    )
    return(NULL)
  }

  scope <- mobj$l2_scope
  if (is.null(scope) || !is.character(scope) || length(scope) != 1L || !nzchar(scope)) {
    scope <- "id_session"
  }
  checkmate::assert_subset(scope, longitudinal_l2_scopes())

  projection_df <- data.frame(
    mobj$metadata[, c("id", "session", "run_number"), drop = FALSE],
    mm,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (identical(l1_session_mode, "pooled")) {
    projection_df <- projection_df[projection_df$id == id, , drop = FALSE]
  } else {
    # Current SPM L1 setup is session-specific.
    projection_df <- projection_df[
      projection_df$id == id & projection_df$session == session,
      ,
      drop = FALSE
    ]
  }

  if (nrow(projection_df) == 0L) {
    lg$warn(
      "No L2 projection rows found for id=%s session=%s in model '%s'. Skipping projection.",
      as.character(id), as.character(session), l2_model_name
    )
    return(NULL)
  }

  if (identical(l1_session_mode, "pooled")) {
    if (anyDuplicated(paste(projection_df$session, projection_df$run_number, sep = "::")) > 0L) {
      stop(
        "Duplicate session/run_number rows detected in L2 projection data for model '",
        l2_model_name, "', id=", as.character(id), ".",
        call. = FALSE
      )
    }
  } else if (anyDuplicated(projection_df$run_number) > 0L) {
    stop(
      "Duplicate run_number rows detected in L2 projection data for model '",
      l2_model_name, "', id=", as.character(id),
      ", session=", as.character(session), ".",
      call. = FALSE
    )
  }

  row_index <- if (identical(l1_session_mode, "pooled")) {
    req_key <- paste(run_sessions, source_run_numbers, sep = "::")
    src_key <- paste(projection_df$session, projection_df$run_number, sep = "::")
    match(req_key, src_key)
  } else {
    match(source_run_numbers, projection_df$run_number)
  }
  if (anyNA(row_index)) {
    missing_runs <- paste0(run_sessions[is.na(row_index)], ":", source_run_numbers[is.na(row_index)])
    stop(
      "Cannot project L2 model '", l2_model_name, "' into SPM L1 for id=",
      as.character(id), ", session=", as.character(session),
      ". Missing session/run_number pairs: ", paste(missing_runs, collapse = ", "),
      call. = FALSE
    )
  }

  projection_df <- projection_df[row_index, , drop = FALSE]

  projected <- projection_df[, regressor_cols, drop = FALSE]
  projected[] <- lapply(projected, as.numeric)

  bad_cols <- regressor_cols[vapply(projected, function(x) any(!is.finite(x)), logical(1))]
  if (length(bad_cols) > 0L) {
    stop(
      "Projected L2 regressor(s) contain non-finite values for model '",
      l2_model_name, "': ", paste(bad_cols, collapse = ", "),
      call. = FALSE
    )
  }

  keep_cols <- regressor_cols[vapply(projected, function(x) length(unique(x)) > 1L, logical(1))]
  dropped_cols <- setdiff(regressor_cols, keep_cols)

  if (length(dropped_cols) > 0L) {
    lg$debug(
      "Dropping constant L2 projection regressors for id=%s session=%s model=%s: %s",
      as.character(id), as.character(session), l2_model_name,
      paste(dropped_cols, collapse = ", ")
    )
  }

  if (length(keep_cols) == 0L) {
    lg$info(
      "All projected L2 regressors were constant for id=%s session=%s model=%s. Skipping projection.",
      as.character(id), as.character(session), l2_model_name
    )
    return(NULL)
  }

  out <- data.frame(
    run_number = run_numbers,
    projected[, keep_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )

  out
}

#' Setup a first-level SPM model using a build_design_matrix object
#'
#' @param d_obj a single-subject design matrix object generated by build_design_matrix
#' @param gpa a gpa (glm_pipeline_arguments) object generated by setup_glm_pipeline
#' @param id subject id.
#' @param session session number.
#' @param model_name a string indicating the model name within \code{gpa} to setup. If
#'   you wish to setup multiple models, this is handled upstream in setup_l1_models.R
#' @param l1_confound_files optional character vector of nuisance regressor files,
#'   one per run, to include in SPM's multi_reg field
#' @param run_nifti optional character vector of NIfTI filenames used in l1 analysis.
#'   Prefer to pass these via the build_design_matrix object in $run_nifti
#' @param run_numbers optional integer vector of run numbers corresponding to \code{run_nifti}
#' @param run_sessions optional session vector corresponding to \code{run_nifti}.
#' @param source_run_numbers optional source run numbers for pooled/session models.
#' @param l1_session_mode optional session handling mode for SPM L1 setup.
#' @param execute_spm whether to execute SPM setup/estimation immediately. Default: FALSE
#'
#' @importFrom checkmate assert_class assert_string assert_character assert_file_exists
#' @importFrom lgr get_logger
#' @importFrom data.table fread
#' @author Michael Hallquist
#' @export
#'
spm_l1_model <- function(
  id = NULL, session = NULL, l1_confound_files = NULL, d_obj, gpa,
  model_name = NULL, run_nifti = NULL, run_numbers = NULL, run_sessions = NULL,
  source_run_numbers = NULL, l1_session_mode = NULL, execute_spm = FALSE
) {

  checkmate::assert_scalar(id, null.ok = FALSE)
  checkmate::assert_scalar(session, null.ok = TRUE)
  checkmate::assert_character(l1_confound_files, null.ok = TRUE)
  checkmate::assert_class(d_obj, "bdm")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_string(model_name) # single string
  checkmate::assert_integerish(run_numbers, null.ok = TRUE, lower = 1L)
  checkmate::assert_integerish(run_sessions, null.ok = TRUE, any.missing = FALSE)
  checkmate::assert_integerish(source_run_numbers, null.ok = TRUE, any.missing = FALSE, lower = 1L)
  checkmate::assert_string(l1_session_mode, null.ok = TRUE)

  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  lg$set_threshold(gpa$lgr_threshold)
  bind_l1_log_appenders(lg = lg, gpa = gpa)

  if (!is.null(d_obj$run_nifti)) {
    lg$debug("Using internal NIfTI files (run_nifti) within d_obj for SPM level 1 setup")
    run_nifti <- d_obj$run_nifti
  } else if (!is.null(d_obj$run_niftis)) {
    lg$debug("Using internal NIfTI files (run_niftis) within d_obj for SPM level 1 setup")
    run_nifti <- d_obj$run_niftis
  }

  if (is.null(run_nifti)) {
    msg <- "spm_l1_model requires run_nifti either in d_obj or passed explicitly."
    lg$error(msg)
    stop(msg)
  }

  if (is.null(run_numbers)) {
    if (!is.null(d_obj$runs_to_output)) {
      run_numbers <- as.integer(d_obj$runs_to_output)
    } else {
      run_numbers <- seq_along(run_nifti)
    }
  }
  run_numbers <- as.integer(run_numbers)
  if (is.null(source_run_numbers)) source_run_numbers <- run_numbers
  source_run_numbers <- as.integer(source_run_numbers)
  if (is.null(run_sessions)) {
    if (is.null(session)) {
      msg <- "session must be provided when run_sessions is NULL."
      lg$error(msg)
      stop(msg)
    }
    run_sessions <- rep.int(as.integer(session), length(run_numbers))
  }
  run_sessions <- as.integer(run_sessions)
  if (length(run_numbers) != length(run_nifti)) {
    msg <- "Length of run_numbers must match length of run_nifti for SPM L1 setup."
    lg$error(msg)
    stop(msg)
  }
  if (length(run_sessions) != length(run_numbers) || length(source_run_numbers) != length(run_numbers)) {
    msg <- "run_numbers, run_sessions, and source_run_numbers must have matching lengths."
    lg$error(msg)
    stop(msg)
  }

  stopifnot(length(run_nifti) == length(d_obj$run_volumes)) # need these to align
  checkmate::assert_character(run_nifti, null.ok = FALSE)
  checkmate::assert_file_exists(run_nifti) # all exist

  stopifnot(model_name %in% names(gpa$l1_models$models))

  spm_defaults <- list(
    matlab_cmd = "matlab",
    matlab_args = "-batch",
    hpf = 100,
    hrf_derivs = "none",
    cvi = "AR(1)",
    estimation_method = "Classical",
    write_residuals = FALSE,
    fmri_t = NULL,
    fmri_t0 = NULL,
    condition_contrasts = TRUE,
    unit_contrasts = TRUE,
    effects_of_interest_F = TRUE,
    l2_projection_model = NULL,
    l2_projection_interactions = NULL,
    l2_projection_interaction_unit_contrasts = TRUE,
    l2_projection_interaction_contrast_modes = NULL,
    l1_session_mode = "separate",
    spm_execute_setup = FALSE,
    spm_execute_glm = FALSE,
    spm_execute_contrasts = FALSE,
    concatenate_runs = FALSE,
    cleanup_tmp = TRUE,
    nifti_tmpdir = NULL,
    spm_path = NULL,
    force_l1_creation = FALSE
  )

  spm_settings <- populate_defaults(gpa$glm_settings$spm, spm_defaults)
  if (!is.null(l1_session_mode)) spm_settings$l1_session_mode <- l1_session_mode
  spm_settings$l1_session_mode <- resolve_spm_l1_session_mode(spm_settings$l1_session_mode, lg = lg)
  session_out <- if (identical(spm_settings$l1_session_mode, "pooled")) 0L else as.integer(session)
  spm_compute_env <- get_compute_environment(gpa, c("spm"))

  lg$debug(
    "SPM L1 setup start: id=%s session=%s model=%s n_runs=%d l1_session_mode=%s",
    id, ifelse(is.null(session), NA, as.character(session)), model_name, length(run_nifti),
    spm_settings$l1_session_mode
  )

  spm_l1_output_dir <- get_output_directory(
    id = id, session = session_out, l1_model = model_name,
    gpa = gpa, glm_software = "spm", what = "l1"
  )

  lg$info("Create l1 spm_l1_output_dir: %s", spm_l1_output_dir)
  dir.create(spm_l1_output_dir, showWarnings = FALSE, recursive = TRUE)
  spm_l1_output_dir <- normalizePath(spm_l1_output_dir)

  projection_model <- resolve_spm_l2_projection_model(
    gpa = gpa,
    l1_model = model_name,
    lg = lg
  )
  interaction_contrast_modes <- resolve_spm_l2_projection_interaction_contrast_modes(
    modes = spm_settings$l2_projection_interaction_contrast_modes,
    legacy_unit_contrasts = isTRUE(spm_settings$l2_projection_interaction_unit_contrasts),
    lg = lg
  )

  if (!is.null(projection_model)) {
    l2_projection_df <- extract_spm_l2_projection(
      gpa = gpa,
      l2_model_name = projection_model,
      id = id,
      session = session,
      run_numbers = run_numbers,
      run_sessions = run_sessions,
      source_run_numbers = source_run_numbers,
      l1_session_mode = spm_settings$l1_session_mode,
      lg = lg
    )
    if (!is.null(l2_projection_df) && ncol(l2_projection_df) > 1L) {
      d_obj$spm_l2_projection <- l2_projection_df
      d_obj$spm_l2_projection_model <- projection_model
      lg$info(
        "Applying %d projected L2 regressor(s) from model '%s' in SPM L1 setup.",
        ncol(l2_projection_df) - 1L,
        projection_model
      )

      interaction_terms <- resolve_spm_l2_projection_interactions(
        interactions = spm_settings$l2_projection_interactions,
        available_terms = setdiff(names(l2_projection_df), "run_number"),
        lg = lg
      )
      if (length(interaction_terms) > 0L) {
        d_obj$spm_l2_projection_interactions <- interaction_terms
        lg$info(
          "Will project %d selected L2-by-L1 interaction term(s): %s",
          length(interaction_terms),
          paste(interaction_terms, collapse = ", ")
        )
        if (length(interaction_contrast_modes) > 0L) {
          lg$info(
            "Projected interaction contrast modes: %s",
            paste(interaction_contrast_modes, collapse = ", ")
          )
        }
      }
    }
  }

  if (is.null(spm_settings$concatenate_runs)) {
    spm_settings$concatenate_runs <- FALSE
  }

  micro <- infer_spm_microtime(
    run_nifti = run_nifti,
    fmri_t = spm_settings$fmri_t,
    fmri_t0 = spm_settings$fmri_t0,
    lg = lg
  )
  spm_settings$fmri_t <- micro$fmri_t
  spm_settings$fmri_t0 <- micro$fmri_t0

  if (!length(l1_confound_files) || all(!is.na(l1_confound_files) & !nzchar(l1_confound_files))) {
    l1_confound_files <- NULL
  }
  if (!is.null(l1_confound_files)) {
    if (isTRUE(spm_settings$concatenate_runs)) {
      if (length(l1_confound_files) != 1L) {
        msg <- "For concatenated runs, l1_confound_files must have length 1 (a single concatenated file)."
        lg$error(msg)
        stop(msg)
      }
    } else {
      if (length(l1_confound_files) != length(run_nifti)) {
        msg <- "Length of l1_confound_files must match length of run_nifti."
        lg$error(msg)
        stop(msg)
      }
    }
  }

  if (is.null(spm_settings$nifti_tmpdir)) {
    spm_settings$nifti_tmpdir <- file.path(spm_l1_output_dir, "nifti_tmp")
  }
  if (isTRUE(execute_spm)) {
    spm_settings$spm_execute_setup <- TRUE
    spm_settings$spm_execute_glm <- TRUE
    spm_settings$spm_execute_contrasts <- TRUE
  }

  lg$debug(
    "SPM L1 settings: spm_path=%s concatenate_runs=%s execute_setup=%s execute_glm=%s execute_contrasts=%s nifti_tmpdir=%s cvi=%s estimation_method=%s write_residuals=%s fmri_t=%s fmri_t0=%s",
    spm_settings$spm_path,
    spm_settings$concatenate_runs,
    spm_settings$spm_execute_setup,
    spm_settings$spm_execute_glm,
    spm_settings$spm_execute_contrasts,
    spm_settings$nifti_tmpdir,
    spm_settings$cvi,
    spm_settings$estimation_method,
    spm_settings$write_residuals,
    spm_settings$fmri_t,
    spm_settings$fmri_t0
  )

  # Determine confounds handling for SPM
  ts_files <- NULL
  if (!is.null(l1_confound_files)) {
    valid_confound <- !is.na(l1_confound_files) & nzchar(l1_confound_files) & file.exists(l1_confound_files)

    if (all(valid_confound)) {
      ts_files <- l1_confound_files
    } else {
      stop("SPM confound files were supplied but are missing or invalid; refusing to omit nuisance adjustment.", call. = FALSE)
    }
  }

  # Ensure SPM uses the run_niftis field
  d_obj$run_niftis <- run_nifti

  spm_status <- get_spm_status(spm_l1_output_dir, lg = lg)

  signature_settings <- spm_settings[!grepl("^force_|^spm_execute_", names(spm_settings))]
  signature <- spm_setup_signature(list(
    bdm = d_obj[c("design", "design_concat", "run_volumes", "drop_volumes", "tr", "runs_to_output",
                  "spm_l2_projection", "spm_l2_projection_interactions")],
    contrasts = gpa$l1_models$models[[model_name]]$contrasts,
    settings = signature_settings, run_sessions = run_sessions, source_runs = source_run_numbers),
    c(run_nifti, ts_files))
  reuse <- !isTRUE(spm_settings$force_l1_creation) && !isTRUE(spm_status$spm_failed) &&
    spm_setup_current(spm_l1_output_dir, signature) &&
    all(file.exists(file.path(spm_l1_output_dir, c("setup_glm_design.m", "run_glm.m", "validate_spm_outputs.m"))))
  spm_syntax <- NULL
  if (reuse) {
    spm_syntax <- readRDS(file.path(spm_l1_output_dir, "spm_setup.rds"))$syntax
  } else {
    spm_invalidate_setup(spm_l1_output_dir)
    lg$info("Generating SPM batch scripts for model: %s", model_name)
    spm_syntax <- tryCatch(
      {
        generate_spm_mat(
          bdm = d_obj, ts_files = ts_files, output_dir = spm_l1_output_dir,
          hpf = spm_settings$hpf, hrf_derivs = spm_settings$hrf_derivs, cvi = spm_settings$cvi,
          estimation_method = spm_settings$estimation_method,
          write_residuals = spm_settings$write_residuals,
          fmri_t = spm_settings$fmri_t,
          fmri_t0 = spm_settings$fmri_t0,
          nifti_tmpdir = spm_settings$nifti_tmpdir, cleanup_tmp = spm_settings$cleanup_tmp,
          condition_contrasts = spm_settings$condition_contrasts,
          unit_contrasts = spm_settings$unit_contrasts,
          effects_of_interest_F = spm_settings$effects_of_interest_F,
          spm_execute_setup = FALSE,
          spm_execute_glm = FALSE,
          spm_execute_contrasts = FALSE,
          concatenate_runs = spm_settings$concatenate_runs,
          spm_path = spm_settings$spm_path,
          matlab_cmd = spm_settings$matlab_cmd,
          matlab_args = spm_settings$matlab_args,
          compute_env = spm_compute_env
        )
      },
      error = function(e) {
        lg$error("Problem running generate_spm_mat for model %s: %s", model_name, as.character(e))
        return(NULL)
      }
    )
  }

  if (is.null(spm_syntax)) {
    msg <- sprintf("SPM L1 setup failed to generate scripts for %s in %s.", model_name, spm_l1_output_dir)
    lg$error(msg)
    stop(msg)
  }

  # Always derive SPM contrasts from the L1 model object (mobj)
  mobj <- gpa$l1_models$models[[model_name]]
  if (is.null(mobj) || is.null(mobj$contrasts) || !inherits(mobj$contrasts, "matrix") || nrow(mobj$contrasts) == 0L) {
    msg <- sprintf("L1 model contrasts are missing or empty for model %s. Cannot generate SPM contrasts.", model_name)
    lg$error(msg)
    stop(msg)
  }

  if (!is.null(spm_syntax$projection_interaction_terms) &&
      length(spm_syntax$projection_interaction_terms) > 0L) {
    mobj <- augment_spm_l1_contrasts_for_projection(
      mobj = mobj,
      interaction_terms = spm_syntax$projection_interaction_terms,
      contrast_modes = interaction_contrast_modes,
      lg = lg
    )
  }

  projection_main_effect_terms <- character(0)
  projection_main_effect_weights <- NULL
  if (identical(spm_settings$l1_session_mode, "pooled") &&
      !isTRUE(spm_settings$concatenate_runs) &&
      !is.null(d_obj$spm_l2_projection) &&
      is.data.frame(d_obj$spm_l2_projection)) {
    projection_main_effect_terms <- setdiff(names(d_obj$spm_l2_projection), "run_number")
    if (length(projection_main_effect_terms) > 0L) {
      projection_main_effect_weights <- d_obj$spm_l2_projection[, projection_main_effect_terms, drop = FALSE]
    }
  }

  spm_contrast_cmds <- generate_spm_contrasts_from_model(
    output_dir = spm_l1_output_dir,
    mobj = mobj,
    spm_path = spm_settings$spm_path,
    execute = FALSE,
    matlab_cmd = spm_settings$matlab_cmd,
    matlab_args = spm_settings$matlab_args,
    projection_interaction_terms = if (!is.null(spm_syntax$projection_interaction_terms)) spm_syntax$projection_interaction_terms else character(0),
    projection_interaction_contrast_modes = interaction_contrast_modes,
    projection_interaction_run_labels = if (identical(spm_settings$l1_session_mode, "pooled")) {
      paste0("s", run_sessions, "r", source_run_numbers)
    } else {
      run_numbers
    },
    projection_main_effect_terms = projection_main_effect_terms,
    projection_main_effect_weights = projection_main_effect_weights
  )

  if (!is.null(spm_contrast_cmds)) {
    cat(
      c("#!/bin/bash", "set -euo pipefail", spm_contrast_cmds$extract_cmd, spm_contrast_cmds$setup_cmd),
      file = file.path(spm_l1_output_dir, "setup_spm_contrasts.sh"),
      sep = "\n"
    )
  }
  spm_syntax$contrast_cmds <- spm_contrast_cmds

  if (!reuse) spm_save_setup(spm_l1_output_dir, signature, spm_syntax,
    isTRUE(spm_settings$write_residuals) && normalize_spm_estimation_method(spm_settings$estimation_method) == "Classical")

  spm_execute_model(spm_l1_output_dir, 1L, spm_settings, spm_compute_env)

  spm_script_files <- list.files(spm_l1_output_dir, pattern = "\\.m$", full.names = FALSE)
  lg$debug(
    "SPM L1 scripts written in %s: %d files (%s)",
    spm_l1_output_dir,
    length(spm_script_files),
    if (length(spm_script_files) > 0L) paste(head(spm_script_files, 5L), collapse = ", ") else "none"
  )

  spm_status <- get_spm_status(spm_l1_output_dir, lg = lg)

  spm_l1_df <- data.frame(
    id = id, session = session_out,
    l1_model = model_name,
    l1_session_mode = spm_settings$l1_session_mode,
    l2_projection_model = ifelse(!is.null(d_obj$spm_l2_projection_model), d_obj$spm_l2_projection_model, NA_character_),
    l2_projection_interaction_modes = if (length(interaction_contrast_modes) > 0L) paste(interaction_contrast_modes, collapse = "|") else NA_character_,
    spm_dir = spm_l1_output_dir, # explicitly named for clarity in run_spm_sepjobs
    n_runs = length(run_nifti),
    run_nifti = I(list(run_nifti)),
    run_sessions = I(list(run_sessions)),
    source_run_numbers = I(list(source_run_numbers)),
    run_volumes = I(list(d_obj$run_volumes)),
    l1_confound_files = I(list(l1_confound_files)),
    concatenate_runs = spm_settings$concatenate_runs,
    gunzip_cmds = I(list(spm_syntax$gunzip_cmds)),
    contrast_cmds = I(list(spm_syntax$contrast_cmds))
  )

  spm_l1_df <- cbind(spm_l1_df, spm_status)
  spm_l1_df$to_run <- !spm_l1_df$spm_complete

  return(spm_l1_df)
}
