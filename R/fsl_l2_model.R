is_l2_passthrough_estimable <- function(dmat, cmat, tol = sqrt(.Machine$double.eps)) {
  if (checkmate::test_data_frame(dmat)) dmat <- as.matrix(dmat)
  if (checkmate::test_data_frame(cmat)) cmat <- as.matrix(cmat)
  if (!is.null(cmat) && !is.matrix(cmat)) {
    cmat <- matrix(as.numeric(cmat), nrow = 1L)
  }

  is.matrix(dmat) &&
    nrow(dmat) == 1L &&
    ncol(dmat) == 1L &&
    isTRUE(abs(as.numeric(dmat[1L, 1L]) - 1) <= tol) &&
    is.matrix(cmat) &&
    nrow(cmat) == 1L &&
    ncol(cmat) == 1L &&
    isTRUE(abs(as.numeric(cmat[1L, 1L]) - 1) <= tol)
}

#' Estimate a level 2 (subject) model using FSL FEAT with fixed effects integration of runs
#'
#' @param l1_df a data.frame containing all runs for a single subject and a single l1 model. This
#'   data.frame defines the inputs for the L2 analysis (i.e., which runs to combine).
#' @param l2_model a model string in gpa$l2_models containing the L2 model to setup
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#'
#' @importFrom dplyr filter
#' @author Michael Hallquist
#' @keywords internal
fsl_l2_model <- function(l1_df=NULL, l2_model, gpa) {
  checkmate::assert_data_frame(l1_df)
  checkmate::assert_subset(c("id", "session", "l1_model", "l1_cope_name", "l1_cope_number", "cope_file"), names(l1_df))
  checkmate::assert_string(l2_model) # single l2 model
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  lg <- lgr::get_logger("glm_pipeline/l2_setup")
  lg$set_threshold(gpa$lgr_threshold)

  if (length(unique(l1_df$id)) > 1L) {
    msg <- "fsl_l2_model is designed for execution on a single id data.frame"
    lg$error(msg)
    stop(msg)
  }

  l2_scope <- gpa$l2_models$models[[l2_model]]$l2_scope
  if (is.null(l2_scope) || !is.character(l2_scope) || length(l2_scope) != 1L || !nzchar(l2_scope)) {
    l2_scope <- "id_session"
  }
  checkmate::assert_subset(l2_scope, longitudinal_l2_scopes())

  if (identical(l2_scope, "id_session") && length(unique(l1_df$session)) > 1L) {
    msg <- "fsl_l2_model with l2_scope='id_session' requires a single-session input data.frame"
    lg$error(msg)
    stop(msg)
  }

  if (length(unique(l1_df$l1_model)) > 1L) {
    msg <- "fsl_l2_model is designed for execution on a single l1 model"
    lg$error(msg)
    stop(msg)
  }

  # elements of metadata for l2
  id <- l1_df$id[1L]
  session <- if (identical(l2_scope, "id")) 0L else l1_df$session[1L]
  l1_model <- l1_df$l1_model[1L]
  l1_cope_name <- l1_df$l1_cope_name[1L]
  l1_cope_number <- l1_df$l1_cope_number[1L]
  l1_cope_files <- l1_df$cope_file

  if (length(unique(l1_df$l1_cope_name)) > 1L || length(unique(l1_df$l1_cope_number)) > 1L) {
    msg <- "fsl_l2_model is designed for execution on a single L1 cope at a time"
    lg$error(msg)
    stop(msg)
  }

  l2_data <- compose_l2_model_data(gpa, lg = lg)
  l2_keys <- unique(l1_df[, c("id", "session", "run_number")])
  l2_data <- dplyr::semi_join(l2_data, l2_keys, by = c("id", "session", "run_number"))

  split_on <- if (identical(l2_scope, "id")) "id" else c("id", "session")
  refit_mobj <- respecify_l2_models_by_subject(
    gpa$l2_models$models[[l2_model]],
    l2_data,
    split_on = split_on,
    aggregated_session = 0L
  )

  if (identical(l2_scope, "id")) {
    ss_df <- refit_mobj$by_subject %>% dplyr::filter(id == !!id)
  } else {
    ss_df <- refit_mobj$by_subject %>% dplyr::filter(id == !!id & session == !!session)
  }

  if (nrow(ss_df) == 0L) {
    lg$warn(
      "No valid L2 design rows remain for id %s, session %s, l1 model %s, l1 cope %s, l2 model %s",
      id, session, l1_model, l1_cope_name, l2_model
    )
    return(NULL)
  }
  if (nrow(ss_df) > 1L) {
    lg$error("More than one subject-specific entry for id %s, session %s", id, session)
    return(NULL)
  }

  dmat <- ss_df$model_matrix[[1L]]
  cmat <- ss_df$contrasts[[1L]]
  if (checkmate::test_data_frame(cmat)) cmat <- as.matrix(cmat)
  if (!is.null(cmat) && !is.matrix(cmat)) {
    cmat <- matrix(as.numeric(cmat), nrow = 1L)
  }

  # Align l1_df to the exact non-missing covariate rows retained in dmat
  model_runs <- attr(dmat, "model_runs")
  if (!is.null(model_runs) && is.data.frame(model_runs)) {
    join_keys <- intersect(c("id", "session", "run_number"), intersect(names(model_runs), names(l1_df)))
    if (length(join_keys) > 0L) {
      l1_df <- dplyr::inner_join(model_runs, as.data.frame(l1_df), by = join_keys)
    }
  }
  l1_cope_files <- l1_df$cope_file

  cope_list <- ss_df$cope_list[[1L]]
  n_l2_copes <- if (is.null(cope_list)) 0L else nrow(cope_list)
  if (nrow(l1_df) == 0L || is.null(cmat) || nrow(cmat) == 0L) {
    lg$warn(
      "No estimable L2 contrasts remain after refitting design for id %s, session %s, l1 model %s, l1 cope %s, l2 model %s",
      id, session, l1_model, l1_cope_name, l2_model
    )
    return(NULL)
  }

  if (nrow(l1_df) == 1L) {
    pt_contrasts <- get_l2_passthrough_contrasts(dmat, cmat, gpa$l2_models$models[[l2_model]])
    if (length(pt_contrasts) > 0L) {
      return(make_l2_passthrough_row(
        l1_df = l1_df,
        l2_model = l2_model,
        gpa = gpa,
        passthrough_contrasts = pt_contrasts,
        lg = lg
      ))
    }

    lg$warn(
      paste(
        "Only one valid run remains for subject %s session %s, L1 model %s, L1 cope %s, L2 model %s,",
        "but the refit L2 design has no intercept/pooling contrast eligible for pass-through. Skipping this L2 analysis."
      ),
      id, session, l1_model, l1_cope_name, l2_model
    )
    return(NULL)
  }

  # tracking data frame for this model
  feat_l2_df <- data.frame(
    id = id, session = session,
    l1_model = l1_model,
    l1_cope_number = l1_cope_number,
    l1_cope_name = l1_cope_name,
    l2_model = l2_model,
    l2_scope = l2_scope,
    l2_input_mode = "cope_files",
    l2_passthrough = FALSE,
    n_l2_copes = n_l2_copes,
    n_input_files = nrow(l1_df),
    passthrough_cope_file = NA_character_,
    stringsAsFactors = FALSE
  )
  feat_l2_df$cope_list <- list(cope_list)

  # generate FSL EV syntax for these regressors
  ev_syntax <- fsl_generate_fsf_ev_syntax(inputs = l1_cope_files, dmat = dmat)

  # generate FSF contrast syntax for this setup
  contrast_syntax <- fsl_generate_fsf_contrast_syntax(cmat)

  l2_fsf_syntax <- readLines(system.file("feat_lvl2_nparam_template.fsf", package = "fmri.pipeline"))
  l2_fsf_syntax <- gsub("set fmri\\(inputtype\\) 1", "set fmri(inputtype) 2", l2_fsf_syntax)

  # Add EVs and contrasts into FSF
  l2_fsf_syntax <- c(l2_fsf_syntax, ev_syntax, contrast_syntax)

  # Get L2 output basename using the configured feat_l2_directory location
  l1_cope_path_name <- fs::path_sanitize(l1_cope_name, replacement = "_")
  fsl_l2_output_base <- get_output_directory(
    id = id, session = session,
    l1_model = l1_model, l2_model = l2_model,
    l1_cope_number = l1_cope_number, l1_contrast = l1_cope_path_name,
    gpa = gpa, glm_software = "fsl", what = "l2"
  )

  feat_l2_fsf <- paste0(sub("\\.fsf$", "", fsl_l2_output_base), ".fsf")
  l2_feat_fsf <- feat_l2_fsf
  l2_feat_dir <- sub("\\.fsf$", ".gfeat", l2_feat_fsf)

  l2_feat_parent_dir <- dirname(l2_feat_fsf)
  if (!dir.exists(l2_feat_parent_dir)) {
    lg$debug("Creating L2 output directory: %s", l2_feat_parent_dir)
    dir.create(l2_feat_parent_dir, recursive = TRUE)
  }

  # add columns regarding whether inputs already exist and FEAT is already complete
  feat_l2_df <- feat_l2_df %>%
    dplyr::bind_cols(get_feat_status(feat_dir = l2_feat_dir, feat_fsf = l2_feat_fsf, lg = lg))

  feat_l2_df$to_run <- !feat_l2_df$feat_complete

  lg$debug("Expected L2 feat directory is: %s", l2_feat_dir)
  lg$debug("Expected L2 feat fsf is: %s", l2_feat_fsf)

  # specify output directory (removing .gfeat suffix)
  # .OUTPUTDIR. : the feat output location
  l2_fsf_syntax <- gsub(".OUTPUTDIR.", sub("\\.gfeat$", "", l2_feat_dir), l2_fsf_syntax, fixed = TRUE)

  # handle custom L1 FSF syntax
  l2_fsf_syntax <- add_custom_feat_syntax(l2_fsf_syntax, gpa$additional$feat_l2_args, lg)

  # Rewrite FSF whenever this model is queued to run or force_l2_creation is TRUE so updated contrasts are never masked by a stale file
  write_l2_fsf <- !file.exists(l2_feat_fsf) ||
    isTRUE(gpa$glm_settings$fsl$force_l2_creation) ||
    isTRUE(feat_l2_df$to_run[1L])

  if (isTRUE(write_l2_fsf)) {
    lg$info("Writing L2 FSF syntax to: %s", l2_feat_fsf)
    cat(l2_fsf_syntax, file = l2_feat_fsf, sep = "\n")
  } else {
    lg$info("Skipping existing L2 FSF syntax: %s", l2_feat_fsf)
  }

  return(feat_l2_df)

}
