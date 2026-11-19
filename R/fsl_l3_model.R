#' Prepare a level 3 FSL FEAT group model from lower-level COPEs
#'
#' @param l3_df a data.frame containing cope inputs for a given l3 model, as well as metadata
#'   that identify the model.
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#'
#' @importFrom dplyr mutate filter select left_join inner_join pull
#' @author Michael Hallquist
#' @keywords internal
fsl_l3_model <- function(l3_df=NULL, gpa) {
  checkmate::assert_data_frame(l3_df)
  checkmate::assert_subset(c("id", "session", "l3_model"), names(l3_df))
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  lg <- lgr::get_logger("glm_pipeline/l3_setup")
  lg$set_threshold(gpa$lgr_threshold)

  if (length(unique(l3_df$l3_model)) > 1L) {
    msg <- "fsl_l3_model is designed for execution on a single l3 model"
    lg$error(msg)
    stop(msg)
  }

  l3_model <- l3_df$l3_model[1L]

  l3_input_mode <- NULL
  if ("l3_input_mode" %in% names(l3_df)) {
    modes <- unique(stats::na.omit(l3_df$l3_input_mode))
    if (length(modes) > 1L) {
      msg <- "fsl_l3_model requires a single l3_input_mode within each input data.frame"
      lg$error(msg)
      stop(msg)
    }
    if (length(modes) == 1L) l3_input_mode <- modes[1L]
  }
  if (is.null(l3_input_mode)) {
    l3_input_mode <- gpa$l3_models$models[[l3_model]]$l3_input_mode
  }
  l3_input_mode <- normalize_l3_input_mode(l3_input_mode)
  checkmate::assert_subset(l3_input_mode, longitudinal_l3_input_modes())

  n_sessions_input <- length(unique(l3_df$session))
  if (identical(l3_input_mode, "per_session") && n_sessions_input > 1L) {
    msg <- "fsl_l3_model with l3_input_mode='per_session' requires a single-session input data.frame"
    lg$error(msg)
    stop(msg)
  }

  if (length(unique(l3_df$l1_model)) > 1L) {
    msg <- "fsl_l3_model is designed for execution on a single l1 model"
    lg$error(msg)
    stop(msg)
  }

  if (isTRUE(gpa$multi_run)) {
    if (length(unique(l3_df$l2_model)) > 1L) {
      msg <- "fsl_l3_model is designed for execution on a single l2 model"
      lg$error(msg)
      stop(msg)
    }
  }

  # elements of metadata for l3
  session <- if (n_sessions_input == 1L) l3_df$session[1L] else 0L
  l1_model <- l3_df$l1_model[1L]
  l2_model <- l3_df$l2_model[1L]
  l1_cope_name <- l1_contrast <- l3_df$l1_cope_name[1L] #use the double assign for synonyms used in glue() output expression
  l2_cope_name <- l2_contrast <- l3_df$l2_cope_name[1L]

  # tracking data frame for this model (column names should follow variable names)
  if (isTRUE(gpa$multi_run)) {
    feat_l3_df <- data.frame(l1_model, l1_cope_name, l2_model, l2_cope_name, l3_model, session, l3_input_mode)
  } else {
    feat_l3_df <- data.frame(l1_model, l1_cope_name, l3_model, session, l3_input_mode)
  }

  # handle model respecification based on available data (e.g., if some subjects failed to run)
  mobj <- respecify_l3_model(gpa$l3_models$models[[l3_model]], new_data = l3_df)

  # align id/session types before inner_join so num2fac = "session" never causes factor/integer join errors
  l3_df$id <- as.character(l3_df$id)
  meta_for_join <- mobj$metadata
  meta_for_join$id <- as.character(meta_for_join$id)
  if (!identical(class(l3_df$session), class(meta_for_join$session))) {
    l3_df$session <- as.character(l3_df$session)
    meta_for_join$session <- as.character(meta_for_join$session)
  }

  # now make sure that the l3_df aligns perfectly with the data being modeled (mostly dropping copes in l3_df that aren't in model)
  l3_df <- l3_df %>% dplyr::inner_join(meta_for_join, by = c("id", "session"))

  dmat_to_write <- if (!is.null(mobj$model_matrix_noalias)) mobj$model_matrix_noalias else mobj$model_matrix
  cmat_to_write <- if (!is.null(mobj$model_matrix_noalias)) mobj$contrasts_noalias else mobj$contrasts
  if (is.null(cmat_to_write) || nrow(cmat_to_write) == 0L) {
    lg$warn("No estimable L3 contrasts remain for L3 model '%s' after respecification. Skipping.", l3_model)
    return(NULL)
  }
  stopifnot(ncol(dmat_to_write) == ncol(cmat_to_write))

  feat_l3_df$cope_list <- list(data.frame(
    l3_cope_number = seq_len(nrow(cmat_to_write)),
    l3_cope_name = rownames(cmat_to_write) %||% paste0("cope", seq_len(nrow(cmat_to_write))),
    stringsAsFactors = FALSE
  ))

  # generate FSL EV syntax for these regressors
  ev_syntax <- fsl_generate_fsf_ev_syntax(inputs = l3_df$cope_file, dmat = dmat_to_write)

  # generate FSF contrast syntax for this setup
  contrast_syntax <- fsl_generate_fsf_contrast_syntax(cmat_to_write)

  l3_fsf_syntax <- readLines(system.file("feat_lvl3_copefiles_template.fsf", package = "fmri.pipeline"))

  # Add EVs and contrasts into FSF
  l3_fsf_syntax <- c(l3_fsf_syntax, ev_syntax, contrast_syntax)

  # need to determine number of copes (contrasts) at level 1, which depends on the model being fit
  # FSL usually reads this from the .feat directories itself, but for batch processing, better to insert into the FSF ourselves
  # Need to put this just after the high pass filter cutoff line for Feat to digest it happily

  #n_l3_models <- length(gpa$l3_models$models)
  #n_l1_models <- length(gpa$l1_models$models)

  l3_outdir <- get_output_directory(
    l1_contrast=l1_cope_name, l1_model=l1_model,
    l2_contrast=l2_cope_name, l2_model=l2_model,
    l3_model=l3_model, gpa = gpa, what = "l3"
  )

  l3_feat_fsf <- file.path(
    l3_outdir,
    glue::glue(gpa$output_locations$feat_l3_fsf) # evaluate glue expression
  )

  n_global_sessions <- n_sessions_input
  if (!is.null(gpa$run_data) && "session" %in% names(gpa$run_data)) {
    n_global_sessions <- length(unique(gpa$run_data$session))
  }
  if (identical(l3_input_mode, "per_session") && n_global_sessions > 1L) {
    if (grepl("\\.fsf$", l3_feat_fsf)) {
      l3_feat_fsf <- sub("\\.fsf$", paste0("_ses-", session, ".fsf"), l3_feat_fsf)
    } else {
      l3_feat_fsf <- paste0(l3_feat_fsf, "_ses-", session, ".fsf")
    }
  }

  l3_feat_dir <- sub(".fsf$", ".gfeat", l3_feat_fsf)

  if (!dir.exists(l3_outdir)) {
    lg$debug("Creating L3 output directory: %s", l3_outdir)
    dir.create(l3_outdir, recursive = TRUE)
  }

  # specify output directory (removing .gfeat suffix)
  # .OUTPUTDIR. : the feat output location
  l3_fsf_syntax <- gsub(".OUTPUTDIR.", sub("\\.gfeat$", "", l3_feat_dir), l3_fsf_syntax, fixed = TRUE)

  # handle custom L3 FSF syntax
  l3_fsf_syntax <- add_custom_feat_syntax(l3_fsf_syntax, gpa$additional$feat_l3_args, lg)

  # handle outlier deweighting flag (ifelse won't keep list names...)
  deweight <- if(isTRUE(mobj$fsl_outlier_deweighting)) list(robust_yn = 1) else list(robust_yn = 0)
  l3_fsf_syntax <- add_custom_feat_syntax(l3_fsf_syntax, deweight, lg)

  feat_save_setup(l3_feat_fsf, l3_fsf_syntax, inputs = feat_higher_inputs(l3_df$cope_file),
    level = 3L, spec = list(design = dmat_to_write, contrasts = cmat_to_write),
    force = isTRUE(gpa$glm_settings$fsl$force_l3_creation))
  feat_l3_df <- dplyr::bind_cols(feat_l3_df,
    get_feat_status(feat_dir = feat_output_directory(l3_fsf_syntax), feat_fsf = l3_feat_fsf, lg = lg))
  feat_l3_df$to_run <- !feat_l3_df$feat_complete

  return(feat_l3_df)
}
