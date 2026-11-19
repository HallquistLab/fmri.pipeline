# TODO: clean up the mess of the 'local' approach, but for now, prefer the glue expression
# consider whether we want a 'what' type argument, like subject directory, analysis directory,
# run directory, model-specific analysis directory, etc.

#' small helper function to return the location of an l1 directory based on
#'   id, session, and run number
#'
#' @param id The id of a participant
#' @param session The session number to lookup
#' @param run_number The run number to lookup
#' @param l1_model,l2_model,l3_model Optional model names used in output templates.
#' @param l1_cope_number,l2_cope_number,l3_cope_number Optional cope numbers used in output templates.
#'   These are combined with the corresponding contrast names into labels such as \code{02_EV_house}.
#' @param l1_contrast,l2_contrast,l3_contrast Optional contrast names used in output templates.
#' @param gpa A \code{glm_pipeline_arguments} object
#' @param glm_software which software is being used for the analysis (since directories may vary)
#' @param what output directory type to return.
#' @param create_if_missing whether to create the directory if it does not exist
get_output_directory <- function(id = NULL, session = NULL, run_number = NULL,
  l1_model = NULL, l2_model = NULL, l3_model = NULL,
  l1_cope_number = NULL, l2_cope_number = NULL, l3_cope_number = NULL,
  l1_contrast = NULL, l2_contrast = NULL, l3_contrast = NULL,
  what="l1", gpa, glm_software = "fsl", create_if_missing = FALSE) {

  checkmate::assert_scalar(id, null.ok=TRUE)
  checkmate::assert_integerish(session, null.ok = TRUE)
  if (is.null(session)) session <- 1
  checkmate::assert_integerish(run_number, lower = 1, null.ok = TRUE)
  checkmate::assert_string(l1_model, null.ok=TRUE)
  checkmate::assert_string(l2_model, null.ok = TRUE)
  checkmate::assert_string(l3_model, null.ok = TRUE)
  checkmate::assert_integerish(l1_cope_number, len = 1L, lower = 1L, null.ok = TRUE)
  checkmate::assert_integerish(l2_cope_number, len = 1L, lower = 1L, null.ok = TRUE)
  checkmate::assert_integerish(l3_cope_number, len = 1L, lower = 1L, null.ok = TRUE)
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  stopifnot("run_data" %in% names(gpa))
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_logical(create_if_missing, len=1L)

  if (!is.null(l1_model)) checkmate::assert_subset(l1_model, names(gpa$l1_models$models))
  if (!is.null(l2_model)) checkmate::assert_subset(l2_model, names(gpa$l2_models$models))
  if (!is.null(l3_model)) checkmate::assert_subset(l3_model, names(gpa$l3_models$models))

  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  lg$set_threshold(gpa$lgr_threshold)

  # create local synonyms for these in case they are used the glue() expression.
  l1_cope_name <- l1_contrast
  l2_cope_name <- l2_contrast
  l3_cope_name <- l3_contrast
  l1_cope_label <- if (!is.null(l1_cope_number) && !is.null(l1_cope_name)) {
    sprintf("%02d_%s", as.integer(l1_cope_number), l1_cope_name)
  } else {
    l1_cope_name
  }
  l2_cope_label <- if (!is.null(l2_cope_number) && !is.null(l2_cope_name)) {
    sprintf("%02d_%s", as.integer(l2_cope_number), l2_cope_name)
  } else {
    l2_cope_name
  }
  l3_cope_label <- if (!is.null(l3_cope_number) && !is.null(l3_cope_name)) {
    sprintf("%02d_%s", as.integer(l3_cope_number), l3_cope_name)
  } else {
    l3_cope_name
  }

  # helper subfunction to select rows where all conditions in filter_list match
  subset_run_data <- function(run_data, filter_list) {
    res <- sapply(seq_along(filter_list), function(xx) {
      run_data[[ names(filter_list)[xx] ]] == filter_list[[xx]]
    })
    all_match <- which(apply(res, 1, all))
    rinfo <- run_data %>% dplyr::slice(all_match)
    if (nrow(rinfo) == 0L) {
      lg$error("Unable to locate a record in gpa$run_data for id %s, session %s, run_number %s.", id, session, run_number)
      return(NULL)
    } else {
      return(rinfo)
    }
  }

  if (what == "l1") {
    # choose software-specific output template if available
    l1_directory <- gpa$output_locations$feat_l1_directory
    if (glm_software == "spm" && "spm_l1_directory" %in% names(gpa$output_locations)) {
      l1_directory <- gpa$output_locations$spm_l1_directory
      spm_mode <- "separate"
      if (!is.null(gpa$glm_settings) &&
          !is.null(gpa$glm_settings$spm) &&
          "l1_session_mode" %in% names(gpa$glm_settings$spm) &&
          !is.null(gpa$glm_settings$spm$l1_session_mode)) {
        spm_mode <- tolower(as.character(gpa$glm_settings$spm$l1_session_mode)[1L])
      }
      if (spm_mode %in% c("pooled", "id", "across_sessions", "subject")) {
        if ("spm_l1_id_scope_directory" %in% names(gpa$output_locations) &&
            checkmate::test_string(gpa$output_locations$spm_l1_id_scope_directory)) {
          l1_directory <- gpa$output_locations$spm_l1_id_scope_directory
        } else if (checkmate::test_string(l1_directory)) {
          # back-compat fallback for older objects lacking spm_l1_id_scope_directory
          l1_directory <- gsub("([/\\\\])ses-\\{session\\}([/\\\\]|$)", "\\1", l1_directory, perl = TRUE)
        }
      }
    } else if (glm_software == "afni") {
      l1_directory <- gpa$output_locations$afni_l1_directory %||% file.path(
        "{gpa$output_directory}", "afni_l1", "sub-{id}", "ses-{session}", "{l1_model}"
      )
    }

    if (l1_directory == "local") {
      checkmate::assert_integerish(run_number, lower = 1, null.ok = FALSE) # enforce run_number for l1 lookup
      rinfo <- subset_run_data(gpa$run_data, list(id=id, session=session, run_number=run_number))
      if ("mr_dir" %in% names(rinfo)) {
        lg$debug("Using mr_dir l1 directory lookup: %s", rinfo$mr_dir[1L])
        if (is.null(l1_model)) {
          lg$debug("Lookup from analysis_name: %s", gpa$analysis_name)
          out_dir <- file.path(normalizePath(file.path(rinfo$mr_dir[1L], "..")), gpa$analysis_name)
        } else {
          lg$debug("Lookup from model outdir: %s", gpa$l1_models$models[[l1_model]]$outdir)
          out_dir <- file.path(normalizePath(file.path(rinfo$mr_dir[1L], "..")), gpa$l1_models$models[[l1_model]]$outdir)
        }
      } else {
        # look in parent folder of relevant run nifti and place l1 model there
        rn <- get_mr_abspath(rinfo[1, , drop = F], "run_nifti")
        lg$debug("Using run_nifti l1 directory lookup: %s", rn)

        if (is.null(l1_model)) {
          lg$debug("Lookup from analysis_name: %s", gpa$analysis_name)
          out_dir <- file.path(
            normalizePath(file.path(dirname(rn), "..")),
            gpa$analysis_name
          )
        } else {
          lg$debug("Lookup from model outdir: %s", gpa$l1_models$models[[l1_model]]$outdir)
          out_dir <- file.path(
            normalizePath(file.path(dirname(rn), "..")),
            gpa$l1_models$models[[l1_model]]$outdir
          )
        }
      }
    } else {
      out_dir <- glue::glue(l1_directory) # evaluate glue expression
    }
  } else if (what == "l2") {
    l2_directory <- gpa$output_locations$feat_l2_directory
    if (!is.null(l2_model)) {
      l2_scope <- gpa$l2_models$models[[l2_model]]$l2_scope
      if (is.null(l2_scope) || !is.character(l2_scope) || length(l2_scope) != 1L || !nzchar(l2_scope)) {
        l2_scope <- "id_session"
      }

      # pooled L2 models are subject-level and should not write to synthetic ses-0 folders
      if (identical(l2_scope, "id")) {
        if ("feat_l2_id_scope_directory" %in% names(gpa$output_locations) &&
            checkmate::test_string(gpa$output_locations$feat_l2_id_scope_directory)) {
          l2_directory <- gpa$output_locations$feat_l2_id_scope_directory
        } else if (checkmate::test_string(l2_directory)) {
          # back-compat fallback for older objects lacking feat_l2_id_scope_directory
          l2_directory <- gsub("([/\\\\])ses-\\{session\\}([/\\\\]|$)", "\\1", l2_directory, perl = TRUE)
        }
      }
    }

    out_dir <- glue::glue(l2_directory) # evaluate glue expression
  } else if (what == "l3") {
    l3_directory <- gpa$output_locations$feat_l3_directory
    if (glm_software == "spm" && "spm_l3_directory" %in% names(gpa$output_locations)) {
      l3_directory <- gpa$output_locations$spm_l3_directory
    } else if (glm_software == "afni") {
      l3_directory <- gpa$output_locations$afni_l3_directory %||% file.path(
        "{gpa$output_directory}", "afni_l3", "L3m-{l3_model}", "L1m-{l1_model}", "l1c-{l1_cope_name}"
      )
    }
    out_dir <- glue::glue(l3_directory) # evaluate glue expression
  } else if (what == "sub") {
    sub_directory <- gpa$output_locations$feat_sub_directory
    if (glm_software == "spm" && "spm_sub_directory" %in% names(gpa$output_locations)) {
      sub_directory <- gpa$output_locations$spm_sub_directory
    } else if (glm_software == "afni" && "afni_sub_directory" %in% names(gpa$output_locations)) {
      sub_directory <- gpa$output_locations$afni_sub_directory
    }

    if (sub_directory == "local") {
      rinfo <- subset_run_data(gpa$run_data, list(id = id, session = session, run_number = run_number))
    } else {
      out_dir <- glue::glue(sub_directory) # evaluate glue expression
    }
  } else if (what == "ses") {
    ses_directory <- gpa$output_locations$feat_ses_directory
    if (glm_software == "spm" && "spm_sub_directory" %in% names(gpa$output_locations)) {
      ses_directory <- gpa$output_locations$spm_sub_directory
    } else if (glm_software == "afni" && "afni_sub_directory" %in% names(gpa$output_locations)) {
      ses_directory <- gpa$output_locations$afni_sub_directory
    }

    if (ses_directory == "local") {
      rinfo <- subset_run_data(gpa$run_data, list(id = id, session = session))
    } else {
      out_dir <- glue::glue(ses_directory) # evaluate glue expression
    }
  }

  #glue returns empty strings when some variables are NULL
  if (identical(glue("{}"), out_dir) || nchar(out_dir) == 0L) {
    lg$error("Unable to sort out %s output directory", what)
    return(NULL)
  } else {
    #out_dir <- make.names(out_dir) #remove weird characters
  }

  if (isTRUE(create_if_missing) && !dir.exists(out_dir)) {
    lg$debug("Create output directory: %s", out_dir)
    dir.create(out_dir, recursive = TRUE)
  }

  return(out_dir)
}
