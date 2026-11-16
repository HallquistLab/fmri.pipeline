#' Output GLM configuration to a '.yaml' or 'json' file
#' @param gpa the `glm_pipeline_arguments` object with relevant fields already populated, including
#'   `$l1_models`, `$l2_models`, and/or `$l3_models`.
#' @param file the filename for the exported configuration file
#' @return invisibly return the YAML/JSON syntax written to `file`
#' @details 
#'   This function exports the configuration of level 1, 2, and 3 GLM models in a `glm_pipeline_arguments`
#'   (aka 'gpa') object to a YAML or JSON file. This allows the user to setup their models from this configuration
#'   file, rather than having to go through the building process manually (with prompts).
#' 
#'   Typically, the user will first run `build_l1_models`, `build_l2_models` (for multi-run data), and `build_l3_models`,
#'   then run `export_glm_config` so that the decisions made during the build process are codified in the resulting file. 
#' 
#' @examples 
#' \dontrun{
#'   gpa <- setup_glm_pipeline(analysis_name="flanker_test",
#'     output_directory="/proj/mnhallqlab/no_backup/flanker_gnomes",
#'     n_expected_runs = 2, tr = 2.0,
#'     run_data = run_df, subject_data = subject_df, trial_data = trial_df,
#'     l1_models=NULL, l2_models=NULL, l3_models=NULL
#'   )
#' 
#'  # setup level 1 models through menu system
#'  gpa <- build_l1_models(gpa)
#' 
#'  # for multi-run data, setup level 2 models through menu system
#'  gpa <- build_l2_models(gpa)
#' 
#'  # setup level 3 models through menu system
#'  gpa <- build_l3_models(gpa)
#' 
#'  # export model configuration to YAML file
#'  export_glm_config(gpa, file="my_glm_config.yaml")
#' 
#'  # in future, we can build the models from the file, rather than menus
#'  gpa <- build_l1_models(gpa, from_spec_file="my_glm_config.yaml")
#'  gpa <- build_l2_models(gpa, from_spec_file="my_glm_config.yaml")
#'  gpa <- build_l3_models(gpa, from_spec_file="my_glm_config.yaml")
#' }
#' @importFrom yaml as.yaml
#' @importFrom jsonlite toJSON
#' @export
export_glm_config <- function(gpa, file="glm_config.yaml") {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_string(file, null.ok = FALSE)

  have_l1 <- !is.null(gpa$l1_models)
  have_l2 <- isTRUE(gpa$multi_run) && !is.null(gpa$l2_models)
  have_l3 <- !is.null(gpa$l3_models)

  if (!any(c(have_l1, have_l2, have_l3))) {
    warning("In export_glm_config, nothing to export. Have you run build_l1_models, build_l2_models, and build_l3_models?")
    return(character(0))
  }
  
  l1_str <- get_l1_config(gpa)
  l2_str <- get_l2_config(gpa)
  l3_str <- get_l3_config(gpa)

  str_out <- c(l1_str, l2_str, l3_str)

  ext <- file_ext(file)
  if (is.na(ext)) { # missing file extension -- default to YAML
    message("Adding .yaml extension to output file: ", file)
    file <- paste0(file, ".yaml")
    ext <- ".yaml"
  }

  if (ext == ".yaml") {
    str_out <- as.yaml(str_out)
  } else if (ext == ".json") {
    str_out <- toJSON(str_out, pretty = TRUE)
  } else {
    stop("Unable to handle file extension: ", ext)
  }

  writeLines(str_out, file)

  return(invisible(str_out))
}


format_custom_contrasts_for_yaml <- function(custom_mat) {
  if (is.null(custom_mat) || length(custom_mat) == 0L) return(NULL)
  if (is.list(custom_mat) && !is.data.frame(custom_mat)) return(custom_mat)
  if (is.data.frame(custom_mat)) custom_mat <- as.matrix(custom_mat)
  if (!is.matrix(custom_mat) || nrow(custom_mat) == 0L) return(NULL)

  rnames <- rownames(custom_mat)
  if (is.null(rnames)) rnames <- paste0("custom_", seq_len(nrow(custom_mat)))
  cnames <- colnames(custom_mat)

  out <- stats::setNames(vector("list", nrow(custom_mat)), rnames)
  for (ii in seq_len(nrow(custom_mat))) {
    row_vec <- as.numeric(custom_mat[ii, ])
    if (!is.null(cnames) && length(cnames) == length(row_vec)) {
      nz <- which(abs(row_vec) > 1e-12)
      if (length(nz) > 0L) {
        out[[ii]] <- as.list(stats::setNames(row_vec[nz], cnames[nz]))
      } else {
        out[[ii]] <- as.list(stats::setNames(row_vec, cnames))
      }
    } else {
      out[[ii]] <- row_vec
    }
  }
  out
}

format_contrast_spec_for_yaml <- function(cspec, fallback_custom = NULL, is_l1_top = FALSE) {
  if (is.null(cspec)) return(list())
  out <- list()
  if (!is.null(cspec$diagonal)) out$diagonal <- cspec$diagonal

  if (!isTRUE(is_l1_top)) {
    if (!is.null(cspec$cell_means)) out$cell_means <- cspec$cell_means
    if (!is.null(cspec$cond_means) && length(cspec$cond_means) > 0L) out$cond_means <- cspec$cond_means
    if (!is.null(cspec$pairwise_diffs) && length(cspec$pairwise_diffs) > 0L) out$pairwise_diffs <- cspec$pairwise_diffs
    if (!is.null(cspec$interaction_diffs) && length(cspec$interaction_diffs) > 0L) out$interaction_diffs <- cspec$interaction_diffs
    if (!is.null(cspec$overall_response)) out$overall_response <- cspec$overall_response
    if (!is.null(cspec$simple_slopes) && length(cspec$simple_slopes) > 0L) out$simple_slopes <- cspec$simple_slopes
    if (!is.null(cspec$slope_diffs) && length(cspec$slope_diffs) > 0L) out$slope_diffs <- cspec$slope_diffs
    if (!is.null(cspec$weights)) out$weights <- cspec$weights
    if (!is.null(cspec$factor_coding) && length(cspec$factor_coding) > 0L) out$factor_coding <- as.list(cspec$factor_coding)
  }

  if (!is.null(cspec$delete) && length(cspec$delete) > 0L) {
    out$delete <- cspec$delete
  }

  custom_out <- format_custom_contrasts_for_yaml(cspec$custom %||% fallback_custom)
  if (!is.null(custom_out) && length(custom_out) > 0L) {
    out$custom <- custom_out
  }

  out
}

#' Function to output l1_models into '.yaml' object
#' @param gpa the \code{gpa} object
#' @return a yaml document of l1_models for an object
#' @details 
#'   This function exports the configuration of level 1 (run-level) GLM settings
#'   in a glm_pipeline_arguments (gpa) object to a YAML file. This allows the user to
#'   then setup their models from this configuration file, rather than having to go through
#'   the `build_l1_models` step manually.
#' 
#' @examples
#' \dontrun{
#'   gpa <- setup_glm_pipeline(analysis_name="flanker_test",
#'     output_directory="/proj/mnhallqlab/no_backup/flanker_gnomes",
#'     n_expected_runs = 2, tr = 2.0,
#'     run_data = run_df, subject_data = subject_df, trial_data = trial_df,
#'     l1_models=NULL, l2_models=NULL, l3_models=NULL
#'   )
#' 
#'  # setup level 1 models through menu system
#'  gpa <- build_l1_models(gpa)
#' 
#'  # get level 1 model settings as a list
#'  l1_settings <- get_l1_config(gpa)
#' 
#'  # we should typically use export_glm_config to write to file, but we could convert manually
#'  l1_yaml <- as.yaml(l1_settings)
#'  writeLines(l1_yaml, "my_l1_config.yaml")
#' 
#'  # in future, we can populate the models using
#'  gpa <- build_l1_models(gpa, from_spec_file="my_l1_config.yaml")
#' }
#' @keywords internal
get_l1_config <- function(gpa) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  l1m <- gpa$l1_models

  if (is.null(l1m)) {
    warning("Cannot export l1 config because $l1_models is NULL. Run build_l1_models first?")
    return(character(0))
  }

  if (!checkmate::test_class(l1m, "l1_model_set")) {
    stop("$l1_models is not of class l1_model_set")
  }

  ## Variable mappings
  str_out <- list(
    onsets = l1m$onsets, # Onset columns
    durations = l1m$durations, # Duration columns
    isis = l1m$isis, # ISI/ITI columns
    wi_factors = l1m$wi_factors, # within-subject factors
    values = l1m$values # parametric modulators
  )

  # export events
  for (ee in l1m$events) {
    eobj <- list()
    checkmate::assert_string(ee$name)
    checkmate::assert_string(ee$onset)
    if (!checkmate::test_string(ee$duration)) {
      checkmate::assert_number(ee$duration, lower = 0)
    }
    checkmate::assert_string(ee$isi, null.ok = TRUE)
    eobj$onset <- ee$onset
    eobj$duration <- ee$duration
    eobj$isi <- ee$isi
    str_out$events[[ee$name]] <- eobj # populate event list
  }

  # export signals
  for (ss in l1m$signals) {
    sobj <- list()
    sobj$event <- ss$event
    sobj$trial_subset_expression <- ss$trial_subset_expression
    if (!is.null(ss$normalization)) {
      checkmate::assert_subset(ss$normalization, c("none", "evtmax_1", "durmax_1"))
      sobj$normalization <- ss$normalization
    }

    if (!is.null(ss$parametric_modulator)) {
      stopifnot(ss$parametric_modulator %in% l1m$values)
      sobj$parametric_modulator <- ss$parametric_modulator
      sobj$value_type <- "parametric"
      sobj$value_fixed <- NULL
    } else if (!is.null(ss$value_fixed)) {
      checkmate::assert_number(ss$value_fixed)
      sobj$value_fixed <- ss$value_fixed
      if (abs(ss$value_fixed - 1) < 1e-5) {
        sobj$value_type <- "unit"
      } else {
        sobj$value_type <- "number"
      }
    }

    # within-subject formula
    if (!is.null(ss$wi_formula)) {
      sobj$wi_formula <- ss$wi_formula
    }

    # within-subject factors
    if (!is.null(ss$wi_factors)) {
      sobj$wi_factors <- ss$wi_factors
    }

    if (!is.null(ss$wi_expansion)) {
      sobj$wi_expansion <- ss$wi_expansion
    }

    if (!is.null(ss$reference_level) && length(ss$reference_level) > 0L) {
      sobj$reference_level <- as.list(ss$reference_level)
    }

    if (!is.null(ss$factor_coding) && length(ss$factor_coding) > 0L) {
      sobj$factor_coding <- as.list(ss$factor_coding)
    }

    # other signal flags
    if (!is.null(ss$demean_convolved)) {
      checkmate::assert_flag(ss$demean_convolved)
      sobj$demean_convolved <- ss$demean_convolved
    }

    if (!is.null(ss$add_deriv)) {
      checkmate::assert_flag(ss$add_deriv)
      sobj$add_deriv <- ss$add_deriv
    }

    if (!is.null(ss$keep_duplicate_occurrences)) {
      checkmate::assert_flag(ss$keep_duplicate_occurrences)
      sobj$keep_duplicate_occurrences <- ss$keep_duplicate_occurrences
    }

    if (!is.null(ss$beta_series)) {
      checkmate::assert_flag(ss$beta_series)
      sobj$beta_series <- ss$beta_series
    }

    if (!is.null(ss$convmax_1)) {
      checkmate::assert_flag(ss$convmax_1)
      sobj$convmax_1 <- ss$convmax_1
    }

    if (!is.null(ss$ts_multiplier) && !checkmate::test_logical(ss$ts_multiplier)) {
      checkmate::assert_string(ss$ts_multiplier)
      sobj$ts_multiplier <- ss$ts_multiplier
    }

    str_out$signals[[ss$name]] <- sobj
  }

  # L1 models
  for (mm in l1m$models) {
    mobj <- list()
    mobj$signals <- mm$signals
    if (!is.null(mm$execution_backend)) mobj$execution_backend <- mm$execution_backend
    mobj$contrasts <- format_contrast_spec_for_yaml(
      mm$contrast_spec,
      fallback_custom = mm$contrast_list$custom,
      is_l1_top = TRUE
    )

    if (!is.null(mm$wi_models)) {
      mobj$contrasts$wi_contrasts <- list()
      for (cc in names(mm$wi_models)) {
        mobj$contrasts$wi_contrasts[[cc]] <- format_contrast_spec_for_yaml(
          mm$wi_models[[cc]]$contrast_spec,
          fallback_custom = mm$wi_models[[cc]]$contrast_list$custom,
          is_l1_top = FALSE
        )
      }
    }

    str_out$l1_models[[mm$name]] <- mobj
  }

  return(invisible(str_out))
}

#' Helper function to output l2_models into '.yaml' object
#' @param gpa the \code{gpa} object
#' @importFrom yaml as.yaml
#' @return a yaml document of l2_models for an object
#' @keywords internal
get_l2_config <- function(gpa) {

  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  # capture whether we are building level 2 (subject) or level 3 (sample) models
  fname <- utils::tail(as.character(match.call()[[1]]), 1L)
  if (fname == "get_l2_config") {
    model_set <- gpa$l2_models
    level <- 2L #model level inside model object
  } else {
    model_set <- gpa$l3_models
    level <- 3L
  }

  if (is.null(model_set)) {
    warning(glue("Cannot export level {level} config because $l{level}_models is NULL. Run build_l{level}_models first?"))
    return(character(0))
  }

  if (!checkmate::test_class(model_set, "hi_model_set")) {
    stop(glue::glue("$l{level}_models is not of class hi_model_set"))
  }

  str_out <- list()

  for (mm in model_set$models) {
    if (!inherits(mm, "hi_model_spec")) {
      warning(glue::glue("In export_l{level}_config, model is not of class 'hi_model_spec'.")) # would be nice if this were more informative, but we don't know the class
      next # skip
    }

    mobj <- list()
    mobj$level <- mm$level
    mobj$model_formula <- as.character(mm$model_formula)
    mobj$num2fac <- mm$num2fac # variables that should be converted from numbers to factors before contrasts are setup
    if (!is.null(mm$execution_backend)) mobj$execution_backend <- mm$execution_backend
    
    # convert named vectors for covariate transformation and reference levels to lists for YAML to capture the key:value pairs
    if (!is.null(mm$covariate_transform)) mobj$covariate_transform <- as.list(mm$covariate_transform) # settings for covariate transformation
    if (!is.null(mm$reference_level)) mobj$reference_level <- as.list(mm$reference_level) # settings for covariate transformation
    if (!is.null(mm$factor_coding)) mobj$factor_coding <- as.list(mm$factor_coding)
    if (level == 2L && !is.null(mm$l2_scope)) mobj$l2_scope <- mm$l2_scope
    if (level == 3L && !is.null(mm$l3_input_mode)) mobj$l3_input_mode <- mm$l3_input_mode
    if (level == 3L && !is.null(mm$producer_backend)) mobj$producer_backend <- mm$producer_backend

    mobj$contrasts <- format_contrast_spec_for_yaml(
      mm$contrast_spec,
      fallback_custom = mm$contrast_list$custom,
      is_l1_top = FALSE
    )

    if (level == 3L) {
      mobj$fsl_outlier_deweighting <- mm$fsl_outlier_deweighting
      
      # 3dLMEr-specific fields
      if (identical(mm$l3_input_mode, "3dlmer")) {
        mobj$random_effects <- mm$random_effects
        mobj$lmer_formula <- mm$lmer_formula
        mobj$lmer_mask <- mm$lmer_mask
        mobj$lmer_njobs <- mm$lmer_njobs
        mobj$lmer_glt_codes <- mm$lmer_glt_codes
      }
    }

    str_out[[mm$name]] <- mobj

  }

  if (level == 2L) {
    str_out <- list(l2_models=str_out) # wrap inside $l2_models top-level element
  } else if (level == 3L) {
    str_out <- list(l3_models=str_out) # wrap inside $l3_models top-level element
  } else {
    stop("Unknown level in get_l2_config")
  }

  return(invisible(str_out))
}

#' Get a list containing configuration settings for l3 models
#' @rdname get_l3_config
#' @name get_l3_config
#' @keywords internal
get_l3_config <- get_l2_config
