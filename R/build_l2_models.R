#' Interactive function to build an l2 model specification for setup_glm_pipeline
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing an analysis pipeline to which $l2_models
#'   should be added. If $l2_models is already present, these will be amended.
#' @param from_spec_file optional YAML or JSON file containing settings to populated into l2 models
#' @param regressor_cols an optional character vector of columns in \code{data} that should be considered
#'   as possible regressors
#'
#' @return a \code{hi_model_set} object containing a list of higher-level regression models
#' @author Michael Hallquist
#' @importFrom checkmate assert_data_frame assert_class assert_subset
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom dplyr select
#' @importFrom stats setNames
#' @importFrom utils menu select.list tail
#' @export
#'
build_l2_models <- function(gpa, from_spec_file = NULL, regressor_cols = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  # capture whether we are building level 2 (subject) or level 3 (sample) models
  fname <- utils::tail(as.character(match.call()[[1]]), 1L)
  if (fname == "build_l2_models") {
    menu_desc <- "second-level (subject)"
    model_set <- gpa$l2_models
    level <- 2L # model level inside model object
    out_field <- "l2_models" # where to put results in gpa
  } else {
    menu_desc <- "third-level (sample)"
    model_set <- gpa$l3_models
    level <- 3L
    out_field <- "l3_models" # where to put results in gpa
  }

  lg <- lgr::get_logger(paste0("glm_pipeline/l", level, "_setup"))
  lg$set_threshold(gpa$lgr_threshold)

  if (level == 2L) {
    data <- compose_l2_model_data(gpa, lg = lg) # L2 pooled modeling frame: run + session predictors
  } else {
    data <- gpa$subject_data # L3
    if (!is.null(gpa$session_data) && is.data.frame(gpa$session_data) &&
        all(c("id", "session") %in% names(data)) && all(c("id", "session") %in% names(gpa$session_data))) {
      extra_sess_cols <- setdiff(names(gpa$session_data), names(data))
      if (length(extra_sess_cols) > 0L) {
        sdat_sub <- gpa$session_data[, c("id", "session", extra_sess_cols), drop = FALSE]
        sdat_sub$session <- if (is.integer(data$session)) as.integer(sdat_sub$session) else as.character(sdat_sub$session)
        data <- dplyr::left_join(data, sdat_sub, by = c("id", "session"))
      }
    }
  }

  # allow deferred model specification upstream
  if (!is.null(model_set) && model_set[1L] == "prompt") model_set <- NULL

  checkmate::assert_class(model_set, "hi_model_set", null.ok = TRUE)
  checkmate::assert_subset(regressor_cols, names(data)) # make sure all parametric regressor columns are in the data frame

  possible_cols <- names(data)
  possible_cols <- possible_cols[!possible_cols %in% "mr_dir"]
  data <- data %>% dplyr::select(all_of(possible_cols))

  if (!is.null(regressor_cols)) { # select only columns requested by user
    data <- data %>% dplyr::select(dplyr::all_of(regressor_cols))
  }

  if (!is.null(from_spec_file)) {
    checkmate::assert_file_exists(from_spec_file)
    lg$info("In build_l%d_models, populated fields from: %s", level, from_spec_file)
    ext <- file_ext(from_spec_file, withdot = FALSE)
    if (ext == "yaml") {
      spec_list <- yaml::read_yaml(from_spec_file)
    } else if (ext == "json") {
      spec_list <- jsonlite::read_json(from_spec_file, simplifyVector = TRUE)
    } else {
      msg <- sprintf("Cannot understand how to input spec file: %s", from_spec_file)
      lg$error(msg)
      stop(msg)
    }

    prompt_models <- FALSE
    if (level == 3L) {
      gpa$fwe_correction <- lapply(spec_list$fwe_correction %||% list(), as_fwe_spec)
      spec_list$fwe_correction <- NULL
    }

    # at present, l2 and l3 model specs don't necessarily have subsections -- each top-level section can be a model
    # given this, standardize the format so that if $l2_models or $l3_models are top-level sections, then
    # spec_list becomes only the relevant part for building l2/l3 models

    # first, for degenerate files, add $l2_models or $l3_models as top-level
    if (is.null(spec_list[[out_field]])) {
      spec_list <- list(spec_list) %>% setNames(out_field)
    }
    # if (level == 2L && is.null(spec_list$l2_models)) {
    #   spec_list <- list(l2_models = spec_list)
    # } else if (level == 3L && is.null(spec_list$l3_models)) {
    #   spec_list <- list(l3_models = spec_list)
    # }

    # next, populate model names
    spec_list <- propagate_spec_names(spec_list) # ensure that $name field is populated for each model

    # finally, subset spec_list down to just hte relevant part
    spec_list <- spec_list[[out_field]]

    nm <- names(spec_list)
    if (any(duplicated(nm))) {
      msg <- sprintf(
        "Duplicate %s model names in %s. Cannot proceed. Resolve duplicates: %s",
        out_field, from_spec_file, paste(unique(nm[duplicated(nm)]), collapse = ", ")
      )
      lg$error(msg)
      stop(msg)
    }

    # if this changes in future, amend propagate_spec_names to work for l2/l3
    # populate the $name field to make it easier for downstream functions to work
    # for (mm in seq_along(spec_list)) {
    #  spec_list[[mm]]$name <- names(spec_list)[mm]
    # }

    ## initialize overall higher-level model set object
    model_set <- list(models = NULL)
    class(model_set) <- c("list", "hi_model_set")

    for (mm in seq_along(spec_list)) {
      # convert covariate_transform and reference_level back to named vectors for object consistency
      # these are converted to lists during the export so that YAML/JSON preserves the key:value pairs
      spec_list[[mm]]$covariate_transform <- unlist(spec_list[[mm]]$covariate_transform)
      spec_list[[mm]]$reference_level <- unlist(spec_list[[mm]]$reference_level)
      spec_list[[mm]]$factor_coding <- unlist(spec_list[[mm]]$factor_coding)
      model_set$models[[names(spec_list)[mm]]] <- create_new_hi_model(data, level = level, cur_model_names = NULL, spec_list = spec_list[[mm]], lg = lg)
    }
  } else {
    prompt_models <- TRUE
  }

  if (is.null(model_set)) {
    ## initialize overall higher-level model set object
    model_set <- list(models = NULL)
    class(model_set) <- c("list", "hi_model_set")
  }

  # regressor manager
  get_regressors <- function(data, regressor_cols = NULL) {
    done_regressors <- FALSE
    while (isFALSE(done_regressors)) {
      cat("Current regressors for this model:\n\n  ", paste(regressor_cols, collapse = ", "), "\n\n")

      # display session-level column annotations when available
      var_origin <- attr(data, "l2_var_origin")
      if (!is.null(var_origin)) {
        session_vars <- names(var_origin)[var_origin == "session"]
        if (length(session_vars) > 0L) {
          cat("  Session-level columns (for between-session effects): ",
            paste(session_vars, collapse = ", "), "\n\n", sep = "")
        }
      }

      action <- menu(c("Add/modify regressors", "Delete regressors", "Done with regressor selection"),
        title = "Would you like to modify the model regressors?"
      )

      if (action == 1L) { # Add/modify
        regressor_cols <- select.list(names(data),
          multiple = TRUE, preselect = regressor_cols,
          title = "Choose all model regressors\n(Command/Control-click to select multiple)"
        )
      } else if (action == 2L) { # Delete
        if (length(regressor_cols) == 0L) {
          message("No regressors yet. Please add at least one.")
        } else {
          which_del <- menu(regressor_cols, title = "Which regressor would you like to remove?")
          if (which_del > 0) {
            proceed <- menu(c("Proceed", "Cancel"),
              title = paste0("Are you sure you want to delete ", regressor_cols[which_del], "?")
            )
            if (proceed == 1) {
              cat("  Deleting ", regressor_cols[which_del], "\n")
              regressor_cols <- regressor_cols[-which_del]
            } else {
              cat("  Not deleting ", regressor_cols[which_del], "\n")
            }
          }
        }
      } else if (action == 3L) { # Done
        if (length(regressor_cols) == 0L) {
          message("No regressors yet. Please add at least one.")
        } else {
          done_regressors <- TRUE
          cat("The following columns were chosen as model regressors.\n\n")
          cat("  ", paste(regressor_cols, collapse = ", "), "\n\n")
        }
      }
    }
    return(regressor_cols)
  }

  ### -- BUILD MODELS ---

  summarize_l2_models <- function(ml) {
    if (length(ml) == 0L) {
      return(invisible(NULL))
    }
    lapply(seq_along(ml), function(ii) {
      this <- ml[[ii]]
      cat("--------\nModel ", ii, "\n\n")
      cat("  Name:", this$name, "\n")
      cat("  Regressors:", paste(this$model_variables, collapse = ", "), "\n")
      if (level == 3L) cat("  Final inference:", if (is.null(this$inference)) "none" else
        paste0(this$inference$method$name, "; FWE alpha ", paste(this$inference$method$fwe_alpha, collapse = ", "),
          "; each map separately"), "\n")
      if (ncol(this$contrasts) < 30) {
        cat("  Contrasts:\n\n")
        print(round(this$contrasts, 3))
        cat("\n--------\n\n")
      } else {
        cat("More than 30 regressor columns. Not printing contrasts\n")
      }
    })
  }

  model_list <- model_set$models
  add_more <- ifelse(isTRUE(prompt_models), 1, 4)
  while (add_more != 4) {
    summarize_l2_models(model_list)

    add_more <- menu(c("Add model", "Modify model", "Delete model", paste("Done with", menu_desc, "model setup")),
      title = paste(sub("^(\\w{1})", "\\U\\1", menu_desc, perl = TRUE), "model setup menu")
    )

    if (add_more == 1L) { # add
      mobj <- create_new_hi_model(data, level = level, cur_model_names = names(model_list), lg = lg)
      if (level == 3L) mobj <- configure_l3_inference(gpa, mobj)
      if (mobj$name %in% names(model_set)) {
        warning("An model with the same name exists: ", mobj$name, ". Overwriting it.")
      }
      model_list[[mobj$name]] <- mobj # add to set
    } else if (add_more == 2L) { # modify
      if (is.null(model_list)) {
        message("No models available to modify. Add at least one model first")
      } else {
        res <- 0L
        while (res == 0L) {
          res <- menu(names(model_list), title = "Which model do you want to modify?")
        }
        model_list[[res]] <- create_new_hi_model(data, to_modify = model_list[[res]], cur_model_names = names(model_list), level = level, lg = lg)
        if (level == 3L) model_list[[res]] <- configure_l3_inference(gpa, model_list[[res]])
      }
    } else if (add_more == 3L) { # delete
      which_del <- menu(names(model_list), title = "Which model would you like to delete?")
      if (which_del > 0) {
        proceed <- menu(c("Proceed", "Cancel"),
          title = paste0("Are you sure you want to delete ", names(model_list)[which_del], "?")
        )
        if (proceed == 1) {
          cat("  Deleting ", names(model_list)[which_del], "\n")
          model_list[[which_del]] <- NULL
        } else {
          cat("  Not deleting ", names(model_list)[which_del], "\n")
        }
      }
    }
  }

  model_set$models <- model_list
  model_set$n_contrasts <- sapply(model_list, function(mobj) { nrow(mobj$contrasts) })

  if (fname == "build_l2_models") {
    gpa$l2_models <- model_set
  } else if (fname == "build_l3_models") {
    gpa$l3_models <- model_set
  }

  has_l2_models <- checkmate::test_class(gpa$l2_models, "hi_model_set") &&
    !is.null(gpa$l2_models$models) &&
    length(gpa$l2_models$models) > 0L
  has_l3_models <- checkmate::test_class(gpa$l3_models, "hi_model_set") &&
    !is.null(gpa$l3_models$models) &&
    length(gpa$l3_models$models) > 0L

  afni_direct_inputs <- FALSE
  if (fname == "build_l3_models" && has_l3_models &&
      identical(normalize_backend_strings(gpa$glm_software), "afni")) {
    requirements <- resolve_model_l3_requirements(gpa)
    afni_direct_inputs <- nrow(requirements) > 0L && all(requirements$producer_backend == "afni")
  }
  if (fname == "build_l3_models" && !isTRUE(has_l2_models) && !afni_direct_inputs) {
    msg <- paste0(
      "L3 models were built before any L2 models were defined. ",
      "L2/L3 signature compatibility cannot be fully validated yet; ",
      "run build_l2_models() and compatibility checks will be applied automatically."
    )
    backend_notes <- character(0)
    if ("spm" %in% gpa$glm_software) {
      backend_notes <- c(
        backend_notes,
        "For SPM, L1 setup also skips projected L2 run/session regressors until L2 models are defined."
      )
    }
    if ("fsl" %in% gpa$glm_software) {
      backend_notes <- c(
        backend_notes,
        "For FSL, L2/L3 pairing will be enforced once L2 models are defined."
      )
    }
    if (length(backend_notes) > 0L) {
      msg <- paste(msg, paste(backend_notes, collapse = " "))
    }
    lg$warn(msg)
    warning(msg)
  }

  if (isTRUE(has_l2_models) && isTRUE(has_l3_models)) {
    pair_df <- enumerate_l2_l3_signature_pairs(gpa, lg = lg)
    incompatible <- pair_df[!pair_df$compatible, , drop = FALSE]
    compatible <- pair_df[pair_df$compatible, , drop = FALSE]

    if (nrow(incompatible) > 0L) {
      detail <- format_l2_l3_incompatibilities(incompatible, max_items = 8L)
      if (nrow(compatible) == 0L) {
        msg <- paste0(
          "All L2/L3 model combinations are currently incompatible by signature. ",
          detail
        )
        lg$warn(msg)
        warning(msg)
      } else {
        msg <- paste0(
          "Some L2/L3 model combinations are incompatible by signature and will be skipped during L3 setup. ",
          detail
        )
        lg$info(msg)
        message(msg)
      }
    }
  }

  return(gpa)
}

#' Interactive function to build an l3 model specification for setup_glm_pipeline
#' @param gpa a \code{glm_pipeline_arguments} object.
#' @param from_spec_file optional YAML or JSON file containing L3 model settings.
#' @param regressor_cols optional character vector of columns to consider as regressors.
#' @rdname build_l3_models
#' @name build_l3_models
#' @export
build_l3_models <- build_l2_models


#' Choose higher-level model regressors interactively
#' @param data data.frame containing candidate regressors
#' @param regressor_cols current selected regressor columns
#' @return character vector of selected regressor columns
#' @importFrom utils menu select.list
#' @noRd
get_regressors <- function(data, regressor_cols = NULL) {
  done_regressors <- FALSE
  while (isFALSE(done_regressors)) {
    cat("Current regressors for this model:\n\n  ", paste(regressor_cols, collapse = ", "), "\n\n")

    var_origin <- attr(data, "l2_var_origin")
    if (!is.null(var_origin)) {
      session_vars <- names(var_origin)[var_origin == "session"]
      if (length(session_vars) > 0L) {
        cat(
          "  Session-level columns (for between-session effects): ",
          paste(session_vars, collapse = ", "), "\n\n", sep = ""
        )
      }
    }

    action <- menu(c("Add/modify regressors", "Delete regressors", "Done with regressor selection"),
      title = "Would you like to modify the model regressors?"
    )

    if (action == 1L) {
      regressor_cols <- select.list(names(data),
        multiple = TRUE, preselect = regressor_cols,
        title = "Choose all model regressors\n(Command/Control-click to select multiple)"
      )
    } else if (action == 2L) {
      if (length(regressor_cols) == 0L) {
        message("No regressors yet. Please add at least one.")
      } else {
        which_del <- menu(regressor_cols, title = "Which regressor would you like to remove?")
        if (which_del > 0) {
          proceed <- menu(c("Proceed", "Cancel"),
            title = paste0("Are you sure you want to delete ", regressor_cols[which_del], "?")
          )
          if (proceed == 1) {
            cat("  Deleting ", regressor_cols[which_del], "\n")
            regressor_cols <- regressor_cols[-which_del]
          } else {
            cat("  Not deleting ", regressor_cols[which_del], "\n")
          }
        }
      }
    } else if (action == 3L) {
      if (length(regressor_cols) == 0L) {
        message("No regressors yet. Please add at least one.")
      } else {
        done_regressors <- TRUE
        cat("The following columns were chosen as model regressors.\n\n")
        cat("  ", paste(regressor_cols, collapse = ", "), "\n\n")
      }
    }
  }
  regressor_cols
}


#' Internal function to walk through l2/l3 model setup, including populated from a specification (YAML) file
#' @param data a data.frame containing data that can be included in the model
#' @param to_modify an existing `hi_model_spec` object to be modified
#' @param level the model level: 2 or 3
#' @param cur_model_names a character vector of the names of all current models in the model list
#' @param spec_list a list containing model specifications read in from a YAML/JSON file
#' @param lg the current logger
#' @details if `spec_list` is passed in, the resulting model object will be populated from the fields
#'   of the specification list, rather than prompting the user for input
#' @return a `hi_model_spec` object containing the L2/L3 model to be added to gpa$l2_models or gpa$l3_models
#' @importFrom checkmate assert_integerish assert_list assert_class test_data_table
#' @importFrom stats alias coef model.matrix relevel terms update.formula
#' @importFrom utils capture.output menu tail
#' @keywords internal
create_new_hi_model <- function(data, to_modify = NULL, level = NULL, cur_model_names = NULL, spec_list = NULL, lg = NULL,
                                retry_name = NULL) {
    if (checkmate::test_data_table(data)) {
      data <- as.data.frame(data) # make subsetting syntax in this function consistent with standard data.frame conventions
    }

    checkmate::assert_class(lg, "Logger", null.ok = TRUE)
    if (is.null(lg)) lg <- lgr::get_logger("fmri.pipeline/hi_model")
    checkmate::assert_class(to_modify, "hi_model_spec", null.ok = TRUE)
    if (is.null(to_modify)) {
      mobj <- list(level = level)
      class(mobj) <- c("list", "hi_model_spec")
      modify <- FALSE
    } else {
      mobj <- to_modify
      modify <- TRUE
    }

    checkmate::assert_integerish(level, lower = 2, upper = 3, len = 1L)
    checkmate::assert_character(cur_model_names, null.ok = TRUE)
    checkmate::assert_list(spec_list, null.ok = TRUE)
    checkmate::assert_string(retry_name, null.ok = TRUE)

    if (level == 2) {
      id_cols <- c("id", "session", "run_number")
    } else {
      id_cols <- c("id", "session")
    }

    choose_l2_scope <- function(default_scope = "id_session") {
      allowed <- longitudinal_l2_scopes()
      checkmate::assert_subset(default_scope, allowed)
      if (!"session" %in% names(data) || length(unique(data$session)) <= 1L) {
        return(default_scope)
      }
      n_sessions <- length(unique(data$session))
      cat(
        "\nYour data contain ", n_sessions, " sessions per subject.\n",
        "This choice determines whether between-session contrasts are available at L2.\n\n",
        sep = ""
      )
      opts <- c(
        "Within-session: combine runs within each session separately (between-session contrasts deferred to L3)",
        "Across-session: stack all runs and sessions together (between-session contrasts available at L2)"
      )
      which_scope <- menu(
        opts,
        title = "How should this L2 model combine sessions?"
      )
      chosen <- if (which_scope == 2L) "id" else "id_session"
      if (identical(chosen, "id_session")) {
        cat(
          "\nNote: With within-session scope, between-session differences cannot be tested at L2.\n",
          "If you need between-session contrasts at L2, re-run model setup and choose across-session scope.\n\n",
          sep = ""
        )
      }
      chosen
    }

    choose_l3_input_mode <- function(default_mode = "per_session") {
      allowed <- longitudinal_l3_input_modes()
      default_mode <- normalize_l3_input_mode(default_mode)
      checkmate::assert_subset(default_mode, allowed)
      if (!"session" %in% names(data) || length(unique(data$session)) <= 1L) {
        return(default_mode)
      }
      opts <- c(
        "Per-session L3 model (requires L2 scope id_session)",
        "Pooled sessions with subject EVs (requires L2 scope id_session)",
        "One row per subject input (requires L2 scope id)",
        "3dLMEr (mixed-effects longitudinal model via AFNI; requires L2 scope id_session)"
      )
      which_mode <- menu(
        opts,
        title = "How should this L3 model use sessions?"
      )
      if (which_mode == 2L) {
        "pooled_sessions_subject_ev"
      } else if (which_mode == 3L) {
        "subject_rows"
      } else if (which_mode == 4L) {
        "3dlmer"
      } else {
        "per_session"
      }
    }

    handle_aliased_model <- function(fit_mobj, model_formula) {
      if (is.null(fit_mobj$aliased_terms)) {
        return(list(mobj = fit_mobj, retry = FALSE))
      }

      alias_msg <- format_hi_model_alias_message(fit_mobj, model_formula = as.character(model_formula))

      if (!is.null(spec_list) || !rlang::is_interactive()) {
        stop(alias_msg, call. = FALSE)
      }

      cat("\n", alias_msg, "\n\n", sep = "")
      retry <- menu(
        c("Respecify model", "Cancel model setup"),
        title = "How would you like to proceed?"
      )
      if (retry == 1L) {
        return(list(mobj = fit_mobj, retry = TRUE))
      }
      stop("Model setup cancelled because the fixed-effects design contains aliased coefficients.", call. = FALSE)
    }

    ### ------ model name ------
    if (isTRUE(modify)) {
      cat("Current model name:", mobj$name, "\n")
      res <- menu(c("No", "Yes"), title = "Change model name?")
      if (res == 2) mobj$name <- NULL # clear out so that it is respecified
    }

    if (!is.null(spec_list$name)) { # populate from specification, if requested
      mobj$name <- spec_list$name
    }
    if (!is.null(retry_name)) {
      mobj$name <- retry_name
    }

    if (!is.null(spec_list$execution_backend)) {
      mobj$execution_backend <- normalize_backend_strings(spec_list$execution_backend)
    }
    if (level == 3L && !is.null(spec_list$producer_backend)) {
      mobj$producer_backend <- normalize_backend_strings(spec_list$producer_backend)
    }

    while (is.null(mobj$name) || mobj$name == "") {
      res <- trimws(readline("Enter the model name: "))
      if (res != "") {
        res <- make.names(res)
        if (res %in% cur_model_names) {
          cat("\nModel name:", res, "already exists. Names must be unique.\n")
          cat("Current models:", paste(cur_model_names, collapse = ", "), "\n")
        } else {
          mobj$name <- res
        }
      }
    }

    if (isTRUE(modify)) {
      cat("Current model regressors:", paste(names(mobj$model_variables), collapse = ", "), "\n")
      res <- menu(c("No", "Yes"), title = "Respecify model regressors (and contrasts)?")
      if (res == 2) { # clear out so that it is respecified
        cat("Okay, resetting model regressors and contrasts\n")
        mobj$model_variables <- mobj$contrasts <- NULL
      }
    }

    # Setup lm formula: let user decide LM formula approach versus walkthrough
    if (!is.null(spec_list$model_formula)) {
      model_approach <- 1L
      model_formula <- as.formula(spec_list$model_formula) # maybe need better error trapping here
    } else {
      model_approach <- 0L
      model_formula <- NULL
    }
    
    while (model_approach == 0L) {
      model_approach <- menu(c("Specify formula", "Model builder"),
        title = "Do you want to specify the model formula or walk through the model builder?"
      )
    }

    if (model_approach == 1L) {
      # formula approach
      while (is.null(model_formula)) {
        cat(c(
          "\nSpecify the right-hand side of the model you wish to fit.",
          "Use variable names in the dataset provided to this function.",
          "Note that this syntax follows standard R model syntax. See ?lm for details.",
          "Example: ~ emotion * wmload + run_number\n",
          "For an intercept-only model, specify: ~1\n"
        ), sep = "\n")

        # annotate columns by origin (run-level vs session-level)
        var_origin <- attr(data, "l2_var_origin")
        if (!is.null(var_origin) && level == 2L) {
          run_vars <- setdiff(names(var_origin)[var_origin == "run"], id_cols)
          session_vars <- names(var_origin)[var_origin == "session"]
          if (length(run_vars) > 0L) {
            cat("Run-level columns:\n")
            cat(strwrap(paste(run_vars, collapse = ", "), 70, exdent = 5), sep = "\n")
          }
          if (length(session_vars) > 0L) {
            cat("Session-level columns (useful for between-session contrasts):\n")
            cat(strwrap(paste(session_vars, collapse = ", "), 70, exdent = 5), sep = "\n")
          }
        } else {
          cat("Available column names:\n")
          cat(strwrap(paste(names(data), collapse = ", "), 70, exdent = 5), sep = "\n")
        }

        res <- trimws(readline("Enter the model formula: "))
        # always trim any LHS specification
        res <- sub("^[^~]*", "", res, perl = TRUE)
        model_formula <- tryCatch(as.formula(res), error = function(e) {
          print(e)
          cat("Problem converting your syntax to formula. Try again\n")
          return(NULL)
        })
      }

      mobj$model_variables <- all.vars(model_formula)

      # if not an intercept-only model, ensure that all variables in formula are present in the data frame
      if (!identical(model_formula, ~1)) {
        checkmate::assert_subset(mobj$model_variables, names(data), empty.ok = FALSE)
      }

    } else if (model_approach == 2L) {
      # walkthrough approach (only support additive model for now)
      mobj$model_variables <- get_regressors(data, regressor_cols = mobj$model_variables)
      model_formula <- as.formula(paste("~", paste(mobj$model_variables, collapse = " + ")))
    }

    # convert integer-like variables to integer
    for (vv in mobj$model_variables) {
      # convert integer-like numbers to integers for possible conversion below
      if (is.numeric(data[[vv]]) && checkmate::test_integerish(data[[vv]])) {
        data[[vv]] <- as.integer(data[[vv]])
      }

      if (is.character(data[[vv]])) {
        cat("Converting character variable to factor for model setup: ", vv, "\n", sep = "")
        data[[vv]] <- as.factor(data[[vv]])
      }

      if (is.integer(data[[vv]]) && length(unique(data[[vv]])) < 20L) {
        if (!is.null(spec_list)) {
          conv <- if (vv %in% spec_list$num2fac) 1L else 2L
        } else if (vv %in% spec_list$num2fac) {
          conv <- 1L
        } else {
          cat(
            "\n",
            vv, " appears to be an integer with < 20 levels.\n",
            "If you convert it to a factor, R will generate dummy code variables in the model.\n",
            "If you leave it as a number, R will include the linear trend in the model.\n\n",
            sep = ""
          )
          conv <- menu(c("Convert to factor", "Leave as number"),
            title = paste("Do you want to convert", vv, "to a factor?")
          )
        }
        if (conv == 1L) {
          cat("Converting", vv, "into a factor.\n")
          data[[vv]] <- as.factor(data[[vv]])
          mobj$num2fac <- unique(c(mobj$num2fac, vv))
        }
      }

      if (is.factor(data[[vv]])) {
        if (is.null(mobj$factor_levels)) mobj$factor_levels <- list()
        mobj$factor_levels[[vv]] <- levels(data[[vv]])
      }
    }

    # Fit a provisional model before asking about centering or reference levels.
    # Those choices do not fix rank deficiency, so aliased formulas should fail
    # or return to model specification before asking follow-up questions.
    if (!is.null(spec_list)) {
      mobj$covariate_transform <- spec_list$covariate_transform
      mobj$reference_level <- spec_list$reference_level
      mobj$factor_coding <- spec_list$factor_coding
    }
    alias_check <- handle_aliased_model(
      mobj_fit_lm(mobj, as.character(model_formula), data, id_cols, lg = lg, report_alias = FALSE),
      model_formula
    )
    if (isTRUE(alias_check$retry)) {
      return(create_new_hi_model(
        data = data,
        to_modify = NULL,
        level = level,
        cur_model_names = setdiff(cur_model_names, alias_check$mobj$name),
        spec_list = NULL,
        lg = lg,
        retry_name = alias_check$mobj$name
      ))
    }
    mobj <- alias_check$mobj

    # need to build a model LM-style
    if (!is.null(spec_list)) {
      mobj$covariate_transform <- spec_list$covariate_transform
      mobj$reference_level <- spec_list$reference_level
      mobj$factor_coding <- spec_list$factor_coding
    } else if (length(mobj$model_variables) > 0L) {
      cat("Summary of variables included in model:\n\n")
      print(summary(data[, mobj$model_variables, drop = FALSE]))
      cat("\n")

      # handle mean centering and reference levels
      cont_vars <- sapply(data[, mobj$model_variables, drop=FALSE], class) %in% c("integer", "numeric")
      if (any(cont_vars)) {
        cat(
          "We will now ask you to indicate whether to transform any of the continuous covariates in the model.",
          "The most common transformation, which is often a good default, is to mean-center the covariate.",
          "Other options are to subtract the minimum (so that zero now represents the lowest covariate value),",
          "subtract the maximum (so that zero represents the highest covariate value),",
          "or to z-score/standardized the covariate to a mean of zero and standard deviation of one.",
          sep = "\n  "
        )
        cont_vars <- mobj$model_variables[cont_vars == TRUE]
        mobj$covariate_transform <- c() # reset any previous centering settings
        for (cc in cont_vars) {
          cat("\n")
          resp <- menu(
            c("Mean center", "Standardize/z-score", "Subtract minimum", "Subtract maximum", "No transformation"),
            title = paste0("Transform ", cc, "?")
          )
          if (resp == 1L) {
            cat(glue("\nMean centering {cc}. Mean centering is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "mean"
          } else if (resp == 2L) {
            cat(glue("\nStandardizing {cc}. Standardization is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "zscore"
          } else if (resp == 3L) {
            cat(glue("\nSubtracting minimum value of {cc}. This is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "min"
          } else if (resp == 4L) {
            cat(glue("\nSubtracting maximum value of {cc}. This is recomputed if runs or subjects are dropped in L2 or L3 analyses.\n", .trim = FALSE))
            mobj$covariate_transform[cc] <- "max"
          }
        }
      }

      # handle reference levels for factors
      cat_vars <- sapply(data[, mobj$model_variables, drop=FALSE], class) %in% c("factor")
      if (any(cat_vars)) {
        cat_vars <- mobj$model_variables[cat_vars == TRUE]
        cat(
          "We will now ask you to choose the reference level for dummy coding factors in the model.",
          "This affects how contrasts are coded by emmeans and it controls how you interpret the Intercept maps",
          "at levels 2 and 3, if you look at these. In general, this choice should not make a big difference,",
          "but we default to the factor level that is most common in the dataset.", sep="\n  "
        )
        for (cc in cat_vars) {
          # default to cell with largest N
          default_lev <- names(which.max(table(data[[cc]])))
          levs <- levels(data[[cc]])
          levs <- levs[c(which(levs == default_lev), which(levs != default_lev))]
          labs <- levs
          labs[1L] <- paste0("Recommended: ", labs[1])
          which_lev <- menu(labs, title = glue("Choose the reference level for the {cc} factor."))
          if (!exists("which_lev") || which_lev == 0L) {
            which_lev <- 1L # choose default on cancel or ctrl+c
          }
          mobj$reference_level[cc] <- levs[which_lev]
        }

      }
    }

    # populate longitudinal signature fields used to derive L2/L3 compatibility
    if (level == 2L) {
      if (!is.null(spec_list$l2_scope)) {
        checkmate::assert_string(spec_list$l2_scope)
        checkmate::assert_subset(spec_list$l2_scope, longitudinal_l2_scopes())
        mobj$l2_scope <- spec_list$l2_scope
      } else if (is.null(mobj$l2_scope)) {
        if (is.null(spec_list)) {
          mobj$l2_scope <- choose_l2_scope(default_scope = "id_session")
        } else {
          mobj$l2_scope <- "id_session"
        }
      } else {
        checkmate::assert_subset(mobj$l2_scope, longitudinal_l2_scopes())
      }
    } else if (level == 3L) {
      if (!is.null(spec_list$l3_input_mode)) {
        checkmate::assert_string(spec_list$l3_input_mode)
        checkmate::assert_subset(spec_list$l3_input_mode, longitudinal_l3_input_modes())
        mobj$l3_input_mode <- spec_list$l3_input_mode
      } else if (is.null(mobj$l3_input_mode)) {
        if (is.null(spec_list)) {
          mobj$l3_input_mode <- choose_l3_input_mode(default_mode = "per_session")
        } else {
          mobj$l3_input_mode <- "per_session"
        }
      } else {
        checkmate::assert_subset(mobj$l3_input_mode, longitudinal_l3_input_modes())
      }

      # Inform user about auto-injection of subject indicator EVs
      if (identical(mobj$l3_input_mode, "pooled_sessions_subject_ev")) {
        lg$info(
          "L3 model '%s' uses pooled_sessions_subject_ev. Per-subject indicator EVs will be auto-injected at setup time if 'id' is not in the model formula.",
          mobj$name
        )
      }
    }

    if (level == 2L) {
      var_origin <- attr(data, "l2_var_origin")
      if (!is.null(var_origin)) {
        session_vars <- names(var_origin)[var_origin == "session"]
        used_session_vars <- intersect(mobj$model_variables, session_vars)
        if (length(used_session_vars) > 0L && identical(mobj$l2_scope, "id")) {
          lg$info(
            "L2 model '%s' includes session-level predictors in pooled id scope. These effects are run-weighted across sessions: %s",
            mobj$name, paste(used_session_vars, collapse = ", ")
          )
        } else if (length(used_session_vars) > 0L && identical(mobj$l2_scope, "id_session")) {
          lg$warn(
            "L2 model '%s' includes session-level predictors with l2_scope='id_session'. These terms are typically aliased within-session and may be dropped: %s",
            mobj$name, paste(used_session_vars, collapse = ", ")
          )
        }
      }
    }

    # fit linear model and populate model object
    # For 3dLMEr, we still fit a standard lm to allow emmeans to work for gltCode generation
    mobj <- mobj_fit_lm(mobj, as.character(model_formula), data, id_cols, lg = lg, report_alias = FALSE)
    alias_check <- handle_aliased_model(mobj, model_formula)
    if (isTRUE(alias_check$retry)) {
      return(create_new_hi_model(
        data = data,
        to_modify = NULL,
        level = level,
        cur_model_names = setdiff(cur_model_names, alias_check$mobj$name),
        spec_list = NULL,
        lg = lg,
        retry_name = alias_check$mobj$name
      ))
    }
    mobj <- alias_check$mobj

    # 3dLMEr-specific: collect lme4-style random-effects specification
    if (level == 3L && identical(mobj$l3_input_mode, "3dlmer")) {
      if (!is.null(spec_list$random_effects)) {
        mobj$random_effects <- spec_list$random_effects
      } else if (!is.null(spec_list$lmer_formula)) {
        mobj$lmer_formula <- spec_list$lmer_formula
      } else {
        cat(
          "\n--- 3dLMEr Random Effects ---\n",
          "Specify the random-effects part of the lme4 formula.\n",
          "Example: (1 | Subj) for random intercepts per subject\n",
          "Example: (1 + Session | Subj) for random intercepts and random slopes\n",
          "Note: Using 'Subj' as the grouping factor is recommended for 3dLMEr compatibility.\n\n",
          sep = ""
        )
        res <- trimws(readline("Enter the random-effects specification: "))
        if (res == "") res <- "(1 | Subj)" # sensible default
        mobj$random_effects <- res
      }

      # Combine fixed and random effects into full lmer formula for AFNI unless lmer_formula was explicitly set without random_effects.
      if (!is.null(mobj$random_effects) || is.null(mobj$lmer_formula)) {
        fixed_rhs <- tail(as.character(model_formula), 1L)
        fixed_rhs <- sub("^~\\s*", "", fixed_rhs, perl = TRUE)
        mobj$lmer_formula <- paste(fixed_rhs, "+", mobj$random_effects)
      }
      mobj$lmer_mask <- spec_list$lmer_mask
      mobj$lmer_njobs <- spec_list$lmer_njobs
      mobj$lmer_glt_codes <- spec_list$lmer_glt_codes

      # Preserve any grouping/random-slope variables referenced only in lmer_formula (e.g., site in (1|site/Subj))
      lmer_all_vars <- tryCatch(
        all.vars(stats::as.formula(paste("~", sub("^~\\s*", "", mobj$lmer_formula)))),
        error = function(e) character(0)
      )
      extra_re_vars <- setdiff(lmer_all_vars, c("Subj", "id", "session", mobj$model_variables))
      if (length(extra_re_vars) > 0L) {
        checkmate::assert_subset(extra_re_vars, names(data))
        mobj$random_effects_vars <- extra_re_vars
        mobj$model_variables <- unique(c(mobj$model_variables, extra_re_vars))
        # Re-fit so model_data / complete-case filtering includes random-effects grouping columns
        mobj <- mobj_fit_lm(mobj, as.character(model_formula), data, id_cols, lg = lg, report_alias = FALSE)
      }

      cat("Full 3dLMEr formula: ", mobj$lmer_formula, "\n")
    }

    # Guided longitudinal prompt: offer between-session pairwise contrasts when applicable
    preselected_pw <- character(0)
    has_multi_session <- "session" %in% names(data) && length(unique(data$session)) > 1L
    if (level == 2L && has_multi_session && identical(mobj$l2_scope, "id") && is.null(spec_list)) {
      var_origin <- attr(data, "l2_var_origin")
      session_vars <- if (!is.null(var_origin)) names(var_origin)[var_origin == "session"] else character(0)
      # identify session-level factor or integer variables in the model that could yield between-session contrasts
      session_contrast_candidates <- intersect(mobj$model_variables, session_vars)
      # also check for session-related factor columns in the fitted model
      if (length(session_contrast_candidates) > 0L) {
        cat(
          "\n--- Between-session contrasts ---\n",
          "Your L2 model includes session-level variables: ",
          paste(session_contrast_candidates, collapse = ", "), "\n",
          "Between-session pairwise differences are commonly desired in longitudinal designs.\n\n",
          sep = ""
        )
        for (sc in session_contrast_candidates) {
          if (is.factor(data[[sc]]) || (is.integer(data[[sc]]) && length(unique(data[[sc]])) < 20L)) {
            add_pw <- menu(
              c("Yes", "No"),
              title = paste0("Add pairwise between-session contrasts for '", sc, "'?")
            )
            if (add_pw == 1L) {
              preselected_pw <- unique(c(preselected_pw, sc))
              cat("  Added pairwise differences for '", sc, "'.\n", sep = "")
            }
          }
        }
      }
    }

    # walk through contrast generation for this model
    mobj <- specify_contrasts(
      mobj,
      spec_list = spec_list,
      preselected_pairwise = if (length(preselected_pw) > 0L) preselected_pw else NULL
    )

    # ask about outlier deweighting
    if (level == 3L) {
      if (!is.null(spec_list)) {
        mobj$fsl_outlier_deweighting <- spec_list$fsl_outlier_deweighting
      } else {
        cat(
          "",
          "Do you want to turn on outlier deweighting in FSL level 3 analysis? This is available for any mixed effects analysis,",
          "most commonly FLAME1 or FLAME1+2. In general, this is a good idea, but we have seen it break down for some legitimate inputs,",
          "generating many complaints about excessive outliers being detected. If you get those messages, you should probably disable this.",
          "",
          sep = "\n  "
        )

        res <- menu(c("No", "Yes"), title = "Turn on outlier deweighting?")
        if (res == 2L) {
          mobj$fsl_outlier_deweighting <- TRUE
        } else {
          mobj$fsl_outlier_deweighting <- FALSE
        }
      }
    }

    if (!is.null(spec_list$inference) && level != 3L) stop("inference is supported only for L3 models.")
    if (level == 3L) {
      inference <- if (!is.null(spec_list)) spec_list$inference else mobj$inference
      mobj$inference <- normalize_l3_inference(inference, mobj$name, rownames(mobj$contrasts))
    }
    return(mobj)
  }
