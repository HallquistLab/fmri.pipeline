#' Interactive function to build an l1 model specification for setup_glm_pipeline
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing an analysis pipeline to which $l1_models
#'   should be added. If $l1_models is already present, these will be amended.
#' @param trial_data a data.frame containing trial-level data for one or more subjects
#' @param ppi_data an optional data.frame containing physiological time-series inputs for PPI signal construction
#' @param l1_model_set optional existing l1_model_set to be modified
#' @param from_spec_file optional YAML or JSON file containing settings to populated into l1 models
#' @param onset_cols an optional character vector of columns in \code{trial_data} that should be
#'   in the set of event onsets
#' @param onset_regex an optional PCRE-compatible regular expression for identifying potential
#'   event onset columns in \code{trial_data}
#' @param duration_cols an optional character vector of columns in \code{trial_data} that should be
#'   in the set of event durations
#' @param duration_regex an optional PCRE-compatible regular expression for identifying potential
#'   event duration columns in \code{trial_data}
#' @param value_cols an optional character vector of columns in \code{trial_data} that should be in the set of signal values
#' @param value_regex an optional PCRE-compatible regular expression for identifying potential
#'   event value columns in \code{trial_data}
#' @param isi_cols an optional character vector of columns in \code{trial_data} that should be in the set of signal isi/iti
#' @param isi_regex an optional PCRE-compatible regular expression for identifying potential
#'   isi/iti columns in \code{trial_data}
#'
#' @details if \code{gpa} is not passed in, then we will work from trial_data and l1_model_set.
#'
#' @return a \code{l1_model_set} object containing events, signals, and models, compatible with build_design_matrix
#' @author Michael Hallquist
#' @importFrom checkmate assert_data_frame assert_class assert_subset
#' @importFrom dplyr bind_cols select mutate
#' @importFrom yaml read_yaml
#' @importFrom jsonlite read_json
#' @importFrom stats coef model.matrix update.formula
#' @importFrom utils head menu select.list
#' @export
build_l1_models <- function(gpa=NULL, trial_data=NULL, ppi_data = NULL, l1_model_set=NULL, from_spec_file=NULL,
                           onset_cols=NULL, onset_regex=".*(onset|time).*",
                           duration_cols=NULL, duration_regex=".*duration.*",
                           value_cols=NULL, value_regex=NULL, isi_cols=NULL, isi_regex="^(iti|isi).*") {

  # Maybe allow glm object to be passed in that would have trial_data and variable_mapping.
  # I guess that would be like "add_l1_model"
  lg <- lgr::get_logger("glm_pipeline/build_l1_models")
  lg$set_threshold(gpa$lgr_threshold)

  # TODO: we should almost definitely eliminate the use case where gpa is null... 
  checkmate::assert_class(gpa, "glm_pipeline_arguments", null.ok = TRUE)
  if (!is.null(gpa)) {
    lg$info("In build_l1_models, using existing gpa object to build l1 models (ignoring trial_data argument etc.)")
    use_gpa <- TRUE
    trial_data <- gpa$trial_data
    l1_model_set <- gpa$l1_models
    ppi_data <- gpa$ppi_data
  } else {
    lg$debug("In build_l1_models, using trial_data passed in, rather than gpa object")
    use_gpa <- FALSE
  }

  checkmate::assert_data_frame(trial_data, null.ok = FALSE)
  checkmate::assert_class(l1_model_set, "l1_model_set", null.ok = TRUE)
  checkmate::assert_subset(onset_cols, names(trial_data)) #make sure that all onset columns are in the data frame
  checkmate::assert_string(onset_regex, null.ok = TRUE)
  checkmate::assert_subset(duration_cols, names(trial_data)) # make sure that all duration columns are in the data frame
  checkmate::assert_string(duration_regex, null.ok = TRUE)
  checkmate::assert_subset(value_cols, names(trial_data)) #make sure all parametric regressor columns are in the data frame
  checkmate::assert_string(value_regex, null.ok = TRUE)
  checkmate::assert_subset(c("id", "session", "run_number", "trial"), names(trial_data)) # required metadata in trial_data

  if (is.null(l1_model_set)) {
    ## initialize overall l1 design object (holds events, signals, and models)
    l1_model_set <- list(onsets=NULL, durations=NULL, values=NULL, wi_factors=NULL, events=NULL, signals=NULL, models=NULL)
    class(l1_model_set) <- c("list", "l1_model_set")
    new_l1 <- TRUE
  } else {
    new_l1 <- FALSE
  }

  if (!is.null(from_spec_file)) {
    checkmate::assert_file_exists(from_spec_file)
    lg$info("In build_l1_models, populated fields from: %s", from_spec_file)
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

    spec_list <- propagate_spec_names(spec_list)

    l1_model_set <- fields_from_spec(l1_model_set, spec_list, trial_data, c("onsets", "durations", "isis", "values", "wi_factors"))
    l1_model_set <- bl1_build_events(l1_model_set, spec_list, trial_data, lg)
    l1_model_set <- signals_from_spec(l1_model_set, spec_list, trial_data, lg)
    l1_model_set <- bl1_build_models(l1_model_set, spec_list, lg)
    new_l1 <- FALSE # always assume that the spec file has enough info not to walk through each stage
    blm1_complete <- TRUE # drop out of model builder when working from a spec file
  } else {
    blm1_complete <- FALSE
  }

  int_vars <- sapply(trial_data, checkmate::test_integerish)
  char_vars <- sapply(trial_data, checkmate::test_character)
  fac_vars <- sapply(trial_data, checkmate::test_factor)
  num_vars <- sapply(trial_data, checkmate::test_numeric)
  possible_factors <- names(which(int_vars | char_vars | fac_vars)) # possible factors must be int, char, or factor
  possible_factors <- setdiff(possible_factors, c("id", "session", "run_number", "trial"))
  possible_pms <- names(which(int_vars | num_vars)) # possible parametric modulators must be numeric
  possible_pms <- setdiff(possible_pms, c("id", "session", "run_number", "trial"))

  # relies on scope of parent function for onset_cols etc.
  take_l1_actions <- function(l1_model_set, actions) {
    checkmate::assert_integerish(actions, lower = 1, upper = 8)
    for (aa in actions) {
      if (aa == 1) {
        # onsets
        l1_model_set <- bl1_get_cols(l1_model_set, trial_data,
          field_name = "onsets", field_desc = "onset",
          select_cols = onset_cols, select_regex = onset_regex
        )
      } else if (aa == 2) {
        # durations
        l1_model_set <- bl1_get_cols(l1_model_set, trial_data,
          field_name = "durations", field_desc = "duration",
          select_cols = duration_cols, select_regex = duration_regex
        )
      } else if (aa == 3) {
        # iti/isi
        l1_model_set <- bl1_get_cols(l1_model_set, trial_data,
          field_name = "isis", field_desc = "ISI/ITI",
          select_cols = isi_cols, select_regex = isi_regex, force_selection=FALSE
        )
      } else if (aa == 4) {
        # parametric values
        l1_model_set <- bl1_get_cols(l1_model_set, trial_data,
          field_name = "values", field_desc = "parametric value",
          select_cols = value_cols, select_regex = value_regex, force_selection=FALSE,
          limit_cols = possible_pms
        )
      } else if (aa == 5) {
        # within-subject factors
        l1_model_set <- bl1_get_cols(l1_model_set, trial_data,
          field_name = "wi_factors", field_desc = "within-subject factor", force_selection=FALSE,
          limit_cols = possible_factors
        )
      } else if (aa == 6) {
        # events
        l1_model_set <- bl1_build_events(l1_model_set, trial_data=trial_data, lg=lg)
      } else if (aa == 7) {
        # signals
        l1_model_set <- bl1_build_signals(l1_model_set, trial_data, ppi_data = ppi_data, lg=lg)
      } else if (aa == 8) {
        # models
        l1_model_set <- bl1_build_models(l1_model_set, lg=lg)
      }
    }

    return(l1_model_set)
  }

  if (isTRUE(new_l1)) {
    cat(
      "Welcome to the l1 model builder. This process will walk you through setting up all level 1 fMRI models for you analyses.",
      "This process consists of eight steps: ",
      "  1. Selection of event onset columns. Event onsets must be in seconds relative to the scan start.",
      "  2. Selection of event duration columns. Event durations must be in seconds relative to the scan start.",
      "  3. Optional selection of interstimulus/intertrial intervals. These must be ISI/ISI durations in seconds.",
      "  4. Selection of all parametric modulator (continuous) event values.",
      "  5. Optional selection of within-subject factor columns that can be used in specifying signals.",
      "  6. Build 'events' which consist of onset times, durations, and optional ITI/ISI.",
      "  7. Build 'signals', which consist of an event, an event value (amplitude), and convolution and regressor settings.",
      "  8. Build 'models', which consist of a set of signals and GLM contrasts for signal-related regressors.",
      sep = "\n"
    )

    l1_model_set <- take_l1_actions(l1_model_set, 1:8) # run user through each step in sequence
  } else {
    cat("Modifying existing l1 model structure\n")
  }

  # drop into l1 model menu, whether new or modification, unless working from a spec file
  while (isFALSE(blm1_complete)) {
    blm1_action <- menu(c(
      "Update event onset columns",                                             # 1
      "Update event duration columns",                                          # 2
      "Update ITI/ISI columns",                                                 # 3
      "Update parametric value columns",                                        # 4
      "Update within-subject factor columns",                                   # 5
      "Update events (each consists of an onset and duration)",                 # 6
      "Update signals (event, factor, parametric value, regressor settings)",   # 7
      "Update models (consisting of signals and contrasts)",                    # 8
      "Done with level 1 model building"                                        # 9
    ), title = "Level 1 model builder")

    if (blm1_action < 9) {
      l1_model_set <- take_l1_actions(l1_model_set, blm1_action)
    } else if (blm1_action == 9) {
      cat("Completing level 1 model building\n")
      blm1_complete <- TRUE
      break
    }
  }

  if (isTRUE(use_gpa)) {
    gpa$l1_models <- l1_model_set
    gpa <- refresh_l1_cope_names(gpa, lg = lg)
    return(gpa)
  } else {
    return(l1_model_set)
  }
}

#' Onset, duration, value column selection helper function
#' @param l1_model_set an l1_model_set object containing onsets etc.
#' @param trial_data a trial + subjects x events data.frame that contains potential onset columns
#' @param field_name the element of \code{l1_model_set} containing columns of a certain purpose (onset, duration, value)
#' @param field_desc the text description of the field being modified (e.g., 'parametric value')
#' @param select_cols a character vector of current columns specified by the user to be added/included
#' @param select_regex a PCRE-compatible regular expression for identifying columns
#' @param limit_cols a vector of variable names in \code{trial_data} that constrains what can be chosen
#' @param force_selection do not allow an empty return for this field
#' @param alpha_sort whether to display eligible columns in alphabetical order
#' @param prompt_input whether to ask user to confirm selections
#' 
#' @return a modified version of l1_model_set that has \code{field_name} updated according to user specification
#' @importFrom glue glue
#' @keywords internal
bl1_get_cols <- function(l1_model_set, trial_data, field_name = NULL, field_desc = NULL, select_cols = NULL, select_regex = NULL,
  limit_cols=NULL, force_selection = TRUE, alpha_sort = TRUE, prompt_input=TRUE) {
  # record of columns before any adjustments.
  current_cols <- l1_model_set[[field_name]]

  # never allow id, session, or run_number as columns
  all_cols <- names(trial_data)
  if (!is.null(limit_cols)) {
    possible_cols <- limit_cols
  } else {
    possible_cols <- setdiff(all_cols, c("id", "session", "run_number"))
  }

  if (isTRUE(alpha_sort)) possible_cols <- sort(possible_cols)

  new_cols <- setdiff(select_cols, current_cols) # any new columns in the argument compared to the l1_model_set?
  chosen_cols <- current_cols # start with current columns
  if (length(new_cols) > 0L) {
    cat(glue("Current {field_desc} columns in the l1 model structure are: {c_string(current_cols)}\n", .trim=FALSE))
    cat(glue("The arguments to build_l1_models also included: {c_string(new_cols)}\n", .trim=FALSE))
    if (isFALSE(prompt_input)) {
      res <- 1L
    } else {
      res <- menu(c("Yes", "No"), title = glue("Do you want to add these columns to possible {field_desc}s?"))
    }

    if (res == 1L) {
      # add chosen_cols to l1 set
      chosen_cols <- c(current_cols, new_cols) # else just keep current_cols
    }
  }

  # handle regex detection of new columns
  if (!is.null(select_regex)) {
    detected_cols <- grep(select_regex, possible_cols, value = TRUE, perl = TRUE)
    uniq_detect <- setdiff(detected_cols, chosen_cols) # only bother the user if the regex reveals new columns
    if (length(detected_cols) > 0L && length(uniq_detect) > 0L) {
      cat(glue("\n\n---\nDetected the following possible event {field_desc} columns:\n\n  {c_string(detected_cols)}\n\n"))
      res <- menu(c("Yes", "No"), title = glue("Do you want to add these columns to possible {field_desc}s?"))
      if (res == 1) {
        # add these to any that were manually specified/current
        chosen_cols <- unique(c(chosen_cols, detected_cols))
      }
    }
  }

  done_cols <- FALSE
  while (isFALSE(done_cols)) {
    if (isFALSE(prompt_input)) {
      action <- 3L # Done
    } else {
      cat(glue("\n\n---\nCurrent {field_desc} columns: {c_string(chosen_cols)}\n\n", .trim=FALSE))
      action <- menu(c(
        glue("Add/modify {field_desc} columns"),
        glue("Delete {field_desc} columns"),
        glue("Done with {field_desc} selection")
      ),
      title = glue("Would you like to modify the event {field_desc} columns?")
      )
    }

    if (action == 1L) { # Add/modify
      # if user selects 0 to cancel, select.list still returns character(0) instead of the current values
      new_cols <- select.list(possible_cols,
        multiple = TRUE, preselect = chosen_cols,
        title = glue("Choose all columns denoting event {field_desc}s\n(Command/Control-click to select multiple)")
      )
      if (!identical(new_cols, character(0))) chosen_cols <- new_cols
    } else if (action == 2L) { # Delete
      if (length(chosen_cols) == 0L && isTRUE(force_selection)) {
        cat(glue("No {field_desc} yet. Please add at least one.\n"))
      } else if (length(chosen_cols) == 1L && isTRUE(force_selection)) {
        cat(glue("Cannot delete the last {field_desc}. Add others before deleting this.\n"))
      } else {
        which_del <- menu(chosen_cols, title = glue("Which {field_desc} column would you like to remove?"))
        if (which_del > 0) {
          proceed <- menu(c("Proceed", "Cancel"),
            title = glue("Are you sure you want to delete {chosen_cols[which_del]}?")
          )

          if (proceed == 1) {
            cat("  Deleting ", chosen_cols[which_del], "\n")
            chosen_cols <- chosen_cols[-which_del]

            # TODO: need mechanism for cascade deleting events, signals, and models that depend on this field!
          } else {
            cat("  Not deleting ", chosen_cols[which_del], "\n")
          }
        }
      }
    } else if (action == 3L) { # Done
      done_cols <- TRUE
      cat(glue("\nThe following columns were chosen as possible event {field_desc}s.\n", .trim=FALSE))
      cat(glue("  {c_string(chosen_cols)}\n", .trim=FALSE))
    }
  }

  l1_model_set[[field_name]] <- chosen_cols
  return(l1_model_set)
}

#' helper function to build events consisting of onsets and durations
#' @param l1_model_set an l1_model_set object that may have extant events in it
#' @param spec_list an optional parsed YAML/JSON model specification used to populate events non-interactively
#' @param trial_data the trial_data object from the \code{gpa} object
#' @param lg an lgr logger used for status and warning messages
#' @return a modified copy of l1_model_set with events added/updated
#' @keywords internal
#' @importFrom dplyr select mutate
#' @importFrom checkmate test_number
#' @importFrom glue glue
bl1_build_events <- function(l1_model_set, spec_list = NULL, trial_data, lg=NULL) {
  cat("Specify all events that can be added to a GLM model. Events consist of an onset time and duration\n")

  summarize_events <- function(l1_model_set) {
    cat("Summary of events available in l1 models:\n--------\n")
    if (!is.null(l1_model_set$events)) {
      lapply(l1_model_set$events, function(x) {
        cat("Event name:     ", x$name, "\n")
        cat("  Event onset:    ", x$onset, "\n")
        cat(
          "    mean [min -- max]: ",
          round(mean(x$data$onset, na.rm = TRUE), 2), "[",
          round(min(x$data$onset, na.rm = TRUE), 2), "--",
          round(max(x$data$onset, na.rm = TRUE), 2), "]\n"
        )
        cat("    First 6 values:", paste(head(x$data$onset), collapse = ", "), "\n\n")
        cat("  Event duration: ", x$duration, "\n")
        cat(
          "    mean [min -- max]: ",
          round(mean(x$data$duration, na.rm = TRUE), 2), "[",
          round(min(x$data$duration, na.rm = TRUE), 2), "--",
          round(max(x$data$duration, na.rm = TRUE), 2), "]\n"
        )
        cat("    First 6 values:", paste(round(head(x$data$duration), 2), collapse = ", "), "\n\n")
        if (!is.null(x$isi)) {
          cat("  Event ISI/ITI:    ", x$isi, "\n")
          cat(
            "    mean [min -- max]: ",
            round(mean(x$data$isi, na.rm = TRUE), 2), "[",
            round(min(x$data$isi, na.rm = TRUE), 2), "--",
            round(max(x$data$isi, na.rm = TRUE), 2), "]\n"
          )
        }
      })
    } else {
      cat("No events in l1 models\n")
    }
  }

  if (!is.null(spec_list)) {
    l1_model_set <- events_from_spec(l1_model_set, spec_list, trial_data)
    prompt_input <- FALSE
  } else {
    prompt_input <- TRUE
  }

  events_complete <- FALSE
  while (isFALSE(events_complete)) {
    summarize_events(l1_model_set)
    if (prompt_input) {
      action <- menu(c("Add event", "Delete event", "Done with event specification"))
    } else {
      action <- 3L
    }

    if (action == 1L) { # add
      eobj <- list()

      # ---- event name ----
      complete <- FALSE
      while (isFALSE(complete)) {
        prompt <- "Enter the event name: "
        nm <- readline(prompt)
        if (nm != "") {
          if (nm %in% names(l1_model_set$events)) {
            cat("Event cannot have the same name as an existing event!\n")
          } else {
            eobj$name <- nm
            complete <- TRUE
          }
        }
      }

      # ---- event onset ----
      complete <- FALSE
      while (isFALSE(complete)) {
        oo <- menu(l1_model_set$onsets, title="Choose event onset")
        if (oo != 0) {
          eobj$onset <- l1_model_set$onsets[oo]
          complete <- TRUE
        }
      }

      # ---- event duration ----
      complete <- FALSE
      while (isFALSE(complete)) {
        choices <- c("Specify fixed duration", l1_model_set$durations)
        oval <- menu(choices, title = "Choose event duration")
        if (oval == 1) {
          duration <- NULL
          while (!checkmate::test_number(duration, lower = 0, upper = 5000)) {
            duration <- as.numeric(readline(paste0("Enter the duration value (in seconds) for ", eobj$name, ": ")))
          }
          if (duration > 50) {
            lg$warn("Duration more than 50s specified. Make sure that your durations are in seconds, not milliseconds!")
          }

          eobj$duration <- duration
          complete <- TRUE
        } else if (oval > 1L) {
          eobj$duration <- choices[oval]
          complete <- TRUE
        }
      }

      # optional selection of ISI if chosen in initial setup
      choose_isi <- menu(c("Yes", "No"), title = "Do you want to specify an ITI or ISI that follows this event?")
      if (choose_isi == 1L) {
        complete <- FALSE
        while (isFALSE(complete)) {
          choices <- c("Specify fixed ISI/ITI", l1_model_set$isis)
          oval <- menu(choices, title = "Choose event ISI/ITI")
          if (oval == 1) {
            isi <- NULL
            while (!checkmate::test_number(isi, lower = 0, upper = 5000)) {
              isi <- as.numeric(readline(paste0("Enter the ISI/ITI value (in seconds) for ", eobj$name, ": ")))
            }
            if (isi > 50) {
              lg$warn("ISI/ITI more than 50s specified. Make sure that your ISI/ITI values are in seconds, not milliseconds!")
            }

            eobj$isi <- isi
            complete <- TRUE
          } else if (oval > 1L) {
            eobj$isi <- choices[oval]
            complete <- TRUE
          }
        }
      }

      # populate data frame for event
      eobj <- populate_event_data(eobj, trial_data)
      class(eobj) <- c("list", "l1_model_set_events") # set 'l1_model_set_events' class
      l1_model_set$events[[nm]] <- eobj
    } else if (action == 2L) { # delete
      event_names <- names(l1_model_set$events)
      which_del <- menu(event_names, title="Which event would you like to delete?")
      if (which_del > 0) {
        proceed <- menu(c("Proceed", "Cancel"), title = glue("Are you sure you want to delete {event_names[which_del]}?"))
        if (proceed==1) {
          cat(glue("  Deleting {event_names[which_del]}\n"))
          l1_model_set$events[[which_del]] <- NULL
        } else {
          cat(glue("  Not deleting {event_names[which_del]}\n"))
        }
      }
    } else if (action == 3L) {
      events_complete <- TRUE
    }
  }

  return(l1_model_set)

}

# helper function to print signal setup
summarize_l1_signals <- function(sl) {
  if (length(sl) == 0L) {
    return(invisible(NULL))
  }
  cat("Summary of signals available in l1 models:\n--------\n")
  lapply(seq_along(sl), function(ii) {
    this <- sl[[ii]]
    cat("--------\nSignal ", ii, "\n\n")
    cat("  Name:", this$name, "\n")
    cat("  Event alignment:", this$event, "\n")
    if (isFALSE(this$trial_subset)) {
      cat("  Trial subset: FALSE", "\n")
    } else if (isTRUE(this$trial_subset)) {
      cat(glue("  Trial subset: {this$trial_subset_expression}"), "\n")
      cat(
        glue("  Proportion of trials included: overall = {round(this$trial_subset_statistics['overall'], 2)}"),
        glue(
          "    By id: mean = {round(this$trial_subset_statistics['mean'], 2)}, ", "sd = {round(this$trial_subset_statistics['sd'], 2)}, ",
          "min = {round(this$trial_subset_statistics['min'], 2)}, ", "max = {round(this$trial_subset_statistics['max'], 2)}"
        ), sep="\n"
      )
    }
    if (this$value_type %in% c("unit", "number")) {
      cat("  Regressor value (constant): ", this$value_fixed[1L], "\n")
    } else {
      cat(
        "  Parametric value: ", this$parametric_modulator, ", mean [min -- max]: ",
        round(mean(this$value$value, na.rm = TRUE), 2), " [ ",
        round(min(this$value$value, na.rm = TRUE), 2), " -- ",
        round(max(this$value$value, na.rm = TRUE), 2), " ] \n", sep=""
      )
    }
    cat("  HRF Normalization:", this$normalization, "\n")
    cat("  Add signal derivative:", as.character(this$add_deriv), "\n")
    cat("  Demean convolved signal:", as.character(this$demean_convolved), "\n")
    cat("  Within-subject factor:", c_string(this$wi_factors), "\n")
    cat("  Generate beta series:", as.character(this$beta_series), "\n")
    cat(
      "  Multiply convolved regressor against time series (PPI):",
      ifelse(is.null(this$ts_multiplier) || isFALSE(this$ts_multiplier),
        "FALSE", this$ts_multiplier
      ), "\n"
    )
    cat("\n--------\n")
  })
}

summarize_l1_models <- function(ml) {
  if (length(ml) == 0L) { return(invisible(NULL)) }
  cat("Summary of l1 models:\n--------\n")
  lapply(seq_along(ml), function(ii) {
    this <- ml[[ii]]
    cat("--------\nModel ", ii, "\n\n")
    cat("  Name:", this$name, "\n")
    cat("  Signals:", paste(this$signals, collapse=", "), "\n")
    if (ncol(this$contrasts) < 20) {
      cat("  Contrasts:\n\n")
      print(round(this$contrasts, 3))
      cat("\n--------\n\n")
    } else { cat("More than 20 regressor columns. Not printing contrasts\n") }
  })
}

#' helper function to build level 1 signals
#' @param l1_model_set An \code{l1_model_set} object whose signals should be created or modified
#' @param trial_data A data.frame containing trial-level signal information
#' @param block_data An optional data.frame containing block-level signal information
#' @param ppi_data An optional data.frame containing physiological signals for PPI analysis
#' @param lg an lgr logger used for status and warning messages
#'
#' @return a modified version of \code{l1_model_set} with updated \code{$signals}
#' @keywords internal
bl1_build_signals <- function(l1_model_set, trial_data, block_data = NULL, ppi_data = NULL, lg=NULL) {
  cat("\nNow, we will build up a set of signals that can be included as regressors in the level 1 model.\n")

  checkmate::assert_class(lg, "Logger")

  signal_list <- l1_model_set$signals
  add_more <- 1
  while (add_more != 4) {
    summarize_l1_signals(signal_list)

    add_more <- menu(c("Add signal", "Modify signal", "Delete signal", "Done with signal setup"),
      title = "Signal setup menu"
    )

    if (add_more == 4) {
      break
    } else if (add_more == 3) {
      if (length(signal_list) == 0L) {
        cat("No signals added yet. Please add one first.\n")
      } else {
        which_del <- menu(names(signal_list), title = "Which signal would you like to delete?")
        if (which_del > 0) {
          proceed <- menu(c("Proceed", "Cancel"),
            title = paste0("Are you sure you want to delete ", names(signal_list)[which_del], "?")
          )
          if (proceed == 1L) {
            cat("  Deleting ", names(signal_list)[which_del], "\n")
            signal_list[[which_del]] <- NULL

            # TODO: cascade delete any model that has this signal
          } else {
            cat("  Not deleting ", names(signal_list)[which_del], "\n")
          }
        }
      }
    } else if (add_more %in% c(1, 2)) {
      modify <- FALSE

      # prompt for signal details
      if (add_more == 2) {
        if (length(signal_list) == 0L) {
          cat("No signals available to modify.\n")
          ss <- list()
        } else {
          res <- 0L
          while (res == 0L) {
            res <- menu(names(signal_list), title = "Which signal do you want to modify?")
          }
          ss <- signal_list[[res]]
          signal_list[[res]] <- NULL # clear out old settings
          modify <- TRUE
        }
      } else {
        ss <- list()
      }

      ### ---- signal name ----
      complete <- FALSE
      while (isFALSE(complete)) {
        prompt <- "Enter the signal name: "
        if (isTRUE(modify)) {
          cat("Current signal name:", ss$name, "\n")
          prompt <- "Enter the signal name (press return to leave unchanged): "
        }

        nm <- readline(prompt)
        if (nm != "") {
          ss$name <- nm
          complete <- TRUE
        } else if (nm == "" && isTRUE(modify)) {
          complete <- TRUE
        }
      }

      ### ---- event alignment ----
      if (isTRUE(modify)) {
        cat(glue("Current signal alignment: {c_string(ss$event)}"), sep = "\n")
        res <- menu(c("No", "Yes"), title = "Change signal alignment?")
        if (res == 2) {
          ss$event <- NULL
        } # clear out event so that it is respecified
      }

      while (is.null(ss$event)) {
        res <- menu(names(l1_model_set$events), title = "With which event is this signal aligned?")
        if (res > 0) {
          ss$event <- names(l1_model_set$events)[res]
        }
      }

      ### ---- trial subsetting ----
      if (isTRUE(modify)) {
        cat(glue("Current trial subset: {c_string(ss$trial_subset)}"), sep = "\n")
        res <- menu(c("No", "Yes"), title = "Change trial subset?")
        if (res == 2) {
          ss$trial_subset <- ss$trial_subset_expression <- ss$trial_subset_statistics <- NULL
        } # clear out subset so that it is respecified
      }

      trial_set <- NULL
      while (is.null(ss$trial_subset)) {
        res <- menu(c("No", "Yes"), title = "Only model this signal for specific trials?")
        if (res == 1L) {
          ss$trial_subset <- FALSE
        } else if (res == 2L) {
          while (is.null(ss$trial_subset_expression)) {
            cat("Specify an R-based expression that will be evaluated against trial_data and return TRUE/FALSE for each row.\n\n")
            cat("For example, you might want to separate out trials on which reaction times were implausibly fast (e.g., less than 100 ms).\n")
            cat("  In this case, if the column in trial_data is called rt and it's in seconds, you would type: rt < .1\n")
            cat("  A compound expression can be used, too, like: rt < .5 & accurate == FALSE\n\n")
            ss_expression <- readline("Enter the trial subsetting expression: ")
            trial_set <- tryCatch(with(trial_data, eval(parse(text = ss_expression))),
              error = function(e) {
                lg$error("Problem evaluating trial subsetting expression: %s. Error: %s", ss_expression, as.character(e))
                return(NULL)
              }
            )
            if (is.null(trial_set)) {
              cat("\nProblem with your subsetting expression. Please try again.\n")
            } else if (sum(trial_set) == 0) {
              cat("\nExpression would retain zero trials! Please try again.\n\n")
            } else {
              ss$trial_subset_expression <- ss_expression
              ss$trial_subset <- TRUE
            }
          }
        }
      }

      ss <- get_trial_subset_stats(ss, trial_data, trial_set)

      ### ---- value of regressor ----
      if (isTRUE(modify)) {
        # repopulate value data.frame in case subset has changed, preserving any existing wi_factors
        ss$value <- get_value_df(
          ss, trial_data, wi_factors = ss$wi_factors, event_data = l1_model_set$events[[ss$event]]$data
        )

        cat(
          "Current signal value:",
          ifelse(length(ss$value) == 1L && is.numeric(ss$value[1L]),
            ss$value[1L],
            paste0(
              ifelse(is.null(ss$parametric_modulator), "", paste0(ss$parametric_modulator, ", ")),
              round(mean(ss$value$value, na.rm = TRUE), 2), " [",
              round(min(ss$value$value, na.rm = TRUE), 2), " -- ",
              round(max(ss$value$value, na.rm = TRUE), 2), "]\n"
            )
          ), "\n"
        )
        res <- menu(c("No", "Yes"), title = "Change signal value?")
        if (res == 2) {
          # clear out event so that it is respecified
          ss$value <- ss$value_fixed <- ss$value_type <- ss$parametric_modulator <- NULL
        }
      }

      while (is.null(ss$value)) {
        regtype <- menu(c(
          "Unit height (1.0)", "Other fixed value (will prompt for details)",
          "Parametric modulator (will prompt for details)"
        ),
        title = "What should be the value of regressor (pre-convolution)?"
        )

        if (regtype == 1L) {
          ss$value_type <- "unit"
          ss$value_fixed <- 1.0
        } else if (regtype == 2L) {
          val <- NULL
          while (!checkmate::test_number(val)) {
            val <- as.numeric(readline("Enter the regressor value/height (pre-convolution): "))
          }
          ss$value_type <- "number"
          ss$value_fixed <- val
        } else if (regtype == 3L) {
          val <- 0L
          while (val == 0L) {
            val <- menu(l1_model_set$values, c("Which value should be used for this signal?"))
            if (val > 0) {
              # TODO: have build_design matrix support a simple value vector, which requires
              # same number of rows as metadata (avoid redundancy)
              # basal data frame for each event

              ss$parametric_modulator <- l1_model_set$values[val] # keep column name
              ss$value_type <- "parametric"
            }
          }
        }

        # populate value data frame for this signal, preserving any existing wi_factors
        ss$value <- get_value_df(
          ss, trial_data, wi_factors = ss$wi_factors, event_data = l1_model_set$events[[ss$event]]$data
        )
      }

      ### ---- within-subject factor modulation ----
      ss <- bl1_specify_wi_factors(ss, l1_model_set, trial_data, modify)

      ### ------ hrf normalization ------
      if (isTRUE(modify)) {
        cat("Current HRF normalization:", ss$normalization, "\n")
        res <- menu(c("No", "Yes"), title = "Change HRF normalization?")
        if (res == 2) {
          ss$normalization <- NULL
        } # clear out so that it is respecified
      }

      while (is.null(ss$normalization)) {
        opt <- c("none", "evtmax_1", "durmax_1")
        res <- menu(c(
          "none",
          "evtmax_1 (aka dmUBLOCK(1); HRF max of 1.0 for each event, regardless of duration)",
          "durmax_1 (aka dmUBLOCK; HRF maxing at 1.0 as events become longer (1.0 around 15 sec)"
        ),
        title = "How should the HRF be normalized in convolution?"
        )
        if (res > 0) ss$normalization <- opt[res]
      }

      prompt_advanced <- FALSE
      if (isTRUE(modify)) {
        cat(
          "Current advanced option settings:\n",
          "  - Temporal derivative:", as.character(ss$add_deriv), "\n",
          "  - Demean convolved signal:", as.character(ss$demean_convolved), "\n",
          "  - Beta series:", as.character(ss$beta_series), "\n",
          "  - Time series multiplier:", as.character(ss$ts_multiplier), "\n\n"
        )

        res <- menu(c("No", "Yes"), title = "Change advanced settings?")
        if (res == 2L) { # clear out so that it is respecified
          ss$add_deriv <- ss$demean_convolved <- ss$beta_series <- ss$ts_multiplier <- NULL
          prompt_advanced <- TRUE
        }
      } else {
        cat("Advanced options defaults:",
          "  - No temporal derivative",
          "  - Demean convolved signal",
          "  - No beta series",
          "  - No time series multiplier (PPI)\n",
          sep = "\n"
        )
        accept_defaults <- menu(c("Yes", "No"), title = "Accept default advanced options for this signal?")
        if (accept_defaults == 2L) {
          prompt_advanced <- TRUE
        }
      }

      if (isFALSE(prompt_advanced)) {
        cat("Using defaults for signal: ", ss$name, "\n")
        ss$add_deriv <- FALSE
        ss$demean_convolved <- TRUE
        ss$beta_series <- FALSE
        ss$ts_multiplier <- FALSE
      } else {
        # derivative
        while (is.null(ss$add_deriv)) {
          res <- menu(c("No", "Yes"), title = "Add temporal derivative?")
          if (res == 1L) {
            ss$add_deriv <- FALSE
          } else if (res == 2L) {
            ss$add_deriv <- TRUE
          }
        }

        # demean
        while (is.null(ss$demean_convolved)) {
          res <- menu(c("No", "Yes"), title = "Demean signal post-convolution?")
          if (res == 1L) {
            ss$demean_convolved <- FALSE
          } else if (res == 2L) {
            ss$demean_convolved <- TRUE
          }
        }

        # beta series
        while (is.null(ss$beta_series)) {
          res <- menu(c("No", "Yes"),
            title = "Generate beta series for this signal (one regressor per trial)?"
          )
          if (res == 1L) {
            ss$beta_series <- FALSE
          } else if (res == 2L) {
            ss$beta_series <- TRUE
          }
        }

        if (is.null(ppi_data)) {
          ss$ts_multiplier <- FALSE
        } else {
          while (is.null(ss$ts_multiplier)) {
            res <- menu(c("No", "Yes"),
              title = "Multiply this signal against a physiological regressor for PPI analysis?"
            )
            if (res == 1L) {
              ss$ts_multiplier <- FALSE
            } else if (res == 2L) {
              possible_cols <- setdiff(names(ppi_data), c("id", "session", "run_number", "volume"))
              ppi_col <- select.list(possible_cols,
                multiple = FALSE,
                title = "Which column in the PPI data should multiply this signal?"
              )
              if (ppi_col[1L] != "") ss$ts_multiplier <- ppi_col
            }
          }
        }        
      }

      signal_list[[ss$name]] <- ss
    }
  }

  class(signal_list) <- c("list", "l1_model_set_signals") # set 'l1_model_set_signals' class

  # populate back into model set
  l1_model_set$signals <- signal_list

  return(l1_model_set)
}

############### BUILD MODELS FROM SIGNALS AND EVENTS
#' helper function to create a set of level 1 models
#' @param l1_model_set a set of level 1 models
#' @param spec_list if creating models from a spec (YAML) file, this is the parsed list
#' @param lg the current logger
#' @keywords internal
bl1_build_models <- function(l1_model_set, spec_list=NULL, lg=NULL) {
  checkmate::assert_class(lg, "Logger")

  create_new_model <- function(signal_list, to_modify=NULL, spec_list=NULL) {
    checkmate::assert_class(to_modify, "l1_model_spec", null.ok=TRUE)
    if (is.null(to_modify)) {
      mobj <- list(level = 1L) #first-level model
      class(mobj) <- c("list", "l1_model_spec")
      modify <- FALSE
    } else {
      mobj <- to_modify
      modify <- TRUE
    }

    ### ------ model name ------
    if (isTRUE(modify)) {
      cat("Current model name:", mobj$name, "\n")
      res <- menu(c("No", "Yes"), title="Change model name?")
      if (res == 2) { mobj$name <- NULL } #clear out so that it is respecified
    }

    if (!is.null(spec_list$name)) { # populate from specification, if requested
      mobj$name <- spec_list$name
    }

    if (!is.null(spec_list$execution_backend)) {
      mobj$execution_backend <- normalize_backend_strings(spec_list$execution_backend)
    }

    while (is.null(mobj$name) || mobj$name == "") {
      res <- trimws(readline("Enter the model name: "))
      if (res != "") {
        res <- make.names(res)
        if (res %in% names(model_list)) {
          cat("\nModel name:", res, "already exists. Names must be unique.\n")
          cat("Current models:", paste(names(model_list), collapse=", "), "\n")
        } else {
          mobj$name <- res
        }
      }
    }

    if (isTRUE(modify)) {
      cat("Current model signals:", paste(names(mobj$signals), collapse=", "), "\n")
      res <- menu(c("No", "Yes"), title="Change model signals (and contrasts)?")
      if (res == 2) { #clear out so that it is respecified
        mobj$signals <- mobj$regressors <- mobj$contrasts <- mobj$contrast_spec <- NULL
      }
    }

    ### ------ model signals ------
    if (!is.null(spec_list$signals)) {
      checkmate::assert_subset(spec_list$signals, names(signal_list))
      mobj$signals <- spec_list$signals
    } else {
      summarize_l1_signals(signal_list) # print summary of signals if user input is needed
    }

    while (is.null(mobj$signals)) {
      signals <- select.list(names(signal_list), multiple = TRUE, preselect = mobj$signals,
        title = "Choose all signals to include in this model\n(Command/Control-click to select multiple)")

      if (length(signals) == 0L) {
        proceed <- menu(c("Yes", "No"), title="Nothing entered. Do you want to cancel model setup?")
        if (proceed == 1L) {
          return(invisible(NULL)) #return nothing from function
        }
      } else {
        mobj$signals <- signals
      }
    }

    #look up what the regressors will be for this.
    mobj$regressors_list <- sapply(signal_list[mobj$signals], get_regressors_from_signal, simplify=FALSE)
    mobj$regressors <- unlist(mobj$regressors_list)

    # NOTE: at level 1, we always essentially have an additive model and rely on the user for contrast specification. Usually, this will
    # just be a diagonal matrix. The only major exception I can think of is the within-subject factor situation where the columns of the design
    # matrix depend on each other and need to be modeled as such. So, in this case, we basically want to populate part of the contrast matrix with
    # contrasts for each wi-factor-modulated signal.

    #contrast editor
    mobj <- specify_contrasts(mobj, signal_list, spec_list)

    return(mobj)
  }

  model_list <- l1_model_set$models
  signal_list <- l1_model_set$signals
  add_more <- 1

  if (!is.null(spec_list)) {
    for (mm in spec_list$l1_models) {
      mobj <- create_new_model(signal_list, spec_list=mm)
      model_list[[mobj$name]] <- mobj # add to set
    }

    # quit out of model builder when working from spec file
    done <- TRUE
  } else {
    done <- FALSE
  }

  while (!done) {
    summarize_l1_models(model_list)

    add_more <- menu(c("Add model", "Modify model", "Delete model", "Done with l1 model setup"),
      title="Level 1 model setup menu")

    if (add_more == 1L) { #add
      mm <- create_new_model(signal_list)
      if (!is.null(mm)) { #NULL indicates that user canceled model setup process
        if (mm$name %in% names(l1_model_set)) {
          lg$warn("A model with the same name exists: ", mm$name, ". Overwriting it.")
        }
        model_list[[mm$name]] <- mm #add to set
      }
    } else if (add_more == 2L) { #modify
      if (is.null(model_list)) {
        message("No models available to modify. Add at least one model first")
      } else {
        res <- 0L
        while (res == 0L) { res <- menu(names(model_list), title="Which model do you want to modify?") }
        model_list[[res]] <- create_new_model(signal_list, to_modify=model_list[[res]])
      }
    } else if (add_more == 3L) { #delete
      which_del <- menu(names(model_list), title="Which model would you like to delete?")
      if (which_del > 0) {
        proceed <- menu(c("Proceed", "Cancel"),
          title=paste0("Are you sure you want to delete ", names(model_list)[which_del], "?"))
        if (proceed==1) {
          cat("  Deleting ", names(model_list)[which_del], "\n")
          model_list[[which_del]] <- NULL
        } else {
          cat("  Not deleting ", names(model_list)[which_del], "\n")
        }
      }
    } else if (add_more == 4L) { # Requesting exit
      if (is.null(model_list)) {
        proceed <- menu(c("Proceed", "Cancel"), title="No l1 models have been setup. Are you sure you want to exit?")          
        if (proceed==1) done <- TRUE
      } else {
        done <- TRUE
      }
    }
  }

  # add class 'l1_model_set_models'
  if (!is.null(model_list))  class(model_list) <- c("list", "l1_model_set_models")
  l1_model_set$models <- model_list
  l1_model_set$n_contrasts <- sapply(model_list, function(mm) { nrow(mm$contrasts) })

  return(l1_model_set)
}

# helper function to specify within-subject factor modulation for a given signal
bl1_specify_wi_factors <- function(ss, l1_model_set, trial_data, modify) {
  if (is.null(l1_model_set$wi_factors)) {
    return(ss) # nothing to do
  }

  if (isTRUE(modify)) {
    cat(glue("Current within-subject factor: {c_string(ss$wi_factors)}"), sep = "\n")
    res <- menu(c("No", "Yes"), title = "Change within-subject factor?")
    if (res == 2L) {
      ss$wi_factors <- NULL # clear out factor so that it is respecified
      query_wi <- TRUE
    } else {
      query_wi <- FALSE
    }
  } else {
    query_wi <- TRUE
  }

  if (isFALSE(query_wi)) {
    return(ss) # nothing to do with specifying/changing within-subject factors
  }

  res <- menu(c("No", "Yes"), title = "Is this signal modulated by one or more within-subject factors?")
  if (res == 1L) {
    ss$wi_formula <- ss$wi_ev_formula <- ss$wi_factors <- ss$wi_factor_levels <- ss$wi_model <- NULL # clear out within-subject specifications
    return(ss)
  }

  # if we arrive here, the user has asked us to specify a within-subject factor
  cat(glue("Available within-subject factors: {c_string(l1_model_set$wi_factors)}\n"))

  wi_formula <- NULL
  while (is.null(wi_formula)) {
    cat(c(
      "\nSpecify the right-hand side of the within-subject model you wish to fit for the factors that modulate this signal.",
      "For example, a one-factor model might look like: ~ trustee",
      "And a two-factor additive model might look like: ~ trustee + share",
      "Note that this syntax follows standard R model syntax. See ?lm for details."
    ), sep = "\n")

    wi_formula <- trimws(readline("Enter the model formula: "))

    # always trim any LHS specification
    wi_formula <- sub("^[^~]*", "", wi_formula, perl = TRUE)
    wi_formula <- tryCatch(as.formula(wi_formula), error = function(e) {
      print(e)
      cat("Problem converting your syntax to formula. Try again\n")
      return(NULL)
    })

    wi_vars <- all.vars(wi_formula)
    invalid_vars <- setdiff(wi_vars, l1_model_set$wi_factors)
    if (length(invalid_vars) > 0L) {
      cat(glue("Invalid within-subject factors specified: {c_string(invalid_vars)}. Please try again.\n"))
      wi_formula <- NULL
    }
  }

  # drop intercept from within-subject model to make contrasts more paradigmatic for L1 (i.e., avoid use of grand intercept)
  ss$wi_formula <- as.character(update.formula(wi_formula, ~ . - 1)) # always store formula as character to avoid storing environment
  ss$wi_factors <- wi_vars

  if (!checkmate::test_data_frame(ss$value)) {
    stop("Signal value element is not a data.frame. Within-subject factors only setup to use data.frames right now.")
  }

  # repopulate trial data for this signal, including the selected wi_factors
  ss$value <- get_value_df(
    ss, trial_data, wi_factors = wi_vars,
    event_data = l1_model_set$events[[ss$event]]$data
  )

  # fit dummy model to populate a set of dummy coefficients, then save those to the object
  ss <- fit_wi_model(ss)

  return(ss)
}


get_value_df <- function(signal, trial_data, wi_factors = NULL, event_data = NULL) {
  trial_set <- get_trial_set_from_signal(signal, trial_data)
  value_df <- trial_data %>%
    dplyr::mutate(occurrence_id = dplyr::row_number()) %>%
    dplyr::select("id", "session", "run_number", "trial", "occurrence_id", dplyr::all_of(wi_factors))

  if (!is.null(trial_set)) {
    stopifnot(length(trial_set) == nrow(trial_data))
    checkmate::assert_logical(trial_set)
    value_df <- value_df[trial_set, , drop = FALSE]
  }

  if (signal$value_type %in% c("unit", "number")) {
    value_df$value <- signal$value_fixed
  } else if (signal$value_type == "parametric") {
    value_df$value <- trial_data[[signal$parametric_modulator]][trial_set]
  } else {
    stop("Failing to populate value column")
  }

  if (!is.null(event_data)) {
    if (!"occurrence_id" %in% names(event_data)) {
      stop("Event data are missing the internal occurrence_id needed for occurrence-level signal alignment.")
    }
    value_df <- value_df %>% dplyr::filter(.data$occurrence_id %in% !!event_data$occurrence_id)
  } else {
    value_df <- value_df %>% dplyr::select(-"occurrence_id")
  }

  return(value_df)
}
