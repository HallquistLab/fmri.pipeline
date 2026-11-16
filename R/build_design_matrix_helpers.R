#' Helper functions for build_design_matrix
#'
#' These functions are internal helpers that implement specific stages of the
#' design matrix construction pipeline. They are split out from the main function
#' to improve testability, maintainability, and clarity.
#'
#' @name build_design_matrix_helpers
#' @keywords internal
NULL

# ==============================================================================
# STAGE 1 HELPERS: Input validation and preprocessing
# ==============================================================================

#' Validate and normalize write_timing_files argument
#'
#' Validates that write_timing_files contains only allowed values and
#' normalizes to lowercase.
#'
#' @param write_timing_files character vector of timing file formats to write,
#'   or NULL if no timing files should be written. Allowed values are:
#'   "convolved", "fsl", "afni", "spm"
#' @return normalized write_timing_files (lowercase), or NULL
#' @keywords internal
#' @noRd
validate_write_timing_files <- function(write_timing_files) {
  checkmate::assert_character(write_timing_files, null.ok = TRUE)
  if (!is.null(write_timing_files)) {
    write_timing_files <- tolower(write_timing_files)
  }
  checkmate::assert_subset(write_timing_files, c("convolved", "fsl", "afni", "spm"))
  return(write_timing_files)
}

#' Validate and prepare run_data data.frame
#'
#' Validates the run_data data.frame structure and adds required columns
#' with defaults if missing.
#'
#' @param run_data data.frame containing run information. Must have columns
#'   for run_nifti (optional), run_volumes (optional), drop_volumes (optional),
#'   and run_number (optional, will be added if missing).
#' @param drop_volumes integer vector of volumes to drop from each run.
#'   If length 1, applied to all runs. Only used if run_data$drop_volumes is NULL.
#' @param lg logger object (from lgr package), or NULL to use default logger
#' @return validated and augmented run_data data.frame with guaranteed
#'   run_number and drop_volumes columns
#' @keywords internal
#' @noRd
validate_run_data <- function(run_data, drop_volumes = 0L, lg = NULL) {
  if (is.null(lg)) lg <- lgr::get_logger()
  
  checkmate::assert_data_frame(run_data)
  
  # verify NIfTIs exist, if specified

  if (!is.null(run_data$run_nifti)) {
    checkmate::assert_file_exists(run_data$run_nifti)
  }
  
  checkmate::assert_integerish(run_data$run_volumes, lower = 1, null.ok = TRUE)
  checkmate::assert_integerish(run_data$drop_volumes, lower = 0, null.ok = TRUE)
  
  # add run_number if missing
  if (is.null(run_data$run_number)) {
    message("No run_number column found in run_data. Assuming that runs numbers are sequential ascending and adding this column.")
    run_data$run_number <- seq_len(nrow(run_data))
  }
  
  checkmate::assert_integerish(run_data$run_number, lower = 1)
  checkmate::assert_integerish(drop_volumes, lower = 0)
  
  # expand drop_volumes to all runs if needed
  if (is.null(run_data$drop_volumes)) {
    if (length(drop_volumes) == nrow(run_data)) {
      run_data$drop_volumes <- drop_volumes
    } else if (length(drop_volumes) == 1L && is.numeric(drop_volumes) && drop_volumes[1L] > 0) {
      lg$info(glue::glue("Using first element of drop_volumes for all runs: {drop_volumes[1L]}"))
      run_data$drop_volumes <- rep(drop_volumes[1L], nrow(run_data))
    } else {
      run_data$drop_volumes <- 0L
    }
  }
  
  return(run_data)
}

#' Validate events data.frame structure and values
#'
#' Validates that the events data.frame has required columns (event, trial,
#' onset, duration) and that onset/duration values are valid (non-NA, non-negative).
#'
#' @param events data.frame containing event information. Required columns:
#'   event (character), trial (integer), onset (numeric), duration (numeric).
#'   Optional: run_number (integer, defaults to 1 if missing).
#' @param lg logger object (from lgr package), or NULL to use default logger
#' @return validated events data.frame with guaranteed run_number column
#' @keywords internal
#' @noRd
validate_events <- function(events, lg = NULL) {
  if (is.null(lg)) lg <- lgr::get_logger()
  
  if (is.null(events)) {
    stop("You must pass in an events data.frame. See ?build_design_matrix for details.")
  }
  
  stopifnot(inherits(events, "data.frame"))
  
  # check required columns
  if (!"event" %in% names(events)) {
    stop("events data.frame must contain event column with the name of the event")
  }
  if (!"trial" %in% names(events)) {
    stop("events data.frame must contain trial column with the trial number for each event")
  }
  if (!"onset" %in% names(events)) {
    stop("events data.frame must contain onset column with the onset time in seconds")
  }
  if (!"duration" %in% names(events)) {
    stop("events data.frame must contain duration column with the event duration in seconds")
  }
  
  # add run_number if missing
  if (!"run_number" %in% names(events)) {
    lg$info("No run_number column found in events. Assuming run_number=1 and adding this column")
    events$run_number <- 1
  }
  
  # validate durations
  if (any(is.na(events$duration))) {
    msg <- "Invalid missing (NA) durations included in events data.frame. Cannot continue."
    lg$error(msg)
    lg$error("%s", capture.output(print(events[is.na(events[["duration"]]), , drop = FALSE])))
    stop(msg)
  } else if (any(events$duration < 0)) {
    msg <- "Invalid negative durations included in events data.frame. Cannot continue."
    lg$error(msg)
    lg$error("%s", capture.output(print(events[events[["duration"]] < 0, , drop = FALSE])))
    stop(msg)
  }
  
  # validate onsets
  if (any(is.na(events$onset))) {
    msg <- "Invalid missing (NA) onsets included in events data.frame. Cannot continue."
    lg$error(msg)
    lg$error("%s", capture.output(print(events[is.na(events[["onset"]]), , drop = FALSE])))
    stop(msg)
  } else if (any(events$onset < 0)) {
    msg <- "Invalid negative onsets included in events data.frame"
    lg$error(msg)
    lg$error("%s", capture.output(print(events[events[["onset"]] < 0, , drop = FALSE])))
    stop(msg)
  }
  
  return(events)
}

#' Validate signals list
#'
#' Validates that signals is a non-empty list.
#'
#' @param signals list of signal specifications
#' @return signals (unchanged if valid)
#' @keywords internal
validate_signals <- function(signals) {
  if (is.null(signals)) {
    stop("You must pass in a signals list. See ?build_design_matrix for details.")
  }
  return(signals)
}

#' Validate TR (repetition time)
#'
#' Validates that TR is a reasonable scalar value.
#'
#' @param tr repetition time in seconds
#' @return tr (unchanged if valid)
#' @keywords internal
#' @noRd
validate_tr <- function(tr) {
  if (is.null(tr)) {
    stop("You must pass in the tr (repetition time) in seconds. See ?build_design_matrix for details.")
  }
  checkmate::assert_number(tr, lower = 0.01, upper = 100)
  return(tr)
}

#' Normalize data.table time-series multipliers to ordinary data.frames
#'
#' @param ts_multipliers NULL, a data.frame, a list of data.frames, or a
#'   character vector of files.
#' @return The input with data.table objects copied to ordinary data.frames.
#' @keywords internal
#' @noRd
normalize_ts_multiplier_tables <- function(ts_multipliers) {
  if (is.data.frame(ts_multipliers)) {
    if (data.table::is.data.table(ts_multipliers)) {
      return(as.data.frame(ts_multipliers))
    }
    return(ts_multipliers)
  }

  if (is.list(ts_multipliers)) {
    return(lapply(ts_multipliers, function(df) {
      if (data.table::is.data.table(df)) as.data.frame(df) else df
    }))
  }

  ts_multipliers
}

#' Validate ts_multipliers structure and numeric content
#'
#' @param ts_multipliers NULL, data.frame, list of data.frames, or character vector of files
#' @param run_data validated run_data data.frame
#' @return validated ts_multipliers, with incoming data.tables copied to data.frames
#' @keywords internal
#' @noRd
validate_ts_multipliers <- function(ts_multipliers, run_data) {
  if (is.null(ts_multipliers)) return(NULL)

  ts_multipliers <- normalize_ts_multiplier_tables(ts_multipliers)

  checkmate::assert_data_frame(run_data)
  n_runs <- nrow(run_data)
  exempt_cols <- c("id", "session", "run_number", "volume")

  validate_numeric_cols <- function(df, label) {
    checkmate::assert_data_frame(df)
    candidate_cols <- setdiff(names(df), exempt_cols)

    if (length(candidate_cols) == 0L) {
      stop("ts_multipliers ", label, " must include at least one multiplier column.", call. = FALSE)
    }

    non_numeric <- candidate_cols[!vapply(df[, candidate_cols, drop = FALSE], is.numeric, logical(1L))]
    if (length(non_numeric) > 0L) {
      stop(
        "ts_multipliers ", label, " has non-numeric multiplier columns: ",
        paste(non_numeric, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (is.data.frame(ts_multipliers)) {
    if (!"run_number" %in% names(ts_multipliers)) {
      stop("If ts_multipliers is a data.frame, it must contain run_number", call. = FALSE)
    }
    validate_numeric_cols(ts_multipliers, "(data.frame)")
    return(ts_multipliers)
  }

  if (is.list(ts_multipliers)) {
    if (length(ts_multipliers) != n_runs) {
      stop(
        "Number of elements in ts_multipliers: ", length(ts_multipliers),
        " does not match number of runs: ", n_runs,
        call. = FALSE
      )
    }
    for (i in seq_along(ts_multipliers)) {
      validate_numeric_cols(ts_multipliers[[i]], paste0("(run ", i, ")"))
    }
    return(ts_multipliers)
  }

  if (is.character(ts_multipliers)) {
    if (length(ts_multipliers) != n_runs) {
      stop(
        "Number of ts_multipliers files: ", length(ts_multipliers),
        " does not match number of runs: ", n_runs,
        call. = FALSE
      )
    }
    checkmate::assert_file_exists(ts_multipliers)
    for (i in seq_along(ts_multipliers)) {
      validate_numeric_cols(utils::read.table(ts_multipliers[i]), paste0("(file run ", i, ")"))
    }
    return(ts_multipliers)
  }

  stop("Don't know how to handle ts_multipliers input", call. = FALSE)
}

#' Get internal processing flags for drop_volumes handling
#'
#' Returns a list of internal flags that control how drop_volumes is applied
#' to different data types. These are not currently exposed to users.
#'
#' @return list with the following logical flags:
#'   \itemize{
#'     \item \code{shift_nifti}: whether to truncate NIfTI files (FALSE, not supported)
#'     \item \code{shift_timing}: whether to shift event timing by drop_volumes*tr (TRUE)
#'     \item \code{shorten_additional}: whether to drop volumes from additional regressors (TRUE)
#'     \item \code{shorten_ts}: whether to drop volumes from ts_multipliers (TRUE)
#'   }
#' @keywords internal
get_drop_volumes_flags <- function() {
  # Internal flags for tracking how drop_volumes is applied. Not exposed to user for now
  # At present, we only support the following:
  #   - NIfTIs can already have drop_volumes applied (externally), or they can be truncated/dropped later (by other programs)
  #   - timing is always shifted according to drop_volumes
  #   - additional regressors are always assumed to be the original (untruncated) length, and volumes are dropped if requested
  #   - ts_multipliers (PPI-style) are always assumed to be the original (untruncated) length, and volumes are dropped if requested
  
  list(
    shift_nifti = FALSE,         # if TRUE, we'd need fslroi x <drop_volumes> approach (not supported)
    shift_timing = TRUE,         # shift drop_volumes*tr from event onset times
    shorten_additional = TRUE,   # apply drop_volumes to additional regressors
    shorten_ts = TRUE            # apply drop_volumes to ts regressors
  )
}

# ==============================================================================
# STAGE 2 HELPERS: Signal expansion and alignment with events
# ==============================================================================

#' Prepare signals for expansion
#'
#' Ensures that beta_series is populated with TRUE or FALSE (default FALSE)
#' for each signal.
#'
#' @param signals list of signal specifications
#' @return signals list with beta_series explicitly set
#' @keywords internal
prepare_signals_for_expansion <- function(signals) {
  if (length(signals) == 0L) {
    stop("signals list is empty. At least one signal must be specified.")
  }
  
  lapply(signals, function(s) {
    s$beta_series <- ifelse(isTRUE(s$beta_series), TRUE, FALSE)
    s
  })
}

#' Expand and flatten signals list
#'
#' Applies expand_signal to each signal and flattens the result into a
#' single list of expanded signals.
#'
#' @param signals list of signal specifications
#' @return named list of expanded signals
#' @keywords internal
#' @noRd
expand_signals_list <- function(signals) {
  signals_expanded <- purrr::list_flatten(unname(lapply(signals, expand_signal)))
  
  # Extract names and check for duplicates
  signal_names <- sapply(signals_expanded, "[[", "name")
  
  if (any(duplicated(signal_names))) {
    dup_names <- unique(signal_names[duplicated(signal_names)])
    stop("Duplicate signal names after expansion: ", paste(dup_names, collapse = ", "),
         ". Each signal must have a unique name.")
  }
  
  names(signals_expanded) <- signal_names
  return(signals_expanded)
}

#' Collapse exact duplicate neural-predictor occurrences
#'
#' Removes repeated occurrence records that would produce the same neural
#' predictor. Set \code{keep_duplicate_occurrences = TRUE} on a signal to
#' preserve repeated records.
#'
#' @param s_aligned data.frame of timing records after signal alignment
#' @param s signal specification list
#' @param signal_name name used in diagnostic logging
#' @param lg logger object for error/info messages
#' @return a data.frame with exact duplicate records removed
#' @keywords internal
#' @noRd
collapse_duplicate_occurrences <- function(s_aligned, s, signal_name, lg = NULL) {
  if (is.null(lg)) lg <- lgr::get_logger()

  keep_duplicates <- if (is.null(s$keep_duplicate_occurrences)) {
    FALSE
  } else {
    checkmate::assert_flag(s$keep_duplicate_occurrences)
    s$keep_duplicate_occurrences
  }
  if (isTRUE(keep_duplicates) || nrow(s_aligned) < 2L) return(s_aligned)

  # A signal's neural predictor is fully determined by its run, optional trial
  # parent, timing, and pre-convolution amplitude.  Collapse exact duplicates
  # here, after signal-specific duration has been applied and before any
  # centering or convolution.  This protects wide trial data that have been
  # repeated to carry additional within-trial timing records.
  key_cols <- intersect(
    c("id", "session", "run_number", "trial", "onset", "duration", "value"),
    names(s_aligned)
  )
  duplicate_rows <- duplicated(s_aligned[, key_cols, drop = FALSE])
  n_duplicates <- sum(duplicate_rows)
  if (n_duplicates == 0L) return(s_aligned)

  s_aligned <- s_aligned[!duplicate_rows, , drop = FALSE]
  lg$info(
    "Signal '%s': collapsed %d exact duplicate occurrence record(s) (%d input records -> %d modeled occurrences). Set keep_duplicate_occurrences: true to retain them.",
    signal_name, n_duplicates, n_duplicates + nrow(s_aligned), nrow(s_aligned)
  )
  s_aligned
}

#' Validate beta-series occurrence cardinality
#'
#' Beta-series expansion creates one regressor for each trial and therefore
#' requires exactly one distinct timing occurrence per trial after duplicate
#' timing records have been collapsed.
#'
#' @param s_aligned data.frame of timing records after duplicate collapse
#' @param s expanded signal specification list
#' @param signal_name name used in error messages
#' @return invisibly TRUE when the signal is not beta-series-expanded or valid
#' @keywords internal
#' @noRd
validate_beta_series_occurrences <- function(s_aligned, s, signal_name) {
  if (is.null(s$beta_series_trial)) return(invisible(TRUE))

  if (nrow(s_aligned) == 0L) {
    stop(
      "Beta-series signal '", signal_name,
      "' requires exactly one distinct timing occurrence per id/session/run/trial after duplicate collapse, but none remained."
    )
  }

  required_cols <- c("run_number", "trial")
  missing_cols <- setdiff(required_cols, names(s_aligned))
  if (length(missing_cols) > 0L || anyNA(s_aligned$trial)) {
    stop(
      "Beta-series signal '", signal_name,
      "' requires a non-missing trial identifier for every modeled occurrence."
    )
  }

  group_cols <- intersect(c("id", "session", "run_number", "trial"), names(s_aligned))
  group_values <- lapply(s_aligned[, group_cols, drop = FALSE], as.character)
  group_key <- do.call(paste, c(group_values, sep = "\r"))
  occurrence_counts <- table(group_key)
  if (any(occurrence_counts != 1L)) {
    bad_keys <- names(occurrence_counts)[occurrence_counts != 1L]
    stop(
      "Beta-series signal '", signal_name,
      "' requires exactly one distinct timing occurrence per id/session/run/trial after duplicate collapse. ",
      "Multiple occurrences were found for: ", paste(bad_keys, collapse = "; "), "."
    )
  }

  invisible(TRUE)
}

#' Align a single signal with events
#'
#' Merges a trial-indexed signal with the time-indexed events data.frame,
#' putting the signal values onto the event timing grid.
#'
#' @param s a signal specification list with event, value, and optionally duration
#' @param events data.frame of events with event, run_number, trial, onset, duration
#' @param lg logger object for error/info messages
#' @return list of data.tables, one per run, with trial, onset, duration, value columns
#' @keywords internal
#' @noRd
align_signal_with_events <- function(s, events, lg = NULL, center_values = TRUE) {
  if (is.null(lg)) lg <- lgr::get_logger()
  
  # Validate signal has required event element
  if (is.null(s$event)) {
    msg <- "Signal does not have event element"
    lg$error(msg)
    lg$error("%s", utils::capture.output(print(s)))
    stop(msg)
  }
  
  # Get signal name for error messages
  signal_name <- if (!is.null(s$name)) s$name else "<unnamed>"
  
  # Default value to 1 if not specified (unit-height regressor)
  if (is.null(s$value)) {
    lg$info("Signal '%s' is missing a 'value' element. Adding 1 for value, assuming a unit-height regressor.", signal_name)
    s$value <- 1
  }
  
  # Determine join columns. Pipeline-generated timing records carry an
  # occurrence_id that identifies their source timing row; this prevents
  # many-to-many joins when a trial has multiple timing occurrences.
  join_cols <- c("run_number", "trial")
  df_events <- dplyr::filter(events, .data$event == s$event)

  # Internal callers normalize aggregated events at their creation boundary,
  # but callers can still supply a data.table directly. Copy only in that
  # defensive case: setDF() would unexpectedly change the caller's object by
  # reference, while ordinary data.frames need no conversion at all.
  if (data.table::is.data.table(df_events)) {
    df_events <- as.data.frame(df_events)
  }
  
  # Check that the event type exists in events
  if (nrow(df_events) == 0L) {
    available_events <- unique(events$event)
    stop("Signal '", signal_name, "' references event '", s$event, 
         "' which has no occurrences in the events data.frame. ",
         "Available event types: ", paste(available_events, collapse = ", "))
  }
  
  event_runs <- factor(sort(unique(df_events$run_number)))
  df_signal <- s$value
  if (data.table::is.data.table(df_signal)) {
    df_signal <- as.data.frame(df_signal)
  }

  use_occurrence_id <- is.data.frame(df_signal) &&
    "occurrence_id" %in% names(df_signal) &&
    "occurrence_id" %in% names(df_events)
  if (use_occurrence_id) join_cols <- c("run_number", "occurrence_id")
  
  # Add additional join columns if present in events
  if ("id" %in% names(df_events)) join_cols <- c(join_cols, "id")
  if ("session" %in% names(df_events)) join_cols <- c(join_cols, "session")
  if ("block_number" %in% names(df_events)) join_cols <- c(join_cols, "block_number")
  
  # Merge signal with events
  if (length(df_signal) == 1L && is.numeric(df_signal)) {
    # Task indicator-type regressor: replicate value for all events
    s_aligned <- df_events
    s_aligned$value <- df_signal
  } else if (is.data.frame(df_signal)) {
    # Validate that signal data.frame has required join columns
    missing_join_cols <- setdiff(join_cols, names(df_signal))
    if (length(missing_join_cols) > 0L) {
      stop("Signal '", signal_name, "' value data.frame is missing required columns: ",
           paste(missing_join_cols, collapse = ", "))
    }
    
    # Validate that signal data.frame has a value column
    if (!"value" %in% names(df_signal)) {
      stop("Signal '", signal_name, "' value data.frame must have a 'value' column")
    }

    if (use_occurrence_id && anyDuplicated(df_signal[, join_cols, drop = FALSE])) {
      stop(
        "Signal '", signal_name,
        "' has duplicate occurrence-level value records. Each id/session/run/occurrence_id combination must be unique."
      )
    }
    if (use_occurrence_id && anyDuplicated(df_events[, join_cols, drop = FALSE])) {
      stop(
        "Signal '", signal_name,
        "' has duplicate occurrence identifiers in its timing records. Each id/session/run/occurrence_id combination must be unique."
      )
    }
    
    # Parametric/factor signals from pipeline timing rows join by occurrence.
    # The trial-key fallback preserves the public legacy build_design_matrix API.
    signal_for_join <- df_signal
    if (use_occurrence_id && "trial" %in% names(signal_for_join)) {
      # The event timing table supplies the parent trial metadata. Keeping the
      # duplicate trial column from the value table would create trial.x/y
      # suffixes even though alignment is intentionally occurrence-based.
      signal_for_join <- signal_for_join %>% dplyr::select(-"trial")
    }
    s_aligned <- signal_for_join %>%
      dplyr::left_join(df_events, by = join_cols) %>%
      dplyr::arrange(.data$run_number, .data$trial)
    
    # Check for join failures (NAs in onset after join)
    na_onsets <- sum(is.na(s_aligned$onset))
    if (na_onsets > 0L) {
      unit_label <- if (use_occurrence_id) "occurrence(s)" else "trial(s)"
      warn_msg <- sprintf("Signal '%s': %d %s in signal did not match events (NA onsets after join). These will be excluded.",
                          signal_name, na_onsets, unit_label)
      lg$warn(warn_msg)
      warning(warn_msg, call. = FALSE)
      s_aligned <- s_aligned[!is.na(s_aligned$onset), , drop = FALSE]
    }
  } else {
    stop("Signal '", signal_name, "': Unknown data type for value. ",
         "Expected numeric scalar or data.frame, got: ", class(df_signal)[1])
  }
  
  # Validate required columns exist before duration application
  validate_aligned_columns(s_aligned, signal_name)
  
  # Handle duration specification
  s_aligned <- apply_signal_duration(s_aligned, s$duration, signal_name)

  na_values <- sum(is.na(s_aligned$value))
  if (na_values > 0L) {
    lg$info(
      "Signal '%s': %d occurrence(s) have missing values and will be omitted before convolution.",
      signal_name, na_values
    )
  }

  # By default, copied trial rows represent one timing occurrence, not repeated
  # neural events. Exact duplicate records are removed before convolution.
  s_aligned <- collapse_duplicate_occurrences(s_aligned, s, signal_name, lg)
  validate_beta_series_occurrences(s_aligned, s, signal_name)

  # For formula-coded parametric signals (e.g., additive ~ A + B or contr.sum), center raw
  # modulator values per run BEFORE multiplying by signed/overlapping contrast weights.
  precentered_param <- FALSE
  if (isTRUE(s$formula_weighted) && isTRUE(s$is_parametric) && "wi_weight" %in% names(s_aligned)) {
    if (isTRUE(center_values) && nrow(s_aligned) > 0L) {
      for (rn in unique(s_aligned$run_number)) {
        r_idx <- which(s_aligned$run_number == rn & !is.na(s_aligned$value))
        if (length(r_idx) <= 1L || stats::sd(s_aligned$value[r_idx]) <= 1e-5) {
          s_aligned$value[r_idx] <- 0
        } else {
          s_aligned$value[r_idx] <- s_aligned$value[r_idx] - mean(s_aligned$value[r_idx])
        }
      }
      precentered_param <- TRUE
    }
    s_aligned$value <- s_aligned$value * s_aligned$wi_weight
  }

  # Split by run and format for dmat
  retsplit <- split_signal_by_run(
    s_aligned,
    event_runs,
    s$event,
    s$physio_only,
    return_deriv = isTRUE(s$return_deriv),
    is_parametric = isTRUE(s$is_parametric),
    precentered_parametric = precentered_param,
    demean_convolved = isTRUE(s$demean_convolved),
    formula_weighted = isTRUE(s$formula_weighted)
  )

  return(retsplit)
}

#' Validate that aligned signal has required columns
#'
#' @param s_aligned data.frame of aligned signal/events
#' @param signal_name character name of signal for error messages
#' @keywords internal
#' @noRd
validate_aligned_columns <- function(s_aligned, signal_name = "<unnamed>") {
  required_cols <- c("run_number", "trial", "onset", "duration", "value")
  missing_cols <- setdiff(required_cols, names(s_aligned))
  
  if (length(missing_cols) > 0L) {
    stop("Signal '", signal_name, "': aligned data is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  invisible(TRUE)
}

#' Apply duration specification to aligned signal
#'
#' Sets the duration column based on the signal's duration specification.
#' Can be a scalar (applied to all rows) or a column name from events.
#'
#' @param s_aligned data.frame of aligned signal/events
#' @param duration_spec the duration specification (NULL, numeric, or column name)
#' @param signal_name character name of signal for error messages
#' @return s_aligned with duration column set appropriately
#' @keywords internal
#' @noRd
apply_signal_duration <- function(s_aligned, duration_spec, signal_name = "<unnamed>") {
  if (length(duration_spec) > 1L) {
    stop("Signal '", signal_name, "': Don't know how to interpret multi-element duration argument: ", 
         paste0(duration_spec, collapse = ", "))
  }
  
  if (!is.null(duration_spec)) {
    if (is.numeric(duration_spec)) {
      s_aligned$duration <- duration_spec
    } else if (is.character(duration_spec)) {
      # Validate that the column exists
      if (!duration_spec %in% names(s_aligned)) {
        stop("Signal '", signal_name, "': duration column '", duration_spec, 
             "' not found in aligned data. Available columns: ",
             paste(names(s_aligned), collapse = ", "))
      }
      s_aligned$duration <- s_aligned[[duration_spec]]
    } else {
      stop("Signal '", signal_name, "': duration must be numeric or a column name string, got: ",
           class(duration_spec)[1])
    }
  }
  
  # Validate durations are valid
  if (any(is.na(s_aligned$duration))) {
    stop("Signal '", signal_name, "': duration contains NA values after applying duration specification")
  }
  
  if (any(s_aligned$duration < 0)) {
    stop("Signal '", signal_name, "': duration contains negative values")
  }
  
  return(s_aligned)
}

#' Split aligned signal by run
#'
#' Converts an aligned signal data.frame into a list of data.tables,
#' one per run, with appropriate attributes.
#'
#' @param s_aligned data.frame with run_number, trial, onset, duration, value
#' @param event_runs factor of run levels
#' @param event_name character name of the event type
#' @param physio_only logical, whether this is a physio-only regressor
#' @return named list of data.tables by run
#' @keywords internal
#' @noRd
split_signal_by_run <- function(s_aligned, event_runs, event_name, physio_only = FALSE,
                                return_deriv = FALSE, is_parametric = FALSE,
                                precentered_parametric = FALSE, demean_convolved = FALSE,
                                formula_weighted = FALSE) {
  retdf <- s_aligned %>%
    dplyr::select("run_number", "trial", "onset", "duration", "value") %>%
    dplyr::mutate(run_number = factor(.data$run_number, levels = event_runs)) %>%
    data.table::setDT()
  
  # Split by run, keeping all factor levels
  retsplit <- split(retdf, by = "run_number", keep.by = FALSE, sorted = TRUE)
  names(retsplit) <- paste0("run_number", names(retsplit))
  
  # Tag with event, physio, derivative, and centering attributes
  retsplit <- lapply(retsplit, function(rr) {
    attr(rr, "event") <- event_name
    if (isTRUE(physio_only)) attr(rr, "physio_only") <- TRUE
    if (isTRUE(return_deriv)) attr(rr, "return_deriv") <- TRUE
    if (isTRUE(is_parametric)) attr(rr, "is_parametric") <- TRUE
    if (isTRUE(precentered_parametric)) attr(rr, "precentered_parametric") <- TRUE
    if (isTRUE(formula_weighted)) attr(rr, "formula_weighted") <- TRUE
    if (isTRUE(demean_convolved)) attr(rr, "demean_convolved") <- TRUE
    return(rr)
  })
  
  return(retsplit)
}

#' Align all signals with events
#'
#' Applies align_signal_with_events to each expanded signal.
#'
#' @param signals_expanded named list of expanded signals
#' @param events data.frame of events
#' @param lg logger object
#' @return named list of aligned signals (each is a list of runs)
#' @keywords internal
#' @noRd
align_all_signals <- function(signals_expanded, events, lg = NULL, center_values = TRUE) {
  lapply(signals_expanded, function(s) {
    align_signal_with_events(s, events, lg, center_values = center_values)
  })
}

#' Extract signal configuration into bdm_args vectors
#'
#' Extracts normalization, beta_series, rm_zeros, convmax_1, and add_deriv
#' settings from each signal into vectors for use during convolution.
#'
#' @param signals_expanded list of expanded signals
#' @return list of configuration vectors:
#'   \itemize{
#'     \item \code{normalizations}: normalization method for each signal
#'     \item \code{beta_series}: whether each signal uses beta series
#'     \item \code{rm_zeros}: whether to remove zeros before convolution
#'     \item \code{convmax_1}: whether to normalize convolved max to 1
#'     \item \code{add_derivs}: whether to add temporal derivatives
#'   }
#' @keywords internal
#' @noRd
extract_signal_config <- function(signals_expanded) {
  config <- list()
  
  # Extract normalization for each regressor
  config$normalizations <- sapply(signals_expanded, function(s) {
    ifelse(is.null(s$normalization), "none", s$normalization)
  })
  
  # Extract beta_series settings
  config$beta_series <- sapply(signals_expanded, function(s) {
    ifelse(isTRUE(s$beta_series), TRUE, FALSE)
  })
  
  # Extract rm_zeros settings (default TRUE)
  config$rm_zeros <- sapply(signals_expanded, function(s) {
    if (is.null(s$rm_zeros)) {
      TRUE
    } else {
      if (!s$rm_zeros %in% c(TRUE, FALSE)) {
        stop("Don't know how to interpret rm_zeros setting of: ", s$rm_zeros)
      } else {
        as.logical(s$rm_zeros)
      }
    }
  })
  
  # Extract convmax_1 settings (default FALSE)
  config$convmax_1 <- sapply(signals_expanded, function(s) {
    if (is.null(s$convmax_1)) {
      FALSE
    } else {
      if (!s$convmax_1 %in% c(TRUE, FALSE)) {
        stop("Don't know how to interpret convmax_1 setting of: ", s$convmax_1)
      } else {
        as.logical(s$convmax_1)
      }
    }
  })
  
  # Extract add_deriv settings (default FALSE)
  config$add_derivs <- sapply(signals_expanded, function(s) {
    if (is.null(s$add_deriv)) {
      FALSE
    } else {
      if (!s$add_deriv %in% c(TRUE, FALSE)) {
        stop("Don't know how to interpret add_deriv setting of: ", s$add_deriv)
      } else {
        as.logical(s$add_deriv)
      }
    }
  })
  
  return(config)
}

#' Perform Stage 2: Signal expansion and alignment
#'
#' Main dispatcher function for Stage 2 that expands signals and aligns
#' them with events timing.
#'
#' @param signals list of signal specifications
#' @param events data.frame of events
#' @param lg logger object
#' @return list with:
#'   \itemize{
#'     \item \code{signals_expanded}: named list of expanded signals
#'     \item \code{signals_aligned}: named list of aligned signals (runs x values)
#'     \item \code{signal_config}: list of configuration vectors for bdm_args
#'   }
#' @keywords internal
#' @noRd
expand_and_align_signals <- function(signals, events, lg = NULL, center_values = TRUE) {
  if (is.null(lg)) lg <- lgr::get_logger()
  
  # Prepare signals (ensure beta_series is set)
  signals <- prepare_signals_for_expansion(signals)
  
  # Expand signals (handle wi_factors, etc.)
  signals_expanded <- expand_signals_list(signals)
  
  # Align each signal with events
  signals_aligned <- align_all_signals(signals_expanded, events, lg, center_values = center_values)
  
  # Extract configuration for convolution
  signal_config <- extract_signal_config(signals_expanded)
  
  return(list(
    signals_expanded = signals_expanded,
    signals_aligned = signals_aligned,
    signal_config = signal_config
  ))
}

# ==============================================================================
# STAGE 6 HELPERS: Write timing files to disk
# ==============================================================================

#' Write timing files to disk for different analysis software
#'
#' Main dispatcher function that writes timing files in various formats
#' (convolved, FSL, AFNI, SPM) based on the requested formats.
#'
#' @param write_timing_files character vector of formats to write. Valid values:
#'   "convolved", "fsl", "afni", "spm"
#' @param dmat 2D array of design matrices (runs x regressors)
#' @param dmat_convolved list of convolved design matrices per run
#' @param output_directory path to output directory for timing files
#' @param runs_to_output numeric vector of run numbers to output
#' @param center_values logical, whether to center regressor values for FSL
#' @return list with file paths for each format written:
#'   \itemize{
#'     \item \code{tf_convolved}: matrix of convolved file paths (runs x regressors)
#'     \item \code{tf_convolved_concat}: named vector of concatenated file paths
#'     \item \code{tf_fsl}: matrix of FSL 3-column file paths
#'     \item \code{tf_afni}: AFNI dmBLOCK file paths (currently not tracked)
#'   }
#' @keywords internal
#' @noRd
write_timing_files_to_disk <- function(write_timing_files, dmat, dmat_convolved,
                                       output_directory, runs_to_output, 
                                       center_values = TRUE) {
  # Initialize return values
  result <- list(
    tf_convolved = NULL,
    tf_convolved_concat = NULL,
    tf_fsl = NULL,
    tf_afni = NULL
  )
  
  if (is.null(write_timing_files)) {
    return(result)
  }
  
  # Create output directory
 dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Write convolved regressors
  if ("convolved" %in% write_timing_files) {
    conv_result <- write_convolved_timing_files(dmat, dmat_convolved, output_directory)
    result$tf_convolved <- conv_result$tf_convolved
    result$tf_convolved_concat <- conv_result$tf_convolved_concat
  }
  
  # Write FSL 3-column format
  if ("fsl" %in% write_timing_files) {
    result$tf_fsl <- write_fsl_timing_files(dmat, output_directory, runs_to_output, center_values)
  }
  
  # Write AFNI dmBLOCK format
  if ("afni" %in% write_timing_files) {
    result$tf_afni <- write_afni_timing_files(dmat, output_directory)
  }
  
  # SPM format not yet implemented
  if ("spm" %in% write_timing_files) {
    warning("SPM timing file format not yet implemented")
  }
  
  return(result)
}

#' Write convolved regressors to disk
#'
#' Writes convolved regressors as 1D files (one per run/regressor) and also
#' creates concatenated versions across runs (for AFNI-style analysis).
#'
#' @param dmat 2D array of design matrices (runs x regressors)
#' @param dmat_convolved list of convolved design matrices per run
#' @param output_directory path to output directory
#' @return list with:
#'   \itemize{
#'     \item \code{tf_convolved}: matrix of file paths (runs x regressors)
#'     \item \code{tf_convolved_concat}: named vector of concatenated file paths
#'   }
#' @keywords internal
#' @noRd
write_convolved_timing_files <- function(dmat, dmat_convolved, output_directory) {
  tf_convolved <- matrix(
    NA_character_,
    nrow = dim(dmat)[1L],
    ncol = dim(dmat)[2L],
    dimnames = dimnames(dmat)
  )
  
  conv_concat <- list()
  
  # Write per-run convolved files and accumulate for concatenation
  for (r in seq_along(dmat_convolved)) {
    for (v in seq_along(dmat_convolved[[r]])) {
      reg_name <- names(dmat_convolved[[r]])[v]
      fname <- paste0(names(dmat_convolved)[r], "_", reg_name, ".1D")
      to_write <- round(dmat_convolved[[r]][[v]], 6)
      
      # Accumulate for concatenated file
      conv_concat[[reg_name]] <- c(conv_concat[[reg_name]], to_write)
      
      # Write individual run file
      ofile <- file.path(output_directory, fs::path_sanitize(fname, replacement = ""))
      tf_convolved[r, v] <- ofile
      write.table(
        to_write,
        file = ofile,
        sep = "\n", eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE
      )
    }
  }
  
  # Write run-concatenated convolved regressors (for AFNI)
  tf_convolved_concat <- write_convolved_concat_files(conv_concat, output_directory)
  
  return(list(
    tf_convolved = tf_convolved,
    tf_convolved_concat = tf_convolved_concat
  ))
}

#' Write concatenated convolved regressor files
#'
#' Creates a single 1D file per regressor containing all runs concatenated.
#'
#' @param conv_concat named list of concatenated regressor vectors
#' @param output_directory path to output directory
#' @return named character vector of file paths
#' @keywords internal
#' @noRd
write_convolved_concat_files <- function(conv_concat, output_directory) {
  if (length(conv_concat) == 0) {
    return(character(0))
  }
  
  tf_convolved_concat <- sapply(seq_along(conv_concat), function(v) {
    fname <- paste0(names(conv_concat)[v], "_concat.1D")
    ofile <- file.path(output_directory, fs::path_sanitize(fname, replacement = ""))
    write.table(
      conv_concat[[v]],
      file = ofile,
      sep = "\n", eol = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE
    )
    return(ofile)
  })
  names(tf_convolved_concat) <- names(conv_concat)
  
  return(tf_convolved_concat)
}

#' Write FSL 3-column timing files
#'
#' Writes timing files in FSL's 3-column format (onset, duration, value).
#' Optionally centers values and removes zero-value events.
#'
#' @param dmat 2D array of design matrices (runs x regressors)
#' @param output_directory path to output directory
#' @param runs_to_output numeric vector of run numbers for naming
#' @param center_values logical, whether to mean-center values
#' @return matrix of file paths (runs x regressors)
#' @keywords internal
write_fsl_timing_files <- function(dmat, output_directory, runs_to_output, 
                                   center_values = TRUE) {
  tf_fsl <- matrix(
    NA_character_,
    nrow = dim(dmat)[1L],
    ncol = dim(dmat)[2L],
    dimnames = dimnames(dmat)
  )
  
  for (i in seq_len(dim(dmat)[1L])) {
    for (reg in seq_len(dim(dmat)[2L])) {
      regout <- dmat[[i, reg]]
      
      # Skip empty regressors
      if (nrow(regout) == 0L) {
        next
      }
      
      # Extract only the columns FSL needs
      regout <- regout[, c("onset", "duration", "value"), drop = FALSE]
      
      # Optionally center values
      regout <- center_fsl_regressor(regout, center_values)
      
      # Write the file
      fname <- paste0("run", runs_to_output[i], "_", dimnames(dmat)[[2L]][reg], "_FSL3col.txt")
      ofile <- file.path(output_directory, fs::path_sanitize(fname, replacement = ""))
      tf_fsl[i, reg] <- ofile
      write.table(regout, file = ofile, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
    }
  }
  
  return(tf_fsl)
}

#' Center FSL regressor values
#'
#' Mean-centers the value column of an FSL regressor, removing zero-value
#' events first. Only centers if there is variation in the values.
#'
#' @param regout data.frame with onset, duration, value columns
#' @param center_values logical, whether to apply centering
#' @return data.frame with potentially centered values
#' @keywords internal
#' @noRd
center_fsl_regressor <- function(regout, center_values = TRUE) {
  if (!center_values) {
    return(regout)
  }
  
  # Skip if all values are zero
  if (all(na.omit(regout[, "value"]) == 0.0)) {
    return(regout)
  }
  
  # Remove zero-value events
  regout <- regout[regout[, "value"] != 0, , drop = FALSE]
  
  # Mean center if there's variation (not just a constant indicator function)
  if (nrow(regout) > 1L && stats::sd(regout[, "value"], na.rm = TRUE) > 0) {
    regout[, "value"] <- regout[, "value"] - mean(regout[, "value"], na.rm = TRUE)
  }
  
  return(regout)
}

#' Write AFNI dmBLOCK timing files
#'
#' Writes timing files in AFNI's dmBLOCK format, which combines timing
#' information across regressors that share onset/duration patterns.
#' Format: time*modulation1,modulation2:duration (one line per run)
#'
#' @param dmat 2D array of design matrices (runs x regressors)
#' @param output_directory path to output directory
#' @return NULL (files are written as side effect)
#' @keywords internal
#' @noRd
write_afni_timing_files <- function(dmat, output_directory) {
  # Extract onsets and durations across all runs for each regressor
  regonsets <- apply(dmat, 2, function(reg) {
    unname(do.call(c, lapply(reg, function(run) run[, "onset"])))
  })
  
  regdurations <- apply(dmat, 2, function(reg) {
    unname(do.call(c, lapply(reg, function(run) run[, "duration"])))
  })
  
  # Find regressors that share onset/duration patterns
  # Uses first row to determine groupings
  first_onset_duration <- paste(regonsets[1, ], regdurations[1, ])
  
  # Group regressors by shared onset/duration patterns
  dt <- data.table::as.data.table(first_onset_duration)[
    , list(comb = list(.I)), 
    by = first_onset_duration
  ]
  
  # Write a file for each group of combined regressors
  lapply(dt$comb, function(comb) {
    write_afni_dmblock_file(dmat, comb, output_directory)
  })
  
  invisible(NULL)
}

#' Write a single AFNI dmBLOCK file for a group of regressors
#'
#' Creates a dmBLOCK file for regressors that share the same onset/duration
#' pattern, combining their amplitude modulations.
#'
#' @param dmat 2D array of design matrices
#' @param comb integer vector of column indices for regressors to combine
#' @param output_directory path to output directory
#' @return NULL (file written as side effect)
#' @keywords internal
#' @noRd
write_afni_dmblock_file <- function(dmat, comb, output_directory) {
  combmat <- dmat[, comb, drop = FALSE]
  
  # Build one line per run in AFNI format
  runvec <- character(dim(combmat)[1L])
  
  for (i in seq_len(dim(combmat)[1L])) {
    # Get onsets, durations from first regressor (shared by definition)
    runonsets <- combmat[[i, 1]][, "onset"]
    rundurations <- combmat[[i, 1]][, "duration"]
    
    # Get values from all regressors in the combination
    runvalues <- do.call(cbind, lapply(combmat[i, ], function(reg) reg[, "value"]))
    
    # AFNI doesn't want the boxcar (all 1s) in dmBLOCK format - remove it
    indicator_func <- apply(runvalues, 2, function(col) all(col == 1.0))
    if (any(indicator_func)) {
      runvalues <- runvalues[, -which(indicator_func), drop = FALSE]
    }
    
    # Format the run's timing string
    runvec[i] <- format_afni_run_string(runonsets, rundurations, runvalues)
  }
  
  # Write the file
  fname <- paste0(paste(dimnames(combmat)[[2L]], collapse = "_"), "_dmBLOCK.txt")
  ofile <- file.path(output_directory, fname)
  writeLines(runvec, ofile)
  
  invisible(NULL)
}

#' Format a single run's AFNI timing string
#'
#' Creates the AFNI dmBLOCK format string for a single run.
#' If no parametric modulation, uses TIME:DURATION format.
#' Otherwise uses TIME*PARAM1,PARAM2:DURATION format.
#'
#' @param onsets numeric vector of event onsets
#' @param durations numeric vector of event durations
#' @param values matrix of modulation values (events x regressors), 
#'   or 0-column matrix if no parametric modulation
#' @return character string in AFNI dmBLOCK format
#' @keywords internal
#' @noRd
format_afni_run_string <- function(onsets, durations, values) {
  if (ncol(values) == 0L) {
    # No parametric modulation: TIME:DURATION format
    event_strings <- sapply(seq_along(onsets), function(j) {
      paste0(round(onsets[j], 6), ":", round(durations[j], 6))
    })
  } else {
    # With parametric modulation: TIME*PARAM1,PARAM2:DURATION format
    event_strings <- sapply(seq_along(onsets), function(j) {
      params <- paste(round(values[j, ], 6), collapse = ",")
      paste0(round(onsets[j], 6), "*", params, ":", round(durations[j], 6))
    })
  }
  
  paste(event_strings, collapse = " ")
}

# ==============================================================================
# STAGE 7 HELPERS: Collinearity diagnostics
# ==============================================================================

#' Compute collinearity diagnostics on convolved design matrices
#'
#' Calculates correlation matrices and variance inflation factors (VIF)
#' for each run's convolved design matrix.
#'
#' @param dmat_convolved list of convolved design matrices per run
#' @return list of collinearity diagnostics per run, each containing:
#'   \itemize{
#'     \item \code{r}: correlation matrix
#'     \item \code{vif}: variance inflation factors
#'   }
#' @keywords internal
#' @noRd
compute_collinearity_convolved <- function(dmat_convolved) {
  lapply(dmat_convolved, function(run) {
    corvals <- suppressWarnings(stats::cor(run, use = "pairwise.complete.obs"))
    col_names <- names(run)
    if (length(col_names) < 2L) {
      return(list(r = corvals, vif = stats::setNames(rep(1, length(col_names)), col_names)))
    }

    sds <- vapply(run, function(col) stats::sd(col, na.rm = TRUE), numeric(1))
    non_const <- !is.na(sds) & sds > 1e-12
    if (sum(non_const) < 2L) {
      var_infl <- stats::setNames(rep(NA_real_, length(col_names)), col_names)
      if (sum(non_const) == 1L) var_infl[non_const] <- 1
      return(list(r = corvals, vif = var_infl))
    }

    active_cols <- col_names[non_const]
    vif_mat <- data.frame(dummy = stats::rnorm(nrow(run)), run[, active_cols, drop = FALSE], check.names = FALSE)
    quoted_cols <- paste0("`", gsub("`", "", active_cols, fixed = TRUE), "`")
    vif_form <- stats::as.formula(paste("dummy ~ 1 +", paste(quoted_cols, collapse = " + ")))

    active_vif <- tryCatch(
      {
        vv <- car::vif(stats::lm(vif_form, data = vif_mat))
        names(vv) <- gsub("^`|`$", "", names(vv))
        vv
      },
      error = function(e) stats::setNames(rep(NA_real_, length(active_cols)), active_cols)
    )

    var_infl <- stats::setNames(rep(NA_real_, length(col_names)), col_names)
    var_infl[names(active_vif)] <- as.numeric(active_vif)

    list(r = corvals, vif = var_infl)
  })
}

#' Merge additional regressors into convolved design matrices
#'
#' Adds confound/nuisance regressors to the convolved and unconvolved
#' design matrices, mean-centering them first.
#'
#' @param dmat_convolved list of convolved design matrices per run
#' @param dmat_unconvolved list of unconvolved design matrices per run
#' @param additional_regressors data.frame of additional regressors with run_number column
#' @return list with updated \code{dmat_convolved} and \code{dmat_unconvolved}
#' @keywords internal
#' @noRd
merge_additional_regressors_to_dmat <- function(dmat_convolved, dmat_unconvolved, 
                                                additional_regressors) {
  if (is.null(additional_regressors)) {
    return(list(
      dmat_convolved = dmat_convolved,
      dmat_unconvolved = dmat_unconvolved
    ))
  }
  
  for (i in seq_along(dmat_convolved)) {
    # Extract run number from names (e.g., "run_number1" -> 1)
    runnum <- as.numeric(sub("run_number(\\d+)", "\\1", names(dmat_convolved)[i], perl = TRUE))
    
    additional_regressors_currun <- additional_regressors %>%
      dplyr::filter(.data$run_number == !!runnum) %>%
      dplyr::select(-"run_number")
    
    # Mean-center additional regressors
    additional_regressors_currun <- as.data.frame(lapply(
      additional_regressors_currun, 
      function(x) x - mean(x, na.rm = TRUE)
    ))
    
    dmat_convolved[[i]] <- dplyr::bind_cols(dmat_convolved[[i]], additional_regressors_currun)
    dmat_unconvolved[[i]] <- dplyr::bind_cols(dmat_unconvolved[[i]], additional_regressors_currun)
  }
  
  return(list(
    dmat_convolved = dmat_convolved,
    dmat_unconvolved = dmat_unconvolved
  ))
}

#' Build concatenated design for AFNI-style analysis
#'
#' Creates a concatenated version of the design where run timing is adjusted
#' to create a single continuous time series across all runs.
#'
#' @param dmat 2D list of signal data (runs x signals)
#' @param run_timing numeric vector of cumulative run timing in seconds
#' @return list of concatenated design matrices, one per signal
#' @keywords internal
#' @noRd
build_design_concat <- function(dmat, run_timing) {
  design_concat <- lapply(seq_len(dim(dmat)[2L]), function(reg) {
    this_reg <- dmat[, reg]
    concat_reg <- do.call(rbind, lapply(seq_along(this_reg), function(run) {
      timing <- this_reg[[run]]
      timing[, "onset"] <- timing[, "onset"] + ifelse(run > 1, run_timing[run - 1], 0)
      return(timing)
    }))
    attr(concat_reg, "event") <- attr(this_reg[[1]], "event")
    concat_reg
  })
  names(design_concat) <- dimnames(dmat)[[2L]]
  
  return(design_concat)
}

#' Extract signal configuration into bdm_args
#'
#' Extracts normalization, beta_series, rm_zeros, convmax_1, and add_deriv
#' settings from each signal into vectors for use during convolution.
#'
#' @param signals_expanded list of expanded signals
#' @return list of configuration vectors:
#'   \itemize{
#'     \item \code{normalizations}: normalization method for each signal
#'     \item \code{beta_series}: whether each signal uses beta series
#'     \item \code{rm_zeros}: whether to remove zeros before convolution
#'     \item \code{convmax_1}: whether to normalize convolved max to 1
#'     \item \code{add_derivs}: whether to add temporal derivatives
#'   }
#' @keywords internal
#' @noRd
extract_signal_config <- function(signals_expanded) {
  bdm_args <- list()
  
  # Extract normalization for each regressor
  bdm_args$normalizations <- sapply(signals_expanded, function(s) {
    ifelse(is.null(s$normalization), "none", s$normalization)
  })
  
  # Extract beta_series settings
  bdm_args$beta_series <- sapply(signals_expanded, function(s) {
    ifelse(isTRUE(s$beta_series), TRUE, FALSE)
  })
  
  # Extract rm_zeros settings (default TRUE)
  bdm_args$rm_zeros <- sapply(signals_expanded, function(s) {
    if (is.null(s$rm_zeros)) {
      TRUE
    } else {
      if (!s$rm_zeros %in% c(TRUE, FALSE)) {
        stop("Don't know how to interpret rm_zeros setting of: ", s$rm_zeros)
      } else {
        as.logical(s$rm_zeros)
      }
    }
  })
  
  # Extract convmax_1 settings (default FALSE)
  bdm_args$convmax_1 <- sapply(signals_expanded, function(s) {
    if (is.null(s$convmax_1)) {
      FALSE
    } else {
      if (!s$convmax_1 %in% c(TRUE, FALSE)) {
        stop("Don't know how to interpret convmax_1 setting of: ", s$convmax_1)
      } else {
        as.logical(s$convmax_1)
      }
    }
  })
  
  # Extract add_deriv settings (default FALSE)
  bdm_args$add_derivs <- sapply(signals_expanded, function(s) {
    if (is.null(s$add_deriv)) {
      FALSE
    } else {
      if (!s$add_deriv %in% c(TRUE, FALSE)) {
        stop("Don't know how to interpret add_deriv setting of: ", s$add_deriv)
      } else {
        as.logical(s$add_deriv)
      }
    }
  })
  
  return(bdm_args)
}
