#' Build concatenated L1 confounds for multi-run analyses
#'
#' Internal helper that uses stored confound column metadata to align shared
#' columns and keep run-specific spike regressors separate with zero padding.
#' @keywords internal
#' @noRd
read_l1_confound_columns <- function(gpa, id, session, run_number, lg = NULL) {
  if (is.null(lg)) lg <- lgr::get_logger("glm_pipeline/l1_setup")
  col_info <- read_df_sqlite(
    gpa = gpa, id = id, session = session, run_number = run_number,
    table = "l1_confound_columns", drop_keys = TRUE, quiet = TRUE
  )
  if (is.null(col_info)) {
    lg$warn("Missing l1_confound_columns for id: %s, session: %s, run_number: %s", id, session, run_number)
    return(NULL)
  }

  if (!all(c("col_index", "col_name") %in% names(col_info))) {
    lg$warn("l1_confound_columns missing required fields for id: %s, session: %s, run_number: %s", id, session, run_number)
    return(NULL)
  }

  col_info <- col_info[order(col_info$col_index), , drop = FALSE]
  if (!"is_spike" %in% names(col_info)) {
    col_info$is_spike <- grepl("_spike_", col_info$col_name)
  } else {
    col_info$is_spike <- as.logical(col_info$is_spike)
  }

  col_info$col_name <- make.unique(col_info$col_name)
  return(col_info)
}

concat_l1_confounds <- function(gpa, id, session, run_numbers, confound_files, output_dir,
                                lg = NULL, file_name = NULL, run_sessions = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_scalar(id, null.ok = FALSE)
  checkmate::assert_integerish(session, null.ok = FALSE)
  checkmate::assert_integerish(run_numbers, null.ok = FALSE)
  checkmate::assert_character(confound_files, null.ok = TRUE)
  checkmate::assert_string(output_dir, null.ok = FALSE)
  checkmate::assert_integerish(run_sessions, null.ok = TRUE, any.missing = FALSE)

  if (is.null(lg)) lg <- lgr::get_logger("glm_pipeline/l1_setup")

  if (length(confound_files) == 0L) return(NULL)
  if (length(confound_files) != length(run_numbers)) {
    lg$warn("Confound files length (%d) does not match run_numbers length (%d).", length(confound_files), length(run_numbers))
    return(NULL)
  }
  if (is.null(run_sessions)) {
    run_sessions <- rep.int(as.integer(session), length(run_numbers))
  }
  run_sessions <- as.integer(run_sessions)
  if (length(run_sessions) != length(run_numbers)) {
    lg$warn("run_sessions length (%d) does not match run_numbers length (%d).", length(run_sessions), length(run_numbers))
    return(NULL)
  }

  valid <- !is.na(confound_files) & nzchar(confound_files) & file.exists(confound_files)
  if (!all(valid)) {
    lg$warn("Some l1_confound_files are missing or invalid; cannot concatenate confounds.")
    lg$warn("Missing: %s", confound_files[!valid])
    return(NULL)
  }

  run_info <- lapply(seq_along(run_numbers), function(ii) {
    run_no <- run_numbers[ii]
    run_session <- run_sessions[ii]
    col_info <- read_l1_confound_columns(gpa, id, run_session, run_no, lg = lg)
    if (is.null(col_info)) return(NULL)
    # Spike indices restart in every run. Qualify them by session and run so
    # unrelated contaminated volumes never share a fitted nuisance coefficient.
    col_info$col_name[col_info$is_spike] <- paste0(
      "ses-", run_session, "_run-", run_no, "_", col_info$col_name[col_info$is_spike]
    )
    list(
      run_number = run_no,
      session = run_session,
      file = confound_files[ii],
      col_names = col_info$col_name,
      is_spike = col_info$is_spike
    )
  })

  if (any(vapply(run_info, is.null, logical(1)))) {
    lg$warn("Cannot concatenate confounds due to missing column metadata.")
    return(NULL)
  }

  confound_dfs <- lapply(run_info, function(info) {
    df <- data.table::fread(info$file, header = FALSE, data.table = FALSE)
    if (ncol(df) != length(info$col_names)) {
      lg$warn("Confound column count mismatch for %s: file has %d, metadata has %d.",
              info$file, ncol(df), length(info$col_names))
      return(NULL)
    }
    names(df) <- info$col_names
    df
  })

  if (any(vapply(confound_dfs, is.null, logical(1)))) {
    lg$warn("Cannot concatenate confounds due to column mismatches.")
    return(NULL)
  }

  ordered_unique <- function(x) x[!duplicated(x)]
  shared_cols <- ordered_unique(unlist(lapply(run_info, function(info) info$col_names[!info$is_spike])))
  spike_cols <- ordered_unique(unlist(lapply(run_info, function(info) info$col_names[info$is_spike])))
  all_cols <- c(shared_cols, spike_cols)

  build_padded <- function(df, col_names) {
    out <- as.data.frame(matrix(0, nrow = nrow(df), ncol = length(all_cols)))
    names(out) <- all_cols
    present <- intersect(col_names, all_cols)
    out[, present] <- df[, present, drop = FALSE]
    out
  }

  padded_runs <- lapply(seq_along(confound_dfs), function(ii) {
    build_padded(confound_dfs[[ii]], run_info[[ii]]$col_names)
  })

  concat_df <- do.call(rbind, padded_runs)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(file_name) || !nzchar(file_name)) {
    file_name <- paste0("l1_confounds_concat_s", session, ".txt")
  }
  out_file <- file.path(output_dir, file_name)
  write.table(concat_df, file = out_file, row.names = FALSE, col.names = FALSE)

  return(out_file)
}
