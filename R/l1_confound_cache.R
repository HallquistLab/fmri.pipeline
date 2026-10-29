# Per-run cache provenance. A cache without a completed, matching manifest is
# rebuilt from source; file existence alone does not establish validity.
l1_confound_input_fingerprint <- function(path, checksum = TRUE) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !nzchar(path)) return(NULL)
  path <- normalizePath(path, mustWork = FALSE)
  info <- file.info(path)
  list(path = path, size = unname(info$size), mtime = as.numeric(info$mtime),
    md5 = if (checksum && file.exists(path)) unname(tools::md5sum(path)) else NULL)
}

prepare_l1_confound_cache <- function(gpa, run_df, demean, analysis_outdir) {
  path <- file.path(analysis_outdir, paste0("run", run_df$run_number, "_l1_confounds_cache.rds"))
  old <- if (file.exists(path)) tryCatch(readRDS(path), error = function(e) NULL) else NULL
  if (!is.list(old) || !identical(old$version, 1L)) old <- NULL

  # finalize_pipeline_configuration passes the previously processed run back on
  # subsequent calls. Restore its original image before applying drops again.
  if (!is.null(old$output) && !identical(old$output$run_nifti, old$source$run_nifti) &&
      identical(as.character(run_df$run_nifti), old$output$run_nifti)) {
    run_df$run_nifti <- old$source$run_nifti
    run_df$run_volumes <- old$source$run_volumes
  } else if (is.null(old) && isTRUE(run_df$orig_volumes > run_df$run_volumes)) {
    stop("Cannot rebuild confounds for a previously truncated run without source provenance. ",
      "Supply the original run_nifti and run_volumes, then rerun configuration.", call. = FALSE)
  }
  inputs <- lapply(c("run_nifti", "confound_input_file", "motion_params_file"), function(column) {
    if (!column %in% names(run_df) || is.na(run_df[[column]][1L])) return(NULL)
    l1_confound_input_fingerprint(get_mr_abspath(run_df, column)[1L], checksum = column != "run_nifti")
  })
  names(inputs) <- c("run_nifti", "confound_input_file", "motion_params_file")
  signature <- list(
    settings = gpa$confound_settings, demean = demean, tr = get_run_tr(run_df, gpa),
    run = as.list(run_df[1L, intersect(c("id", "session", "run_number", "run_volumes", "drop_volumes",
      "confound_input_file_present", "motion_params_present",
      setdiff(grep("^last_", names(run_df), value = TRUE), "last_volume")),
      names(run_df)), drop = FALSE]),
    inputs = inputs
  )
  list(path = path, signature = signature, previous = old, run_df = run_df,
    current = !is.null(old) && identical(old$signature, signature),
    refresh_image = is.null(old) || !identical(old$signature$inputs$run_nifti, inputs$run_nifti))
}

write_l1_confound_cache <- function(cache, output = NULL) {
  value <- list(version = 1L,
    source = list(run_nifti = as.character(cache$run_df$run_nifti), run_volumes = cache$run_df$run_volumes),
    output = if (is.null(output)) cache$previous$output else list(run_nifti = as.character(output$run_nifti)),
    signature = if (is.null(output)) NULL else cache$signature)
  temporary <- tempfile(".confound-cache-", tmpdir = dirname(cache$path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary)
  # Windows cannot rename over an existing file. Removing the old manifest is
  # safe: an interrupted write leaves an invalid cache, never a false cache hit.
  if (file.exists(cache$path)) unlink(cache$path)
  if (!file.rename(temporary, cache$path)) stop("Could not write confound cache manifest: ", cache$path)
}

invalidate_l1_confound_cache <- function(gpa, id, session, run_number) {
  con <- DBI::dbConnect(RSQLite::SQLite(), gpa$output_locations$sqlite_db)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  RSQLite::sqliteSetBusyHandler(con, get_default_sqlite_busy_timeout())
  DBI::dbExecute(con, "BEGIN IMMEDIATE")
  committed <- FALSE
  on.exit(if (!committed) DBI::dbRollback(con), add = TRUE, after = FALSE)
  for (table in c("l1_run_calculations", "l1_confounds", "l1_confound_columns",
                  "l1_exclusion_data", "l1_truncation_data", "l1_confound_inputs", "l1_motion_parameters")) {
    if (DBI::dbExistsTable(con, table)) {
      query <- glue::glue_sql("DELETE FROM {`table`} WHERE id = {id} AND session = {session} AND run_number = {run_number}", .con = con)
      DBI::dbExecute(con, query)
    }
  }
  DBI::dbCommit(con)
  committed <- TRUE
  invisible(NULL)
}
