#' Lookup FSL FEAT output images by model and contrast
#'
#' Build a tidy lookup table that maps the human-readable pipeline model and
#' contrast names onto the FEAT folder structure for level 1, 2, and/or 3 FSL
#' outputs. The default result separates the structurally different FEAT
#' levels into \code{$l1}, \code{$l2}, and \code{$l3} tables with stable
#' level-specific schemas; levels that were not requested are present as empty
#' tables.
#'
#' L2 passthroughs remain logical L2 rows (\code{output_level = 2}) but are
#' marked with \code{is_passthrough = TRUE},
#' \code{materialized_level = 1}, and
#' \code{analysis_status = "passthrough"}. Because no L2 FEAT analysis is
#' materialized, \code{analysis_fsf} and \code{analysis_dir} are missing and
#' \code{image_feat_dir} identifies the originating L1 FEAT directory. A
#' passthrough has only a \code{"cope"} image row.
#'
#' @param gpa a \code{glm_pipeline_arguments} object. The function is most
#'   informative when populated FSL setup tables are present, but \code{source =
#'   "auto"} can also search scheduler caches and crawl FEAT folders.
#' @param level integer vector selecting FEAT levels to lookup. Valid values are
#'   \code{1}, \code{2}, and \code{3}. Defaults to all levels.
#' @param what statistic image types to include. Valid values are \code{"cope"},
#'   \code{"varcope"}, \code{"zstat"}, and \code{"tstat"}.
#' @param include_missing if \code{TRUE}, include expected output paths even when
#'   the image does not exist. If \code{FALSE}, return only existing images.
#' @param include_internal if \code{TRUE}, retain setup/debug columns such as
#'   FEAT execution timestamps and raw FEAT status fields. The default \code{FALSE}
#'   returns compact, level-specific user-facing tables.
#' @param format return format. \code{"by_level"} returns a
#'   \code{feat_output_lookup} object with \code{$l1}, \code{$l2}, and
#'   \code{$l3} tables. \code{"long"} returns their combined data.frame form.
#' @param source where to look for output metadata. \code{"setup"} uses only the
#'   \code{gpa$l*_model_setup$fsl} tables, \code{"cache"} searches scheduler
#'   \code{run_pipeline_cache*.RData} files, \code{"filesystem"} crawls FEAT
#'   folders, and \code{"auto"} tries these in that order.
#' @param cache_dir optional directory containing scheduler batch caches. If
#'   \code{NULL}, caches are searched below \code{gpa$output_directory} and
#'   \code{gpa$output_locations$scheduler_scripts}.
#' @param refresh_status if \code{TRUE}, refresh the FEAT status columns in the
#'   setup table before building the lookup.
#' @param lg optional \code{lgr::Logger} object.
#'
#' @return With \code{format = "by_level"}, a \code{feat_output_lookup}
#'   object containing level-specific \code{$l1}, \code{$l2}, and \code{$l3}
#'   data.frames. With \code{format = "long"}, a data.frame with one row per
#'   analysis instance, contrast, and statistic image. The long form adds
#'   \code{output_model}, \code{output_cope_number}, and
#'   \code{output_cope_name} as level-independent convenience columns.
#' @export
lookup_feat_outputs <- function(gpa, level = c(1L, 2L, 3L),
                                what = c("cope", "varcope", "zstat", "tstat"),
                                include_missing = TRUE,
                                include_internal = FALSE,
                                source = c("auto", "setup", "cache", "filesystem"),
                                cache_dir = NULL,
                                refresh_status = FALSE,
                                lg = NULL,
                                format = c("by_level", "long")) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, any.missing = FALSE)
  checkmate::assert_subset(what, c("cope", "varcope", "zstat", "tstat"))
  checkmate::assert_logical(include_missing, len = 1L)
  checkmate::assert_logical(include_internal, len = 1L)
  format <- match.arg(format)
  source <- match.arg(source)
  checkmate::assert_string(source)
  checkmate::assert_string(cache_dir, null.ok = TRUE)
  checkmate::assert_logical(refresh_status, len = 1L)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (is.null(lg)) lg <- lgr::get_logger("glm_pipeline/lookup_feat_outputs")

  level <- unique(as.integer(level))
  what <- unique(what)

  if (isTRUE(refresh_status)) {
    for (ll in level) {
      gpa <- refresh_glm_status(gpa, level = ll, glm_software = "fsl", lg = lg)
    }
  }

  out <- list()
  remaining <- level

  if (source %in% c("auto", "setup")) {
    setup_out <- lookup_feat_outputs_from_setup(
      gpa, level = remaining, what = what, lg = lg,
      warn_missing = identical(source, "setup")
    )
    if (nrow(setup_out) > 0L) {
      setup_out$lookup_source <- "setup"
      out[[length(out) + 1L]] <- setup_out
      remaining <- setdiff(remaining, unique(setup_out$level))
    }
    if (identical(source, "setup")) remaining <- integer(0)
  }

  if (length(remaining) > 0L && source %in% c("auto", "cache")) {
    cache_out <- lookup_feat_outputs_from_cache(
      gpa = gpa, level = remaining, what = what,
      cache_dir = cache_dir, lg = lg
    )
    if (nrow(cache_out) > 0L) {
      cache_out$lookup_source <- "cache"
      out[[length(out) + 1L]] <- cache_out
      remaining <- setdiff(remaining, unique(cache_out$level))
    }
    if (identical(source, "cache")) remaining <- integer(0)
  }

  if (length(remaining) > 0L && source %in% c("auto", "filesystem")) {
    fs_out <- lookup_feat_outputs_from_filesystem(gpa, level = remaining, what = what, lg = lg)
    if (nrow(fs_out) > 0L) {
      fs_out$lookup_source <- "filesystem"
      out[[length(out) + 1L]] <- fs_out
      remaining <- setdiff(remaining, unique(fs_out$level))
    }
  }

  out <- dplyr::bind_rows(out)
  if (nrow(out) > 0L) {
    out$image_exists <- !is.na(out$image_file) & file.exists(out$image_file)
    if (isFALSE(include_missing)) {
      out <- out %>% dplyr::filter(.data$image_exists == TRUE)
    }
    out <- lookup_feat_add_public_fields(out)
  }

  result <- new_feat_output_lookup(
    out,
    requested_levels = level,
    include_internal = include_internal
  )
  if (identical(format, "long")) return(as.data.frame(result))
  result
}

lookup_feat_add_public_fields <- function(out) {
  out <- as.data.frame(out)
  out$output_level <- as.integer(out$level)

  is_passthrough <- out$output_level == 2L &
    out$l2_input_mode %in% "l1_cope_file_passthrough"
  normal_l2 <- out$output_level == 2L & !is_passthrough
  missing_l2_mode <- normal_l2 &
    (is.na(out$l2_input_mode) | !nzchar(as.character(out$l2_input_mode)))
  out$l2_input_mode[missing_l2_mode] <- "cope_files"

  out$is_passthrough <- is_passthrough
  out$materialized_level <- ifelse(is_passthrough, 1L, out$output_level)

  out$output_model <- NA_character_
  out$output_cope_number <- NA_integer_
  out$output_cope_name <- NA_character_
  for (ll in 1:3) {
    use <- out$output_level == ll
    out$output_model[use] <- as.character(out[[paste0("l", ll, "_model")]][use])
    out$output_cope_number[use] <- as.integer(out[[paste0("l", ll, "_cope_number")]][use])
    out$output_cope_name[use] <- as.character(out[[paste0("l", ll, "_cope_name")]][use])
  }

  out$analysis_fsf <- as.character(out$feat_fsf)
  out$analysis_dir <- as.character(out$feat_dir)
  out$analysis_fsf[is_passthrough] <- NA_character_
  out$analysis_dir[is_passthrough] <- NA_character_

  out$image_feat_dir <- as.character(out$feat_dir)
  out$image_feat_dir[normal_l2] <- file.path(out$feat_dir[normal_l2], "cope1.feat")
  is_l3 <- out$output_level == 3L
  out$image_feat_dir[is_l3] <- file.path(out$feat_dir[is_l3], "cope1.feat")
  out$image_stats_dir <- ifelse(
    is.na(out$image_file),
    NA_character_,
    dirname(out$image_file)
  )

  out$image_cope_number <- NA_integer_
  image_pattern <- "^.*(?:varcope|cope|zstat|tstat)([0-9]+)\\.nii(?:\\.gz)?$"
  has_image_number <- !is.na(out$image_file) &
    grepl(image_pattern, out$image_file, perl = TRUE)
  out$image_cope_number[has_image_number] <- as.integer(sub(
    image_pattern,
    "\\1",
    out$image_file[has_image_number],
    perl = TRUE
  ))
  missing_image_number <- is.na(out$image_cope_number)
  out$image_cope_number[missing_image_number] <-
    out$output_cope_number[missing_image_number]

  failed <- out$feat_failed %in% TRUE
  complete <- out$feat_complete %in% TRUE
  dir_exists <- out$feat_dir_exists %in% TRUE |
    (!is.na(out$feat_dir) & dir.exists(out$feat_dir))
  out$analysis_status <- rep("not_started", nrow(out))
  out$analysis_status[dir_exists] <- "incomplete"
  out$analysis_status[failed] <- "failed"
  out$analysis_status[complete] <- "complete"
  out$analysis_status[is_passthrough] <- "passthrough"

  out
}

lookup_feat_level_user_columns <- function(level) {
  image_cols <- c(
    "materialized_level", "statistic", "image_cope_number",
    "image_file", "image_exists",
    "analysis_fsf", "analysis_dir", "image_feat_dir", "image_stats_dir",
    "analysis_status", "lookup_source", "cache_file"
  )

  switch(as.character(level),
    "1" = c(
      "output_level", "id", "session", "run_number",
      "l1_model", "l1_cope_number", "l1_cope_name",
      image_cols
    ),
    "2" = c(
      "output_level", "id", "session", "l2_scope",
      "l1_model", "l1_cope_number", "l1_cope_name",
      "l2_model", "l2_cope_number", "l2_cope_name",
      "l2_input_mode", "is_passthrough",
      image_cols
    ),
    "3" = c(
      "output_level", "session",
      "l1_model", "l1_cope_number", "l1_cope_name",
      "l2_model", "l2_cope_number", "l2_cope_name",
      "l3_model", "l3_cope_number", "l3_cope_name",
      "l3_input_mode",
      image_cols
    )
  )
}

lookup_feat_level_internal_columns <- function(level) {
  common <- c(
    "feat_fsf", "feat_dir", "stats_dir",
    "feat_complete", "feat_failed", "feat_dir_exists", "feat_fsf_exists",
    "to_run", "feat_fsf_modified_date",
    "feat_execution_start", "feat_execution_end", "feat_execution_min",
    "feat_auto_retried"
  )
  level_specific <- switch(as.character(level),
    "1" = c("run_volumes", "run_nifti", "l1_confound_file"),
    "2" = c(
      "l2_passthrough", "passthrough_cope_file", "l2_cope_dir",
      "n_l2_copes", "n_input_files"
    ),
    "3" = character(0)
  )
  c(common, level_specific)
}

lookup_feat_column_type <- function(column) {
  if (column %in% c(
    "output_level", "session", "run_number",
    "l1_cope_number", "l2_cope_number", "l3_cope_number",
    "materialized_level", "image_cope_number",
    "n_l2_copes", "n_input_files", "run_volumes"
  )) return("integer")
  if (column %in% c(
    "is_passthrough", "image_exists",
    "feat_complete", "feat_failed", "feat_dir_exists", "feat_fsf_exists",
    "to_run", "feat_auto_retried", "l2_passthrough"
  )) return("logical")
  if (column %in% "feat_execution_min") return("numeric")
  if (column %in% c(
    "feat_fsf_modified_date", "feat_execution_start", "feat_execution_end"
  )) return("POSIXct")
  "character"
}

lookup_feat_typed_column <- function(column, n = 0L) {
  switch(lookup_feat_column_type(column),
    integer = rep(NA_integer_, n),
    logical = rep(NA, n),
    numeric = rep(NA_real_, n),
    POSIXct = rep(as.POSIXct(NA), n),
    character = rep(NA_character_, n)
  )
}

lookup_feat_cast_column <- function(x, column) {
  switch(lookup_feat_column_type(column),
    integer = suppressWarnings(as.integer(x)),
    logical = as.logical(x),
    numeric = as.numeric(x),
    POSIXct = {
      if (inherits(x, "POSIXct")) x else as.POSIXct(x)
    },
    character = as.character(x)
  )
}

lookup_feat_level_table <- function(out, level, include_internal = FALSE) {
  user_cols <- lookup_feat_level_user_columns(level)
  columns <- user_cols
  if (isTRUE(include_internal)) {
    columns <- unique(c(columns, lookup_feat_level_internal_columns(level)))
  }

  if (!is.data.frame(out) || nrow(out) == 0L) {
    empty <- stats::setNames(
      lapply(columns, lookup_feat_typed_column),
      columns
    )
    return(as.data.frame(empty, optional = TRUE, stringsAsFactors = FALSE))
  }

  if (level == 2L) {
    id_scope <- out$l2_scope %in% "id"
    out$session[id_scope] <- NA_integer_
  } else if (level == 3L) {
    out$session[out$session %in% 0L] <- NA_integer_
  }

  for (column in columns) {
    if (!column %in% names(out)) {
      out[[column]] <- lookup_feat_typed_column(column, nrow(out))
    } else {
      out[[column]] <- lookup_feat_cast_column(out[[column]], column)
    }
  }

  out <- out[, columns, drop = FALSE]
  sort_cols <- intersect(
    c(
      "l3_model", "l2_model", "l1_model", "id", "session", "run_number",
      "l1_cope_number", "l2_cope_number", "l3_cope_number", "statistic"
    ),
    names(out)
  )
  if (length(sort_cols) > 0L && nrow(out) > 1L) {
    ord <- do.call(order, c(out[sort_cols], list(na.last = TRUE)))
    out <- out[ord, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

new_feat_output_lookup <- function(out, requested_levels = 1:3,
                                   include_internal = FALSE) {
  tables <- lapply(1:3, function(level) {
    level_rows <- if (is.data.frame(out) && nrow(out) > 0L) {
      out[out$output_level == level, , drop = FALSE]
    } else {
      data.frame()
    }
    lookup_feat_level_table(
      level_rows,
      level = level,
      include_internal = include_internal
    )
  })
  names(tables) <- c("l1", "l2", "l3")
  structure(
    tables,
    class = c("feat_output_lookup", "list"),
    requested_levels = as.integer(requested_levels),
    include_internal = isTRUE(include_internal)
  )
}

#' @rdname lookup_feat_outputs
#' @param x a \code{feat_output_lookup} object.
#' @param row.names,optional standard \code{as.data.frame} arguments.
#' @param ... unused.
#' @export
as.data.frame.feat_output_lookup <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
  out <- dplyr::bind_rows(lapply(x, as.data.frame))
  n <- nrow(out)
  out$output_model <- rep(NA_character_, n)
  out$output_cope_number <- rep(NA_integer_, n)
  out$output_cope_name <- rep(NA_character_, n)

  for (ll in 1:3) {
    use <- out$output_level == ll
    model_col <- paste0("l", ll, "_model")
    number_col <- paste0("l", ll, "_cope_number")
    name_col <- paste0("l", ll, "_cope_name")
    if (model_col %in% names(out)) {
      out$output_model[use] <- out[[model_col]][use]
      out$output_cope_number[use] <- out[[number_col]][use]
      out$output_cope_name[use] <- out[[name_col]][use]
    }
  }

  if ("is_passthrough" %in% names(out)) {
    out$is_passthrough[is.na(out$is_passthrough)] <- FALSE
  }

  front <- c(
    "output_level", "output_model", "output_cope_number", "output_cope_name"
  )
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]
  if (n > 1L) {
    sort_cols <- intersect(
      c(
        "output_level", "l3_model", "l2_model", "l1_model",
        "id", "session", "run_number",
        "l1_cope_number", "l2_cope_number", "l3_cope_number", "statistic"
      ),
      names(out)
    )
    ord <- do.call(order, c(out[sort_cols], list(na.last = TRUE)))
    out <- out[ord, , drop = FALSE]
  }
  rownames(out) <- row.names
  as.data.frame(out, optional = optional)
}

#' @rdname lookup_feat_outputs
#' @export
print.feat_output_lookup <- function(x, ...) {
  requested <- attr(x, "requested_levels")
  cat("<feat_output_lookup>\n")
  for (ll in 1:3) {
    label <- paste0("  L", ll, ": ")
    if (!ll %in% requested) {
      cat(label, "not requested\n", sep = "")
      next
    }
    tab <- x[[paste0("l", ll)]]
    existing <- if ("image_exists" %in% names(tab)) {
      sum(tab$image_exists %in% TRUE)
    } else {
      0L
    }
    cat(label, nrow(tab), " image row(s), ", existing, " existing\n", sep = "")
  }
  invisible(x)
}

lookup_feat_outputs_from_setup <- function(gpa, level, what, lg, warn_missing = TRUE) {
  out <- lapply(level, function(ll) {
    switch(as.character(ll),
      "1" = lookup_feat_outputs_l1(gpa, what = what, lg = lg, warn_missing = warn_missing),
      "2" = lookup_feat_outputs_l2(gpa, what = what, lg = lg, warn_missing = warn_missing),
      "3" = lookup_feat_outputs_l3(gpa, what = what, lg = lg, warn_missing = warn_missing)
    )
  })

  dplyr::bind_rows(out)
}

lookup_feat_outputs_l1 <- function(gpa, what, lg, warn_missing = TRUE) {
  setup <- get_fsl_setup_table(gpa, level = 1L, lg = lg, warn_missing = warn_missing)
  if (nrow(setup) == 0L) return(lookup_feat_outputs_empty())

  subj_df <- setup %>%
    dplyr::select("id", "session") %>%
    dplyr::distinct()

  model_set <- setup %>%
    dplyr::select("l1_model") %>%
    dplyr::distinct()

  cope_df <- get_l1_cope_df(gpa, model_set = model_set, subj_df = subj_df)

  base <- setup %>%
    dplyr::left_join(cope_df, by = c("id", "session", "l1_model")) %>%
    dplyr::mutate(
      level = 1L,
      l2_model = NA_character_,
      l2_cope_number = NA_integer_,
      l2_cope_name = NA_character_,
      l3_model = NA_character_,
      l3_cope_number = NA_integer_,
      l3_cope_name = NA_character_,
      stats_dir = file.path(.data$feat_dir, "stats")
    )

  lookup_feat_expand_stats(
    base = base,
    what = what,
    number_col = "l1_cope_number",
    name_col = "l1_cope_name"
  )
}

lookup_feat_outputs_l2 <- function(gpa, what, lg, warn_missing = TRUE) {
  setup <- get_fsl_setup_table(gpa, level = 2L, lg = lg, warn_missing = warn_missing)
  if (nrow(setup) == 0L) return(lookup_feat_outputs_empty())

  model_set <- setup %>%
    dplyr::select("l1_model", "l2_model") %>%
    dplyr::distinct()

  l2_df <- get_l2_cope_df(gpa, model_set = model_set)
  join_cols <- intersect(
    c("id", "session", "l1_model", "l1_cope_number", "l1_cope_name", "l2_model", "l2_input_mode"),
    intersect(names(setup), names(l2_df))
  )

  base <- setup %>%
    dplyr::inner_join(l2_df, by = join_cols) %>%
    dplyr::mutate(
      level = 2L,
      l3_model = NA_character_,
      l3_cope_number = NA_integer_,
      l3_cope_name = NA_character_,
      l2_cope_dir = dplyr::if_else(.data$l2_input_mode == "cope_files", "cope1.feat", NA_character_),
      stats_dir = dplyr::if_else(
        .data$l2_input_mode == "cope_files",
        file.path(.data$feat_dir, .data$l2_cope_dir, "stats"),
        NA_character_
      )
    )

  lookup <- lookup_feat_expand_stats(
    base = base,
    what = what,
    number_col = "l2_cope_number",
    name_col = "l2_cope_name"
  )

  is_passthrough <- lookup$l2_input_mode %in% "l1_cope_file_passthrough"
  lookup <- lookup[!is_passthrough | lookup$statistic == "cope", , drop = FALSE]
  is_passthrough <- lookup$l2_input_mode %in% "l1_cope_file_passthrough"
  if ("passthrough_cope_file" %in% names(lookup)) {
    lookup$image_file[is_passthrough] <- lookup$passthrough_cope_file[is_passthrough]
  }

  lookup
}

lookup_feat_outputs_l3 <- function(gpa, what, lg, warn_missing = TRUE) {
  setup <- get_fsl_setup_table(gpa, level = 3L, lg = lg, warn_missing = warn_missing)
  if (nrow(setup) == 0L) return(lookup_feat_outputs_empty())

  meta <- lookup_feat_l3_metadata(gpa, setup)
  join_cols <- intersect(
    c("l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model"),
    intersect(names(setup), names(meta))
  )

  base <- setup %>%
    dplyr::left_join(meta, by = join_cols, relationship = "many-to-many") %>%
    dplyr::mutate(
      level = 3L,
      id = NA_character_,
      run_number = NA_integer_,
      stats_dir = file.path(.data$feat_dir, "cope1.feat", "stats")
    )

  lookup_feat_expand_stats(
    base = base,
    what = what,
    number_col = "l3_cope_number",
    name_col = "l3_cope_name"
  )
}

get_fsl_setup_table <- function(gpa, level, lg, warn_missing = TRUE) {
  setup_name <- paste0("l", level, "_model_setup")
  if (!setup_name %in% names(gpa) || is.null(gpa[[setup_name]]) || is.null(gpa[[setup_name]]$fsl)) {
    if (isTRUE(warn_missing)) lg$warn("No FSL setup table found at gpa$%s$fsl.", setup_name)
    return(lookup_feat_outputs_empty())
  }

  setup <- gpa[[setup_name]]$fsl
  if (!is.data.frame(setup) || nrow(setup) == 0L) {
    if (isTRUE(warn_missing)) lg$warn("FSL setup table at gpa$%s$fsl is empty.", setup_name)
    return(lookup_feat_outputs_empty())
  }

  setup
}

lookup_feat_outputs_from_cache <- function(gpa, level, what, cache_dir = NULL, lg) {
  cache_files <- lookup_feat_cache_files(gpa, cache_dir = cache_dir)
  if (length(cache_files) == 0L) {
    lg$warn("No scheduler cache files found for lookup_feat_outputs cache fallback.")
    return(lookup_feat_outputs_empty())
  }

  out <- list()
  remaining <- level
  for (ff in cache_files) {
    cached_gpa <- lookup_feat_load_cache_gpa(ff, lg = lg)
    if (is.null(cached_gpa)) next

    cache_out <- lookup_feat_outputs_from_setup(
      gpa = cached_gpa, level = remaining, what = what, lg = lg,
      warn_missing = FALSE
    )
    if (nrow(cache_out) == 0L) next

    cache_out$cache_file <- ff
    out[[length(out) + 1L]] <- cache_out
    remaining <- setdiff(remaining, unique(cache_out$level))
    if (length(remaining) == 0L) break
  }

  dplyr::bind_rows(out)
}

lookup_feat_cache_files <- function(gpa, cache_dir = NULL) {
  roots <- cache_dir
  if (is.null(roots)) {
    roots <- c(
      gpa$output_locations$active_batch_directory,
      gpa$output_locations$scheduler_scripts,
      file.path(gpa$output_directory, "scheduler_scripts")
    )
  }
  roots <- unique(roots[!is.na(roots) & nzchar(roots) & dir.exists(roots)])
  if (length(roots) == 0L) return(character(0))

  cache_files <- unique(unlist(lapply(roots, function(root) {
    c(
      Sys.glob(file.path(root, "run_pipeline_cache*.RData")),
      Sys.glob(file.path(root, "batch_*", "run_pipeline_cache*.RData"))
    )
  }), use.names = FALSE))

  cache_files <- cache_files[file.exists(cache_files)]
  if (length(cache_files) == 0L) return(character(0))

  # Prefer backend-specific FSL caches, then shared caches, newest first.
  info <- file.info(cache_files)
  priority <- ifelse(grepl("_fsl\\.RData$", cache_files), 1L,
    ifelse(grepl("run_pipeline_cache\\.RData$", cache_files), 2L, 3L)
  )
  cache_files[order(priority, -as.numeric(info$mtime))]
}

lookup_feat_load_cache_gpa <- function(cache_file, lg) {
  env <- new.env(parent = emptyenv())
  ok <- tryCatch({
    load(cache_file, envir = env)
    TRUE
  }, error = function(e) {
    lg$warn("Could not load scheduler cache for lookup fallback: %s", cache_file)
    FALSE
  })
  if (!isTRUE(ok) || !"gpa" %in% ls(env)) return(NULL)
  gpa <- env$gpa
  if (!inherits(gpa, "glm_pipeline_arguments")) return(NULL)
  gpa
}

lookup_feat_outputs_from_filesystem <- function(gpa, level, what, lg) {
  out <- lapply(level, function(ll) {
    switch(as.character(ll),
      "1" = lookup_feat_filesystem_l1(gpa, what = what, lg = lg),
      "2" = lookup_feat_filesystem_l2(gpa, what = what, lg = lg),
      "3" = lookup_feat_filesystem_l3(gpa, what = what, lg = lg)
    )
  })

  dplyr::bind_rows(out)
}

lookup_feat_filesystem_l1 <- function(gpa, what, lg) {
  root <- file.path(gpa$output_directory, "feat_l1")
  feat_dirs <- lookup_feat_list_dirs(root, "\\.feat$")
  feat_dirs <- feat_dirs[grepl("FEAT_LVL1_run[0-9]+\\.feat$", feat_dirs)]
  if (length(feat_dirs) == 0L) return(lookup_feat_outputs_empty())

  out <- lapply(feat_dirs, function(feat_dir) {
    rel <- lookup_feat_relative_parts(root, feat_dir)
    id <- lookup_feat_extract_part(rel, "^sub-(.+)$")
    session <- lookup_feat_extract_part(rel, "^ses-(.+)$")
    run_number <- as.integer(sub(".*FEAT_LVL1_run([0-9]+)\\.feat$", "\\1", feat_dir))
    model_idx <- max(which(!grepl("^(sub-|ses-|FEAT_LVL1_run)", rel)))
    l1_model <- if (is.finite(model_idx)) rel[model_idx] else NA_character_
    cope_df <- scan_feat_stat_files(feat_dir, statistics = what, include_missing = TRUE)
    lookup_feat_rows_from_cope_df(
      cope_df = cope_df, level = 1L, what = what,
      feat_dir = feat_dir, feat_fsf = sub("\\.feat$", ".fsf", feat_dir),
      id = id, session = lookup_feat_as_integer(session), run_number = run_number,
      l1_model = l1_model,
      contrast_level = 1L
    )
  })

  dplyr::bind_rows(out)
}

lookup_feat_filesystem_l2 <- function(gpa, what, lg) {
  root <- file.path(gpa$output_directory, "feat_l2")
  gfeat_dirs <- lookup_feat_list_dirs(root, "\\.gfeat$")
  if (length(gfeat_dirs) == 0L) return(lookup_feat_outputs_empty())

  out <- lapply(gfeat_dirs, function(feat_dir) {
    rel <- lookup_feat_relative_parts(root, feat_dir)
    id <- lookup_feat_extract_part(rel, "^sub-(.+)$")
    session <- lookup_feat_extract_part(rel, "^ses-(.+)$")
    l1_model <- lookup_feat_strip_prefix(lookup_feat_extract_part(rel, "^L1m-(.+)$"), "L1m-")
    l1c <- lookup_feat_extract_part(rel, "^l1c-(.+)$")
    l1_info <- lookup_feat_parse_numbered_label(l1c)
    l2_model <- lookup_feat_strip_prefix(lookup_feat_extract_part(rel, "^L2m-(.+)\\.gfeat$|^L2m-(.+)$"), "L2m-")
    cope_df <- scan_feat_stat_files(
      file.path(feat_dir, "cope1.feat"),
      statistics = what,
      include_missing = TRUE
    )
    rows <- lookup_feat_rows_from_cope_df(
      cope_df = cope_df, level = 2L, what = what,
      feat_dir = feat_dir, feat_fsf = sub("\\.gfeat$", ".fsf", feat_dir),
      id = id, session = lookup_feat_as_integer(session), run_number = NA_integer_,
      l1_model = l1_model,
      l1_cope_number = l1_info$number, l1_cope_name = l1_info$name,
      l2_model = l2_model,
      stats_dir = file.path(feat_dir, "cope1.feat", "stats"),
      contrast_level = 2L
    )
    rows$l2_input_mode <- "cope_files"
    rows$l2_passthrough <- FALSE
    rows
  })

  dplyr::bind_rows(out)
}

lookup_feat_filesystem_l3 <- function(gpa, what, lg) {
  root <- file.path(gpa$output_directory, "feat_l3")
  gfeat_dirs <- lookup_feat_list_dirs(root, "\\.gfeat$")
  if (length(gfeat_dirs) == 0L) return(lookup_feat_outputs_empty())

  out <- lapply(gfeat_dirs, function(feat_dir) {
    rel <- lookup_feat_relative_parts(root, feat_dir)
    l3_model <- lookup_feat_strip_prefix(lookup_feat_extract_part(rel, "^L3m-(.+)$"), "L3m-")
    l1_model <- lookup_feat_strip_prefix(lookup_feat_extract_part(rel, "^L1m-(.+)$"), "L1m-")
    l2_model <- lookup_feat_strip_prefix(lookup_feat_extract_part(rel, "^L2m-(.+)$"), "L2m-")
    l2_cope_name <- lookup_feat_strip_prefix(lookup_feat_extract_part(rel, "^l2c-(.+)$"), "l2c-")
    l1_cope_name <- sub("^FEAT_l1c-(.+)\\.gfeat$", "\\1", basename(feat_dir))
    if (identical(l1_cope_name, basename(feat_dir))) l1_cope_name <- NA_character_
    cope_df <- scan_feat_stat_files(
      file.path(feat_dir, "cope1.feat"),
      statistics = what,
      include_missing = TRUE
    )
    lookup_feat_rows_from_cope_df(
      cope_df = cope_df, level = 3L, what = what,
      feat_dir = feat_dir, feat_fsf = sub("\\.gfeat$", ".fsf", feat_dir),
      id = NA_character_, session = NA_integer_, run_number = NA_integer_,
      l1_model = l1_model, l1_cope_name = l1_cope_name,
      l2_model = l2_model, l2_cope_name = l2_cope_name,
      l3_model = l3_model,
      stats_dir = file.path(feat_dir, "cope1.feat", "stats"),
      contrast_level = 3L
    )
  })

  dplyr::bind_rows(out)
}

lookup_feat_list_dirs <- function(root, pattern) {
  if (!dir.exists(root)) return(character(0))
  dirs <- list.dirs(root, full.names = TRUE, recursive = TRUE)
  dirs[grepl(pattern, dirs)]
}

lookup_feat_relative_parts <- function(root, path) {
  rel <- sub(paste0("^", gsub("([][{}()+*?.^$|\\\\])", "\\\\\\1", normalizePath(root, mustWork = FALSE)), "/?"), "", normalizePath(path, mustWork = FALSE))
  strsplit(rel, .Platform$file.sep, fixed = TRUE)[[1]]
}

lookup_feat_extract_part <- function(parts, pattern) {
  idx <- grep(pattern, parts, perl = TRUE)
  if (length(idx) == 0L) return(NA_character_)
  val <- parts[idx[1L]]
  m <- regexec(pattern, val, perl = TRUE)
  reg <- regmatches(val, m)[[1]]
  reg <- reg[-1L]
  reg <- reg[!is.na(reg) & nzchar(reg)]
  if (length(reg) == 0L) return(val)
  reg[1L]
}

lookup_feat_strip_prefix <- function(x, prefix) {
  if (is.na(x)) return(NA_character_)
  sub(paste0("^", prefix), "", x)
}

lookup_feat_parse_numbered_label <- function(x) {
  if (is.na(x) || !nzchar(x)) return(list(number = NA_integer_, name = NA_character_))
  m <- regexec("^([0-9]+)_(.+)$", x)
  reg <- regmatches(x, m)[[1]]
  if (length(reg) == 3L) {
    return(list(number = as.integer(reg[2L]), name = reg[3L]))
  }
  list(number = NA_integer_, name = x)
}

lookup_feat_as_integer <- function(x) {
  if (is.na(x) || !nzchar(x)) return(NA_integer_)
  suppressWarnings(as.integer(x))
}

lookup_feat_rows_from_cope_df <- function(cope_df, level, what, feat_dir, feat_fsf,
                                          id = NA_character_, session = NA_integer_,
                                          run_number = NA_integer_,
                                          l1_model = NA_character_,
                                          l1_cope_number = NA_integer_,
                                          l1_cope_name = NA_character_,
                                          l2_model = NA_character_,
                                          l2_cope_number = NA_integer_,
                                          l2_cope_name = NA_character_,
                                          l3_model = NA_character_,
                                          l3_cope_number = NA_integer_,
                                          l3_cope_name = NA_character_,
                                          stats_dir = NULL,
                                          contrast_level) {
  if (is.null(cope_df) || nrow(cope_df) == 0L) return(lookup_feat_outputs_empty())
  if (is.null(stats_dir)) stats_dir <- file.path(feat_dir, "stats")

  if ("zstat" %in% what && "z" %in% names(cope_df) && !"zstat" %in% names(cope_df)) {
    cope_df$zstat <- cope_df$z
  }
  if ("tstat" %in% what && "t" %in% names(cope_df) && !"tstat" %in% names(cope_df)) {
    cope_df$tstat <- cope_df$t
  }

  rows <- cope_df
  rows <- rows[, intersect(c("cope_number", "contrast_name", what), names(rows)), drop = FALSE]
  if (!any(what %in% names(rows))) return(lookup_feat_outputs_empty())
  rows <- rows %>%
    tidyr::pivot_longer(
      cols = tidyselect::any_of(what),
      names_to = "statistic",
      values_to = "image_file"
    ) %>%
    dplyr::mutate(
      level = level,
      id = id,
      session = session,
      run_number = run_number,
      l1_model = l1_model,
      l1_cope_number = if (contrast_level == 1L) .data$cope_number else l1_cope_number,
      l1_cope_name = if (contrast_level == 1L) .data$contrast_name else l1_cope_name,
      l2_model = l2_model,
      l2_cope_number = if (contrast_level == 2L) .data$cope_number else l2_cope_number,
      l2_cope_name = if (contrast_level == 2L) .data$contrast_name else l2_cope_name,
      l3_model = l3_model,
      l3_cope_number = if (contrast_level == 3L) .data$cope_number else l3_cope_number,
      l3_cope_name = if (contrast_level == 3L) .data$contrast_name else l3_cope_name,
      stat_number = .data$cope_number,
      feat_fsf = feat_fsf,
      feat_dir = feat_dir,
      stats_dir = stats_dir,
      feat_complete = file.exists(file.path(feat_dir, ".feat_complete")),
      feat_failed = ifelse(file.exists(file.path(feat_dir, ".feat_fail")), TRUE, NA)
    )

  lookup_feat_select_columns(rows)
}

lookup_feat_l3_metadata <- function(gpa, setup) {
  if (!is.null(gpa$l3_model_setup$metadata) &&
      is.data.frame(gpa$l3_model_setup$metadata) &&
      nrow(gpa$l3_model_setup$metadata) > 0L) {
    meta_cols <- intersect(
      c(
        "l1_model", "l1_cope_number", "l1_cope_name",
        "l2_model", "l2_cope_number", "l2_cope_name",
        "l3_model", "l3_cope_number", "l3_cope_name"
      ),
      names(gpa$l3_model_setup$metadata)
    )
    return(gpa$l3_model_setup$metadata %>%
      dplyr::select(dplyr::all_of(meta_cols)) %>%
      dplyr::distinct())
  }

  model_set <- setup %>%
    dplyr::select(dplyr::any_of(c("l1_model", "l2_model", "l3_model"))) %>%
    dplyr::distinct()

  l3_df <- get_l3_cope_df(gpa, model_set = model_set) %>%
    dplyr::select("l3_model", "l3_cope_number", "l3_cope_name") %>%
    dplyr::distinct()

  l1_df <- NULL
  if (all(c("l1_model", "l1_cope_name") %in% names(setup))) {
    l1_df <- get_l1_cope_df(gpa, model_set = model_set) %>%
      dplyr::select("l1_model", "l1_cope_number", "l1_cope_name") %>%
      dplyr::distinct()
  }

  l2_df <- NULL
  if (all(c("l2_model", "l2_cope_name") %in% names(setup)) &&
      !is.null(gpa$l2_model_setup) &&
      !is.null(gpa$l2_model_setup$fsl)) {
    l2_df <- get_l2_cope_df(gpa, model_set = model_set) %>%
      dplyr::select(
        dplyr::any_of(c(
          "l1_model", "l1_cope_number", "l1_cope_name",
          "l2_model", "l2_cope_number", "l2_cope_name"
        ))
      ) %>%
      dplyr::distinct()
  }

  meta <- setup %>%
    dplyr::select(dplyr::any_of(c("l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model"))) %>%
    dplyr::distinct() %>%
    dplyr::left_join(l3_df, by = "l3_model", relationship = "many-to-many")

  if (!is.null(l1_df)) {
    meta <- meta %>% dplyr::left_join(l1_df, by = c("l1_model", "l1_cope_name"))
  }

  if (!is.null(l2_df)) {
    meta <- meta %>% dplyr::left_join(
      l2_df,
      by = intersect(
        c("l1_model", "l1_cope_number", "l1_cope_name", "l2_model", "l2_cope_name"),
        intersect(names(meta), names(l2_df))
      )
    )
  }

  meta
}

lookup_feat_expand_stats <- function(base, what, number_col, name_col) {
  if (nrow(base) == 0L) return(lookup_feat_outputs_empty())

  expanded <- lapply(what, function(stat) {
    out <- base
    out$statistic <- stat
    out$stat_number <- out[[number_col]]
    out$contrast_name <- out[[name_col]]
    out$image_file <- ifelse(
      is.na(out$stats_dir) | is.na(out$stat_number),
      NA_character_,
      file.path(out$stats_dir, paste0(stat, out$stat_number, ".nii.gz"))
    )
    out
  })

  out <- dplyr::bind_rows(expanded)
  lookup_feat_select_columns(out)
}

lookup_feat_select_columns <- function(out, include_internal = TRUE) {
  out <- as.data.frame(out)

  user_cols <- c(
    "level", "id", "session", "run_number",
    "l1_model", "l1_cope_number", "l1_cope_name",
    "l2_model", "l2_cope_number", "l2_cope_name",
    "l3_model", "l3_cope_number", "l3_cope_name",
    "contrast_name", "statistic", "stat_number",
    "feat_fsf", "feat_dir", "stats_dir", "image_file", "image_exists",
    "feat_complete", "feat_failed", "feat_dir_exists", "feat_fsf_exists",
    "lookup_source", "cache_file"
  )

  internal_cols <- c(
    "l2_scope", "l2_input_mode", "l2_passthrough", "passthrough_cope_file",
    "l2_cope_dir", "l3_input_mode", "to_run", "run_volumes", "run_nifti",
    "l1_confound_file", "n_l2_copes", "feat_fsf_modified_date",
    "feat_execution_start", "feat_execution_end", "feat_execution_min",
    "feat_auto_retried"
  )
  canonical <- if (isTRUE(include_internal)) c(user_cols, internal_cols) else user_cols

  missing_cols <- setdiff(canonical, names(out))
  for (cc in missing_cols) {
    out[[cc]] <- NA
  }

  list_cols <- names(out)[vapply(out, is.list, logical(1L))]
  out <- out[, setdiff(names(out), list_cols), drop = FALSE]

  extra_cols <- character(0)
  if (isTRUE(include_internal)) {
    extra_cols <- setdiff(names(out), canonical)
  }

  out %>%
    dplyr::select(dplyr::all_of(intersect(canonical, names(out))), dplyr::all_of(extra_cols))
}

lookup_feat_outputs_empty <- function() {
  data.frame(
    level = integer(0),
    id = character(0),
    session = integer(0),
    run_number = integer(0),
    l1_model = character(0),
    l1_cope_number = integer(0),
    l1_cope_name = character(0),
    l2_model = character(0),
    l2_cope_number = integer(0),
    l2_cope_name = character(0),
    l3_model = character(0),
    l3_cope_number = integer(0),
    l3_cope_name = character(0),
    contrast_name = character(0),
    statistic = character(0),
    stat_number = integer(0),
    feat_fsf = character(0),
    feat_dir = character(0),
    stats_dir = character(0),
    image_file = character(0),
    image_exists = logical(0),
    stringsAsFactors = FALSE
  )
}
