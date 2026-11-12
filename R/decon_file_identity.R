# Legacy filenames and BIDS-style source paths remain a fallback. Explicit
# metadata is needed when the naming convention does not identify a run.
parse_decon_identity <- function(path, fields = c("id", "run", "session")) {
  labels <- c(id = "sub", run = "run", session = "ses")
  ignore_words <- c(
    "subject", "subjects", "subcortical", "subset", "subsets", "subspace",
    "subanalysis", "subproject", "subprojects",
    "session", "sessions",
    "run", "runs", "running", "runtime", "runtimes"
  )
  base_tokens <- unlist(strsplit(basename(path), "[/\\\\._]", perl = TRUE), use.names = FALSE)
  dir_part <- dirname(path)
  dir_tokens <- if (!identical(dir_part, ".") && nzchar(dir_part)) {
    unlist(strsplit(dir_part, "[/\\\\._]", perl = TRUE), use.names = FALSE)
  } else {
    character(0)
  }

  extract_values <- function(toks, label, allow_legacy_alpha) {
    if (!length(toks)) return(character(0))
    toks <- toks[nzchar(toks) & !(tolower(toks) %in% ignore_words)]
    if (!length(toks)) return(character(0))
    explicit_pat <- paste0("^", label, "(-[[:alnum:]]+|[0-9]+[[:alnum:]]*)$")
    explicit_hits <- toks[grepl(explicit_pat, toks)]
    if (length(explicit_hits)) {
      return(unique(sub(paste0("^", label, "-?"), "", explicit_hits)))
    }
    if (!allow_legacy_alpha) return(character(0))
    legacy_pat <- paste0("^", label, "[[:alpha:]][[:alnum:]]*$")
    legacy_hits <- toks[grepl(legacy_pat, toks)]
    unique(sub(paste0("^", label), "", legacy_hits))
  }

  result <- list()
  for (field in fields) {
    label <- labels[[field]]
    values <- extract_values(base_tokens, label, allow_legacy_alpha = TRUE)
    if (!length(values) && length(dir_tokens)) {
      values <- extract_values(dir_tokens, label, allow_legacy_alpha = FALSE)
    }
    if (length(values) > 1L) {
      stop("Ambiguous ", field, " in path; supply file_manifest: ", path)
    }
    if (length(values)) result[[field]] <- values
  }
  result
}

validate_decon_manifest <- function(file_manifest, d_files) {
  if (is.null(file_manifest)) return(NULL)
  checkmate::assert_data_frame(file_manifest)
  checkmate::assert_subset(c("file", "id", "run"), names(file_manifest))
  file_manifest <- as.data.frame(file_manifest)
  file_manifest$file <- as.character(file_manifest$file)
  if (anyNA(file_manifest$file) || any(!nzchar(file_manifest$file))) {
    stop("file_manifest must have nonmissing, nonempty file paths")
  }
  file_manifest$file <- normalizePath(file_manifest$file, winslash = "/", mustWork = FALSE)
  if (anyDuplicated(file_manifest$file)) {
    stop("file_manifest must have one unique, nonmissing file path per row")
  }
  for (column in c("id", "run")) {
    if (anyNA(file_manifest[[column]]) || any(!nzchar(as.character(file_manifest[[column]])))) {
      stop("file_manifest must contain nonmissing id and run values")
    }
  }
  if (any(!normalizePath(d_files, winslash = "/", mustWork = FALSE) %in% file_manifest$file)) {
    stop("file_manifest is missing one or more requested deconvolution files")
  }
  file_manifest
}

resolve_decon_identity <- function(fname, data, file_manifest = NULL) {
  if (!is.null(file_manifest)) {
    row <- file_manifest[file_manifest$file == normalizePath(fname, winslash = "/"), , drop = FALSE]
  } else if (file.exists(paste0(fname, ".metadata.csv"))) {
    row <- as.data.frame(readr::read_csv(paste0(fname, ".metadata.csv"),
      col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE))
  } else if (all(c("id", "run") %in% names(data))) {
    row <- unique(data[intersect(c("id", "run", "session"), names(data))])
  } else row <- data.frame()
  if (nrow(row) > 1L) stop("Multiple identities found for file: ", fname)
  identity <- as.list(row[intersect(c("id", "run", "session"), names(row))])
  identity <- identity[vapply(identity, function(x) length(x) == 1L && !is.na(x) && nzchar(as.character(x)), logical(1L))]
  source <- if (".nifti" %in% names(row) && nrow(row) == 1L && !is.na(row$.nifti)) row$.nifti else basename(fname)
  missing <- setdiff(c("id", "run", "session"), names(identity))
  if (length(missing)) identity <- c(identity, parse_decon_identity(source, missing))
  if (!all(c("id", "run") %in% names(identity))) {
    stop("Cannot determine id and run for ", fname,
      "; supply file_manifest with file, id, run, and optional session columns")
  }
  identity
}

match_decon_trials <- function(identity, trial_df) {
  checkmate::assert_subset(c("id", "run"), names(trial_df))
  rows <- rep(TRUE, nrow(trial_df))
  for (column in intersect(c("id", "run", "session"), names(identity))) {
    if (!column %in% names(trial_df)) stop("Trial data lack identity column: ", column)
    values <- as.character(trial_df[[column]])
    target <- as.character(identity[[column]])
    equal <- !is.na(values) & values == target
    # Preserve exact string identities, including leading zeros. Numeric-style
    # legacy labels can match numeric trial data only when no exact match exists.
    if (!any(rows & equal) && grepl("^[0-9]+$", target)) {
      digits <- !is.na(values) & grepl("^[0-9]+$", values)
      equal <- digits & sub("^0+(?=[0-9])", "", values, perl = TRUE) ==
        sub("^0+(?=[0-9])", "", target, perl = TRUE)
    }
    rows <- rows & equal
  }
  selected <- trial_df[which(rows), , drop = FALSE]
  if (nrow(selected) == 0L) stop("No trial data match the file's id/run/session")
  keys <- intersect(c("id", "run", "session"), names(selected))
  if (nrow(unique(selected[keys])) != 1L) {
    stop("Ambiguous trial identity; supply exact id, run, and session in file_manifest")
  }
  selected
}
