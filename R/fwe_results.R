# Collected FWE results and persistent artifact manifests

fwe_result_schema_version <- function() 1L

fwe_artifact_manifest_schema_version <- function() 1L

fwe_artifact_manifest_column_classes <- function() {
  columns <- c(
    "manifest_schema_version", "correction_name", "method", "artifact_id",
    "task_id", "target_id", "target_key", "analysis_id", "artifact_role",
    "fwe_alpha", "artifact_file", "required", "exists", "artifact_status",
    "task_status", "execution_status", "job_id", "can_use_as_mask",
    "source_statistic_file", "correction_mask_file", "output_directory",
    "analysis_dir", "image_feat_dir", "collected_at", "session",
    "l1_model", "l1_cope_number", "l1_cope_name", "l2_model",
    "l2_cope_number", "l2_cope_name", "l3_model", "l3_cope_number",
    "l3_cope_name", "log_file"
  )
  classes <- stats::setNames(rep("character", length(columns)), columns)
  classes[c(
    "manifest_schema_version", "session", "l1_cope_number",
    "l2_cope_number", "l3_cope_number"
  )] <- "integer"
  classes["fwe_alpha"] <- "numeric"
  classes[c("required", "exists", "can_use_as_mask")] <- "logical"
  classes
}

fwe_result_paths <- function(output_directory) {
  manifest_directory <- file.path(output_directory, "manifest")
  list(
    directory = manifest_directory,
    manifest_file = file.path(manifest_directory, "artifacts.csv"),
    result_file = file.path(manifest_directory, "result.rds")
  )
}

fwe_atomic_write <- function(file, writer) {
  parent <- dirname(file)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
  temporary <- tempfile(
    pattern = paste0(".", basename(file), "-"),
    tmpdir = parent
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writer(temporary)
  if (!file.rename(temporary, file)) {
    copied <- file.copy(temporary, file, overwrite = TRUE)
    if (!isTRUE(copied)) {
      stop("Could not persist FWE result file: ", file, call. = FALSE)
    }
    unlink(temporary)
  }
  normalizePath(file, mustWork = TRUE)
}

fwe_collection_context <- function(x) {
  if (inherits(x, "fwe_run")) {
    return(list(
      plan = x$plan,
      execution = x$execution,
      job_ids = x$job_ids
    ))
  }
  if (inherits(x, "fwe_result")) {
    return(list(
      plan = x$plan,
      execution = x$execution,
      job_ids = x$job_ids
    ))
  }
  if (inherits(x, "fwe_plan")) {
    return(list(plan = x, execution = NULL, job_ids = character(0)))
  }
  if (checkmate::test_string(x) && dir.exists(x)) {
    return(fwe_collection_context(read_fwe_result(x)))
  }
  if (checkmate::test_string(x) && file.exists(x)) {
    value <- readRDS(x)
    return(fwe_collection_context(value))
  }
  stop(
    "x must be an fwe_plan, fwe_run, fwe_result, or saved RDS result.",
    call. = FALSE
  )
}

fwe_manifest_value <- function(source, rows, field, default = NA_character_) {
  if (is.null(source) || !field %in% names(source)) {
    return(rep(default, length(rows)))
  }
  source[[field]][rows]
}

build_fwe_artifact_manifest <- function(plan, execution = NULL,
                                        collected_at = Sys.time()) {
  plan <- refresh_fwe_plan(plan)
  artifacts <- plan$artifacts
  task_rows <- match(artifacts$task_id, plan$tasks$task_id)
  target_rows <- match(plan$tasks$target_id[task_rows], plan$targets$target_id)
  analysis_rows <- match(
    plan$tasks$analysis_id[task_rows], plan$analyses$analysis_id
  )

  alpha_id <- ifelse(
    is.na(artifacts$fwe_alpha), "none",
    vapply(artifacts$fwe_alpha, function(value) {
      if (is.na(value)) "none" else format_fwe_probability(value)
    }, character(1L))
  )
  target_id <- plan$tasks$target_id[task_rows]
  artifact_id <- paste(
    target_id, artifacts$artifact_role, alpha_id, sep = "::"
  )
  artifact_status <- ifelse(
    artifacts$exists, ifelse(plan$tasks$status[task_rows] == "complete", "available", "unverified"),
    ifelse(artifacts$required, "missing_required", "missing_optional")
  )

  manifest <- data.frame(
    manifest_schema_version = fwe_artifact_manifest_schema_version(),
    correction_name = plan$name,
    method = plan$method,
    artifact_id = artifact_id,
    task_id = artifacts$task_id,
    target_id = target_id,
    target_key = plan$tasks$target_key[task_rows],
    analysis_id = plan$tasks$analysis_id[task_rows],
    artifact_role = artifacts$artifact_role,
    fwe_alpha = artifacts$fwe_alpha,
    artifact_file = artifacts$artifact_file,
    required = artifacts$required,
    exists = artifacts$exists,
    artifact_status = artifact_status,
    task_status = plan$tasks$status[task_rows],
    can_use_as_mask = artifacts$artifact_role == "thresholded_statistic" & artifact_status == "available",
    source_statistic_file = fwe_manifest_value(
      plan$tasks, task_rows, "zstat_file"
    ),
    correction_mask_file = fwe_manifest_value(
      plan$tasks, task_rows, "correction_mask_file"
    ),
    output_directory = fwe_manifest_value(
      plan$tasks, task_rows, "output_directory"
    ),
    analysis_dir = fwe_manifest_value(
      plan$analyses, analysis_rows, "analysis_dir"
    ),
    image_feat_dir = fwe_manifest_value(
      plan$analyses, analysis_rows, "image_feat_dir"
    ),
    collected_at = format(
      collected_at, "%Y-%m-%dT%H:%M:%OS6%z", tz = "UTC"
    ),
    stringsAsFactors = FALSE
  )

  selector_fields <- c(
    "session",
    "l1_model", "l1_cope_number", "l1_cope_name",
    "l2_model", "l2_cope_number", "l2_cope_name",
    "l3_model", "l3_cope_number", "l3_cope_name"
  )
  for (field in selector_fields) {
    manifest[[field]] <- fwe_manifest_value(
      plan$targets, target_rows, field
    )
  }

  execution_rows <- if (is.null(execution)) {
    rep(NA_integer_, nrow(manifest))
  } else {
    match(manifest$task_id, execution$task_id)
  }
  manifest$execution_status <- fwe_manifest_value(
    execution, execution_rows, "execution_status"
  )
  manifest$job_id <- fwe_manifest_value(
    execution, execution_rows, "job_id"
  )
  manifest$log_file <- fwe_manifest_value(
    execution, execution_rows, "log_file"
  )

  front <- c(
    "manifest_schema_version", "correction_name", "method", "artifact_id",
    "task_id", "target_id", "target_key", "analysis_id",
    "artifact_role", "fwe_alpha", "artifact_file", "required", "exists",
    "artifact_status", "task_status", "execution_status", "job_id",
    "can_use_as_mask"
  )
  manifest <- manifest[, c(front, setdiff(names(manifest), front)), drop = FALSE]
  rownames(manifest) <- NULL
  validate_fwe_artifact_manifest(manifest)
  list(plan = plan, manifest = manifest)
}

#' Validate a persistent FWE artifact manifest
#'
#' @param x a data frame created by \code{collect_fwe_results()}.
#' @return \code{x}, invisibly, when valid.
#' @keywords internal
validate_fwe_artifact_manifest <- function(x) {
  checkmate::assert_data_frame(x, min.rows = 1L)
  required <- c(
    "manifest_schema_version", "correction_name", "method", "artifact_id",
    "task_id", "target_id", "analysis_id", "artifact_role", "fwe_alpha",
    "artifact_file", "required", "exists", "artifact_status", "task_status",
    "can_use_as_mask", "collected_at"
  )
  checkmate::assert_names(names(x), must.include = required)
  if (!all(x$manifest_schema_version ==
           fwe_artifact_manifest_schema_version())) {
    stop("Unsupported FWE artifact manifest schema version.", call. = FALSE)
  }
  if (anyDuplicated(x$artifact_id)) {
    stop("FWE artifact manifest artifact_id values must be unique.", call. = FALSE)
  }
  checkmate::assert_character(x$artifact_file, any.missing = FALSE)
  checkmate::assert_logical(x$required, any.missing = FALSE)
  checkmate::assert_logical(x$exists, any.missing = FALSE)
  checkmate::assert_logical(x$can_use_as_mask, any.missing = FALSE)
  checkmate::assert_subset(
    x$artifact_status,
    c("available", "unverified", "missing_required", "missing_optional"),
    empty.ok = FALSE
  )
  invisible(x)
}

#' Validate a collected FWE result
#'
#' @param x an \code{fwe_result} or equivalent list.
#' @return \code{x}, invisibly, when valid.
#' @keywords internal
validate_fwe_result <- function(x) {
  if (!is.list(x)) stop("An FWE result must be a list.", call. = FALSE)
  required <- c(
    "schema_version", "name", "method", "collected_at", "complete",
    "plan", "manifest", "execution", "job_ids", "manifest_file",
    "result_file"
  )
  checkmate::assert_names(names(x), identical.to = required)
  if (!identical(as.integer(x$schema_version), fwe_result_schema_version())) {
    stop("Unsupported FWE result schema version: ", x$schema_version,
         call. = FALSE)
  }
  checkmate::assert_string(x$name)
  checkmate::assert_string(x$method)
  checkmate::assert_logical(x$complete, len = 1L, any.missing = FALSE)
  validate_fwe_plan(x$plan)
  validate_fwe_artifact_manifest(x$manifest)
  if (!identical(x$name, x$plan$name) || !identical(x$method, x$plan$method)) {
    stop("FWE result identity does not match its plan.", call. = FALSE)
  }
  manifest_complete <- all(
    !x$manifest$required | (x$manifest$exists & x$manifest$task_status == "complete")
  )
  if (!identical(x$complete, manifest_complete)) {
    stop("FWE result completeness disagrees with its artifact manifest.",
         call. = FALSE)
  }
  if (!is.null(x$execution)) checkmate::assert_data_frame(x$execution)
  checkmate::assert_character(x$job_ids, any.missing = FALSE)
  checkmate::assert_string(x$manifest_file)
  checkmate::assert_string(x$result_file)
  invisible(x)
}

#' Collect FWE results and persist the artifact manifest
#'
#' Refresh expected artifacts, join them to task and semantic target metadata,
#' and optionally write an atomic CSV manifest plus an RDS result snapshot.
#' Repeated collection updates the same canonical files, making this suitable
#' for collecting scheduler results after jobs finish.
#'
#' @param x an \code{fwe_plan}, \code{fwe_run}, \code{fwe_result}, saved RDS
#'   result, or result directory.
#' @param output_directory correction output directory. Defaults to the plan's
#'   output directory.
#' @param require_complete if \code{TRUE}, fail unless every required artifact
#'   exists.
#' @param persist if \code{TRUE}, atomically write \code{manifest/artifacts.csv}
#'   and \code{manifest/result.rds}.
#'
#' @return an \code{fwe_result} containing the refreshed plan and normalized
#'   artifact manifest.
#' @export
collect_fwe_results <- function(x, output_directory = NULL,
                                require_complete = FALSE, persist = TRUE) {
  checkmate::assert_logical(require_complete, len = 1L, any.missing = FALSE)
  checkmate::assert_logical(persist, len = 1L, any.missing = FALSE)
  context <- fwe_collection_context(x)
  plan <- context$plan
  validate_fwe_plan(plan)
  if (is.null(output_directory)) output_directory <- plan$output_directory
  checkmate::assert_string(output_directory, min.chars = 1L)
  output_directory <- as.character(fs::path_abs(path.expand(output_directory)))
  paths <- fwe_result_paths(output_directory)
  collected_at <- Sys.time()
  built <- build_fwe_artifact_manifest(
    plan, execution = context$execution, collected_at = collected_at
  )
  complete <- all(built$plan$tasks$status == "complete")
  if (isTRUE(require_complete) && !complete) {
    missing <- built$manifest$artifact_id[
      built$manifest$required & built$manifest$artifact_status != "available"
    ]
    stop(
      "Required FWE artifacts are missing or unverified: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  result <- list(
    schema_version = fwe_result_schema_version(),
    name = built$plan$name,
    method = built$plan$method,
    collected_at = collected_at,
    complete = complete,
    plan = built$plan,
    manifest = built$manifest,
    execution = context$execution,
    job_ids = as.character(context$job_ids),
    manifest_file = as.character(fs::path_abs(paths$manifest_file)),
    result_file = as.character(fs::path_abs(paths$result_file))
  )
  class(result) <- c("fwe_result", "list")
  validate_fwe_result(result)

  if (isTRUE(persist)) {
    result$manifest_file <- fwe_atomic_write(
      paths$manifest_file,
      function(file) utils::write.csv(
        result$manifest, file, row.names = FALSE, na = ""
      )
    )
    result$result_file <- as.character(fs::path_abs(paths$result_file))
    result$result_file <- fwe_atomic_write(
      paths$result_file,
      function(file) saveRDS(result, file)
    )
  }
  result
}

#' Read a collected FWE result or artifact manifest
#'
#' @param path an RDS/CSV file or correction output directory.
#' @return \code{read_fwe_result()} returns an \code{fwe_result};
#'   \code{read_fwe_artifact_manifest()} returns a validated data frame.
#' @export
read_fwe_result <- function(path) {
  checkmate::assert_string(path, min.chars = 1L)
  if (dir.exists(path)) path <- fwe_result_paths(path)$result_file
  checkmate::assert_file_exists(path)
  value <- readRDS(path)
  if (!inherits(value, "fwe_result")) class(value) <- c("fwe_result", "list")
  validate_fwe_result(value)
  value
}

#' @rdname read_fwe_result
#' @export
read_fwe_artifact_manifest <- function(path) {
  checkmate::assert_string(path, min.chars = 1L)
  if (dir.exists(path)) path <- fwe_result_paths(path)$manifest_file
  checkmate::assert_file_exists(path)
  value <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    na.strings = "",
    colClasses = fwe_artifact_manifest_column_classes()
  )
  validate_fwe_artifact_manifest(value)
  value
}

as_fwe_artifact_manifest <- function(x) {
  if (is.data.frame(x)) {
    validate_fwe_artifact_manifest(x)
    return(x)
  }
  if (inherits(x, "fwe_result")) return(x$manifest)
  if (inherits(x, "fwe_run") && !is.null(x$result)) return(x$result$manifest)
  if (inherits(x, "fwe_run") || inherits(x, "fwe_plan")) {
    return(collect_fwe_results(x, persist = FALSE)$manifest)
  }
  if (checkmate::test_string(x) && dir.exists(x)) {
    return(read_fwe_artifact_manifest(x))
  }
  if (checkmate::test_string(x) && file.exists(x)) {
    if (grepl("\\.csv$", x, ignore.case = TRUE)) {
      return(read_fwe_artifact_manifest(x))
    }
    return(read_fwe_result(x)$manifest)
  }
  stop("Cannot obtain an FWE artifact manifest from x.", call. = FALSE)
}

#' Select artifacts from a collected FWE result
#'
#' @param x an FWE plan, run, result, artifact manifest, or persisted path.
#' @param artifact_role optional artifact roles to retain.
#' @param fwe_alpha optional FWE alpha values to retain.
#' @param target_id,l3_cope_name optional semantic target selectors.
#' @param require_existing if \code{TRUE}, retain only existing artifacts.
#' @return a filtered artifact-manifest data frame.
#' @export
fwe_result_artifacts <- function(
    x, artifact_role = NULL, fwe_alpha = NULL, target_id = NULL,
    l3_cope_name = NULL, require_existing = TRUE) {
  checkmate::assert_logical(require_existing, len = 1L, any.missing = FALSE)
  manifest <- as_fwe_artifact_manifest(x)
  keep <- rep(TRUE, nrow(manifest))
  if (!is.null(artifact_role)) {
    checkmate::assert_character(artifact_role, any.missing = FALSE)
    keep <- keep & manifest$artifact_role %in% artifact_role
  }
  if (!is.null(fwe_alpha)) {
    checkmate::assert_numeric(fwe_alpha, any.missing = FALSE)
    alpha_match <- vapply(manifest$fwe_alpha, function(value) {
      !is.na(value) && any(abs(value - fwe_alpha) < sqrt(.Machine$double.eps))
    }, logical(1L))
    keep <- keep & alpha_match
  }
  if (!is.null(target_id)) {
    checkmate::assert_character(target_id, any.missing = FALSE)
    keep <- keep & manifest$target_id %in% target_id
  }
  if (!is.null(l3_cope_name)) {
    checkmate::assert_character(l3_cope_name, any.missing = FALSE)
    keep <- keep & manifest$l3_cope_name %in% l3_cope_name
  }
  if (isTRUE(require_existing)) keep <- keep & manifest$exists & manifest$artifact_status == "available"
  manifest[keep, , drop = FALSE]
}

#' @keywords internal
#' @export
print.fwe_result <- function(x, ...) {
  available <- sum(x$manifest$exists)
  cat(sprintf("<fwe_result> %s\n", x$name))
  cat(sprintf("  method: %s\n", x$method))
  cat(sprintf("  complete: %s\n", x$complete))
  cat(sprintf("  artifacts: %d/%d available\n", available, nrow(x$manifest)))
  cat(sprintf("  manifest: %s\n", x$manifest_file))
  invisible(x)
}
