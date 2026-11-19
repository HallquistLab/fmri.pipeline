# Compiled familywise-error correction plans

fwe_plan_schema_version <- function() 1L

as_fwe_method_adapter <- function(method) {
  method <- normalize_fwe_method(method)
  class(method) <- c(
    paste0("fwe_method_", method$name),
    "fwe_method",
    "list"
  )
  method
}

compile_fwe_method <- function(method, targets, gpa, output_directory,
                               require_inputs = TRUE) {
  method <- as_fwe_method_adapter(method)
  UseMethod("compile_fwe_method", method)
}

#' @keywords internal
#' @export
compile_fwe_method.default <- function(method, targets, gpa, output_directory,
                                       require_inputs = TRUE) {
  stop(
    "No FWE planning adapter is implemented for method: ", method$name,
    call. = FALSE
  )
}

resolve_fwe_scheduler <- function(spec, gpa) {
  scheduler <- spec$compute$scheduler
  if (identical(scheduler, "inherit")) {
    scheduler <- gpa$scheduler
    if (is.null(scheduler)) scheduler <- "local"
  }
  scheduler <- tolower(scheduler)
  switch(scheduler,
    sbatch = "slurm",
    qsub = "torque",
    sh = "local",
    scheduler
  )
}

fwe_analysis_key <- function(targets) {
  fields <- intersect(
    c(
      "session", "l1_model", "l1_cope_number", "l1_cope_name",
      "l2_model", "l2_cope_number", "l2_cope_name", "l3_model",
      "analysis_dir", "image_feat_dir"
    ),
    names(targets)
  )
  apply(targets[, fields, drop = FALSE], 1L, function(row) {
    values <- ifelse(is.na(row), "<NA>", as.character(row))
    paste(paste0(fields, "=", utils::URLencode(values, reserved = TRUE)), collapse = "|")
  })
}

fwe_analysis_id <- function(analyses) {
  fields <- intersect(
    c("session", "l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model"),
    names(analyses)
  )
  ids <- apply(analyses[, fields, drop = FALSE], 1L, function(row) {
    values <- ifelse(is.na(row), "none", as.character(row))
    fs::path_sanitize(paste(paste0(fields, "-", values), collapse = "__"), replacement = "_")
  })
  make.unique(ids, sep = "__")
}

fwe_nifti_parts <- function(path) {
  filename <- basename(path)
  extension <- if (grepl("\\.nii\\.gz$", filename, ignore.case = TRUE)) {
    ".nii.gz"
  } else if (grepl("\\.nii$", filename, ignore.case = TRUE)) {
    ".nii"
  } else {
    stop("Expected a NIfTI filename ending in .nii or .nii.gz: ", path, call. = FALSE)
  }
  list(
    stem = substr(filename, 1L, nchar(filename) - nchar(extension)),
    extension = extension
  )
}

format_fwe_probability <- function(x) {
  format(x, digits = 15L, scientific = FALSE, trim = TRUE)
}

ptfce_artifact_table <- function(task_id, zstat_file, task_output_directory,
                                 fwe_alpha, write_threshold_images) {
  parts <- fwe_nifti_parts(zstat_file)
  base <- file.path(task_output_directory, paste0(parts$stem, "_ptfce"))
  artifacts <- data.frame(
    task_id = rep(task_id, 2L),
    artifact_role = c("ptfce_statistic", "threshold_table"),
    fwe_alpha = c(NA_real_, NA_real_),
    artifact_file = c(
      paste0(base, parts$extension),
      paste0(base, "_zthresh.csv")
    ),
    required = TRUE,
    stringsAsFactors = FALSE
  )

  if (isTRUE(write_threshold_images)) {
    thresholded <- data.frame(
      task_id = rep(task_id, length(fwe_alpha)),
      artifact_role = rep("thresholded_statistic", length(fwe_alpha)),
      fwe_alpha = as.numeric(fwe_alpha),
      artifact_file = vapply(fwe_alpha, function(alpha) {
        paste0(
          base, "_fwep_", format_fwe_probability(alpha), parts$extension
        )
      }, character(1L)),
      required = TRUE,
      stringsAsFactors = FALSE
    )
    artifacts <- dplyr::bind_rows(artifacts, thresholded)
  }
  artifacts$exists <- file.exists(artifacts$artifact_file)
  artifacts
}

#' @keywords internal
#' @export
compile_fwe_method.fwe_method_ptfce <- function(
    method, targets, gpa, output_directory, require_inputs = TRUE) {
  spec <- attr(targets, "fwe_spec")
  if (is.null(spec)) {
    stop("Resolved FWE targets are missing their fwe_spec attribute.", call. = FALSE)
  }
  targets <- as.data.frame(targets)
  targets <- targets[order(targets$target_key), , drop = FALSE]
  rownames(targets) <- NULL
  targets$analysis_key <- fwe_analysis_key(targets)

  analysis_fields <- intersect(
    c(
      "analysis_key", "session",
      "l1_model", "l1_cope_number", "l1_cope_name",
      "l2_model", "l2_cope_number", "l2_cope_name",
      "l3_model", "l3_input_mode", "analysis_dir", "image_feat_dir",
      "image_stats_dir", "analysis_status", "lookup_source"
    ),
    names(targets)
  )
  analyses <- unique(targets[, analysis_fields, drop = FALSE])
  analyses <- analyses[order(analyses$analysis_key), , drop = FALSE]
  rownames(analyses) <- NULL
  analyses$analysis_id <- fwe_analysis_id(analyses)
  analyses$model_mask_file <- file.path(analyses$image_feat_dir, "mask.nii.gz")
  analyses$smoothness_file <- file.path(analyses$image_stats_dir, "smoothness")
  analyses$model_mask_exists <- file.exists(analyses$model_mask_file)
  analyses$smoothness_exists <- file.exists(analyses$smoothness_file)

  analysis_map <- stats::setNames(analyses$analysis_id, analyses$analysis_key)
  tasks <- data.frame(
    task_id = sprintf("ptfce-%03d", seq_len(nrow(targets))),
    target_id = targets$target_id,
    target_key = targets$target_key,
    analysis_id = unname(analysis_map[targets$analysis_key]),
    method = "ptfce",
    l3_cope_number = targets$l3_cope_number,
    l3_cope_name = targets$l3_cope_name,
    zstat_file = targets$image_file,
    stringsAsFactors = FALSE
  )

  analysis_rows <- match(tasks$analysis_id, analyses$analysis_id)
  if (identical(spec$correction_mask$source, "file")) {
    correction_mask <- fs::path_abs(path.expand(spec$correction_mask$path))
    tasks$correction_mask_file <- rep(correction_mask, nrow(tasks))
  } else {
    tasks$correction_mask_file <- analyses$model_mask_file[analysis_rows]
  }
  tasks$smoothness_file <- analyses$smoothness_file[analysis_rows]
  tasks$model_mask_file <- analyses$model_mask_file[analysis_rows]
  tasks$two_sided <- method$two_sided
  tasks$fwe_alpha <- rep(list(as.numeric(method$fwe_alpha)), nrow(tasks))
  tasks$write_threshold_images <- method$write_threshold_images
  tasks$output_directory <- file.path(output_directory, "targets", tasks$target_id)

  artifacts <- dplyr::bind_rows(lapply(seq_len(nrow(tasks)), function(ii) {
    ptfce_artifact_table(
      task_id = tasks$task_id[ii],
      zstat_file = tasks$zstat_file[ii],
      task_output_directory = tasks$output_directory[ii],
      fwe_alpha = method$fwe_alpha,
      write_threshold_images = method$write_threshold_images
    )
  }))

  tasks$missing_inputs <- ptfce_missing_inputs(tasks)
  tasks$status <- ifelse(nzchar(tasks$missing_inputs), "blocked", "ready")
  tasks$input_fingerprint <- fwe_task_fingerprints(tasks)
  tasks$status[fwe_complete_tasks(tasks, artifacts)] <- "complete"

  if (isTRUE(require_inputs) && any(tasks$status == "blocked")) {
    blocked <- tasks[tasks$status == "blocked", c("task_id", "l3_cope_name", "missing_inputs"), drop = FALSE]
    details <- apply(blocked, 1L, function(row) {
      sprintf("%s (%s): %s", row[["task_id"]], row[["l3_cope_name"]], row[["missing_inputs"]])
    })
    stop(
      "Missing required pTFCE input(s): ", paste(details, collapse = "; "),
      ". Use require_inputs = FALSE to inspect a blocked plan.",
      call. = FALSE
    )
  }

  front <- c("analysis_id", "analysis_key")
  analyses <- analyses[, c(front, setdiff(names(analyses), front)), drop = FALSE]
  list(analyses = analyses, tasks = tasks, artifacts = artifacts)
}

#' Compile an FWE specification into an executable plan
#'
#' Resolve an FWE specification against completed GLM outputs and compile it
#' through a method adapter. Planning is side-effect free: output directories
#' and result files are described but not created. The pTFCE adapter creates one
#' task per selected z-statistic while retaining shared GFEAT-level inputs in an
#' analysis table.
#'
#' @param gpa a \code{glm_pipeline_arguments} containing group-analysis output
#'   metadata.
#' @param spec an \code{fwe_spec} or equivalent list.
#' @param output_directory root directory for this correction. The default is
#'   \code{<gpa$output_directory>/fwe/<spec$name>}.
#' @param require_inputs if \code{TRUE}, error when a required statistic, mask,
#'   or FSL smoothness file is missing. If \code{FALSE}, return a plan with
#'   blocked tasks for inspection.
#' @param source,cache_dir,refresh_status passed to
#'   \code{resolve_fwe_targets()}.
#' @param lg optional \code{lgr::Logger}.
#'
#' @return an object of class \code{fwe_plan} containing \code{analyses},
#'   \code{tasks}, and expected \code{artifacts} tables.
#' @export
plan_fwe_correction <- function(
    gpa, spec, output_directory = NULL, require_inputs = TRUE,
    source = c("auto", "setup", "cache", "filesystem"),
    cache_dir = NULL, refresh_status = FALSE, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  spec <- as_fwe_spec(spec)
  assert_fwe_executable(spec)
  checkmate::assert_logical(require_inputs, len = 1L, any.missing = FALSE)
  source <- match.arg(source)

  if (is.null(output_directory)) {
    checkmate::assert_string(gpa$output_directory, min.chars = 1L)
    output_directory <- file.path(gpa$output_directory, "fwe", spec$name)
  }
  checkmate::assert_string(output_directory, min.chars = 1L)
  output_directory <- fs::path_abs(path.expand(output_directory))

  targets <- resolve_fwe_targets(
    gpa = gpa,
    spec = spec,
    require_existing = require_inputs,
    require_complete = require_inputs,
    source = source,
    cache_dir = cache_dir,
    refresh_status = refresh_status,
    lg = lg
  )

  compiled <- compile_fwe_method(
    method = spec$method,
    targets = targets,
    gpa = gpa,
    output_directory = output_directory,
    require_inputs = require_inputs
  )

  plan <- list(
    schema_version = fwe_plan_schema_version(),
    name = spec$name,
    method = spec$method$name,
    output_directory = output_directory,
    scheduler = resolve_fwe_scheduler(spec, gpa),
    spec = spec,
    targets = targets,
    analyses = compiled$analyses,
    tasks = compiled$tasks,
    artifacts = compiled$artifacts
  )
  class(plan) <- c("fwe_plan", "list")
  validate_fwe_plan(plan, require_inputs = require_inputs)
  plan
}

#' Validate a compiled FWE plan
#'
#' @param x an \code{fwe_plan} or equivalent list.
#' @param require_inputs if \code{TRUE}, require every task to have all inputs.
#' @return \code{x}, invisibly, when valid.
#' @keywords internal
validate_fwe_plan <- function(x, require_inputs = FALSE) {
  if (!is.list(x)) stop("An FWE plan must be a list.", call. = FALSE)
  required <- c(
    "schema_version", "name", "method", "output_directory", "scheduler",
    "spec", "targets", "analyses", "tasks", "artifacts"
  )
  missing <- setdiff(required, names(x))
  unknown <- setdiff(names(x), required)
  if (length(missing) > 0L) {
    stop("Missing FWE plan field(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(unknown) > 0L) {
    stop("Unknown FWE plan field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  checkmate::assert_integerish(x$schema_version, len = 1L)
  if (!identical(as.integer(x$schema_version), fwe_plan_schema_version())) {
    stop("Unsupported FWE plan schema version: ", x$schema_version, call. = FALSE)
  }
  checkmate::assert_string(x$name)
  checkmate::assert_string(x$method)
  checkmate::assert_string(x$output_directory)
  checkmate::assert_string(x$scheduler)
  validate_fwe_spec(x$spec)
  if (!identical(x$name, x$spec$name)) {
    stop("FWE plan name does not match its specification.", call. = FALSE)
  }
  if (!identical(x$method, x$spec$method$name)) {
    stop("FWE plan method does not match its specification.", call. = FALSE)
  }
  if (!x$scheduler %in% c("local", "slurm", "torque")) {
    stop("Unsupported FWE plan scheduler: ", x$scheduler, call. = FALSE)
  }
  checkmate::assert_data_frame(x$targets, min.rows = 1L)
  checkmate::assert_data_frame(x$analyses, min.rows = 1L)
  checkmate::assert_data_frame(x$tasks, min.rows = 1L)
  checkmate::assert_data_frame(x$artifacts, min.rows = 1L)
  checkmate::assert_names(names(x$targets), must.include = "target_id")
  checkmate::assert_names(names(x$analyses), must.include = "analysis_id")
  checkmate::assert_names(
    names(x$tasks),
    must.include = c("task_id", "analysis_id", "target_id", "method", "status")
  )
  checkmate::assert_names(
    names(x$artifacts),
    must.include = c("task_id", "artifact_role", "artifact_file", "required")
  )

  if (anyDuplicated(x$tasks$task_id)) stop("FWE plan task_id values must be unique.", call. = FALSE)
  if (anyDuplicated(x$analyses$analysis_id)) stop("FWE plan analysis_id values must be unique.", call. = FALSE)
  if (!all(x$tasks$analysis_id %in% x$analyses$analysis_id)) {
    stop("At least one FWE task references an unknown analysis_id.", call. = FALSE)
  }
  if (!all(x$tasks$target_id %in% x$targets$target_id)) {
    stop("At least one FWE task references an unknown target_id.", call. = FALSE)
  }
  if (!all(x$tasks$method == x$method)) {
    stop("At least one FWE task uses a method that differs from its plan.", call. = FALSE)
  }
  if (!all(x$tasks$status %in% c("blocked", "ready", "complete"))) {
    stop("FWE plan tasks contain an unsupported status.", call. = FALSE)
  }
  if (!all(x$artifacts$task_id %in% x$tasks$task_id)) {
    stop("At least one FWE artifact references an unknown task_id.", call. = FALSE)
  }
  if (anyDuplicated(x$artifacts$artifact_file)) {
    stop("FWE plan artifact_file values must be unique.", call. = FALSE)
  }
  if (!is.logical(x$artifacts$required) || anyNA(x$artifacts$required)) {
    stop("FWE plan artifact required flags must be non-missing logical values.", call. = FALSE)
  }
  if (isTRUE(require_inputs) && any(x$tasks$status == "blocked")) {
    stop("FWE plan contains blocked tasks with missing inputs.", call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
#' @export
print.fwe_plan <- function(x, ...) {
  statuses <- table(x$tasks$status)
  status_text <- paste(paste0(names(statuses), "=", as.integer(statuses)), collapse = ", ")
  cat(sprintf("<fwe_plan> %s\n", x$name))
  cat(sprintf("  method: %s\n", x$method))
  cat(sprintf("  analyses: %d\n", nrow(x$analyses)))
  cat(sprintf("  tasks: %d (%s)\n", nrow(x$tasks), status_text))
  cat(sprintf("  output: %s\n", x$output_directory))
  invisible(x)
}

#' Write or read a compiled FWE plan
#'
#' @param x an \code{fwe_plan}.
#' @param file an \code{.rds} file path.
#' @param overwrite whether an existing file may be replaced.
#' @return \code{write_fwe_plan()} invisibly returns the normalized file path;
#'   \code{read_fwe_plan()} returns a validated \code{fwe_plan}.
#' @export
write_fwe_plan <- function(x, file, overwrite = FALSE) {
  validate_fwe_plan(x)
  checkmate::assert_string(file, min.chars = 1L)
  checkmate::assert_logical(overwrite, len = 1L, any.missing = FALSE)
  if (!grepl("\\.rds$", file, ignore.case = TRUE)) {
    stop("Compiled FWE plans must use an .rds extension.", call. = FALSE)
  }
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite existing FWE plan: ", file, call. = FALSE)
  }
  parent <- dirname(file)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
  saveRDS(x, file)
  invisible(normalizePath(file, mustWork = FALSE))
}

#' @rdname write_fwe_plan
#' @export
read_fwe_plan <- function(file) {
  checkmate::assert_file_exists(file)
  value <- readRDS(file)
  validate_fwe_plan(value)
  if (!inherits(value, "fwe_plan")) class(value) <- c("fwe_plan", "list")
  value
}
