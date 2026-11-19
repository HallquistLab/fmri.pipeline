# Familywise-error correction specifications

fwe_spec_schema_version <- function() 1L

fwe_target_fields <- function(level = 3L) {
  if (!identical(as.integer(level), 3L)) {
    stop("FWE specifications currently support level = 3 only.", call. = FALSE)
  }

  c(
    "session",
    "l1_model", "l1_cope_number", "l1_cope_name",
    "l2_model", "l2_cope_number", "l2_cope_name",
    "l3_model", "l3_cope_number", "l3_cope_name"
  )
}

fwe_method_defaults <- function(method_name) {
  switch(method_name,
    ptfce = list(
      name = "ptfce",
      fwe_alpha = 0.05,
      two_sided = TRUE,
      write_threshold_images = TRUE
    ),
    afni_3dclustsim_permutation = list(
      name = "afni_3dclustsim_permutation",
      iter = 30000L,
      residual_njobs = 32L,
      ncpus = 8L,
      pthr = c(0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002, 0.0001),
      athr = c(0.05, 0.02, 0.01, 0.005, 0.002, 0.001),
      voxel_p = 0.001,
      cluster_alpha = 0.05,
      NN = 1L,
      sided = "bi"
    ),
    stop("Unsupported FWE method: ", method_name, call. = FALSE)
  )
}

normalize_fwe_method <- function(method, options = list()) {
  if (checkmate::test_string(method)) {
    method <- list(name = method)
  }
  if (!is.list(method) || is.null(names(method)) || !"name" %in% names(method)) {
    stop("method must be a method name or a named list containing `name`.", call. = FALSE)
  }
  checkmate::assert_string(method$name)

  aliases <- c(
    "3dclustsim_permutation" = "afni_3dclustsim_permutation",
    "clustsim_permutation" = "afni_3dclustsim_permutation"
  )
  method_name <- tolower(method$name)
  if (method_name %in% names(aliases)) method_name <- unname(aliases[[method_name]])
  method$name <- method_name

  if (length(options) > 0L) {
    if (is.null(names(options)) || any(!nzchar(names(options)))) {
      stop("All method options supplied through `...` must be named.", call. = FALSE)
    }
    method[names(options)] <- options
  }

  if (identical(method_name, "ptfce") && "min_cluster_voxels" %in% names(method)) {
    warning("min_cluster_voxels never filtered pTFCE inference outputs and is deprecated. ",
      "Use ptfce_spec$get_clusters(clust_nvox = ...) for separate cluster reporting.", call. = FALSE)
    method$min_cluster_voxels <- NULL
  }

  defaults <- fwe_method_defaults(method_name)
  unknown <- setdiff(names(method), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      "Unknown option(s) for FWE method ", method_name, ": ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  method <- utils::modifyList(defaults, method)
  validate_fwe_method(method)
  method
}

validate_fwe_method <- function(method) {
  probability_vector <- function(x, field) {
    if (!checkmate::test_numeric(
      x, lower = 1e-10, upper = 1 - 1e-10,
      any.missing = FALSE, unique = TRUE, min.len = 1L
    )) {
      stop("method$", field, " must contain unique probabilities strictly between 0 and 1.", call. = FALSE)
    }
  }

  if (identical(method$name, "ptfce")) {
    probability_vector(method$fwe_alpha, "fwe_alpha")
    checkmate::assert_logical(method$two_sided, len = 1L)
    checkmate::assert_logical(method$write_threshold_images, len = 1L)
  } else if (identical(method$name, "afni_3dclustsim_permutation")) {
    checkmate::assert_integerish(method$iter, len = 1L, lower = 10L)
    checkmate::assert_integerish(method$residual_njobs, len = 1L, lower = 1L)
    checkmate::assert_integerish(method$ncpus, len = 1L, lower = 1L)
    probability_vector(method$pthr, "pthr")
    probability_vector(method$athr, "athr")
    checkmate::assert_number(method$voxel_p, lower = 1e-10, upper = 1 - 1e-10)
    checkmate::assert_number(method$cluster_alpha, lower = 1e-10, upper = 1 - 1e-10)
    checkmate::assert_integerish(method$NN, len = 1L, lower = 1L, upper = 3L)
    checkmate::assert_string(method$sided)
    checkmate::assert_subset(method$sided, c("1", "2", "bi"), empty.ok = FALSE)

    if (!any(abs(method$pthr - method$voxel_p) < 1e-12)) {
      stop("method$voxel_p must be included in method$pthr.", call. = FALSE)
    }
    if (!any(abs(method$athr - method$cluster_alpha) < 1e-12)) {
      stop("method$cluster_alpha must be included in method$athr.", call. = FALSE)
    }
  } else {
    stop("Unsupported FWE method: ", method$name, call. = FALSE)
  }
  invisible(method)
}

normalize_fwe_targets <- function(targets, level = 3L) {
  if (!is.list(targets) || length(targets) == 0L || is.null(names(targets))) {
    stop("targets must be a non-empty named list of model/contrast selectors.", call. = FALSE)
  }
  if (any(!nzchar(names(targets))) || anyDuplicated(names(targets))) {
    stop("Every target selector must have a unique, non-empty name.", call. = FALSE)
  }

  allowed <- fwe_target_fields(level)
  unknown <- setdiff(names(targets), allowed)
  if (length(unknown) > 0L) {
    stop("Unknown FWE target selector(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  number_fields <- c("session", "l1_cope_number", "l2_cope_number", "l3_cope_number")
  for (field in names(targets)) {
    value <- targets[[field]]
    if (length(value) == 0L || is.list(value) || anyNA(value)) {
      stop("Target selector `", field, "` must contain one or more non-missing values.", call. = FALSE)
    }
    if (field %in% number_fields) {
      checkmate::assert_integerish(value, any.missing = FALSE)
      value <- as.integer(value)
    } else {
      checkmate::assert_character(value, any.missing = FALSE, min.chars = 1L)
      value <- as.character(value)
    }
    targets[[field]] <- unique(value)
  }
  targets
}

normalize_fwe_correction_mask <- function(correction_mask) {
  if (checkmate::test_string(correction_mask)) {
    if (identical(correction_mask, "model")) {
      correction_mask <- list(source = "model")
    } else {
      correction_mask <- list(source = "file", path = correction_mask)
    }
  }
  if (!is.list(correction_mask) || is.null(names(correction_mask))) {
    stop("correction_mask must be `model`, a file path, or a named list.", call. = FALSE)
  }

  allowed <- c("source", "path")
  unknown <- setdiff(names(correction_mask), allowed)
  if (length(unknown) > 0L) {
    stop("Unknown correction_mask field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  checkmate::assert_string(correction_mask$source)
  checkmate::assert_subset(correction_mask$source, c("model", "file"), empty.ok = FALSE)

  if (identical(correction_mask$source, "file")) {
    checkmate::assert_string(correction_mask$path, min.chars = 1L)
  } else {
    correction_mask$path <- NULL
  }
  correction_mask
}

normalize_fwe_compute <- function(compute) {
  if (is.null(compute)) compute <- list()
  if (!is.list(compute) || (length(compute) > 0L && is.null(names(compute)))) {
    stop("compute must be a named list.", call. = FALSE)
  }
  defaults <- list(scheduler = "inherit", memgb_per_command = 1,
    time_per_command = "00:00:45", commands_per_cpu = 1L, cpus_per_job = 1L)
  unknown <- setdiff(names(compute), names(defaults))
  if (length(unknown) > 0L) {
    stop("Unknown compute field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  compute <- utils::modifyList(defaults, compute)
  checkmate::assert_string(compute$scheduler)
  checkmate::assert_subset(
    compute$scheduler,
    c("inherit", "local", "sh", "slurm", "sbatch", "torque", "qsub"),
    empty.ok = FALSE
  )
  checkmate::assert_number(compute$memgb_per_command, lower = 0.1, finite = TRUE)
  checkmate::assert_string(compute$time_per_command)
  validate_dhms(compute$time_per_command)
  checkmate::assert_integerish(compute$commands_per_cpu, lower = 1, len = 1L)
  checkmate::assert_integerish(compute$cpus_per_job, lower = 1, len = 1L)
  compute
}

#' Create a familywise-error correction specification
#'
#' Create a declarative, serializable specification describing which completed
#' group-analysis contrasts should receive a familywise-error correction and
#' which correction method should be used. The specification contains no live
#' execution objects and can therefore be stored safely in a pipeline object or
#' YAML file.
#'
#' @param name unique name for this correction.
#' @param targets non-empty named list selecting level-3 model and contrast
#'   fields. Valid fields are \code{session}, \code{l1_model},
#'   \code{l1_cope_number}, \code{l1_cope_name}, \code{l2_model},
#'   \code{l2_cope_number}, \code{l2_cope_name}, \code{l3_model},
#'   \code{l3_cope_number}, and \code{l3_cope_name}. Values may be vectors.
#' @param method correction method name or named method configuration. Currently
#'   supported specifications are \code{"ptfce"} and
#'   \code{"afni_3dclustsim_permutation"}. Only pTFCE currently has a planning
#'   and execution adapter; the AFNI specification can be saved for future use.
#' @param level GLM level to correct. Only level 3 is currently supported.
#' @param correction_mask \code{"model"} to use each analysis mask, a NIfTI
#'   path, or a list with \code{source} equal to \code{"model"} or
#'   \code{"file"} (and \code{path} for the latter).
#' @param compute named list of compute settings. By default the scheduler is
#'   inherited from the pipeline.
#' @param schema_version specification schema version. Normally left at its
#'   default.
#' @param ... named method-specific options overriding method defaults.
#'
#' @return an object of class \code{fwe_spec}.
#' @export
fwe_spec <- function(name, targets, method = "ptfce", level = 3L,
                     correction_mask = "model",
                     compute = list(scheduler = "inherit"),
                     schema_version = fwe_spec_schema_version(), ...) {
  spec <- list(
    schema_version = as.integer(schema_version),
    name = name,
    level = as.integer(level),
    targets = normalize_fwe_targets(targets, level = level),
    method = normalize_fwe_method(method, options = list(...)),
    correction_mask = normalize_fwe_correction_mask(correction_mask),
    compute = normalize_fwe_compute(compute)
  )
  class(spec) <- c("fwe_spec", "list")
  validate_fwe_spec(spec)
  spec
}

#' Validate an FWE correction specification
#'
#' @param x an \code{fwe_spec} or an equivalent plain list.
#' @return \code{x}, invisibly, if validation succeeds. Otherwise an error is
#'   raised with the invalid field.
#' @export
validate_fwe_spec <- function(x) {
  if (!is.list(x)) stop("An FWE specification must be a list.", call. = FALSE)
  required <- c(
    "schema_version", "name", "level", "targets", "method",
    "correction_mask", "compute"
  )
  missing <- setdiff(required, names(x))
  unknown <- setdiff(names(x), required)
  if (length(missing) > 0L) {
    stop("Missing FWE specification field(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(unknown) > 0L) {
    stop("Unknown FWE specification field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  checkmate::assert_integerish(x$schema_version, len = 1L, lower = 1L)
  if (!identical(as.integer(x$schema_version), fwe_spec_schema_version())) {
    stop(
      "Unsupported FWE specification schema version: ", x$schema_version,
      ". This package supports version ", fwe_spec_schema_version(), ".",
      call. = FALSE
    )
  }
  checkmate::assert_string(x$name, min.chars = 1L, pattern = "\\S")
  if (!identical(x$name, trimws(x$name)) || x$name %in% c(".", "..") || grepl("[/\\\\]", x$name)) {
    stop("FWE specification names cannot have surrounding whitespace or contain path separators.", call. = FALSE)
  }
  checkmate::assert_integerish(x$level, len = 1L, lower = 1L, upper = 3L)
  if (!identical(as.integer(x$level), 3L)) {
    stop("FWE specifications currently support level = 3 only.", call. = FALSE)
  }

  normalize_fwe_targets(x$targets, level = x$level)
  normalize_fwe_method(x$method)
  normalize_fwe_correction_mask(x$correction_mask)
  normalize_fwe_compute(x$compute)
  invisible(x)
}

as_fwe_spec <- function(x) {
  validate_fwe_spec(x)
  fwe_spec(
    name = x$name,
    targets = x$targets,
    method = x$method,
    level = x$level,
    correction_mask = x$correction_mask,
    compute = x$compute,
    schema_version = x$schema_version
  )
}

#' @keywords internal
#' @export
print.fwe_spec <- function(x, ...) {
  cat(sprintf("<fwe_spec> %s\n", x$name))
  cat(sprintf("  level: %d\n", x$level))
  cat(sprintf("  method: %s\n", x$method$name))
  cat("  targets:\n")
  for (field in names(x$targets)) {
    cat(sprintf("    %s: %s\n", field, paste(x$targets[[field]], collapse = ", ")))
  }
  invisible(x)
}

format_fwe_selectors <- function(targets) {
  paste(
    vapply(names(targets), function(field) {
      sprintf("%s={%s}", field, paste(targets[[field]], collapse = ", "))
    }, character(1L)),
    collapse = "; "
  )
}

fwe_target_key <- function(targets) {
  key_fields <- intersect(
    c(
      "session",
      "l1_model", "l1_cope_number", "l1_cope_name",
      "l2_model", "l2_cope_number", "l2_cope_name",
      "l3_model", "l3_cope_number", "l3_cope_name"
    ),
    names(targets)
  )
  apply(targets[, key_fields, drop = FALSE], 1L, function(row) {
    values <- ifelse(is.na(row), "<NA>", as.character(row))
    paste(paste0(key_fields, "=", utils::URLencode(values, reserved = TRUE)), collapse = "|")
  })
}

fwe_target_id <- function(targets) {
  id_fields <- intersect(
    c("session", "l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model", "l3_cope_name"),
    names(targets)
  )
  ids <- apply(targets[, id_fields, drop = FALSE], 1L, function(row) {
    values <- ifelse(is.na(row), "none", as.character(row))
    fs::path_sanitize(paste(paste0(id_fields, "-", values), collapse = "__"), replacement = "_")
  })
  make.unique(ids, sep = "__")
}

#' Resolve an FWE specification to concrete group-analysis targets
#'
#' Resolve semantic model and contrast selectors against the canonical FSL
#' output inventory returned by \code{lookup_feat_outputs()}. Resolution never
#' constructs FEAT paths independently.
#'
#' @param gpa a \code{glm_pipeline_arguments} with level-3 FSL setup metadata.
#' @param spec an \code{fwe_spec} or equivalent list.
#' @param require_existing if \code{TRUE}, fail when a selected statistic image
#'   is absent. If \code{FALSE}, unresolved output paths are retained for
#'   planning and preview.
#' @param require_complete if \code{TRUE}, fail when a selected group analysis
#'   is not complete.
#' @param source,cache_dir,refresh_status passed to \code{lookup_feat_outputs()}.
#' @param lg optional \code{lgr::Logger}.
#'
#' @return a data.frame of class \code{fwe_target_set}, with one row per selected
#'   group contrast and stable \code{target_key} and \code{target_id} columns.
#' @export
resolve_fwe_targets <- function(gpa, spec, require_existing = TRUE,
                                require_complete = TRUE,
                                source = c("auto", "setup", "cache", "filesystem"),
                                cache_dir = NULL, refresh_status = FALSE,
                                lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  spec <- as_fwe_spec(spec)
  checkmate::assert_logical(require_existing, len = 1L)
  checkmate::assert_logical(require_complete, len = 1L)
  source <- match.arg(source)

  inventory <- lookup_feat_outputs(
    gpa = gpa,
    level = spec$level,
    what = "zstat",
    include_missing = TRUE,
    include_internal = TRUE,
    source = source,
    cache_dir = cache_dir,
    refresh_status = refresh_status,
    lg = lg,
    format = "by_level"
  )$l3

  if (!is.data.frame(inventory) || nrow(inventory) == 0L) {
    stop("No level-3 FSL statistic outputs are available for FWE target resolution.", call. = FALSE)
  }

  keep <- rep(TRUE, nrow(inventory))
  for (field in names(spec$targets)) {
    if (!field %in% names(inventory)) {
      stop("Target field `", field, "` is unavailable in the level-3 output inventory.", call. = FALSE)
    }
    keep <- keep & !is.na(inventory[[field]]) & inventory[[field]] %in% spec$targets[[field]]
  }
  targets <- inventory[keep, , drop = FALSE]

  if (nrow(targets) == 0L) {
    stop(
      "No level-3 outputs match FWE selectors: ",
      format_fwe_selectors(spec$targets),
      call. = FALSE
    )
  }

  incomplete <- !targets$analysis_status %in% "complete"
  if (isTRUE(require_complete) && any(incomplete)) {
    stop(
      sum(incomplete), " selected FWE target(s) are not complete. ",
      "Use require_complete = FALSE to preview them.",
      call. = FALSE
    )
  }
  missing_images <- is.na(targets$image_file) | !targets$image_exists
  if (isTRUE(require_existing) && any(missing_images)) {
    stop(
      sum(missing_images), " selected FWE target image(s) do not exist. ",
      "Use require_existing = FALSE to preview expected paths.",
      call. = FALSE
    )
  }

  targets$correction_name <- spec$name
  targets$target_key <- fwe_target_key(targets)
  targets$target_id <- fwe_target_id(targets)
  targets$target_ready <- !incomplete & !missing_images

  front <- c("correction_name", "target_id", "target_key", "target_ready")
  targets <- targets[, c(front, setdiff(names(targets), front)), drop = FALSE]
  rownames(targets) <- NULL
  structure(
    targets,
    class = c("fwe_target_set", "data.frame"),
    fwe_spec = spec
  )
}

#' Write or read an FWE correction specification as YAML
#'
#' @param x an \code{fwe_spec} or equivalent list.
#' @param file YAML file path.
#' @param overwrite whether an existing file may be replaced.
#'
#' @return \code{write_fwe_spec()} invisibly returns the normalized output path;
#'   \code{read_fwe_spec()} returns an \code{fwe_spec}.
#' @export
write_fwe_spec <- function(x, file, overwrite = FALSE) {
  x <- as_fwe_spec(x)
  checkmate::assert_string(file, min.chars = 1L)
  checkmate::assert_logical(overwrite, len = 1L)
  if (!grepl("\\.ya?ml$", file, ignore.case = TRUE)) {
    stop("FWE specification files must use a .yaml or .yml extension.", call. = FALSE)
  }
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite existing FWE specification: ", file, call. = FALSE)
  }
  parent <- dirname(file)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)

  yaml::write_yaml(unclass(x), file = file)
  invisible(normalizePath(file, mustWork = FALSE))
}

#' @rdname write_fwe_spec
#' @export
read_fwe_spec <- function(file) {
  checkmate::assert_file_exists(file)
  if (!grepl("\\.ya?ml$", file, ignore.case = TRUE)) {
    stop("FWE specification files must use a .yaml or .yml extension.", call. = FALSE)
  }
  value <- yaml::read_yaml(file)
  as_fwe_spec(value)
}
