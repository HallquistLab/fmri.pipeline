# GLM backend specification and resolution helpers

# Canonical empty data.frame constructors for the requirement and preflight schemas.
# Using constructors keeps column names/types defined in one place.
empty_requirement_df <- function() {
  data.frame(
    artifact = character(0),
    producer_backend = character(0),
    producer_level = integer(0),
    stringsAsFactors = FALSE
  )
}

empty_model_requirement_df <- function() {
  data.frame(
    model_name = character(0),
    execution_backend = character(0),
    artifact = character(0),
    producer_backend = character(0),
    producer_level = integer(0),
    stringsAsFactors = FALSE
  )
}

glm_backend_not_implemented <- function(backend, feature) {
  force(backend)
  force(feature)
  fn <- function(...) {
    stop(sprintf("GLM backend '%s' does not implement '%s' yet.", backend, feature), call. = FALSE)
  }
  attr(fn, "glm_backend_not_implemented") <- TRUE
  fn
}

resolve_backend_fn <- function(fn_name, pkg = "fmri.pipeline") {
  if (is.null(fn_name)) return(NULL)
  checkmate::assert_string(fn_name)
  checkmate::assert_string(pkg)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required namespace '%s' not available for backend resolution.", pkg), call. = FALSE)
  }
  get(fn_name, envir = asNamespace(pkg))
}

default_glm_backend_specs <- function() {
  list(
    fsl = list(
      name = "fsl",
      runs_l1 = TRUE,
      runs_l2 = TRUE,
      runs_l3 = TRUE,
      multi_run_strategy = "explicit_l2",
      produced_artifacts = c("run_level_contrasts", "subject_session_contrasts", "group_level_stats"),
      l3_required_artifacts = c("subject_session_contrasts"),
      l3_input_provider_backends = c("fsl"),
      l1_setup = "fsl_l1_model",
      l2_setup = "fsl_l2_model",
      l3_setup = "fsl_l3_model",
      l1_status = "get_feat_status",
      l2_status = "get_feat_status",
      l3_status = "get_feat_status",
      l1_status_inputs = c("feat_dir", "feat_fsf"),
      l2_status_inputs = c("feat_dir", "feat_fsf"),
      l3_status_inputs = c("feat_dir", "feat_fsf"),
      output_dir = "get_output_directory",
      l1_run = "run_feat_sepjobs",
      l2_run = "run_feat_sepjobs",
      l3_run = "run_feat_sepjobs"
    ),
    spm = list(
      name = "spm",
      runs_l1 = TRUE,
      runs_l2 = FALSE,
      runs_l3 = TRUE,
      multi_run_strategy = "concat_in_l1",
      produced_artifacts = c("run_level_contrasts", "subject_session_contrasts", "group_level_stats"),
      l3_required_artifacts = c("subject_session_contrasts"),
      l3_input_provider_backends = c("spm"),
      l1_setup = "spm_l1_model",
      l2_setup = NULL,
      l3_setup = "spm_l3_model",
      l1_status = "get_spm_status",
      l2_status = NULL,
      l3_status = "get_spm_status",
      l1_status_inputs = c("spm_dir"),
      l2_status_inputs = character(0),
      l3_status_inputs = c("spm_dir"),
      output_dir = "get_output_directory",
      l1_run = "run_spm_sepjobs",
      l2_run = NULL,
      l3_run = "run_spm_sepjobs"
    ),
    afni = list(
      name = "afni",
      runs_l1 = TRUE,
      runs_l2 = FALSE,
      runs_l3 = TRUE,
      multi_run_strategy = "concat_in_l1",
      produced_artifacts = c("subject_session_contrasts", "group_level_stats"),
      l3_required_artifacts = c("subject_session_contrasts"),
      l3_input_provider_backends = c("afni"),
      l1_setup = "afni_l1_model",
      l2_setup = NULL,
      l3_setup = "setup_l3_backend_afni",
      l1_status = "get_afni_glm_status",
      l2_status = NULL,
      l3_status = "get_afni_l3_status",
      l1_status_inputs = c("afni_dir"),
      l2_status_inputs = character(0),
      l3_status_inputs = c("output_file"),
      output_dir = "get_output_directory",
      l1_run = "run_afni_sepjobs",
      l2_run = NULL,
      l3_run = "run_afni_sepjobs"
    )
  )
}

# Upgrade the previous built-in AFNI placeholder in serialized pipeline objects.
# Custom backend implementations stay intact.
upgrade_glm_backend_specs <- function(specs) {
  if (is.null(specs)) return(default_glm_backend_specs())
  if (identical(specs$afni$l1_setup, "__not_implemented__") &&
      identical(specs$afni$l3_setup, "afni_3dlmer_setup")) {
    specs$afni <- default_glm_backend_specs()$afni
  }
  specs
}

backend_name <- function(backend) {
  if (!is.list(backend)) return("")
  if (!is.null(backend$name) && is.character(backend$name) && length(backend$name) == 1L) {
    return(tolower(backend$name))
  }
  ""
}

normalize_backend_strings <- function(x) {
  if (is.null(x) || length(x) == 0L) return(character(0))
  if (!is.character(x)) x <- as.character(x)
  x <- tolower(trimws(x))
  unique(x[nzchar(x) & !is.na(x)])
}

normalize_model_backend_override_map <- function(x) {
  if (is.null(x)) return(list())
  if (is.character(x)) {
    if (is.null(names(x)) || any(!nzchar(names(x)))) {
      stop("Model backend overrides must be a named character vector keyed by model name.", call. = FALSE)
    }
    out <- lapply(seq_along(x), function(ii) normalize_backend_strings(x[ii]))
    names(out) <- names(x)
    return(out[lengths(out) > 0L])
  }
  if (!is.list(x)) {
    stop("Model backend overrides must be a named character vector or named list.", call. = FALSE)
  }
  if (length(x) == 0L) return(list())
  if (is.null(names(x)) || any(!nzchar(names(x)))) {
    stop("Model backend overrides must be keyed by model name.", call. = FALSE)
  }
  out <- lapply(x, normalize_backend_strings)
  out[lengths(out) > 0L]
}

normalize_backend_override_config <- function(backend_overrides = NULL) {
  empty <- list(
    execution = list(l1 = list(), l2 = list(), l3 = list()),
    producer = list(l1 = list(), l2 = list(), l3 = list())
  )
  if (is.null(backend_overrides)) return(empty)
  if (!is.list(backend_overrides)) {
    stop("backend_overrides must be a list.", call. = FALSE)
  }

  out <- empty

  # Shorthand: list(l1 = c(model = "fsl"), l3 = c(model = "afni"))
  level_keys <- intersect(names(backend_overrides), c("l1", "l2", "l3"))
  for (level_key in level_keys) {
    entry <- backend_overrides[[level_key, exact = TRUE]]
    if (is.list(entry) && any(names(entry) %in% c("execution", "producer"))) {
      if (!is.null(entry$execution)) out$execution[[level_key]] <- normalize_model_backend_override_map(entry$execution)
      if (!is.null(entry$producer)) out$producer[[level_key]] <- normalize_model_backend_override_map(entry$producer)
    } else {
      out$execution[[level_key]] <- normalize_model_backend_override_map(entry)
    }
  }

  # Explicit shape: list(execution = list(l1 = ...), producer = list(l3 = ...))
  for (type in c("execution", "producer")) {
    type_entry <- backend_overrides[[type, exact = TRUE]]
    if (is.null(type_entry)) next
    if (!is.list(type_entry)) {
      stop(sprintf("backend_overrides$%s must be a list keyed by level.", type), call. = FALSE)
    }
    for (level_key in intersect(names(type_entry), c("l1", "l2", "l3"))) {
      out[[type]][[level_key]] <- normalize_model_backend_override_map(type_entry[[level_key]])
    }
  }

  out
}

backend_runs_level <- function(backend, level) {
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L)
  if (!is.list(backend)) return(FALSE)
  level_field <- switch(as.integer(level),
    "1" = "runs_l1",
    "2" = "runs_l2",
    "3" = "runs_l3"
  )
  isTRUE(backend[[level_field]])
}

backend_multi_run_strategy <- function(backend) {
  if (!is.list(backend)) return("not_applicable")
  val <- backend$multi_run_strategy
  if (!is.null(val) && is.character(val) && length(val) == 1L && nzchar(val)) {
    return(tolower(val))
  }
  "not_applicable"
}

backend_produced_artifacts <- function(backend) {
  if (!is.list(backend)) return(character(0))
  normalize_backend_strings(backend$produced_artifacts)
}

backend_produces_artifact <- function(backend, artifact) {
  checkmate::assert_string(artifact)
  artifact <- tolower(trimws(artifact))
  artifact %in% backend_produced_artifacts(backend)
}

backend_artifact_production_level <- function(backend, artifact, multi_run = TRUE) {
  checkmate::assert_string(artifact)
  checkmate::assert_flag(multi_run)
  if (!is.list(backend)) return(NA_integer_)

  artifact <- tolower(trimws(artifact))
  if (!backend_produces_artifact(backend, artifact)) return(NA_integer_)

  if (identical(artifact, "run_level_contrasts")) {
    return(if (backend_runs_level(backend, 1L)) 1L else NA_integer_)
  }

  if (identical(artifact, "subject_session_contrasts")) {
    strategy <- backend_multi_run_strategy(backend)
    if (identical(strategy, "explicit_l2")) {
      if (isTRUE(multi_run) && backend_runs_level(backend, 2L)) return(2L)
      if (backend_runs_level(backend, 1L)) return(1L)
      return(NA_integer_)
    }
    if (identical(strategy, "concat_in_l1")) {
      return(if (backend_runs_level(backend, 1L)) 1L else NA_integer_)
    }
  }

  if (identical(artifact, "group_level_stats")) {
    return(if (backend_runs_level(backend, 3L)) 3L else NA_integer_)
  }

  if (backend_runs_level(backend, 1L)) return(1L)
  if (backend_runs_level(backend, 2L)) return(2L)
  if (backend_runs_level(backend, 3L)) return(3L)
  NA_integer_
}

backend_l3_required_artifacts <- function(backend) {
  if (!is.list(backend)) return(character(0))
  val <- normalize_backend_strings(backend$l3_required_artifacts)
  if (length(val) > 0L) return(val)

  bname <- backend_name(backend)
  if (bname %in% c("fsl", "spm", "afni")) return("subject_session_contrasts")
  character(0)
}

backend_l3_input_provider_backends <- function(backend) {
  if (!is.list(backend)) return(character(0))
  val <- normalize_backend_strings(backend$l3_input_provider_backends)
  if (length(val) > 0L) return(val)

  legacy <- backend$l3_l2_source_backend
  val <- normalize_backend_strings(legacy)
  if (length(val) > 0L) return(val)

  bname <- backend_name(backend)
  if (bname == "fsl") return("fsl")
  if (bname == "spm") return("spm")
  if (bname == "afni") return("afni")
  character(0)
}

get_backends_producing_artifact <- function(backends, artifact) {
  checkmate::assert_list(backends)
  checkmate::assert_string(artifact)
  out <- names(backends)[vapply(backends, backend_produces_artifact, logical(1), artifact = artifact)]
  normalize_backend_strings(out)
}

resolve_l3_input_producers <- function(backends, l3_backend) {
  checkmate::assert_list(backends)
  if (!is.list(l3_backend)) return(character(0))

  explicit_providers <- backend_l3_input_provider_backends(l3_backend)
  if (length(explicit_providers) > 0L) return(explicit_providers)

  required_artifacts <- backend_l3_required_artifacts(l3_backend)
  if (length(required_artifacts) == 0L) return(character(0))

  producers <- names(backends)
  for (artifact in required_artifacts) {
    producers <- intersect(producers, get_backends_producing_artifact(backends, artifact))
  }

  normalize_backend_strings(producers)
}

resolve_l3_producer_requirements <- function(backends, l3_backend, multi_run = TRUE,
                                             producer_backends = NULL) {
  checkmate::assert_list(backends)
  checkmate::assert_flag(multi_run)
  checkmate::assert_character(producer_backends, null.ok = TRUE)
  if (!is.list(l3_backend)) {
    return(empty_requirement_df())
  }

  required_artifacts <- backend_l3_required_artifacts(l3_backend)
  if (length(required_artifacts) == 0L) {
    return(empty_requirement_df())
  }

  explicit_providers <- normalize_backend_strings(producer_backends)
  if (length(explicit_providers) == 0L) {
    explicit_providers <- backend_l3_input_provider_backends(l3_backend)
  }

  out <- list()
  for (artifact in required_artifacts) {
    providers <- explicit_providers
    if (length(providers) == 0L) {
      providers <- get_backends_producing_artifact(backends, artifact)
    }
    if (length(providers) == 0L) next

    providers <- providers[vapply(
      providers,
      function(provider_name) {
        provider <- backends[[provider_name]]
        !is.null(provider) && backend_produces_artifact(provider, artifact)
      },
      logical(1)
    )]

    if (length(providers) == 0L) next

    for (provider_name in providers) {
      provider <- backends[[provider_name]]
      producer_level <- backend_artifact_production_level(
        provider,
        artifact = artifact,
        multi_run = multi_run
      )
      if (is.na(producer_level)) next
      out[[length(out) + 1L]] <- data.frame(
        artifact = artifact,
        producer_backend = provider_name,
        producer_level = producer_level,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out) == 0L) {
    return(empty_requirement_df())
  }

  unique(do.call(rbind, out))
}

resolve_model_l3_requirements <- function(gpa, l3_model_names = NULL, execution_backend_map = NULL,
                                          producer_backend_map = NULL, specs = NULL,
                                          multi_run = isTRUE(gpa$multi_run)) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_list(execution_backend_map, null.ok = TRUE)
  checkmate::assert_list(producer_backend_map, null.ok = TRUE)
  checkmate::assert_list(specs, null.ok = TRUE)
  checkmate::assert_flag(multi_run)

  if (is.null(specs)) specs <- upgrade_glm_backend_specs(gpa$glm_backend_specs)
  if (is.null(specs)) specs <- default_glm_backend_specs()
  resolved_backends <- resolve_glm_backends(specs)

  if (is.null(l3_model_names)) l3_model_names <- get_model_names_for_level(gpa, 3L)
  if (is.null(execution_backend_map)) {
    execution_backend_map <- get_effective_model_backends(gpa, level = 3L, model_names = l3_model_names)
  }
  if (is.null(producer_backend_map)) {
    producer_backend_map <- get_effective_model_backends(
      gpa,
      level = 3L,
      model_names = l3_model_names,
      type = "producer"
    )
  }

  out <- list()
  for (model_name in l3_model_names) {
    execution_backends <- normalize_backend_strings(execution_backend_map[[model_name]])
    explicit_producers <- normalize_backend_strings(producer_backend_map[[model_name]])
    if (length(execution_backends) == 0L) next

    for (execution_backend in execution_backends) {
      l3_backend <- resolved_backends[[execution_backend]]
      if (is.null(l3_backend)) next

      req <- resolve_l3_producer_requirements(
        backends = resolved_backends,
        l3_backend = l3_backend,
        multi_run = multi_run,
        producer_backends = explicit_producers
      )
      if (nrow(req) == 0L) next

      req$model_name <- model_name
      req$execution_backend <- execution_backend
      out[[length(out) + 1L]] <- req[, c(
        "model_name", "execution_backend", "artifact", "producer_backend", "producer_level"
      )]
    }
  }

  if (length(out) == 0L) {
    return(empty_model_requirement_df())
  }

  unique(do.call(rbind, out))
}

get_requirement_producer_backends <- function(requirements, producer_level = NULL) {
  checkmate::assert_data_frame(requirements)
  checkmate::assert_integerish(producer_level, lower = 1L, upper = 3L, len = 1L, null.ok = TRUE)
  if (nrow(requirements) == 0L) return(character(0))
  req <- requirements
  if (!is.null(producer_level)) {
    req <- req[req$producer_level == as.integer(producer_level), , drop = FALSE]
  }
  normalize_backend_strings(req$producer_backend)
}

build_backend_preflight_report <- function(gpa, l1_model_backend_map = list(), l2_model_backend_map = list(),
                                           l3_model_backend_map = list(), l3_requirement_df = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_list(l1_model_backend_map)
  checkmate::assert_list(l2_model_backend_map)
  checkmate::assert_list(l3_model_backend_map)
  checkmate::assert_data_frame(l3_requirement_df, null.ok = TRUE)

  empty_report <- data.frame(
    level = integer(0),
    model_name = character(0),
    execution_backend = character(0),
    producer_backend = character(0),
    producer_level = integer(0),
    producer_artifact = character(0),
    stringsAsFactors = FALSE
  )

  add_execution_rows <- function(level, backend_map) {
    if (length(backend_map) == 0L) return(empty_report)
    rows <- lapply(names(backend_map), function(model_name) {
      exec_backends <- normalize_backend_strings(backend_map[[model_name]])
      if (length(exec_backends) == 0L) return(NULL)
      data.frame(
        level = rep.int(as.integer(level), length(exec_backends)),
        model_name = rep.int(model_name, length(exec_backends)),
        execution_backend = exec_backends,
        producer_backend = NA_character_,
        producer_level = NA_integer_,
        producer_artifact = NA_character_,
        stringsAsFactors = FALSE
      )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) return(empty_report)
    do.call(rbind, rows)
  }

  l1_rows <- add_execution_rows(1L, l1_model_backend_map)
  l2_rows <- add_execution_rows(2L, l2_model_backend_map)

  l3_rows <- if (length(l3_model_backend_map) == 0L) {
    empty_report
  } else {
    rows <- lapply(names(l3_model_backend_map), function(model_name) {
      exec_backends <- normalize_backend_strings(l3_model_backend_map[[model_name]])
      if (length(exec_backends) == 0L) return(NULL)

      lapply(exec_backends, function(exec_backend) {
        req <- if (!is.null(l3_requirement_df) && nrow(l3_requirement_df) > 0L) {
          l3_requirement_df[
            l3_requirement_df$model_name == model_name &
              l3_requirement_df$execution_backend == exec_backend,
            ,
            drop = FALSE
          ]
        } else {
          NULL
        }

        producer_backend <- if (!is.null(req) && nrow(req) > 0L) {
          paste(normalize_backend_strings(req$producer_backend), collapse = ",")
        } else {
          NA_character_
        }
        producer_level <- if (!is.null(req) && nrow(req) > 0L) {
          as.integer(min(unique(req$producer_level)))
        } else {
          NA_integer_
        }
        producer_artifact <- if (!is.null(req) && nrow(req) > 0L) {
          paste(normalize_backend_strings(req$artifact), collapse = ",")
        } else {
          NA_character_
        }

        data.frame(
          level = 3L,
          model_name = model_name,
          execution_backend = exec_backend,
          producer_backend = producer_backend,
          producer_level = producer_level,
          producer_artifact = producer_artifact,
          stringsAsFactors = FALSE
        )
      })
    })
    rows <- unlist(rows, recursive = FALSE)
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) == 0L) empty_report else do.call(rbind, rows)
  }

  report <- rbind(l1_rows, l2_rows, l3_rows)
  if (nrow(report) == 0L) return(empty_report)
  report[order(report$level, report$model_name, report$execution_backend), , drop = FALSE]
}

log_backend_preflight_report <- function(report, lg) {
  checkmate::assert_data_frame(report)
  checkmate::assert_class(lg, "Logger")

  if (nrow(report) == 0L) {
    lg$info("No backend preflight report entries were generated.")
    return(invisible(report))
  }

  lg$info("Resolved backend plan before submission:")
  for (ii in seq_len(nrow(report))) {
    row <- report[ii, , drop = FALSE]
    if (isTRUE(row$level == 3L) && !is.na(row$producer_backend)) {
      lg$info(
        "L%d model '%s': execute with %s; producer=%s at level %d for %s",
        row$level,
        row$model_name,
        row$execution_backend,
        row$producer_backend,
        row$producer_level,
        row$producer_artifact
      )
    } else {
      lg$info(
        "L%d model '%s': execute with %s",
        row$level,
        row$model_name,
        row$execution_backend
      )
    }
  }

  invisible(report)
}

any_l3_backend_requires_l2 <- function(backends, multi_run = TRUE) {
  checkmate::assert_list(backends)
  checkmate::assert_flag(multi_run)
  if (length(backends) == 0L) return(FALSE)
  if (isFALSE(multi_run)) return(FALSE)
  any(vapply(backends, function(backend) {
    if (!is.list(backend)) return(FALSE)
    req <- resolve_l3_producer_requirements(
      backends = backends,
      l3_backend = backend,
      multi_run = multi_run
    )
    nrow(req) > 0L && any(req$producer_level == 2L)
  }, logical(1)))
}

get_l3_dependency_backends <- function(backends, multi_run = TRUE) {
  checkmate::assert_list(backends)
  checkmate::assert_flag(multi_run)
  if (isFALSE(multi_run)) return(character(0))
  deps <- unique(unlist(lapply(backends, function(backend) {
    if (!is.list(backend)) return(NULL)
    req <- resolve_l3_producer_requirements(
      backends = backends,
      l3_backend = backend,
      multi_run = multi_run
    )
    if (nrow(req) > 0L && any(req$producer_level == 2L)) {
      return(unique(req$producer_backend[req$producer_level == 2L]))
    }
    NULL
  })))
  deps[!is.na(deps) & nzchar(deps)]
}

normalize_level_backend_config <- function(gpa, specs = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (is.null(specs)) specs <- upgrade_glm_backend_specs(gpa$glm_backend_specs)
  if (is.null(specs)) specs <- default_glm_backend_specs()

  legacy_backends <- normalize_backend_strings(gpa$glm_software)
  level_backends <- gpa$level_backends
  if (is.null(level_backends) || !is.list(level_backends)) level_backends <- list()

  out <- list()
  for (level in 1:3) {
    level_name <- paste0("l", level)
    configured <- level_backends[[level_name]]
    if (is.null(configured)) configured <- legacy_backends
    out[[level_name]] <- normalize_backend_strings(configured)
  }

  if (length(legacy_backends) == 0L && all(lengths(out) == 0L)) {
    out$l1 <- "fsl"
    out$l2 <- "fsl"
    out$l3 <- "fsl"
  }

  out
}

merge_level_backend_overrides <- function(gpa, level_backends = NULL, specs = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (is.null(specs)) specs <- upgrade_glm_backend_specs(gpa$glm_backend_specs)
  if (is.null(specs)) specs <- default_glm_backend_specs()

  out <- normalize_level_backend_config(gpa, specs = specs)
  if (is.null(level_backends)) return(out)
  if (!is.list(level_backends)) {
    stop("level_backends must be a list keyed by l1/l2/l3.", call. = FALSE)
  }
  for (level_key in intersect(names(level_backends), c("l1", "l2", "l3"))) {
    out[[level_key]] <- normalize_backend_strings(level_backends[[level_key]])
  }
  out
}

get_model_set_for_level <- function(gpa, level) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L)
  gpa[[paste0("l", level, "_models")]]
}

get_model_names_for_level <- function(gpa, level) {
  model_set <- get_model_set_for_level(gpa, level)
  if (is.null(model_set) || is.null(model_set$models)) return(character(0))
  names(model_set$models)
}

get_model_spec_for_level <- function(gpa, level, model_name) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L)
  checkmate::assert_string(model_name)
  model_set <- get_model_set_for_level(gpa, level)
  if (is.null(model_set) || is.null(model_set$models) || !model_name %in% names(model_set$models)) return(NULL)
  model_set$models[[model_name]]
}

get_model_backend_override <- function(gpa, level, model_name, type = c("execution", "producer"),
                                       backend_overrides = NULL) {
  type <- match.arg(type)
  level_key <- paste0("l", as.integer(level))
  overrides <- if (is.null(backend_overrides)) gpa$backend_overrides else normalize_backend_override_config(backend_overrides)
  if (is.null(overrides[[type]][[level_key]][[model_name]])) return(character(0))
  normalize_backend_strings(overrides[[type]][[level_key]][[model_name]])
}

validate_l3_backend_resolution <- function(gpa, l3_model_names, execution_backend_map,
                                           producer_backend_map = NULL, requested_backends = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_list(execution_backend_map)
  checkmate::assert_list(producer_backend_map, null.ok = TRUE)
  checkmate::assert_character(requested_backends, null.ok = TRUE)

  if (is.null(l3_model_names) || length(l3_model_names) == 0L) {
    return(invisible(TRUE))
  }

  requested_backends <- normalize_backend_strings(requested_backends)

  for (model_name in l3_model_names) {
    model_spec <- get_model_spec_for_level(gpa, level = 3L, model_name = model_name)
    l3_input_mode <- if (!is.null(model_spec)) normalize_l3_input_mode(model_spec$l3_input_mode) else "per_session"
    execution_backends <- normalize_backend_strings(execution_backend_map[[model_name]])
    producer_backends <- normalize_backend_strings(producer_backend_map[[model_name]])

    if (length(requested_backends) > 0L &&
        length(intersect(execution_backends, requested_backends)) == 0L) {
      stop(
        sprintf(
          "Requested backend filter '%s' excludes resolved execution backend(s) for L3 model '%s': %s",
          paste(requested_backends, collapse = ", "),
          model_name,
          paste(execution_backends, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    for (execution in execution_backends) {
      allowed <- if (execution == "afni") c("fsl", "afni") else execution
      if (length(producer_backends) && any(!producer_backends %in% allowed)) {
        stop(sprintf("L3 backend '%s' supports producer_backend=%s; requested: %s",
          execution, paste(sprintf("'%s'", allowed), collapse = " or "),
          paste(producer_backends, collapse = ", ")), call. = FALSE)
      }
    }
    if ("afni" %in% execution_backends && length(producer_backends) > 1L) {
      stop("Each AFNI L3 model must select one producer backend.", call. = FALSE)
    }

    if (identical(l3_input_mode, "3dlmer")) {
      if (length(execution_backends) == 0L || any(execution_backends != "afni")) {
        stop(
          sprintf(
            "L3 model '%s' uses l3_input_mode='3dlmer' and must execute only with backend 'afni'. Resolved execution backend(s): %s",
            model_name,
            if (length(execution_backends) > 0L) paste(execution_backends, collapse = ", ") else "<none>"
          ),
          call. = FALSE
        )
      }


    }
  }

  invisible(TRUE)
}

get_effective_model_backends <- function(gpa, level, model_names = NULL, type = c("execution", "producer"),
                                         level_backends = NULL, backend_overrides = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L)
  checkmate::assert_character(model_names, null.ok = TRUE)
  type <- match.arg(type)

  if (is.null(model_names)) model_names <- get_model_names_for_level(gpa, level)
  level_key <- paste0("l", as.integer(level))
  merged_level_backends <- merge_level_backend_overrides(gpa, level_backends = level_backends)
  default_backends <- if (type == "execution") merged_level_backends[[level_key]] else character(0)
  overrides <- if (is.null(backend_overrides)) gpa$backend_overrides else normalize_backend_override_config(backend_overrides)

  out <- setNames(vector("list", length(model_names)), model_names)
  for (model_name in model_names) {
    runtime_override <- get_model_backend_override(
      gpa = gpa, level = level, model_name = model_name, type = type, backend_overrides = overrides
    )
    if (length(runtime_override) > 0L) {
      out[[model_name]] <- runtime_override
      next
    }

    model_spec <- get_model_spec_for_level(gpa, level, model_name)
    spec_field <- if (type == "execution") "execution_backend" else "producer_backend"
    model_override <- if (!is.null(model_spec)) normalize_backend_strings(model_spec[[spec_field]]) else character(0)
    if (length(model_override) > 0L) {
      out[[model_name]] <- model_override
    } else if (isTRUE(level == 3L) && !is.null(model_spec) &&
               identical(normalize_l3_input_mode(model_spec$l3_input_mode), "3dlmer")) {
      out[[model_name]] <- if (identical(type, "execution")) "afni" else {
        # Preserve the established FSL -> 3dLMEr default for mixed pipelines.
        l1_backends <- merged_level_backends$l1
        if ("fsl" %in% l1_backends) "fsl" else "afni"
      }
    } else {
      out[[model_name]] <- default_backends
    }
  }
  out
}

group_models_by_backend <- function(model_backend_map) {
  checkmate::assert_list(model_backend_map)
  out <- list()
  for (model_name in names(model_backend_map)) {
    backends <- normalize_backend_strings(model_backend_map[[model_name]])
    for (backend_name in backends) {
      out[[backend_name]] <- unique(c(out[[backend_name]], model_name))
    }
  }
  out
}

get_level_backend_names <- function(gpa, level = NULL, must_exist = TRUE, specs = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_flag(must_exist)
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L, null.ok = TRUE)
  if (is.null(specs)) specs <- upgrade_glm_backend_specs(gpa$glm_backend_specs)
  if (is.null(specs)) specs <- default_glm_backend_specs()

  level_backends <- normalize_level_backend_config(gpa, specs = specs)
  backend_names <- if (is.null(level)) {
    normalize_backend_strings(unlist(level_backends, use.names = FALSE))
  } else {
    level_backends[[paste0("l", level)]]
  }

  missing_backends <- setdiff(backend_names, names(specs))
  if (length(missing_backends) > 0L) {
    if (isTRUE(must_exist)) {
      stop(sprintf("Unknown GLM backend '%s'. Available: %s", missing_backends[1L], paste(names(specs), collapse = ", ")), call. = FALSE)
    }
    backend_names <- setdiff(backend_names, missing_backends)
  }

  backend_names
}

resolve_glm_backends <- function(specs, pkg = "fmri.pipeline") {
  checkmate::assert_list(specs)

  # Declarative fields that every backend spec must provide
  required_fields <- c(
    "runs_l1", "runs_l2", "runs_l3",
    "multi_run_strategy", "produced_artifacts"
  )

  out <- list()
  for (backend_name in names(specs)) {
    spec <- specs[[backend_name]]
    if (is.null(spec) || !is.list(spec)) next
    resolved <- spec
    if (is.null(resolved$name) || !is.character(resolved$name) || length(resolved$name) != 1L || !nzchar(resolved$name)) {
      resolved$name <- backend_name
    }

    missing <- setdiff(required_fields, names(resolved))
    if (length(missing) > 0L) {
      stop(sprintf(
        "Backend spec '%s' is missing required fields: %s",
        backend_name, paste(missing, collapse = ", ")
      ), call. = FALSE)
    }

    for (field in c("l1_setup", "l2_setup", "l3_setup", "l1_status", "l2_status", "l3_status", "output_dir", "l1_run", "l2_run", "l3_run")) {
      val <- resolved[[field]]
      if (is.null(val)) next
      if (identical(val, "__not_implemented__")) {
        resolved[[field]] <- glm_backend_not_implemented(backend_name, field)
      } else if (is.character(val)) {
        resolved[[field]] <- resolve_backend_fn(val, pkg = pkg)
      }
    }
    out[[backend_name]] <- resolved
  }
  out
}

get_glm_backends <- function(gpa, must_exist = TRUE, level = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_flag(must_exist)
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L, null.ok = TRUE)

  specs <- upgrade_glm_backend_specs(gpa$glm_backend_specs)
  if (is.null(specs)) specs <- default_glm_backend_specs()
  backend_names <- get_level_backend_names(gpa, level = level, must_exist = must_exist, specs = specs)

  # Always resolve lazily from specs to avoid persisting stale function handles
  resolved <- resolve_glm_backends(specs)

  out <- list()
  for (name in backend_names) {
    if (!name %in% names(resolved)) {
      if (isTRUE(must_exist)) {
        stop(sprintf("Unknown GLM backend '%s'. Available: %s", name, paste(names(resolved), collapse = ", ")), call. = FALSE)
      }
      next
    }
    out[[name]] <- resolved[[name]]
  }

  out
}

initialize_glm_backends <- function(gpa) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  gpa$glm_backend_specs <- upgrade_glm_backend_specs(gpa$glm_backend_specs)
  gpa$level_backends <- normalize_level_backend_config(gpa, specs = gpa$glm_backend_specs)
  gpa$glm_software <- get_level_backend_names(gpa, must_exist = FALSE, specs = gpa$glm_backend_specs)
  # Do not persist function handles in the gpa object
  gpa$glm_backends <- NULL
  gpa
}
