# A single capability table drives configuration and execution preflight.
fwe_method_capabilities <- function() {
  data.frame(method = c("ptfce", "afni_3dclustsim_permutation"),
    executable = c(TRUE, FALSE), backend = c("fsl", "fsl"),
    label = c("pTFCE (each map corrected separately)", "AFNI residual simulation (standalone legacy wrapper)"))
}

assert_fwe_executable <- function(spec) {
  methods <- fwe_method_capabilities()
  if (!spec$method$name %in% methods$method[methods$executable]) {
    stop("No FWE planning adapter is implemented for method: ", spec$method$name,
      ". Automatic L3 correction currently supports pTFCE for FSL outputs.", call. = FALSE)
  }
  invisible(spec)
}

normalize_l3_inference <- function(inference, model_name, contrasts = NULL) {
  if (is.null(inference) || identical(inference, FALSE)) return(NULL)
  if (!is.list(inference)) stop("L3 inference must be a named list or FALSE.", call. = FALSE)
  if (!is.null(inference$schema_version)) {
    spec <- as_fwe_spec(inference)
  } else {
    inference$name <- inference$name %||% paste0(model_name, "_ptfce")
    inference$targets <- inference$targets %||% list(l3_model = model_name)
    inference$targets$l3_model <- inference$targets$l3_model %||% model_name
    spec <- do.call(fwe_spec, inference)
  }
  if (!identical(spec$targets$l3_model, model_name)) {
    stop("Model-level inference must target its own L3 model: ", model_name, call. = FALSE)
  }
  if (!is.null(contrasts) && !is.null(spec$targets$l3_cope_name)) {
    checkmate::assert_subset(spec$targets$l3_cope_name, contrasts)
  }
  assert_fwe_executable(spec)
  spec
}

pipeline_fwe_specs <- function(gpa, model_names = names(gpa$l3_models$models)) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_subset(model_names, names(gpa$l3_models$models))
  specs <- list()
  if (!length(model_names)) return(specs)
  for (name in model_names) {
    model <- gpa$l3_models$models[[name]]
    spec <- normalize_l3_inference(model$inference, name, rownames(model$contrasts))
    if (!is.null(spec)) specs[[length(specs) + 1L]] <- spec
  }
  for (value in gpa$fwe_correction) {
    if (!is.list(value) || is.null(value$method)) {
      stop("Legacy live FWE objects cannot run automatically. Reconfigure with build_fwe_correction() or fwe_spec().", call. = FALSE)
    }
    spec <- as_fwe_spec(value)
    if (!is.null(spec$targets$l3_model)) checkmate::assert_subset(spec$targets$l3_model, names(gpa$l3_models$models))
    selected <- if (is.null(spec$targets$l3_model)) model_names else intersect(spec$targets$l3_model, model_names)
    if (!length(selected)) next
    spec$targets$l3_model <- selected
    assert_fwe_executable(spec)
    specs[[length(specs) + 1L]] <- spec
  }
  names(specs) <- vapply(specs, `[[`, character(1L), "name")
  if (anyDuplicated(names(specs))) stop("FWE correction names must be unique across model and pipeline settings.")
  specs
}

validate_pipeline_fwe <- function(gpa, specs, automatic = FALSE) {
  for (spec in specs) {
    assert_fwe_executable(spec)
    model_names <- spec$targets$l3_model
    checkmate::assert_subset(model_names, names(gpa$l3_models$models))
    backends <- get_effective_model_backends(gpa, level = 3L, model_names = model_names)
    if (any(vapply(backends, function(x) !identical(x, "fsl"), logical(1L)))) {
      stop("FWE correction '", spec$name, "' requires an FSL L3 model; AFNI/SPM adapters are not available yet.", call. = FALSE)
    }
    if (automatic && spec$compute$scheduler != "inherit") {
      scheduler <- normalize_fwe_run_scheduler(spec$compute$scheduler)
      if (!scheduler %in% c("local", normalize_fwe_run_scheduler(gpa$scheduler))) {
        stop("Automatic FWE correction must use the pipeline scheduler or local execution. ",
          "Use run_fwe_correction() to submit independently to another scheduler.", call. = FALSE)
      }
    }
  }
  invisible(specs)
}

prompt_ptfce_spec <- function(model, existing = NULL) {
  methods <- fwe_method_capabilities()
  choices <- methods[methods$executable, , drop = FALSE]
  action <- utils::menu(c("No automatic correction", choices$label,
    if (!is.null(existing)) "Keep current correction"), title = "Final L3 group inference")
  if (action == 0L || action == nrow(choices) + 2L) return(existing)
  if (action == 1L) return(NULL)
  targets <- existing$targets %||% list(l3_model = model$name)
  selected <- utils::select.list(rownames(model$contrasts), multiple = TRUE,
    preselect = existing$targets$l3_cope_name %||% rownames(model$contrasts),
    title = "Contrasts to correct (FWE controlled separately within each map)")
  if (!length(selected) || !any(nzchar(selected))) return(existing)
  targets$l3_cope_name <- selected
  default_alpha <- existing$method$fwe_alpha %||% 0.05
  alpha <- readline(paste0("FWE alpha values, separated by spaces [", paste(default_alpha, collapse = " "), "]: "))
  alpha <- if (!nzchar(trimws(alpha))) default_alpha else scan(text = alpha, quiet = TRUE)
  sided <- utils::menu(c("Two-sided (positive and negative)", "One-sided (positive)"), title = "Direction of inference")
  default_region <- if (identical(existing$correction_mask$source, "file")) existing$correction_mask$path else "model"
  region <- readline(paste0("Correction mask: model or path to a prespecified binary ROI [", default_region, "]: "))
  if (!nzchar(trimws(region))) region <- default_region
  fwe_spec(name = existing$name %||% paste0(model$name, "_ptfce"),
    targets = targets, method = choices$method[action - 1L],
    fwe_alpha = alpha, two_sided = if (sided == 0L) existing$method$two_sided %||% TRUE else sided != 2L,
    correction_mask = region, compute = existing$compute %||% list(),
    write_threshold_images = existing$method$write_threshold_images %||% TRUE)
}

configure_l3_inference <- function(gpa, model) {
  model$inference <- normalize_l3_inference(model$inference, model$name, rownames(model$contrasts))
  temporary <- gpa
  temporary$l3_models$models[[model$name]] <- model
  backends <- get_effective_model_backends(temporary, 3L, model$name)[[model$name]]
  if (!identical(backends, "fsl")) {
    if (!is.null(model$inference)) stop("Automatic pTFCE requires an FSL L3 model.")
    message("Automatic FWE correction is currently available for FSL L3 models only.")
    return(model)
  }
  model$inference <- prompt_ptfce_spec(model, model$inference)
  model
}

save_fwe_pipeline_state <- function(gpa) {
  path <- gpa$output_locations$object_cache
  if (!is.null(path)) fwe_atomic_write(path, function(file) saveRDS(gpa, file))
  invisible(gpa)
}

#' Run configured group inference without refitting GLMs
#'
#' Runs the declarative corrections stored in L3 models' \code{inference}
#' fields and in \code{gpa$fwe_correction}. Scheduler runs return after
#' submission; use \code{refresh_fwe_correction()} to collect them.
#' @param gpa a \code{glm_pipeline_arguments} object.
#' @param l3_model_names models whose configured corrections should run;
#'   defaults to all L3 models.
#' @param force rerun corrections even when verified outputs exist.
#' @param dry_run preview plans and commands without writing or submitting.
#' @param ... additional arguments to \code{run_fwe_plan()}.
#' @return Updated \code{gpa}, with \code{fwe_runs} and \code{fwe_status}.
#' @export
run_fwe_correction <- function(gpa, l3_model_names = names(gpa$l3_models$models),
                               force = FALSE, dry_run = FALSE, ...) {
  checkmate::assert_flag(force)
  checkmate::assert_flag(dry_run)
  specs <- pipeline_fwe_specs(gpa, l3_model_names)
  validate_pipeline_fwe(gpa, specs)
  run_pipeline_fwe_specs(gpa, specs, force = force, dry_run = dry_run, ...)
}

run_pipeline_fwe_specs <- function(gpa, specs, force = FALSE, dry_run = FALSE, ...) {
  arguments <- utils::modifyList(list(sched_args = gpa$parallel$sched_args,
    tracking_sqlite_db = gpa$output_locations$sqlite_db), list(...))
  for (name in names(specs)) {
    run <- tryCatch({
      plan <- plan_fwe_correction(gpa, specs[[name]])
      do.call(run_fwe_plan, c(list(plan = plan, force = force, dry_run = dry_run), arguments))
    }, error = function(e) {
        if (inherits(e, "fwe_execution_error")) gpa$fwe_runs[[name]] <- e$run
        gpa$fwe_status[[name]] <- "failed"
        if (!dry_run) save_fwe_pipeline_state(gpa)
        stop(e)
      })
    gpa$fwe_runs[[name]] <- run
    gpa$fwe_status[[name]] <- if (dry_run) "planned" else if (all(run$plan$tasks$status == "complete")) {
      "complete"
    } else if (any(run$execution$execution_status %in% c("failed", "incomplete", "submission_failed"))) {
      "failed"
    } else if (length(run$job_ids)) "submitted" else "incomplete"
    if (!dry_run) save_fwe_pipeline_state(gpa)
  }
  gpa
}

#' Refresh configured group inference results
#' @param gpa a \code{glm_pipeline_arguments} with \code{fwe_runs}.
#' @param correction_names corrections to refresh; defaults to all recorded runs.
#' @param require_complete error after persisting state if any correction is incomplete.
#' @return Updated \code{gpa}, including result manifests and separate FWE status.
#' @export
refresh_fwe_correction <- function(gpa, correction_names = names(gpa$fwe_runs), require_complete = FALSE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_subset(correction_names, names(gpa$fwe_runs))
  checkmate::assert_flag(require_complete)
  for (name in correction_names) {
    run <- gpa$fwe_runs[[name]]
    result <- collect_fwe_results(run)
    run$plan <- result$plan
    run$result <- result
    gpa$fwe_runs[[name]] <- run
    gpa$fwe_status[[name]] <- if (result$complete) "complete" else "incomplete"
  }
  save_fwe_pipeline_state(gpa)
  if (require_complete && any(unlist(gpa$fwe_status[correction_names]) != "complete")) {
    stop("Group estimation completed, but requested FWE correction is incomplete. Inspect gpa$fwe_runs and retry run_fwe_correction().")
  }
  gpa
}
