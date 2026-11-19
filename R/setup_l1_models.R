#' Function for setting up level 1 model specifications and corresponding files
#'
#' @param gpa A glm_model_arguments function setup by the \code{setup_glm_pipeline} function
#' @param l1_model_names A character vector of model names within \code{gpa$l1_models} whose l1
#'    inputs should be generated. If omitted, all models within \code{gpa$l1_models} will be generated.
#'
#' @details The \code{l1_model_names} argument allows the creation of multiple l1 models to be parallelized at
#'   a superordinate level or, even better, to be spread across independent jobs on a cluster. This function
#'   already provides the option to parallelize over subjects for a single model if \code{gpa$l1_setup_cpus}
#'   is greater than 1.
#'
#' @author Michael Hallquist
#'
#' @importFrom checkmate assert_class assert_character
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ
#' @importFrom iterators iter
#' @importFrom parallel makeCluster stopCluster
#' @importFrom RNifti niftiHeader
#' @importFrom dplyr filter pull
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @export
setup_l1_models <- function(gpa, l1_model_names=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_names(
    names(gpa$run_data),
    must.include = c("id", "session", "run_number", "run_nifti", "exclude_run")
  )

  ctx <- init_l1_setup_context(gpa = gpa, l1_model_names = l1_model_names)
  gpa <- ctx$gpa

  # TODO: This is a mess if use the l1_model_names since we will always be overwriting what's in the l1_model_setup
  # data.frame. We need more of an append/update approach, perhaps like the sqlite setup for the overall pipeline.

  cl <- register_l1_parallel_backend(gpa = gpa, lg = ctx$lg)
  if (!is.null(cl)) {
    on.exit(try(parallel::stopCluster(cl)), add = TRUE)
  }

  all_subj_l1_list <- foreach(
    subj_df = iter(gpa$subject_data, by = "row"), .inorder = FALSE,
    .packages = c("dplyr", "fmri.pipeline"), .errorhandling = "remove",
    .export = c(
      "truncate_runs", "fsl_l1_model", "spm_l1_model", "get_spm_status",
      "get_mr_abspath", "get_output_directory", "run_fsl_command",
      "get_feat_status", "add_custom_feat_syntax"
    )
  ) %dopar% {
    subj_df <- subj_df # avoid complaints about visible global binding in R CMD check
    tryCatch(
      setup_l1_subject(subj_df = subj_df, gpa = gpa, ctx = ctx),
      error = function(e) {
        ctx$lg$error(
          "Error in setup_l1_subject for id=%s, session=%s: %s",
          as.character(subj_df$id[1L]),
          if ("session" %in% names(subj_df)) as.character(subj_df$session[1L]) else "NA",
          conditionMessage(e)
        )
        NULL
      }
    )
  }

  gpa$l1_model_setup <- combine_l1_subject_results(all_subj_l1_list)
  gpa <- refresh_glm_status(gpa, level=1L, lg=ctx$lg)

  return(gpa)
}

bind_l1_log_appenders <- function(lg, gpa) {
  checkmate::assert_class(lg, "Logger")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  sync_appender <- function(name, enabled, factory, target_path) {
    if (!isTRUE(enabled)) return(invisible(NULL))
    existing <- lg$appenders[[name]]
    existing_path <- if (!is.null(existing)) existing$file %||% NULL else NULL
    if (!is.null(existing) && !identical(existing_path, target_path)) {
      lg$remove_appender(name)
      existing <- NULL
    }
    if (is.null(existing)) {
      lg$add_appender(factory(target_path), name = name)
    }
    invisible(NULL)
  }

  sync_appender(
    name = "setup_l1_log_json",
    enabled = gpa$log_json,
    factory = function(path) lgr::AppenderJson$new(path),
    target_path = gpa$output_locations$setup_l1_log_json
  )

  sync_appender(
    name = "setup_l1_log_txt",
    enabled = gpa$log_txt,
    factory = function(path) lgr::AppenderFile$new(path),
    target_path = gpa$output_locations$setup_l1_log_txt
  )
}


# Internal helper: resolve backends, logger state, and SPM session-mode context once.
init_l1_setup_context <- function(gpa, l1_model_names) {
  if (is.null(l1_model_names)) {
    l1_model_names <- names(gpa$l1_models$models)
  }

  backend_specs <- gpa$glm_backend_specs
  if (is.null(backend_specs)) backend_specs <- default_glm_backend_specs()
  gpa$glm_backend_specs <- backend_specs
  l1_model_backend_map <- get_effective_model_backends(gpa, level = 1L, model_names = l1_model_names)
  requested_l1_backends <- normalize_backend_strings(unlist(l1_model_backend_map, use.names = FALSE))
  resolved_backends <- resolve_glm_backends(backend_specs)
  glm_backends <- resolved_backends[intersect(requested_l1_backends, names(resolved_backends))]

  lg <- lgr::get_logger("glm_pipeline/l1_setup")
  lg$set_threshold(gpa$lgr_threshold)
  bind_l1_log_appenders(lg = lg, gpa = gpa)

  gpa <- refresh_l1_cope_names(gpa, lg = lg)
  gpa <- refresh_glm_status(gpa, level = 1L, lg = lg)

  l1_backend_names <- names(glm_backends)
  spm_l1_session_mode <- "separate"
  spm_anchor_by_id <- NULL
  if ("spm" %in% l1_backend_names) {
    spm_mode_raw <- NULL
    if (!is.null(gpa$glm_settings) && !is.null(gpa$glm_settings$spm) &&
      "l1_session_mode" %in% names(gpa$glm_settings$spm)) {
      spm_mode_raw <- gpa$glm_settings$spm$l1_session_mode
    }
    spm_l1_session_mode <- resolve_spm_l1_session_mode(spm_mode_raw, lg = lg)
    if (identical(spm_l1_session_mode, "pooled")) {
      subj_ids <- unique(as.character(gpa$subject_data$id))
      spm_anchor_by_id <- setNames(
        vapply(subj_ids, function(sid) {
          min(as.integer(gpa$subject_data$session[gpa$subject_data$id == sid]), na.rm = TRUE)
        }, integer(1)),
        subj_ids
      )
      lg$info("SPM L1 session mode is 'pooled'; SPM L1 models will be created once per id across all sessions.")
    } else {
      lg$info("SPM L1 session mode is 'separate'; SPM L1 models will be created per id/session.")
    }
  }

  list(
    gpa = gpa,
    l1_model_names = l1_model_names,
    l1_model_backend_map = l1_model_backend_map,
    glm_backends = glm_backends,
    l1_backend_names = l1_backend_names,
    spm_l1_session_mode = spm_l1_session_mode,
    spm_anchor_by_id = spm_anchor_by_id,
    lg = lg
  )
}


# Internal helper: register the requested parallel backend and return the cluster if one was created.
register_l1_parallel_backend <- function(gpa, lg) {
  if (!is.null(gpa$parallel$l1_setup_cores) && gpa$parallel$l1_setup_cores[1L] > 1L) {
    lg$info("Initializing l1 setup cluster with %d cores", gpa$parallel$l1_setup_cores[1L])
    cl <- parallel::makeCluster(gpa$parallel$l1_setup_cores[1L])
    doParallel::registerDoParallel(cl)
    return(cl)
  }

  lg$info("Initializing l1 setup with serial execution")
  foreach::registerDoSEQ()
  NULL
}


# Internal helper: process one id/session row and return the same per-subject result shape used previously.
setup_l1_subject <- function(subj_df, gpa, ctx) {
  subject_ctx <- prepare_subject_run_context(subj_df = subj_df, gpa = gpa, ctx = ctx)
  if (is.null(subject_ctx)) return(NULL)
  slg <- subject_ctx$lg

  l1_file_setup <- list(fsl = list(), spm = list(), afni = list(), metadata = subject_ctx$mr_df)

  for (this_model in ctx$l1_model_names) {
    model_ctx <- prepare_l1_model_inputs(
      gpa = gpa,
      model_name = this_model,
      subj_id = subject_ctx$subj_id,
      subj_session = subject_ctx$subj_session,
      lg = slg
    )

    output_ctx <- ensure_l1_output_dir(
      gpa = gpa,
      subj_id = subject_ctx$subj_id,
      subj_session = subject_ctx$subj_session,
      model_name = this_model,
      lg = slg
    )

    d_obj <- load_or_build_l1_bdm(
      subject_ctx = subject_ctx,
      model_ctx = model_ctx,
      gpa = gpa,
      model_name = this_model,
      l1_output_dir = output_ctx$l1_output_dir,
      bdm_out_file = output_ctx$bdm_out_file,
      lg = slg
    )
    if (is.null(d_obj)) next

    feat_l1_df <- setup_l1_model_fsl(
      subject_ctx = subject_ctx,
      model_ctx = model_ctx,
      d_obj = d_obj,
      gpa = gpa,
      model_name = this_model,
      backend_fsl = ctx$glm_backends[["fsl"]],
      lg = slg
    )
    if (!is.null(feat_l1_df)) {
      l1_file_setup$fsl <- rbind(l1_file_setup$fsl, feat_l1_df)
    }

    spm_l1_df <- setup_l1_model_spm(
      subject_ctx = subject_ctx,
      model_ctx = model_ctx,
      d_obj = d_obj,
      gpa = gpa,
      model_name = this_model,
      ctx = ctx
    )
    if (!is.null(spm_l1_df)) {
      l1_file_setup$spm <- rbind(l1_file_setup$spm, spm_l1_df)
    }
  }

  slg$info("Completed processing of subject: %s", subject_ctx$subj_id)
  return(l1_file_setup)
}


# Internal helper: collect, validate, and prune run-level inputs for a single id/session row.
prepare_subject_run_context <- function(subj_df, gpa, ctx) {
  subj_id <- subj_df$id
  subj_session <- subj_df$session
  slg <- get_subject_logger(
    base_logger = "glm_pipeline/l1_setup",
    id = subj_id,
    gpa = gpa,
    log_prefix = "setup_l1_models"
  )
  is_spm_anchor <- !identical(ctx$spm_l1_session_mode, "pooled") ||
    (!is.null(ctx$spm_anchor_by_id) &&
      as.character(subj_id) %in% names(ctx$spm_anchor_by_id) &&
      as.integer(subj_session) == as.integer(ctx$spm_anchor_by_id[[as.character(subj_id)]]))

  rdata <- gpa$run_data %>%
    dplyr::filter(id == !!subj_id & session == !!subj_session & run_nifti_present == TRUE)

  if (!"tr" %in% names(rdata)) {
    slg$error("No tr column in run data for subject: %s, session: %d", subj_id, subj_session)
    return(NULL)
  }

  if (length(unique(rdata$tr)) > 1L) {
    slg$error(
      "More than one TR value for runs within a subject. This is not currently supported! subject: %s, session: %d",
      subj_id, subj_session
    )
    return(NULL)
  }

  l1_confound_files <- if ("l1_confound_file" %in% names(rdata)) {
    rdata$l1_confound_file
  } else {
    rep("", nrow(rdata))
  }

  run_df <- data.frame(
    id = subj_id,
    session = subj_session,
    run_number = rdata$run_number,
    run_nifti = get_mr_abspath(rdata, "run_nifti"),
    l1_confound_file = l1_confound_files,
    run_volumes = rdata$run_volumes,
    nvoxels = rdata$nvoxels,
    exclude_run = rdata$exclude_run,
    stringsAsFactors = FALSE
  )

  subj_mr_dir <- subj_df$mr_dir
  if (nrow(run_df) == 0L) {
    slg$warn("Unable to find any preprocessed fMRI files in dir: %s", subj_mr_dir)
    return(NULL)
  }

  nii_present <- file.exists(run_df$run_nifti)
  if (any(!nii_present)) {
    slg$warn("Could not find some of the expected preprocessed fMRI files. These will be dropped.")
    slg$warn("Missing: %s", run_df$run_nifti[!nii_present])
    run_df <- run_df[nii_present, , drop = FALSE]
  } else {
    slg$debug(paste("MR files to analyze:", run_df$run_nifti))
  }

  if (nrow(run_df) == 0L) {
    slg$warn("Unable to find any preprocessed fMRI files in dir: %s", subj_mr_dir)
    return(NULL)
  }

  slg$debug("Volumes in run_nifti: %s", paste(run_df$run_volumes, collapse = ", "))

  # This is a fresh aggregation, so normalize it by reference before passing it
  # to data.frame-oriented design helpers.
  m_events <- data.table::rbindlist(
    lapply(gpa$l1_models$events, function(this_event) {
      this_event$data %>% dplyr::filter(id == !!subj_id & session == !!subj_session)
    })
  ) |>
    data.table::setDF()
  dropped <- drop_runs_without_events(
    run_df = run_df,
    events = m_events,
    subj_id = subj_id,
    subj_session = subj_session,
    lg = slg
  )
  run_df <- dropped$run_df
  m_events <- dropped$events

  if (nrow(run_df) == 0L) {
    slg$warn(
      "No analyzable runs remain after dropping runs with missing events for id=%s, session=%s. Skipping subject.",
      subj_id, subj_session
    )
    return(NULL)
  }

  mr_df <- data.frame(
    id = subj_id,
    session = subj_session,
    run_number = run_df$run_number,
    run_nifti = run_df$run_nifti,
    l1_confound_file = run_df$l1_confound_file,
    run_volumes = run_df$run_volumes,
    exclude_run = run_df$exclude_run,
    row.names = NULL
  )

  l1_confound_files_spm <- run_df$l1_confound_file
  if ("spm" %in% ctx$l1_backend_names && nrow(run_df) > 1L) {
    confound_outdir <- if (length(run_df$l1_confound_file) > 0L && file.exists(run_df$l1_confound_file[1L])) {
      dirname(run_df$l1_confound_file[1L])
    } else {
      get_output_directory(
        id = subj_id, session = subj_session, gpa = gpa,
        what = "sub", create_if_missing = TRUE
      )
    }
    l1_confound_files_spm <- resolve_spm_confound_inputs(
      gpa = gpa,
      id = subj_id,
      session = subj_session,
      run_numbers = run_df$run_number,
      confound_files = run_df$l1_confound_file,
      output_dir = confound_outdir,
      lg = slg,
      pooled = FALSE
    )
  }

  list(
    subj_id = subj_id,
    subj_session = subj_session,
    subj_mr_dir = subj_mr_dir,
    lg = slg,
    is_spm_anchor = is_spm_anchor,
    rdata = rdata,
    tr = rdata$tr[1L],
    run_df = run_df,
    mr_df = mr_df,
    m_events = m_events,
    l1_confound_files_spm = l1_confound_files_spm
  )
}


# Internal helper: drop runs that have imaging data but no matching event rows.
drop_runs_without_events <- function(run_df, events, subj_id, subj_session = NULL, model_name = NULL, lg,
                                     pooled = FALSE) {
  event_run_nums <- if (is.data.frame(events) && nrow(events) > 0L && "run_number" %in% names(events)) {
    sort(unique(events$run_number))
  } else {
    integer(0)
  }
  missing_event_runs <- setdiff(run_df$run_number, event_run_nums)
  if (length(missing_event_runs) == 0L) {
    return(list(run_df = run_df, events = events))
  }

  keep_runs <- run_df$run_number %in% event_run_nums
  if (isTRUE(pooled)) {
    lg$warn(
      paste(
        "Dropping %d pooled SPM run(s) with missing event rows for id=%s, model=%s:",
        "pooled runs=%s.",
        "These runs have imaging inputs but no matching event/trial data."
      ),
      length(missing_event_runs), subj_id, model_name,
      paste(missing_event_runs, collapse = ", ")
    )
  } else {
    lg$warn(
      paste(
        "Dropping %d run(s) with missing event rows for id=%s, session=%s:",
        "runs=%s.",
        "These runs have imaging inputs but no matching event/trial data."
      ),
      length(missing_event_runs), subj_id, subj_session,
      paste(missing_event_runs, collapse = ", ")
    )
  }

  run_df <- run_df[keep_runs, , drop = FALSE]
  if (is.data.frame(events) && nrow(events) > 0L && "run_number" %in% names(events)) {
    events <- events %>% dplyr::filter(.data$run_number %in% !!run_df$run_number)
  }

  list(run_df = run_df, events = events)
}


# Internal helper: filter model-specific signals and ts-multiplier data for one subject context.
prepare_l1_model_inputs <- function(gpa, model_name, subj_id, subj_session = NULL, lg,
                                    remap_run_numbers = NULL, warn_on_empty = TRUE) {
  model_signal_names <- gpa$l1_models$models[[model_name]]$signals
  model_backend_names <- normalize_backend_strings(get_effective_model_backends(gpa, level = 1L, model_names = model_name)[[model_name]])

  ts_multiplier_cols <- unlist(lapply(gpa$l1_models$signals[model_signal_names], function(this_signal) {
    if (isFALSE(this_signal$ts_multiplier)) {
      return(NULL)
    }
    this_signal$ts_multiplier
  }))

  ts_multiplier_data <- if (!is.null(ts_multiplier_cols)) {
    if (is.null(subj_session)) {
      gpa$ppi_data %>% dplyr::filter(id == !!subj_id)
    } else {
      gpa$ppi_data %>% dplyr::filter(id == !!subj_id & session == !!subj_session)
    }
  } else {
    NULL
  }
  if (is.function(remap_run_numbers) && is.data.frame(ts_multiplier_data) && nrow(ts_multiplier_data) > 0L) {
    ts_multiplier_data <- remap_run_numbers(ts_multiplier_data, what = "ts_multipliers")
  }

  m_signals <- lapply(gpa$l1_models$signals[model_signal_names], function(this_signal) {
    if (!inherits(this_signal$value, "data.frame")) {
      msg <- "Unable to sort out how to refit wi_model with signal that doesn't have a data.frame in the $value slot"
      lg$error(msg)
      stop(msg)
    }

    if (is.null(subj_session)) {
      this_signal$value <- this_signal$value %>% dplyr::filter(id == !!subj_id)
    } else {
      this_signal$value <- this_signal$value %>% dplyr::filter(id == !!subj_id & session == !!subj_session)
    }

    if (is.function(remap_run_numbers)) {
      this_signal$value <- remap_run_numbers(this_signal$value, what = paste0("signal ", this_signal$name))
    }

    if (isTRUE(warn_on_empty) && nrow(this_signal$value) == 0L) {
      msg <- if (is.null(subj_session)) {
        glue(
          "In L1 model setup, failed to find any rows in in gpa$l1_models$signals${this_signal$name}$value",
          " for id: {subj_id}, model: {model_name}.\n  We will create an empty regressor.",
          .trim = FALSE
        )
      } else {
        glue(
          "In L1 model setup, failed to find any rows in in gpa$l1_models$signals${this_signal$name}$value",
          " for id: {subj_id}, session: {subj_session}, model: {model_name}.\n  We will create an empty regressor.",
          .trim = FALSE
        )
      }
      lg$warn(msg)
      warning(msg)
    }

    if (!is.null(this_signal$wi_factors)) {
      this_signal <- fit_wi_model(this_signal)
    }
    this_signal
  })

  list(
    model_backend_names = model_backend_names,
    ts_multiplier_cols = ts_multiplier_cols,
    ts_multiplier_data = ts_multiplier_data,
    m_signals = m_signals
  )
}


# Internal helper: ensure the per-model output directory exists and return the cache path.
ensure_l1_output_dir <- function(gpa, subj_id, subj_session, model_name, lg) {
  l1_output_dir <- get_output_directory(
    id = subj_id, session = subj_session,
    l1_model = model_name, gpa = gpa, what = "l1"
  )
  if (!dir.exists(l1_output_dir)) {
    lg$info("Creating subject output directory: %s", l1_output_dir)
    dir.create(l1_output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  list(
    l1_output_dir = l1_output_dir,
    bdm_out_file = file.path(l1_output_dir, paste0(gpa$l1_models$models[[model_name]]$name, "_bdm_setup.RData"))
  )
}


# Internal helper: load an existing BDM cache without rebuilding if it is still valid.
load_cached_l1_bdm <- function(bdm_out_file, lg, expected_args = NULL) {
  if (!file.exists(bdm_out_file)) {
    return(list(run_bdm = TRUE, d_obj = NULL))
  }

  lg$info("Loading BDM info from extant file: %s", bdm_out_file)
  out <- tryCatch(
    {
      env <- new.env(parent = baseenv())
      load(bdm_out_file, envir = env)
      if (!is.null(expected_args)) {
        saved_args <- env$bdm_args
        if (is.null(saved_args) || !identical(saved_args[setdiff(names(saved_args), "lg")],
                                             expected_args[setdiff(names(expected_args), "lg")])) {
          return(list(run_bdm = TRUE, d_obj = NULL))
        }
      }
      d_obj <- if (exists("d_obj", envir = env, inherits = FALSE)) env$d_obj else NULL
      if (is.null(d_obj) || !inherits(d_obj, "bdm")) {
        lg$warn("Regenerating design because d_obj is missing or not a valid bdm object in extant BDM file: %s.", bdm_out_file)
        return(list(run_bdm = TRUE, d_obj = NULL))
      }
      if ("run_4d_files" %in% names(d_obj)) {
        d_obj$run_nifti <- d_obj$run_4d_files
      }
      list(run_bdm = FALSE, d_obj = d_obj)
    },
    error = function(e) {
      lg$error("Failed to load BDM file: %s with error: %s. I will regenerate the design matrix.", bdm_out_file, e)
      list(run_bdm = TRUE, d_obj = NULL)
    }
  )

  out
}


# Internal helper: assemble the argument list used for build_design_matrix.
make_l1_bdm_args <- function(gpa, events, signals, tr, write_timing_files, run_data,
                             runs_to_output, output_directory, lg, ts_multipliers = NULL) {
  bdm_args <- gpa$additional$bdm_args
  bdm_args$events <- events
  bdm_args$signals <- signals
  bdm_args$tr <- tr
  bdm_args$write_timing_files <- write_timing_files
  bdm_args$drop_volumes <- gpa$drop_volumes
  bdm_args$run_data <- run_data
  bdm_args$runs_to_output <- runs_to_output
  bdm_args$output_directory <- output_directory
  bdm_args$lg <- lg
  bdm_args$ts_multipliers <- ts_multipliers
  bdm_args
}


# Internal helper: run build_design_matrix with consistent error logging.
run_l1_build_design_matrix <- function(bdm_args, subj_id, model_name, lg, subj_session = NULL, pooled = FALSE) {
  tryCatch(
    do.call(build_design_matrix, bdm_args),
    error = function(e) {
      if (isTRUE(pooled)) {
        lg$error("Failed pooled SPM build_design_matrix for id: %s, model: %s", subj_id, model_name)
      } else {
        lg$error("Failed build_design_matrix for id: %s, session: %s, model: %s", subj_id, subj_session, model_name)
      }
      lg$error("Error message: %s", as.character(e))
      return(NULL)
    }
  )
}


# Internal helper: share the same inputs between BDM cache checks and creation.
l1_bdm_setup_args <- function(subject_ctx, model_ctx, gpa, l1_output_dir, lg) {
  t_out <- model_ctx$model_backend_names
  if (isTRUE(gpa$use_preconvolve)) t_out <- c("convolved", t_out)

  make_l1_bdm_args(
    gpa = gpa,
    events = subject_ctx$m_events,
    signals = model_ctx$m_signals,
    tr = subject_ctx$tr,
    write_timing_files = t_out,
    run_data = subject_ctx$mr_df,
    runs_to_output = subject_ctx$run_df$run_number,
    output_directory = file.path(l1_output_dir, "timing_files"),
    lg = lg,
    ts_multipliers = model_ctx$ts_multiplier_data
  )
}

build_and_cache_l1_bdm <- function(subject_ctx, model_ctx, gpa, model_name, l1_output_dir, bdm_out_file, lg) {
  bdm_args <- l1_bdm_setup_args(subject_ctx, model_ctx, gpa, l1_output_dir, lg)

  d_obj <- run_l1_build_design_matrix(
    bdm_args = bdm_args,
    subj_id = subject_ctx$subj_id,
    subj_session = subject_ctx$subj_session,
    model_name = model_name,
    lg = lg
  )

  mr_df <- subject_ctx$mr_df
  mr_run_nums <- subject_ctx$run_df$run_number
  subj_mr_dir <- subject_ctx$subj_mr_dir
  run_nifti <- subject_ctx$run_df$run_nifti
  run_volumes <- subject_ctx$run_df$run_volumes
  exclude_run <- subject_ctx$run_df$exclude_run
  nvoxels <- subject_ctx$run_df$nvoxels
  l1_confound_files <- subject_ctx$run_df$l1_confound_file
  subj_id <- subject_ctx$subj_id
  subj_session <- subject_ctx$subj_session
  this_model <- model_name

  save(
    d_obj, bdm_args, mr_df, mr_run_nums, subj_mr_dir, run_nifti, run_volumes, exclude_run, nvoxels,
    l1_confound_files, subj_id, subj_session, this_model,
    file = bdm_out_file
  )

  d_obj
}


# Internal helper: prefer a cached BDM if it exists, otherwise rebuild and save it.
load_or_build_l1_bdm <- function(subject_ctx, model_ctx, gpa, model_name, l1_output_dir, bdm_out_file, lg) {
  force_rebuild <- isTRUE(gpa$glm_settings$fsl$force_l1_creation) ||
    isTRUE(gpa$glm_settings$spm$force_l1_creation) ||
    isTRUE(gpa$glm_settings$afni$force_l1_creation)

  if (!isTRUE(force_rebuild)) {
    expected_args <- if ("spm" %in% model_ctx$model_backend_names) {
      l1_bdm_setup_args(subject_ctx, model_ctx, gpa, l1_output_dir, lg)
    } else NULL
    cache <- load_cached_l1_bdm(bdm_out_file = bdm_out_file, lg = lg, expected_args = expected_args)
    if (!isTRUE(cache$run_bdm)) {
      return(cache$d_obj)
    }
  } else if (file.exists(bdm_out_file)) {
    lg$info("force_l1_creation is TRUE; rebuilding BDM design instead of loading cached %s", bdm_out_file)
  }

  build_and_cache_l1_bdm(
    subject_ctx = subject_ctx,
    model_ctx = model_ctx,
    gpa = gpa,
    model_name = model_name,
    l1_output_dir = l1_output_dir,
    bdm_out_file = bdm_out_file,
    lg = lg
  )
}


# Internal helper: run the FSL backend branch for one subject/model combination.
setup_l1_model_fsl <- function(subject_ctx, model_ctx, d_obj, gpa, model_name, backend_fsl, lg) {
  if (is.null(backend_fsl) || !"fsl" %in% model_ctx$model_backend_names) {
    return(NULL)
  }

  tryCatch(
    backend_fsl$l1_setup(
      id = subject_ctx$subj_id,
      session = subject_ctx$subj_session,
      l1_confound_files = subject_ctx$run_df$l1_confound_file,
      d_obj = d_obj,
      gpa = gpa,
      model_name,
      nvoxels = subject_ctx$run_df$nvoxels
    ),
    error = function(e) {
      lg$error(
        "Problem with fsl_l1_model. Model: %s, Subject: %s, Session: %s",
        model_name, subject_ctx$subj_id, subject_ctx$subj_session
      )
      lg$error("Error message: %s", as.character(e))
      return(NULL)
    }
  )
}


# Internal helper: dispatch SPM setup while keeping separate and pooled session modes isolated.
setup_l1_model_spm <- function(subject_ctx, model_ctx, d_obj, gpa, model_name, ctx) {
  lg <- subject_ctx$lg
  if (is.null(lg)) lg <- ctx$lg

  if (!"spm" %in% model_ctx$model_backend_names) {
    lg$debug(
      "Skipping SPM L1 setup for model %s because its effective execution backend excludes SPM.",
      model_name
    )
    return(NULL)
  }

  backend_spm <- ctx$glm_backends[["spm"]]
  if (is.null(backend_spm)) {
    lg$warn(
      "SPM backend is not configured; skipping SPM L1 setup for subject %s session %s.",
      subject_ctx$subj_id, subject_ctx$subj_session
    )
    return(NULL)
  }

  if (identical(ctx$spm_l1_session_mode, "pooled") && !isTRUE(subject_ctx$is_spm_anchor)) {
    lg$debug(
      "Skipping SPM pooled L1 setup for non-anchor row id=%s session=%s model=%s (anchor session=%s).",
      subject_ctx$subj_id, subject_ctx$subj_session, model_name,
      ctx$spm_anchor_by_id[[as.character(subject_ctx$subj_id)]]
    )
    return(NULL)
  }

  if (!identical(ctx$spm_l1_session_mode, "pooled")) {
    return(setup_l1_model_spm_separate(
      subject_ctx = subject_ctx,
      d_obj = d_obj,
      gpa = gpa,
      model_name = model_name,
      backend_spm = backend_spm,
      lg = lg
    ))
  }

  setup_l1_model_spm_pooled(
    subject_ctx = subject_ctx,
    gpa = gpa,
    model_name = model_name,
    backend_spm = backend_spm,
    lg = lg
  )
}


# Internal helper: run the standard per-session SPM setup branch.
setup_l1_model_spm_separate <- function(subject_ctx, d_obj, gpa, model_name, backend_spm, lg) {
  lg$info(
    "Setting up SPM L1 model %s for subject %s session %s.",
    model_name, subject_ctx$subj_id, subject_ctx$subj_session
  )
  spm_l1_df <- tryCatch(
    backend_spm$l1_setup(
      id = subject_ctx$subj_id,
      session = subject_ctx$subj_session,
      l1_confound_files = subject_ctx$l1_confound_files_spm,
      d_obj = d_obj,
      gpa = gpa,
      model_name = model_name,
      run_nifti = subject_ctx$run_df$run_nifti,
      run_numbers = subject_ctx$run_df$run_number,
      run_sessions = rep.int(as.integer(subject_ctx$subj_session), nrow(subject_ctx$run_df)),
      source_run_numbers = subject_ctx$run_df$run_number,
      l1_session_mode = "separate"
    ),
    error = function(e) {
      lg$error(
        "Problem running spm_l1_model. Model: %s, Subject: %s, Session: %s",
        model_name, subject_ctx$subj_id, subject_ctx$subj_session
      )
      lg$error("Error message: %s", as.character(e))
      return(NULL)
    }
  )

  if (is.null(spm_l1_df)) {
    lg$warn(
      "SPM L1 model setup returned NULL for model %s, subject %s, session %s.",
      model_name, subject_ctx$subj_id, subject_ctx$subj_session
    )
  }
  spm_l1_df
}


# Internal helper: construct the pooled-session SPM inputs, including run remapping.
prepare_pooled_spm_context <- function(subject_ctx, gpa, model_name, lg) {
  pooled_rdata <- gpa$run_data %>% dplyr::filter(id == !!subject_ctx$subj_id)
  if ("run_nifti_present" %in% names(pooled_rdata)) {
    pooled_rdata <- pooled_rdata %>% dplyr::filter(run_nifti_present == TRUE)
  }
  pooled_rdata <- pooled_rdata[order(pooled_rdata$session, pooled_rdata$run_number), , drop = FALSE]

  if (nrow(pooled_rdata) == 0L) {
    lg$warn("No run data found for pooled SPM L1 setup: id=%s model=%s", subject_ctx$subj_id, model_name)
    return(NULL)
  }
  if (!"tr" %in% names(pooled_rdata) || length(unique(pooled_rdata$tr)) > 1L) {
    lg$error("SPM pooled L1 requires one TR per subject across sessions. id=%s model=%s", subject_ctx$subj_id, model_name)
    return(NULL)
  }

  pooled_run_df <- data.frame(
    id = subject_ctx$subj_id,
    session = as.integer(pooled_rdata$session),
    run_number = seq_len(nrow(pooled_rdata)),
    source_run_number = as.integer(pooled_rdata$run_number),
    run_nifti = get_mr_abspath(pooled_rdata, "run_nifti"),
    l1_confound_file = if ("l1_confound_file" %in% names(pooled_rdata)) pooled_rdata$l1_confound_file else rep("", nrow(pooled_rdata)),
    run_volumes = pooled_rdata$run_volumes,
    exclude_run = pooled_rdata$exclude_run,
    stringsAsFactors = FALSE
  )

  run_key <- paste(pooled_run_df$session, pooled_run_df$source_run_number, sep = "::")
  remap_run_numbers <- function(df, what = "data") {
    if (!is.data.frame(df) || nrow(df) == 0L) return(df)
    if (!all(c("session", "run_number") %in% names(df))) return(df)
    idx <- match(paste(df$session, df$run_number, sep = "::"), run_key)
    if (anyNA(idx)) {
      lg$warn(
        "Dropping %d row(s) with unmatched session/run_number while remapping pooled SPM %s for id=%s model=%s.",
        sum(is.na(idx)), what, subject_ctx$subj_id, model_name
      )
      df <- df[!is.na(idx), , drop = FALSE]
      idx <- idx[!is.na(idx)]
    }
    df$run_number <- pooled_run_df$run_number[idx]
    df
  }

  m_events_spm <- data.table::rbindlist(
    lapply(gpa$l1_models$events, function(this_event) {
      this_event$data %>% dplyr::filter(id == !!subject_ctx$subj_id)
    })
  ) |>
    data.table::setDF()
  m_events_spm <- remap_run_numbers(m_events_spm, what = "events")

  dropped <- drop_runs_without_events(
    run_df = pooled_run_df,
    events = m_events_spm,
    subj_id = subject_ctx$subj_id,
    model_name = model_name,
    lg = lg,
    pooled = TRUE
  )
  pooled_run_df <- dropped$run_df
  m_events_spm <- dropped$events

  if (nrow(pooled_run_df) == 0L) {
    lg$warn(
      "No analyzable pooled SPM runs remain after dropping runs with missing events for id=%s, model=%s. Skipping pooled SPM setup.",
      subject_ctx$subj_id, model_name
    )
    return(NULL)
  }

  model_ctx_spm <- prepare_l1_model_inputs(
    gpa = gpa,
    model_name = model_name,
    subj_id = subject_ctx$subj_id,
    subj_session = NULL,
    lg = lg,
    remap_run_numbers = remap_run_numbers,
    warn_on_empty = FALSE
  )

  mr_df_spm <- data.frame(
    id = subject_ctx$subj_id,
    session = pooled_run_df$session,
    run_number = pooled_run_df$run_number,
    run_nifti = pooled_run_df$run_nifti,
    l1_confound_file = pooled_run_df$l1_confound_file,
    run_volumes = pooled_run_df$run_volumes,
    exclude_run = pooled_run_df$exclude_run,
    row.names = NULL
  )

  spm_outdir <- get_output_directory(
    id = subject_ctx$subj_id, session = 0L, l1_model = model_name,
    gpa = gpa, glm_software = "spm", what = "l1", create_if_missing = TRUE
  )
  t_out_spm <- "spm"
  if (isTRUE(gpa$use_preconvolve)) t_out_spm <- c("convolved", t_out_spm)
  bdm_args_spm <- make_l1_bdm_args(
    gpa = gpa,
    events = m_events_spm,
    signals = model_ctx_spm$m_signals,
    tr = pooled_rdata$tr[1L],
    write_timing_files = t_out_spm,
    run_data = mr_df_spm,
    runs_to_output = pooled_run_df$run_number,
    output_directory = file.path(spm_outdir, "timing_files"),
    lg = lg,
    ts_multipliers = model_ctx_spm$ts_multiplier_data
  )

  list(
    pooled_rdata = pooled_rdata,
    pooled_run_df = pooled_run_df,
    mr_df_spm = mr_df_spm,
    m_events_spm = m_events_spm,
    model_ctx_spm = model_ctx_spm,
    spm_outdir = spm_outdir,
    bdm_args_spm = bdm_args_spm
  )
}


# Internal helper: resolve whether SPM should use per-run or concatenated confound files.
resolve_spm_confound_inputs <- function(gpa, id, session, run_numbers, confound_files, output_dir, lg,
                                        run_sessions = NULL, file_name = NULL, model_name = NULL, pooled = FALSE) {
  spm_settings <- populate_defaults(gpa$glm_settings$spm, list(concatenate_runs = FALSE))
  if (is.null(spm_settings$concatenate_runs)) {
    spm_settings$concatenate_runs <- FALSE
  }
  if (!length(confound_files) || all(!is.na(confound_files) & !nzchar(confound_files))) return(NULL)
  if (anyNA(confound_files) || any(!nzchar(confound_files)) || any(!file.exists(confound_files))) {
    stop("SPM confounds are missing or invalid; refusing to omit nuisance adjustment.")
  }
  if (!isTRUE(spm_settings$concatenate_runs)) return(confound_files)

  concat_args <- list(gpa = gpa, id = id, session = session, run_numbers = run_numbers,
    confound_files = confound_files, output_dir = output_dir, lg = lg)
  if (!is.null(run_sessions)) concat_args$run_sessions <- run_sessions
  if (!is.null(file_name)) concat_args$file_name <- file_name
  concat_file <- do.call(concat_l1_confounds, concat_args)
  if (is.null(concat_file)) stop("Unable to concatenate SPM confounds; refusing to omit nuisance adjustment.")
  concat_file
}


# Internal helper: run the pooled-session SPM branch after the anchor session has been selected.
setup_l1_model_spm_pooled <- function(subject_ctx, gpa, model_name, backend_spm, lg) {
  lg$info("Setting up pooled-session SPM L1 model %s for subject %s.", model_name, subject_ctx$subj_id)
  pooled_ctx <- prepare_pooled_spm_context(
    subject_ctx = subject_ctx,
    gpa = gpa,
    model_name = model_name,
    lg = lg
  )
  if (is.null(pooled_ctx)) {
    return(NULL)
  }

  d_obj_spm <- run_l1_build_design_matrix(
    bdm_args = pooled_ctx$bdm_args_spm,
    subj_id = subject_ctx$subj_id,
    model_name = model_name,
    lg = lg,
    pooled = TRUE
  )
  if (is.null(d_obj_spm)) {
    return(NULL)
  }

  l1_confound_files_spm_use <- resolve_spm_confound_inputs(
    gpa = gpa,
    id = subject_ctx$subj_id,
    session = 0L,
    run_numbers = pooled_ctx$pooled_run_df$source_run_number,
    run_sessions = pooled_ctx$pooled_run_df$session,
    confound_files = pooled_ctx$pooled_run_df$l1_confound_file,
    output_dir = pooled_ctx$spm_outdir,
    file_name = paste0("l1_confounds_concat_id", subject_ctx$subj_id, ".txt"),
    model_name = model_name,
    lg = lg,
    pooled = TRUE
  )

  spm_l1_df <- tryCatch(
    backend_spm$l1_setup(
      id = subject_ctx$subj_id,
      session = 0L,
      l1_confound_files = l1_confound_files_spm_use,
      d_obj = d_obj_spm,
      gpa = gpa,
      model_name = model_name,
      run_nifti = pooled_ctx$pooled_run_df$run_nifti,
      run_numbers = pooled_ctx$pooled_run_df$run_number,
      run_sessions = pooled_ctx$pooled_run_df$session,
      source_run_numbers = pooled_ctx$pooled_run_df$source_run_number,
      l1_session_mode = "pooled"
    ),
    error = function(e) {
      lg$error("Problem running pooled spm_l1_model. Model: %s, Subject: %s", model_name, subject_ctx$subj_id)
      lg$error("Error message: %s", as.character(e))
      return(NULL)
    }
  )

  if (is.null(spm_l1_df)) {
    lg$warn("SPM pooled L1 model setup returned NULL for model %s, subject %s.", model_name, subject_ctx$subj_id)
  }
  spm_l1_df
}


# Internal helper: combine per-subject results into the historical l1_setup object shape.
combine_l1_subject_results <- function(all_subj_l1_list) {
  spm_list <- lapply(all_subj_l1_list, "[[", "spm")
  spm_combined <- NULL
  if (!all(vapply(spm_list, is.null, logical(1)))) {
    spm_combined <- rbindlist(spm_list, fill = TRUE) |> data.table::setDF()
  }

  fsl_combined <- rbindlist(lapply(all_subj_l1_list, "[[", "fsl"), fill = TRUE) |> data.table::setDF()
  metadata_combined <- rbindlist(lapply(all_subj_l1_list, "[[", "metadata"), fill = TRUE) |> data.table::setDF()

  all_subj_l1_combined <- list(
    fsl = fsl_combined,
    spm = spm_combined,
    metadata = metadata_combined
  )
  class(all_subj_l1_combined) <- c("l1_setup", "list")

  all_subj_l1_combined
}
