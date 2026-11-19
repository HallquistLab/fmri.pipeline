# TODO: Support ability not to run models at a certain level (e.g., level 3) if only lower-level models are needed, or also
# run just higher-level models if lower levels are already complete.

#' primary function for running a GLM analysis pipeline
#' 
#' @param gpa a glm_pipeline_arguments object containing a model specification (created by setup_glm_pipeline)
#' @param l1_model_names a character vector of level 1 model names (specified during build_l1_models) that should be executed
#' @param l2_model_names a character vector of level 2 model names (specified during build_l2_models) that should be executed
#' @param l3_model_names a character vector of level 3 model names (specified during build_l3_models) that should be executed
#' @param glm_software which GLM software should be used for model estimation (fsl, spm, or afni)
#' @param level_backends optional per-level backend override list keyed by `l1`, `l2`, and/or `l3`
#' @param backend_overrides optional model-specific backend override list
#' @importFrom checkmate assert_string assert_class assert_subset assert_integerish
#' @importFrom stats setNames
#' @importFrom utils menu modifyList
#' @export
run_glm_pipeline <- function(gpa, l1_model_names = "prompt", l2_model_names = "prompt",
l3_model_names = "prompt", glm_software = NULL, level_backends = NULL, backend_overrides = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_character(glm_software, null.ok = TRUE)

  lg <- lgr::get_logger("glm_pipeline/run_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  if (!is.null(glm_software)) {
    glm_software <- normalize_backend_strings(glm_software)
    gpa$glm_software <- glm_software
    gpa$level_backends <- list(l1 = glm_software, l2 = glm_software, l3 = glm_software)
  }

  backend_specs <- upgrade_glm_backend_specs(gpa$glm_backend_specs)
  if (is.null(backend_specs)) backend_specs <- default_glm_backend_specs()
  gpa$glm_backend_specs <- backend_specs
  gpa$level_backends <- merge_level_backend_overrides(gpa, level_backends = level_backends, specs = backend_specs)
  gpa$backend_overrides <- normalize_backend_override_config(backend_overrides)
  gpa <- initialize_glm_backends(gpa)

  resolved_backends <- resolve_glm_backends(backend_specs)

  model_list <- choose_glm_set(gpa, l1_model_names, l2_model_names, l3_model_names, lg)
  if (is.null(model_list)) { return(invisible(NULL)) } # user canceled

  l1_model_backend_map <- get_effective_model_backends(gpa, level = 1L, model_names = model_list$l1_model_names)
  l2_model_backend_map <- get_effective_model_backends(gpa, level = 2L, model_names = model_list$l2_model_names)
  l3_model_backend_map <- get_effective_model_backends(gpa, level = 3L, model_names = model_list$l3_model_names)
  l3_producer_backend_map <- get_effective_model_backends(gpa, level = 3L, model_names = model_list$l3_model_names, type = "producer")
  validate_l3_backend_resolution(
    gpa = gpa,
    l3_model_names = model_list$l3_model_names,
    execution_backend_map = l3_model_backend_map,
    producer_backend_map = l3_producer_backend_map
  )
  l3_requirement_df <- if (!is.null(model_list$l3_model_names)) {
    resolve_model_l3_requirements(
      gpa = gpa,
      l3_model_names = model_list$l3_model_names,
      execution_backend_map = l3_model_backend_map,
      producer_backend_map = l3_producer_backend_map,
      specs = backend_specs,
      multi_run = isTRUE(gpa$multi_run)
    )
  } else {
    empty_model_requirement_df()
  }

  unsupported_producers <- setdiff(
    normalize_backend_strings(unlist(l3_producer_backend_map, use.names = FALSE)),
    c("", "fsl", "afni")
  )
  if (length(unsupported_producers) > 0L) {
    stop(
      sprintf(
        "L3 producer backends are currently implemented for FSL and AFNI. Unsupported producer override(s): %s",
        paste(unsupported_producers, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  l1_models_by_backend <- group_models_by_backend(l1_model_backend_map)
  l2_models_by_backend <- group_models_by_backend(l2_model_backend_map)
  l3_models_by_backend <- group_models_by_backend(l3_model_backend_map)

  requested_l1_backend_names <- names(l1_models_by_backend)
  requested_l2_backend_names <- names(l2_models_by_backend)
  requested_l3_backend_names <- names(l3_models_by_backend)

  l3_dependency_backends <- if (!is.null(model_list$l3_model_names)) {
    get_requirement_producer_backends(l3_requirement_df)
  } else {
    character(0)
  }

  reassign_unrunnable_models <- function(model_backend_map, level, fallback_backends = character(0)) {
    if (length(model_backend_map) == 0L) return(model_backend_map)
    fallback_backends <- fallback_backends[vapply(
      fallback_backends,
      function(backend_name) {
        backend <- resolved_backends[[backend_name]]
        !is.null(backend) && backend_runs_level(backend, level)
      },
      logical(1)
    )]

    out <- model_backend_map
    for (model_name in names(out)) {
      current_backends <- normalize_backend_strings(out[[model_name]])
      runnable_backends <- current_backends[vapply(
        current_backends,
        function(backend_name) {
          backend <- resolved_backends[[backend_name]]
          !is.null(backend) && backend_runs_level(backend, level)
        },
        logical(1)
      )]
      if (length(runnable_backends) > 0L) {
        out[[model_name]] <- runnable_backends
      } else if (length(fallback_backends) > 0L) {
        out[[model_name]] <- fallback_backends
      } else {
        out[[model_name]] <- current_backends
      }
    }
    out
  }

  l1_model_backend_map <- reassign_unrunnable_models(l1_model_backend_map, level = 1L, fallback_backends = l3_dependency_backends)
  l2_model_backend_map <- reassign_unrunnable_models(l2_model_backend_map, level = 2L, fallback_backends = l3_dependency_backends)
  # A runnable L1 backend may differ from an explicitly requested L3 producer.
  # In that case estimate both, so the producer's L2/L3 jobs have actual parents.
  for (producer in l3_dependency_backends) {
    if (backend_runs_level(resolved_backends[[producer]], 1L)) {
      for (model in names(l1_model_backend_map)) {
        l1_model_backend_map[[model]] <- unique(c(l1_model_backend_map[[model]], producer))
      }
    }
  }
  l1_models_by_backend <- group_models_by_backend(l1_model_backend_map)
  l2_models_by_backend <- group_models_by_backend(l2_model_backend_map)
  requested_l1_backend_names <- names(l1_models_by_backend)
  requested_l2_backend_names <- names(l2_models_by_backend)

  execution_level_backends <- list(
    l1 = requested_l1_backend_names,
    l2 = requested_l2_backend_names,
    l3 = requested_l3_backend_names
  )
  for (backend_name in l3_dependency_backends) {
    backend <- resolved_backends[[backend_name]]
    if (is.null(backend)) next
    if (backend_runs_level(backend, 1L)) {
      execution_level_backends$l1 <- normalize_backend_strings(c(execution_level_backends$l1, backend_name))
    }
    if (backend_runs_level(backend, 2L)) {
      execution_level_backends$l2 <- normalize_backend_strings(c(execution_level_backends$l2, backend_name))
    }
  }
  execution_backend_names <- normalize_backend_strings(unlist(execution_level_backends, use.names = FALSE))
  missing_execution_backends <- setdiff(execution_backend_names, names(resolved_backends))
  if (length(missing_execution_backends) > 0L) {
    stop(
      sprintf(
        "Missing backend specs for required execution backends: %s",
        paste(missing_execution_backends, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  gpa$level_backends <- execution_level_backends
  gpa$glm_software <- execution_backend_names
  gpa$backend_preflight_report <- build_backend_preflight_report(
    gpa = gpa,
    l1_model_backend_map = l1_model_backend_map,
    l2_model_backend_map = l2_model_backend_map,
    l3_model_backend_map = l3_model_backend_map,
    l3_requirement_df = l3_requirement_df
  )
  log_backend_preflight_report(gpa$backend_preflight_report, lg)

  batch_id <- uuid::UUIDgenerate()
  if (is.null(gpa$batch_run)) gpa$batch_run <- list()

  # setup a unique directory for this GLM batch submission
  # this will contain all scripts for this run and contain a snapshot of the gpa object
  # that is amended as each batch job runs and completes.
  batch_directory <- file.path(gpa$output_locations$scheduler_scripts, paste0("batch_", batch_id))
  if (!dir.exists(batch_directory)) dir.create(batch_directory, recursive = TRUE)
  gpa$batch_run <- modifyList(gpa$batch_run %||% list(), list(
    sequence_id = batch_id,
    batch_directory = batch_directory
  ))
  gpa$output_locations$active_batch_directory <- batch_directory
  gpa_cache <- file.path(batch_directory, "run_pipeline_cache.RData")
  save(gpa, file=gpa_cache)

  # batch job to finalize pipeline configuration
  f_batch <- R_batch_job$new(
    job_name = "finalize_configuration", batch_directory = batch_directory, scheduler = gpa$scheduler,
    input_rdata_file = gpa_cache, output_rdata_file = gpa_cache,
    n_nodes = 1, n_cpus = 1, wall_time = gpa$parallel$finalize_time,
    mem_total = "16G",
    r_code = "gpa <- finalize_pipeline_configuration(gpa)", r_packages = "fmri.pipeline",
    batch_code = get_compute_environment(gpa),
    scheduler_options = gpa$parallel$sched_args,
    sqlite_db = gpa$output_locations$sqlite_db
  )

  if (is.null(gpa$finalize_complete) || isFALSE(gpa$finalize_complete)) {
    lg$info("finalize_pipeline_configuration has not been run on this object. We will start with this step.")
    run_finalize <- TRUE
  } else {
    run_finalize <- FALSE
  }

  # initialize batch sequence elements as NULL so that they are ignored in R_batch_sequence if not used
  l1_setup_batch <- split_cache_batch <- NULL
  l1_execute_batches <- list()
  l3_execute_batches <- list()
  l1_backend_names <- execution_level_backends$l1
  l2_backend_names <- execution_level_backends$l2
  l3_backend_names <- execution_level_backends$l3
  l3_glm_backends <- resolved_backends[intersect(l3_backend_names, names(resolved_backends))]
  execution_glm_backends <- resolved_backends[execution_backend_names]
  l3_level2_requirements <- if (nrow(l3_requirement_df) > 0L) {
    unique(l3_requirement_df[l3_requirement_df$producer_level == 2L, c("execution_backend", "producer_backend", "producer_level"), drop = FALSE])
  } else {
    data.frame(
      execution_backend = character(0),
      producer_backend = character(0),
      producer_level = integer(0),
      stringsAsFactors = FALSE
    )
  }
  l3_cross_backend_requirements <- if (nrow(l3_requirement_df) > 0L) {
    unique(l3_requirement_df[
      l3_requirement_df$producer_backend != l3_requirement_df$execution_backend,
      c("execution_backend", "producer_backend", "producer_level"),
      drop = FALSE
    ])
  } else {
    data.frame(
      execution_backend = character(0),
      producer_backend = character(0),
      producer_level = integer(0),
      stringsAsFactors = FALSE
    )
  }
  backend_cache_map <- setNames(
    file.path(batch_directory, paste0("run_pipeline_cache_", execution_backend_names, ".RData")),
    execution_backend_names
  )

  # Note: one tricky job dependency problem occurs when a job spawns multiple child jobs.
  # In this case, downstream jobs cannot know the job ids to wait for because these have not been generated yet.
  # This occurs in the case of L1 job submissions, where run_feat_sepjobs submit many jobs that need to complete before
  # L2 jobs should start. The new approach is to keep the parent job active until all children complete. This is
  # implemented by wait_for_children in R_batch_job. Here, we run the l1 model setup, the launch all feat runs in batch
  # jobs and wait for all of these to complete before the l1_setup_batch (parent) job completes.

  if (!is.null(model_list$l1_model_names)) {
    # batch job for setting up l1 models -- calls setup_l1_models to create relevant FSFs
    l1_setup_batch <- f_batch$copy(
      job_name = "setup_l1", n_cpus = gpa$parallel$l1_setup_cores,
      wall_time = gpa$parallel$l1_setup_time, mem_total = gpa$parallel$l1_setup_memgb,
      r_code = sprintf(
        "gpa <- setup_l1_models(gpa, l1_model_names=%s)", paste(deparse(model_list$l1_model_names), collapse = "")
      )
    )

    l1_setup_batch$depends_on_parents <- "finalize_configuration"

    if (length(execution_backend_names) > 0L) {
      split_cache_batch <- f_batch$copy(
        job_name = "split_backend_caches", n_cpus = 1,
        wall_time = "0:10:00",
        r_code = c(
          sprintf("shared_cache <- '%s'", gpa_cache),
          "if (!file.exists(shared_cache)) {",
          "  stop('Shared GPA cache not found: ', shared_cache, '. The setup_l1 job may have failed.')",
          "}",
          sprintf(
            "backend_cache_paths <- c(%s)",
            paste(shQuote(unname(backend_cache_map)), collapse = ", ")
          ),
          "for (path in backend_cache_paths) {",
          "  success <- file.copy(shared_cache, path, overwrite = TRUE)",
          "  if (!isTRUE(success)) {",
          "    stop('Failed to copy shared cache to backend cache: ', path)",
          "  }",
          "}"
        )
      )

      split_cache_batch$depends_on_parents <- "setup_l1"
    }

    # batch jobs for executing l1 jobs (and waiting) after setup, per backend
    backend_l1_exec_time <- function(backend_name) {
      if (backend_name == "fsl") {
        return(gpa$parallel$fsl$l1_feat_alljobs_time)
      }
      if (backend_name == "spm") {
        return(gpa$parallel$spm$l1_spm_alljobs_time)
      }
      if (backend_name == "afni") {
        return(gpa$parallel$afni$l1_glm_alljobs_time %||% "24:00:00")
      }
      stop(sprintf("No L1 execution time configured for backend '%s'.", backend_name))
    }

    for (backend_name in names(l1_models_by_backend)) {
      spec <- backend_specs[[backend_name]]
      run_fn_name <- if (is.null(spec)) NULL else spec$l1_run
      if (is.null(run_fn_name) || identical(run_fn_name, "__not_implemented__")) next
      if (!is.character(run_fn_name)) {
        lg$warn("Skipping backend '%s' L1 runner because l1_run is not a character function name.", backend_name)
        next
      }
      if (!is.null(execution_glm_backends[[backend_name]]$l1_run) &&
        isTRUE(attr(execution_glm_backends[[backend_name]]$l1_run, "glm_backend_not_implemented"))) {
        next
      }

      l1_exec <- f_batch$copy(
        job_name = paste0("run_l1_", backend_name), n_cpus = 1,
        wall_time = backend_l1_exec_time(backend_name),
        r_code = c(
          "child_job_ids <- c()",
          sprintf(
            "child_job_ids <- c(child_job_ids, %s(gpa, level = 1L, model_names=%s))",
            run_fn_name, paste(deparse(l1_models_by_backend[[backend_name]]), collapse = "")
          )
        ),
        post_children_r_code = c(
          sprintf("gpa <- refresh_glm_status(gpa, level = 1L, glm_software = '%s')", backend_name)
        ) # refresh l1 status after all jobs complete
      )

      l1_exec$depends_on_parents <- if (!is.null(split_cache_batch)) "split_backend_caches" else "setup_l1"
      l1_exec$wait_for_children <- TRUE # need to wait for l1 jobs to complete before moving to l2/l3
      if (!is.null(backend_cache_map[[backend_name]])) {
        l1_exec$input_rdata_file <- backend_cache_map[[backend_name]]
        l1_exec$output_rdata_file <- backend_cache_map[[backend_name]]
      }
      l1_execute_batches[[backend_name]] <- l1_exec
    }

    run_l1 <- TRUE
  } else {
    run_l1 <- FALSE
  }

  # If L1 is skipped but later levels will run, split the shared cache for each backend.
  if (isFALSE(run_l1) && is.null(split_cache_batch) && length(execution_backend_names) > 0L &&
      (!is.null(model_list$l2_model_names) || !is.null(model_list$l3_model_names))) {
    split_cache_batch <- f_batch$copy(
      job_name = "split_backend_caches", n_cpus = 1,
      wall_time = "0:10:00",
      r_code = c(
        sprintf("shared_cache <- '%s'", gpa_cache),
        "if (!file.exists(shared_cache)) {",
        "  stop('Shared GPA cache not found: ', shared_cache, '. The finalize job may have failed.')",
        "}",
        sprintf(
          "backend_cache_paths <- c(%s)",
          paste(shQuote(unname(backend_cache_map)), collapse = ", ")
        ),
        "for (path in backend_cache_paths) {",
        "  success <- file.copy(shared_cache, path, overwrite = TRUE)",
        "  if (!isTRUE(success)) {",
        "    stop('Failed to copy shared cache to backend cache: ', path)",
        "  }",
        "}"
      )
    )
    split_cache_batch$depends_on_parents <- if (isTRUE(run_finalize)) "finalize_configuration" else NULL
  }


  # todo
  # gpa <- verify_lv1_runs(gpa)

  # Standalone L2 execution is only needed when a resolved producer emits the required
  # subject/session contrasts at level 2.
  l2_execution_backends <- l2_backend_names[vapply(
    l2_backend_names,
    function(backend_name) {
      backend <- resolved_backends[[backend_name]]
      !is.null(backend) && backend_runs_level(backend, 2L)
    },
    logical(1)
  )]
  can_run_l2 <- isTRUE(gpa$multi_run) && length(l2_execution_backends) > 0L
  if (!is.null(model_list$l2_model_names) && isFALSE(can_run_l2)) {
    reasons <- c()
    if (!isTRUE(gpa$multi_run)) reasons <- c(reasons, "dataset is not multi-run")
    if (length(l2_execution_backends) == 0L) reasons <- c(reasons, "no configured backend provides standalone level 2 execution")
    spm_note <- if ("spm" %in% execution_backend_names) {
      "SPM does not run standalone L2 jobs; its L2 model spec is projected during SPM L1 setup."
    } else {
      NULL
    }
    lg$warn(
      "L2 models were requested but cannot be run as standalone jobs: %s%s",
      if (length(reasons) > 0L) paste(reasons, collapse = "; ") else "unknown reason",
      if (!is.null(spm_note)) paste0(". ", spm_note) else "."
    )
  }
  l2_execute_batches <- list()
  if (isTRUE(can_run_l2) && (!is.null(model_list$l2_model_names) || nrow(l3_level2_requirements) > 0L)) {

    for (backend_name in l2_execution_backends) {
      spec <- backend_specs[[backend_name]]
      run_fn_name <- if (is.null(spec)) NULL else spec$l2_run
      if (is.null(run_fn_name) || identical(run_fn_name, "__not_implemented__")) next
      if (!is.character(run_fn_name)) {
        lg$warn("Skipping backend '%s' L2 runner because l2_run is not a character function name.", backend_name)
        next
      }

      l2_model_subset <- l2_models_by_backend[[backend_name]] %||% model_list$l2_model_names

      l2_batch <- f_batch$copy(
        job_name = paste0("setup_run_l2_", backend_name), n_cpus = gpa$parallel$l2_setup_cores,
        wall_time = gpa$parallel$l2_setup_run_time,
        r_code = c(
          sprintf(
            "gpa <- setup_l2_models(gpa, l1_model_names=%s, l2_model_names=%s, backend = '%s')",
            paste(deparse(model_list$l1_model_names), collapse = ""),
            paste(deparse(l2_model_subset), collapse = ""),
            backend_name
          ),
          sprintf(
            "child_job_ids <- %s(gpa, level = 2L, model_names=%s)",
            run_fn_name, paste(deparse(l2_model_subset), collapse = "")
          )
        )
      )

      if (isTRUE(run_l1)) {
        l2_batch$depends_on_parents <- paste0("run_l1_", backend_name)
      } else if (!is.null(split_cache_batch)) {
        l2_batch$depends_on_parents <- "split_backend_caches"
      } else if (isTRUE(run_finalize)) {
        l2_batch$depends_on_parents <- "finalize_configuration"
      } else {
        l2_batch$depends_on_parents <- NULL
      }
      l2_batch$wait_for_children <- TRUE # need to wait for l2 jobs to complete before moving to l3
      if (!is.null(backend_cache_map[[backend_name]])) {
        l2_batch$input_rdata_file <- backend_cache_map[[backend_name]]
        l2_batch$output_rdata_file <- backend_cache_map[[backend_name]]
      }
      l2_execute_batches[[backend_name]] <- l2_batch
    }

    run_l2 <- TRUE
  } else {
    run_l2 <- FALSE
  }

  sync_backend_caches_batches <- list()
  if (nrow(l3_cross_backend_requirements) > 0L) {
    for (producer_level in sort(unique(l3_cross_backend_requirements$producer_level))) {
      sync_req <- unique(l3_cross_backend_requirements[
        l3_cross_backend_requirements$producer_level == producer_level,
        ,
        drop = FALSE
      ])
      if (nrow(sync_req) == 0L) next

      sync_lines <- c()
      for (ii in seq_len(nrow(sync_req))) {
        source_backend <- sync_req$producer_backend[ii]
        dest_backend <- sync_req$execution_backend[ii]
        if (is.null(backend_cache_map[[source_backend]]) || is.null(backend_cache_map[[dest_backend]])) next
        sync_lines <- c(
          sync_lines,
          sprintf("src <- %s", shQuote(backend_cache_map[[source_backend]])),
          sprintf("dest <- %s", shQuote(backend_cache_map[[dest_backend]])),
          "if (!file.exists(src)) stop('Required backend cache not found: ', src)",
          "success <- file.copy(src, dest, overwrite = TRUE)",
          "if (!isTRUE(success)) stop('Failed to sync backend cache from ', src, ' to ', dest)"
        )
      }
      if (length(sync_lines) == 0L) next

      job_name <- paste0("sync_l", producer_level, "_backend_caches")
      sync_job <- f_batch$copy(
        job_name = job_name,
        n_cpus = 1,
        wall_time = "0:10:00",
        r_code = sync_lines
      )

      if (producer_level == 1L) {
        parent_jobs <- unique(vapply(sync_req$producer_backend, function(source_backend) {
          if (!is.null(l1_execute_batches[[source_backend]])) {
            paste0("run_l1_", source_backend)
          } else if (!is.null(split_cache_batch)) {
            "split_backend_caches"
          } else if (!is.null(l1_setup_batch)) {
            "setup_l1"
          } else if (isTRUE(run_finalize)) {
            "finalize_configuration"
          } else {
            NA_character_
          }
        }, character(1)))
        parent_jobs <- parent_jobs[!is.na(parent_jobs) & nzchar(parent_jobs)]
        sync_job$depends_on_parents <- if (length(parent_jobs) > 0L) parent_jobs else NULL
      } else if (producer_level == 2L) {
        l2_job_names <- paste0("setup_run_l2_", names(l2_execute_batches))
        sync_job$depends_on_parents <- if (length(l2_job_names) > 0L) l2_job_names else NULL
      }

      sync_backend_caches_batches[[as.character(producer_level)]] <- sync_job
    }
  }

  
  if (!is.null(model_list$l3_model_names)) {
    for (backend_name in l3_backend_names) {
      spec <- backend_specs[[backend_name]]
      run_fn_name <- if (is.null(spec)) NULL else spec$l3_run
      if (is.null(run_fn_name) || identical(run_fn_name, "__not_implemented__")) next
      if (!is.character(run_fn_name)) {
        lg$warn("Skipping backend '%s' L3 runner because l3_run is not a character function name.", backend_name)
        next
      }
      if (!is.null(l3_glm_backends[[backend_name]]$l3_run) &&
        isTRUE(attr(l3_glm_backends[[backend_name]]$l3_run, "glm_backend_not_implemented"))) {
        next
      }

      l1_parent <- if (!is.null(l1_execute_batches[[backend_name]])) {
        paste0("run_l1_", backend_name)
      } else if (!is.null(split_cache_batch)) {
        "split_backend_caches"
      } else if (!is.null(l1_setup_batch)) {
        "setup_l1"
      } else if (isTRUE(run_finalize)) {
        "finalize_configuration"
      } else {
        NULL
      }

      producer_reqs <- if (nrow(l3_requirement_df) > 0L) {
        unique(l3_requirement_df[
          l3_requirement_df$execution_backend == backend_name,
          c("producer_backend", "producer_level"),
          drop = FALSE
        ])
      } else {
        data.frame(
          producer_backend = character(0),
          producer_level = integer(0),
          stringsAsFactors = FALSE
        )
      }

      l3_exec <- f_batch$copy(
        job_name = paste0("setup_run_l3_", backend_name), n_cpus = gpa$parallel$l2_setup_cores,
        wall_time = gpa$parallel$l3_setup_run_time,
        r_code = c(
          sprintf(
            "gpa <- setup_l3_models(gpa, l1_model_names=%s, l2_model_names=%s, l3_model_names=%s, backend = '%s')",
            paste(deparse(model_list$l1_model_names), collapse = ""),
            paste(deparse(model_list$l2_model_names), collapse = ""),
            paste(deparse(l3_models_by_backend[[backend_name]]), collapse = ""),
            backend_name
          ),
          "child_job_ids <- c()",
          sprintf(
            "child_job_ids <- c(child_job_ids, %s(gpa, level = 3L, model_names=%s))",
            run_fn_name, paste(deparse(l3_models_by_backend[[backend_name]]), collapse = "")
          )
        ),
        post_children_r_code = c(
          sprintf("gpa <- refresh_glm_status(gpa, level = 3L, glm_software = '%s')", backend_name)
        )
      )

      if (nrow(producer_reqs) > 0L) {
        producer_parents <- c()
        for (ii in seq_len(nrow(producer_reqs))) {
          producer_backend <- producer_reqs$producer_backend[ii]
          producer_level <- producer_reqs$producer_level[ii]
          sync_job <- sync_backend_caches_batches[[as.character(producer_level)]]

          if (!is.null(sync_job) && !identical(producer_backend, backend_name)) {
            producer_parents <- c(producer_parents, sync_job$job_name)
            next
          }

          if (producer_level == 2L) {
            l2_parent <- paste0("setup_run_l2_", producer_backend)
            if (!is.null(l2_execute_batches[[producer_backend]])) {
              producer_parents <- c(producer_parents, l2_parent)
            } else if (isTRUE(run_l2)) {
              producer_parents <- c(
                producer_parents,
                paste0("setup_run_l2_", names(l2_execute_batches))
              )
            }
          } else if (producer_level == 1L) {
            if (!is.null(l1_execute_batches[[producer_backend]])) {
              producer_parents <- c(producer_parents, paste0("run_l1_", producer_backend))
            } else {
              producer_parents <- c(producer_parents, l1_parent)
            }
          }
        }
        producer_parents <- unique(producer_parents[!is.na(producer_parents) & nzchar(producer_parents)])
        l3_exec$depends_on_parents <- if (length(producer_parents) > 0L) producer_parents else l1_parent
      } else {
        l3_exec$depends_on_parents <- l1_parent
      }
      l3_exec$wait_for_children <- TRUE
      if (!is.null(backend_cache_map[[backend_name]])) {
        l3_exec$input_rdata_file <- backend_cache_map[[backend_name]]
        l3_exec$output_rdata_file <- backend_cache_map[[backend_name]]
      }
      l3_execute_batches[[backend_name]] <- l3_exec
    }

    run_l3 <- length(l3_execute_batches) > 0L
  } else {
    run_l3 <- FALSE
  }

  # cleanup step: refresh l3 feat status and copy gpa back to main directory
  backend_cache_vec <- NULL
  if (length(backend_cache_map) > 0L) {
    backend_cache_vec <- paste(
      sprintf("%s=%s", names(backend_cache_map), shQuote(unname(backend_cache_map))),
      collapse = ", "
    )
  }

  cleanup_code <- c()
  if (!is.null(backend_cache_vec) && nzchar(backend_cache_vec)) {
    cleanup_code <- c(cleanup_code, sprintf("backend_cache_paths <- c(%s)", backend_cache_vec))
    cleanup_code <- c(cleanup_code, "gpa <- cleanup_glm_pipeline(gpa, backend_cache_paths = backend_cache_paths)")
  } else {
    cleanup_code <- c(cleanup_code, "gpa <- cleanup_glm_pipeline(gpa)")
  }

  cleanup_batch <- f_batch$copy(
    job_name = "cleanup_glm", n_cpus = 1,
    wall_time = "30:00", # 30 minutes should be plenty
    r_code = cleanup_code
  )

  if (length(l3_execute_batches) > 0L) {
    cleanup_batch$depends_on_parents <- vapply(l3_execute_batches, function(x) x$job_name, character(1))
  } else if (length(l2_execute_batches) > 0L) {
    cleanup_batch$depends_on_parents <- vapply(l2_execute_batches, function(x) x$job_name, character(1))
  } else if (length(l1_execute_batches) > 0L) {
    cleanup_batch$depends_on_parents <- vapply(l1_execute_batches, function(x) x$job_name, character(1))
  } else {
    cleanup_batch$depends_on_parents <- if (isTRUE(run_finalize)) "finalize_configuration" else NULL
  }
  
  run_cleanup <- TRUE # always TRUE for now

  if (isTRUE(run_finalize)) {
    glm_batch <- R_batch_sequence$new(
      joblist = c(
        list(f_batch, l1_setup_batch, split_cache_batch),
        l1_execute_batches,
        l2_execute_batches,
        sync_backend_caches_batches,
        l3_execute_batches,
        list(cleanup_batch)
      ),
      sequence_id = batch_id
    )
  } else {
    glm_batch <- R_batch_sequence$new(
      joblist = c(
        list(l1_setup_batch, split_cache_batch),
        l1_execute_batches,
        l2_execute_batches,
        sync_backend_caches_batches,
        l3_execute_batches,
        list(cleanup_batch)
      ),
      sequence_id = batch_id
    )
  }
  glm_batch$submit()
  
  l1_exec_flags <- if (length(l1_execute_batches) > 0L) {
    setNames(lapply(l1_execute_batches, function(x) !is.null(x)), paste0("l1_", names(l1_execute_batches)))
  } else {
    list()
  }
  l3_exec_flags <- if (length(l3_execute_batches) > 0L) {
    setNames(lapply(l3_execute_batches, function(x) !is.null(x)), paste0("l3_", names(l3_execute_batches)))
  } else {
    list()
  }

  job_sequence <- c(
    list(
      finalize = run_finalize,
      l1 = run_l1,
      split_backend_caches = !is.null(split_cache_batch),
      l2 = run_l2,
      l3 = run_l3,
      cleanup = run_cleanup
    ),
    l1_exec_flags,
    l3_exec_flags
  ) # logicals of which jobs were submitted
  update_project_config(gpa = gpa, job_sequence = job_sequence, 
                        sequence_id = batch_id, batch_directory) # update config file with this run's details
}

#' helper function to guide user through process of choosing which models to run in GLM pipeline
#'
#' @param l1_model_names a character vector of level 1 model names
#' @param l2_model_names a character vector of level 2 model names
#' @param l3_model_names a character vector of level 3 model names
#' 
#' @return a named list containing all models that were selected along with additional information
#'   about whether to rerun existing models
#' @keywords internal
choose_glm_set <- function(gpa, l1_model_names=NULL, l2_model_names=NULL, l3_model_names=NULL, lg=NULL) {
  checkmate::assert_class(lg, "Logger")

  if (is.null(l1_model_names)) {
    lg$debug("l1_model_names was NULL. Defaulting to running all l1 models")
    l1_model_names <- "all"
  }

  if (is.null(l2_model_names)) {
    lg$debug("l2_model_names was NULL. Defaulting to running all l2 models")
    l2_model_names <- "all"
  }

  if (is.null(l3_model_names)) {
    lg$debug("l3_model_names was NULL. Defaulting to running all l3 models")
    l3_model_names <- "all"
  }

  m_list <- named_list(l1_model_names, l2_model_names, l3_model_names)

  for (nn in names(m_list)) {
    # enforce that all, none, and prompt must be singleton arguments
    if (any(m_list[[nn]] %in% c("all", "none", "prompt")) && length(m_list[[nn]]) > 1L) {
      msg <- sprintf(
        "Argument %s has value 'all', 'none', or 'prompt'. These must be passed alone, not with other model names.", 
        m_list[[nn]]
      )
      lg$error(msg)
      stop(msg)
    }
  }

  m_string <- function(str) { if (is.null(str)) "none" else str }
  is_prompt <- function(x) {
    !is.null(x) && length(x) > 0L && identical(x[1L], "prompt")
  }

  prompt_selection_requested <- is_prompt(l1_model_names) ||
    is_prompt(l3_model_names) ||
    (isTRUE(gpa$multi_run) && is_prompt(l2_model_names))

  models_specified <- FALSE
  while (isFALSE(models_specified)) {
    l1_model_names <- choose_glm_models(gpa, l1_model_names, level = 1)
    if (isTRUE(gpa$multi_run)) l2_model_names <- choose_glm_models(gpa, l2_model_names, level = 2)
    l3_model_names <- choose_glm_models(gpa, l3_model_names, level = 3)

    cat("\nGLM models to run:\n------------------\n\n")
    cat("Level 1: ", paste(m_string(l1_model_names), collapse = ", "), "\n")
    if (isTRUE(gpa$multi_run)) cat("Level 2: ", paste(m_string(l2_model_names), collapse = ", "), "\n")
    cat("Level 3: ", paste(m_string(l3_model_names), collapse = ", "), "\n")
    cat("\n------------------\n\n")

    if (isTRUE(gpa$multi_run)) {
      menu_options <- c(
        "Yes (run)", "No, respecify level 1 models",
        "No, respecify level 2 models", "No, respecify level 3 models", "Cancel"
      )
    } else {
      menu_options <- c(
        "Yes (run)", "No, respecify level 1 models",
        "No, respecify level 3 models", "Cancel"
      )
    }

    if (!interactive()) {
      lg$info("Non-interactive session detected. Proceeding with selected models without prompt.")
      models_specified <- TRUE
      next
    }

    if (isFALSE(prompt_selection_requested)) {
      lg$info("Model names were provided as input arguments. Proceeding without confirmation prompt.")
      models_specified <- TRUE
      next
    }

    cat("\n")
    respecify <- menu(menu_options, title = "Do you want to continue (this submits the models for execution)?")
    if (respecify == 0L) {
      respecify <- "Cancel"
    } else {
      respecify <- menu_options[respecify]
    }

    if (respecify == "Yes (run)") {
      models_specified <- TRUE
    } else if (respecify == "No, respecify level 1 models") {
      l1_model_names <- "prompt"
    } else if (respecify == "No, respecify level 2 models") {
      l2_model_names <- "prompt"
    } else if (respecify == "No, respecify level 3 models") {
      l3_model_names <- "prompt"
    } else if (respecify == "Cancel") {
      return(invisible(NULL))
    }
  }

  return(named_list(l1_model_names, l2_model_names, l3_model_names))

}
