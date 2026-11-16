#' This function generates the inputs for level 3 analyses, where multi-subject data are analyzed
#'   in group analyses.
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing analysis speecification
#' @param l3_model_names a subset of L3 models to be setup by this function. If not specified,
#'   all models in gpa$l2_models will be included
#' @param l2_model_names a subset of L2 models to be setup by this function. If not specified,
#'   all models in gpa$l2_models will be included
#' @param l1_model_names a subset of L1 models to be passed to L2 by this function. If not
#'   specified, all models in gpa$l1_models will be included
#' @param backend optional backend filter (e.g., "spm"). If supplied, only those backends
#'   will be processed.
#'
#' @details
#'   This function will setup FSL level 2 (subject) .fsf files for all combinations of
#'   \code{l2_model_names} and \code{l1_model_names}.
#'
#' FSL 2-level versus 3-level setup
#'
#' 2-level setup (one run per subject)
#'   - Pass L1 .feat folders as input to L3 .fsf setup
#'   - In this approach, the copes in the .fsf pertain to the L1 cope numbers
#'   - Requires one .fsf per L3 model
#'
#' 3-level setup (multiple runs per subject, combined at L2)
#'   - Pass individual cope*.feat folders within subject .gfeat folders
#'   - The folder cope numbers pertain to L1 copes
#'   - The cope*.nii.gz in the cope*.feat subfolders pertain to the L2 contrasts
#'   - Requires one .fsf per L1 cope x L3 model combination
#'   - Example: FSL_L2.gfeat/cope3.feat/stats/cope1.nii.gz
#'      ==> cope3 is the third contrast in the L1 feat model
#'      ==> cope1 is the first contrast in the L2 feat model
#'
#' @author Michael Hallquist
#' @importFrom checkmate assert_class assert_character assert_data_frame
#' @importFrom data.table is.data.table
#' @importFrom lgr get_logger
#' @importFrom stats setNames
#' @importFrom utils capture.output
#' @export
setup_l3_models <- function(gpa, l3_model_names = NULL, l2_model_names = NULL, l1_model_names = NULL, backend = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_character(backend, null.ok = TRUE)
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_class(gpa$l3_models, "hi_model_set")
  checkmate::assert_class(gpa$l1_models, "l1_model_set")

  if (is.null(l3_model_names)) l3_model_names <- names(gpa$l3_models$models)
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)

  backend_specs <- gpa$glm_backend_specs
  if (is.null(backend_specs)) backend_specs <- default_glm_backend_specs()
  gpa$glm_backend_specs <- backend_specs
  l3_model_backend_map <- get_effective_model_backends(gpa, level = 3L, model_names = l3_model_names)
  l3_producer_backend_map <- get_effective_model_backends(gpa, level = 3L, model_names = l3_model_names, type = "producer")
  validate_l3_backend_resolution(
    gpa = gpa,
    l3_model_names = l3_model_names,
    execution_backend_map = l3_model_backend_map,
    producer_backend_map = l3_producer_backend_map,
    requested_backends = backend
  )
  requested_l3_backends <- normalize_backend_strings(unlist(l3_model_backend_map, use.names = FALSE))
  resolved_backends <- resolve_glm_backends(backend_specs)
  glm_backends <- resolved_backends[intersect(requested_l3_backends, names(resolved_backends))]
  backend_names <- names(glm_backends)
  if (!is.null(backend)) {
    backend <- tolower(backend)
    checkmate::assert_subset(backend, c("fsl", "spm", "afni"))
    missing_backends <- setdiff(backend, backend_names)
    if (length(missing_backends) > 0L) {
      warning(
        sprintf("Backends not available for setup_l3_models: %s", paste(missing_backends, collapse = ", "))
      )
    }
    backend_names <- intersect(backend_names, backend)
  }

  level_backends_orig <- gpa$level_backends
  if (!is.null(backend)) {
    level_backends <- normalize_level_backend_config(gpa)
    level_backends$l3 <- backend_names
    gpa$level_backends <- level_backends
    on.exit({ gpa$level_backends <- level_backends_orig }, add = TRUE)
  }

  use_fsl <- "fsl" %in% backend_names
  use_spm <- "spm" %in% backend_names
  use_afni <- "afni" %in% backend_names
  l3_models_in_scope <- l3_model_names[vapply(
    l3_model_names,
    function(model_name) any(backend_names %in% normalize_backend_strings(l3_model_backend_map[[model_name]])),
    logical(1)
  )]
  l3_requirement_df <- resolve_model_l3_requirements(
    gpa = gpa,
    l3_model_names = l3_models_in_scope,
    execution_backend_map = l3_model_backend_map,
    producer_backend_map = l3_producer_backend_map,
    specs = backend_specs,
    multi_run = isTRUE(gpa$multi_run)
  )
  producer_stage_requirements <- if (nrow(l3_requirement_df) > 0L) {
    unique(l3_requirement_df[, c("producer_backend", "producer_level"), drop = FALSE])
  } else {
    data.frame(
      producer_backend = character(0),
      producer_level = integer(0),
      stringsAsFactors = FALSE
    )
  }
  uses_level2_inputs <- nrow(producer_stage_requirements) > 0L &&
    any(producer_stage_requirements$producer_level == 2L)

  # Validate model name subsets
  checkmate::assert_subset(l3_model_names, names(gpa$l3_models$models))
  checkmate::assert_subset(l1_model_names, names(gpa$l1_models$models))

  # full 3-level analysis (runs, subject, sample)
  if (isTRUE(gpa$multi_run) && isTRUE(uses_level2_inputs)) {
    checkmate::assert_class(gpa$l2_models, "hi_model_set")
    checkmate::assert_subset(l2_model_names, names(gpa$l2_models$models))

    # if no l2 model subset is requested, output all models
    if (is.null(l2_model_names)) l2_model_names <- names(gpa$l2_models$models)
  } else if (isTRUE(gpa$multi_run)) {
    l2_model_names <- NULL
  }

  lg <- lgr::get_logger("glm_pipeline/l3_setup")
  lg$set_threshold(gpa$lgr_threshold)

  # normalize model signatures now in case this function is called before full finalization
  gpa <- normalize_longitudinal_model_signatures(gpa, lg)
  gpa <- refresh_l1_cope_names(gpa, lg = lg)

  add_log_suffix <- function(path, suffix) {
    if (is.null(path) || is.null(suffix) || !nzchar(suffix)) return(path)
    ext <- tools::file_ext(path)
    if (nzchar(ext)) {
      base <- sub(paste0("\\.", ext, "$"), "", path)
      return(paste0(base, "_", suffix, ".", ext))
    }
    paste0(path, "_", suffix)
  }
  log_suffix <- if (!is.null(backend)) paste(backend, collapse = "_") else NULL
  setup_l3_log_txt <- add_log_suffix(gpa$output_locations$setup_l3_log_txt, log_suffix)
  setup_l3_log_json <- add_log_suffix(gpa$output_locations$setup_l3_log_json, log_suffix)

  add_base_logger_appenders(
    lg = lg,
    gpa = gpa,
    log_txt_path = setup_l3_log_txt,
    log_json_path = setup_l3_log_json,
    txt_appender_name = "setup_l3_log_txt",
    json_appender_name = "setup_l3_log_json"
  )

  lg$debug("In setup_l3_models, setting up the following L3 models:")
  lg$debug("L3 model: %s", l3_model_names)
  if (isTRUE(uses_level2_inputs)) {
    lg$debug("In setup_l3_models, passing the following L2 models to L3:")
    lg$debug("L2 model: %s", l2_model_names)
  }
  lg$debug("In setup_l3_models, passing the following L1 models to L3:")
  lg$debug("L1 model: %s", l1_model_names)
  if (nrow(producer_stage_requirements) > 0L) {
    for (ii in seq_len(nrow(producer_stage_requirements))) {
      lg$info(
        "L3 prerequisites resolved to backend '%s' outputs at level %d.",
        producer_stage_requirements$producer_backend[ii],
        producer_stage_requirements$producer_level[ii]
      )
    }
  }

  l2_l3_pairs <- NULL
  if (isTRUE(uses_level2_inputs)) {
    l2_l3_pairs <- resolve_l2_l3_compatible_pairs(
      gpa,
      l2_model_names = l2_model_names,
      l3_model_names = l3_model_names,
      lg = lg
    )
    if (nrow(l2_l3_pairs) == 0L) {
      pair_df <- enumerate_l2_l3_signature_pairs(
        gpa,
        l2_model_names = l2_model_names,
        l3_model_names = l3_model_names,
        lg = lg
      )
      detail <- format_l2_l3_incompatibilities(pair_df, max_items = 8L)
      msg <- "No L2+L3 model combinations remain after applying L2/L3 signature compatibility rules."
      if (nzchar(detail)) msg <- paste(msg, detail)
      lg$warn(msg)
      warning(msg)
      gpa$l3_setup_status <- list(
        success = FALSE,
        reason = msg,
        timestamp = Sys.time()
      )
      return(gpa)
    }
  }

  if (!isTRUE(use_fsl) && !isTRUE(use_spm) && !isTRUE(use_afni)) {
    lg$info("No supported L3 backends requested in setup_l3_models. Skipping.")
    return(gpa)
  }

  # handle refresh of feat status for lower-level models
  # Resolve and refresh the concrete producer stage(s) that feed the requested L3 models.
  if (nrow(producer_stage_requirements) == 0L) {
    gpa <- refresh_glm_status(gpa, level = 1L, lg = lg)
  } else {
    for (ii in seq_len(nrow(producer_stage_requirements))) {
      gpa <- refresh_glm_status(
        gpa,
        level = producer_stage_requirements$producer_level[ii],
        lg = lg,
        glm_software = producer_stage_requirements$producer_backend[ii]
      )
    }
  }

  excluded_runs <- gpa$run_data %>%
    dplyr::select("id", "session", "run_number", "exclude_run", "exclude_subject") %>%
    dplyr::filter(.data$exclude_run == TRUE | .data$exclude_subject == TRUE)

  if (nrow(excluded_runs) > 1L) {
    lg$info("In setup_l3_models, the following runs will be excluded from L3 modeling: ")
    lg$info(
      "  subject: %s, session: %s, run_number: %s",
      excluded_runs$id, excluded_runs$session, excluded_runs$run_number
    )
  }

  if (nrow(producer_stage_requirements) == 0L) {
    enforce_glms_complete(gpa, level = 1L, lg)
  } else {
    for (ii in seq_len(nrow(producer_stage_requirements))) {
      enforce_glms_complete(
        gpa,
        level = producer_stage_requirements$producer_level[ii],
        lg = lg,
        glm_software = producer_stage_requirements$producer_backend[ii]
      )
    }
  }

  if (isTRUE(uses_level2_inputs)) {
    lg$info("In setup_l3_models, using producer outputs from level 2 to supply subject/session contrasts to L3.")
  } else {
    if (isTRUE(gpa$multi_run)) {
      lg$info("In setup_l3_models, using multi-run inputs produced without a standalone L2 step.")
    } else {
      lg$info("In setup_l3_models, using a single run 2-level setup with subjects (l1), sample (l3)")
    }
    # only retain good runs and subjects
    run_data <- gpa$run_data %>%
      dplyr::filter(.data$exclude_run == FALSE & .data$exclude_subject == FALSE)

    if (nrow(run_data) == 0L) {
      msg <- "In setup_l3_models, no runs survived the exclude_subject and exclude_run step."
      lg$warn(msg)
      warning(msg)
      gpa$l3_setup_status <- list(
        success = FALSE,
        reason = msg,
        n_excluded = nrow(excluded_runs),
        timestamp = Sys.time()
      )
      return(gpa)
    }
  }

  # subjects and sessions to run at l3
  subj_df <- gpa$subject_data %>%
    dplyr::filter(.data$exclude_subject == FALSE) %>%
    dplyr::select("id", "session")

  # loop over requested backends and setup all requested combinations of L1/L2/L3 models
  backend_results <- list()
  metadata <- NULL
  for (backend_name in backend_names) {
    backend_l3_model_names <- l3_model_names[vapply(
      l3_model_names,
      function(model_name) backend_name %in% normalize_backend_strings(l3_model_backend_map[[model_name]]),
      logical(1)
    )]
    if (length(backend_l3_model_names) == 0L) next

    helper <- get_l3_backend_helper(backend_name)
    if (is.null(helper)) {
      lg$warn("No L3 helper registered for backend '%s'. Skipping.", backend_name)
      next
    }

    res <- helper(
      gpa = gpa,
      backend = glm_backends[[backend_name]],
      lg = lg,
      l1_model_names = l1_model_names,
      l2_model_names = l2_model_names,
      l3_model_names = backend_l3_model_names,
      l2_l3_pairs = l2_l3_pairs,
      subj_df = subj_df,
      requires_l2 = uses_level2_inputs
    )
    backend_results[[backend_name]] <- res
    if (!is.null(res$metadata)) {
      metadata <- if (is.null(metadata)) res$metadata else dplyr::bind_rows(metadata, res$metadata)
    }
  }

  all_subj_l3_combined <- c(
    list(metadata = metadata),
    setNames(lapply(backend_results, `[[`, "data"), names(backend_results))
  )

  class(all_subj_l3_combined) <- c("l3_setup", "list")

  # Combine the generated l3_model_setup with any existing l3_model_setup information
  # This is important so that if a user requests a model subset (using input arguments, and from run_glm_pipeline),
  # we don't clear out information about many other models that may have already completed
  existing_l3 <- if (!is.null(gpa$l3_model_setup) && inherits(gpa$l3_model_setup, "l3_setup")) gpa$l3_model_setup else NULL

  if (!is.null(existing_l3)) {
    if (is.null(all_subj_l3_combined$metadata) && !is.null(existing_l3$metadata)) {
      all_subj_l3_combined$metadata <- existing_l3$metadata
    } else if (!is.null(all_subj_l3_combined$metadata) && !is.null(existing_l3$metadata)) {
      meta_keys <- intersect(
        c("session", "l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model"),
        intersect(names(existing_l3$metadata), names(all_subj_l3_combined$metadata))
      )
      if (length(meta_keys) > 0L) {
        kept_old <- dplyr::anti_join(existing_l3$metadata, all_subj_l3_combined$metadata, by = meta_keys)
        all_subj_l3_combined$metadata <- dplyr::bind_rows(kept_old, all_subj_l3_combined$metadata)
      }
    }

    backend_union <- union(setdiff(names(existing_l3), "metadata"), names(backend_results))
    for (backend_name in backend_union) {
      new_df <- all_subj_l3_combined[[backend_name]]
      existing_df <- existing_l3[[backend_name]]
      res <- backend_results[[backend_name]]
      id_cols <- if (!is.null(res)) res$id_cols else NULL

      if (!is.null(existing_df) && !is.null(new_df) && is.data.frame(new_df) &&
        nrow(new_df) > 0L && !is.null(id_cols)) {
        all_subj_l3_combined[[backend_name]] <- update_df(
          current = existing_df, new = new_df, id_cols = id_cols
        )
      } else if (is.null(new_df) || (is.data.frame(new_df) && nrow(new_df) == 0L)) {
        all_subj_l3_combined[[backend_name]] <- existing_df
      }
    }
  }

  # append l3 setup to gpa
  gpa$l3_model_setup <- all_subj_l3_combined

  # refresh l3 model status in $l3_model_setup
  gpa <- refresh_glm_status(gpa, level = 3L, lg = lg)

  # the expressions for l3 output locations should not generate any duplicate fsfs
  if (isTRUE(use_fsl) && !is.null(gpa$l3_model_setup$fsl)) {
    dupe_fsfs <- duplicated(gpa$l3_model_setup$fsl$feat_fsf)
    if (any(dupe_fsfs, na.rm = TRUE)) {
      lg$warn("There are duplicate fsfs in gpa$l3_model_setup$fsl.")
      lg$warn("This suggests a problem with gpa$output_locations$feat_l3_fsf and gpa$output_locations$feat_l3_directory settings.")
      lg$warn("gpa$output_locations$feat_l3_fsf: %s", gpa$output_locations$feat_l3_fsf)
      lg$warn("gpa$output_locations$feat_l3_directory: %s", gpa$output_locations$feat_l3_directory)
      lg$debug("%s", capture.output(print(gpa$l3_model_setup$fsl[dupe_fsfs, ])))
    }
  }

  n_models <- 0L
  for (backend_name in names(backend_results)) {
    df <- all_subj_l3_combined[[backend_name]]
    if (!is.null(df)) n_models <- n_models + nrow(df)
  }

  # record successful completion status
  gpa$l3_setup_status <- list(
    success = TRUE,
    n_models = n_models,
    backend = backend,
    timestamp = Sys.time()
  )

  return(gpa)
}


get_l3_backend_helper <- function(backend_name) {
  backend_name <- tolower(backend_name)
  switch(
    backend_name,
    fsl = setup_l3_backend_fsl,
    spm = setup_l3_backend_spm,
    afni = afni_3dlmer_setup,
    NULL
  )
}

setup_l3_backend_fsl <- function(gpa, backend, lg, l1_model_names, l2_model_names, l3_model_names, l2_l3_pairs, subj_df, requires_l2) {
  if (is.null(backend)) {
    return(list(metadata = NULL, data = NULL, id_cols = NULL))
  }

  if (isTRUE(requires_l2)) {
    if (is.null(l2_l3_pairs)) {
      model_set <- expand.grid(
        l1_model = l1_model_names, l2_model = l2_model_names,
        l3_model = l3_model_names, stringsAsFactors = FALSE
      )
    } else {
      checkmate::assert_data_frame(l2_l3_pairs)
      checkmate::assert_subset(c("l2_model", "l3_model"), names(l2_l3_pairs))
      pair_cols <- c("l2_model", "l3_model")
      if ("l2_scope" %in% names(l2_l3_pairs)) pair_cols <- c(pair_cols, "l2_scope")
      if ("l3_input_mode" %in% names(l2_l3_pairs)) pair_cols <- c(pair_cols, "l3_input_mode")
      model_set <- merge(
        data.frame(l1_model = l1_model_names, stringsAsFactors = FALSE),
        l2_l3_pairs[, pair_cols, drop = FALSE],
        by = NULL
      )
    }
  } else {
    model_set <- expand.grid(l1_model = l1_model_names, l3_model = l3_model_names, stringsAsFactors = FALSE)
  }

  default_fsl_id_cols <- if (isTRUE(gpa$multi_run)) {
    c("session", "l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model")
  } else {
    c("session", "l1_model", "l1_cope_name", "l3_model")
  }

  if (nrow(model_set) == 0L) {
    lg$warn("No model combinations were generated for FSL L3 setup.")
    return(list(metadata = NULL, data = NULL, id_cols = default_fsl_id_cols))
  }

  l3_cope_config <- get_fsl_l3_model_df(gpa, model_set, subj_df)

  l3_cope_input_df <- l3_cope_config %>%
    dplyr::filter(.data$l3_cope_number == 1L) %>%
    dplyr::select(-c("l3_cope_number", "l3_cope_name"))

  to_run <- get_feat_l3_inputs(gpa, l3_cope_input_df, lg)

  if (length(to_run) == 0L) {
    lg$warn("No feat l3 inputs returned from get_feat_l3_inputs. Is it possible no lower-level analyses are complete?")
    return(list(metadata = l3_cope_config, data = NULL,
                id_cols = default_fsl_id_cols))
  }

  all_l3_list <- foreach(
    model_info = iter(to_run), .inorder = FALSE,
    .packages = c("fmri.pipeline", "dplyr", "data.table"), .export = c("lg", "gpa", "backend")
  ) %dopar% {
    model_info <- model_info # to avoid complaints about global variable binding in R CMD check

    l3_file_setup <- list(fsl = list(), spm = list(), afni = list())

    if (nrow(model_info) <= 3) {
      lg$warn(
        "Fewer than 4 complete feat input directories for l1 model %s, l2 model %s, l3 model %s",
        model_info$l1_model[1L],
        if ("l2_model" %in% names(model_info)) model_info$l2_model[1L] else "none",
        model_info$l3_model[1L]
      )
      l3_file_setup$fsl <- NULL
      return(l3_file_setup)
    }

    l3_file_setup$fsl <- tryCatch(backend$l3_setup(model_info, gpa = gpa),
      error = function(e) {
        lg$error(
          "Problem with fsl_l3_model. L1 Model: %s, L2 Model: %s, L3 model %s",
          model_info$l1_model[1L],
          if ("l2_model" %in% names(model_info)) model_info$l2_model[1L] else "none",
          model_info$l3_model[1L]
        )
        lg$error("Error message: %s", as.character(e))
        return(NULL)
      }
    )

    return(l3_file_setup)
  }

  fsl_df <- dplyr::bind_rows(lapply(all_l3_list, "[[", "fsl"))

  # Synchronize l3_cope_config metadata with actual post-refit contrasts (e.g., when pooled_sessions_subject_ev
  # suppresses EV_id* diagonals or prunes aliased between-subject contrasts)
  if (nrow(fsl_df) > 0L && "cope_list" %in% names(fsl_df)) {
    analysis_keys <- intersect(default_fsl_id_cols, names(fsl_df))
    refit_configs <- lapply(seq_len(nrow(fsl_df)), function(ii) {
      analysis <- fsl_df[ii, analysis_keys, drop = FALSE]
      # Pooled outputs use session=0, while input metadata retains real sessions.
      # Match those inputs across sessions; partitioned outputs match one session.
      join_keys <- analysis_keys
      if ("session" %in% join_keys && identical(as.character(analysis$session), "0")) {
        join_keys <- setdiff(join_keys, "session")
      }
      base <- l3_cope_config %>%
        dplyr::semi_join(analysis, by = join_keys) %>%
        dplyr::select(-c("l3_cope_number", "l3_cope_name")) %>%
        dplyr::distinct()
      cl <- fsl_df$cope_list[[ii]]
      if (is.null(cl)) {
        cl <- data.frame(l3_cope_number = integer(0), l3_cope_name = character(0))
      }
      list(base = base, metadata = dplyr::cross_join(base, cl))
    })
    updated_bases <- dplyr::bind_rows(lapply(refit_configs, "[[", "base"))
    keep_cfg <- dplyr::anti_join(l3_cope_config, updated_bases, by = names(updated_bases))
    l3_cope_config <- dplyr::bind_rows(
      keep_cfg, dplyr::bind_rows(lapply(refit_configs, "[[", "metadata"))
    )
    fsl_df$cope_list <- NULL
  }

  resolved_id_cols <- if (nrow(fsl_df) > 0L) intersect(default_fsl_id_cols, names(fsl_df)) else default_fsl_id_cols
  list(
    metadata = l3_cope_config,
    data = fsl_df,
    id_cols = resolved_id_cols
  )
}

setup_l3_backend_spm <- function(gpa, backend, lg, l1_model_names, l2_model_names, l3_model_names, l2_l3_pairs, subj_df, requires_l2) {
  if (is.null(backend)) {
    return(list(metadata = NULL, data = NULL, id_cols = NULL))
  }

  if (isTRUE(gpa$multi_run)) {
    lg$info(
      "In setup_l3_models, using SPM L1->L3 setup (no standalone L2 step; L2 specs are projected during SPM L1 setup)."
    )
  }

  spm_model_set <- expand.grid(
    l1_model = l1_model_names, l3_model = l3_model_names,
    stringsAsFactors = FALSE
  )

  spm_cope_config <- get_spm_l3_model_df(gpa, spm_model_set, subj_df)

  spm_cope_input_df <- spm_cope_config %>%
    dplyr::filter(.data$l3_cope_number == 1L) %>%
    dplyr::select(-c("l3_cope_number", "l3_cope_name"))

  spm_to_run <- get_spm_l3_inputs(gpa, spm_cope_input_df, lg)

  if (length(spm_to_run) == 0L) {
    lg$warn("No SPM L3 inputs returned. Is it possible no lower-level SPM analyses are complete?")
    return(list(metadata = spm_cope_config, data = NULL,
                id_cols = c("l1_model", "l1_cope_name", "l3_model")))
  }

  all_spm_list <- foreach(
    model_info = iter(spm_to_run), .inorder = FALSE,
    .packages = c("fmri.pipeline", "dplyr", "data.table"), .export = c("lg", "gpa", "backend")
  ) %dopar% {
    model_info <- model_info # to avoid complaints about global variable binding in R CMD check

    l3_file_setup <- list(fsl = list(), spm = list(), afni = list())

    # Validate minimum number of inputs (same check as FSL backend)
    if (nrow(model_info) <= 3) {
      lg$warn(
        "Fewer than 4 complete SPM input directories for l1 model %s, l3 model %s",
        model_info$l1_model[1L], model_info$l3_model[1L]
      )
      l3_file_setup$spm <- NULL
      return(l3_file_setup)
    }

    l3_file_setup$spm <- tryCatch(backend$l3_setup(model_info, gpa = gpa),
      error = function(e) {
        lg$error(
          "Problem with spm_l3_model. L1 Model: %s, L3 model %s",
          model_info$l1_model[1L], model_info$l3_model[1L]
        )
        lg$error("Error message: %s", as.character(e))
        return(NULL)
      }
    )

    return(l3_file_setup)
  }

  spm_df <- dplyr::bind_rows(lapply(all_spm_list, "[[", "spm"))
  list(
    metadata = spm_cope_config,
    data = spm_df,
    id_cols = c("l1_model", "l1_cope_name", "l3_model")
  )
}

#' Internal function to setup AFNI 3dLMEr longitudinal models
#'
#' @param gpa a glm_pipeline_arguments object
#' @param backend the afni backend specification
#' @param lg the current logger
#' @param l1_model_names L1 models to process
#' @param l2_model_names L2 models to process
#' @param l3_model_names L3 models to process
#' @param l2_l3_pairs data.frame of compatible L2/L3 model pairs
#' @param subj_df filtered subject data.frame
#' @param requires_l2 whether the resolved upstream producer emits required inputs at level 2
#'
#' @return a setup list with metadata and status data.frame
#' @keywords internal
afni_3dlmer_setup <- function(gpa, backend, lg, l1_model_names, l2_model_names, l3_model_names, l2_l3_pairs, subj_df, requires_l2) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_class(lg, "Logger")

  if (is.null(backend)) backend <- list(name = "afni")

  harvested <- list()
  producer_backend_map <- get_effective_model_backends(gpa, level = 3L, model_names = l3_model_names, type = "producer")
  for (l3_name in l3_model_names) {
    provider_backends <- normalize_backend_strings(producer_backend_map[[l3_name]])
    if (length(provider_backends) > 0L && any(provider_backends != "fsl")) {
      stop(
        sprintf(
          "AFNI 3dLMEr currently supports only FSL as an upstream producer backend. Model '%s' requested: %s",
          l3_name, paste(provider_backends, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    this_backend <- backend
    if (length(provider_backends) > 0L) {
      this_backend$l3_input_provider_backends <- provider_backends
    }

    harvested_one <- harvest_l3_inputs(
      gpa = gpa,
      l3_backend = this_backend,
      l1_model_names = l1_model_names,
      l2_model_names = l2_model_names,
      l3_model_names = l3_name,
      lg = lg
    )
    if (!is.null(harvested_one)) {
      harvested[[l3_name]] <- harvested_one[[l3_name]]
    }
  }

  if (is.null(harvested) || length(harvested) == 0) {
    lg$warn("No harvested inputs for AFNI 3dLMEr setup.")
    return(list(metadata = NULL, data = NULL, id_cols = NULL))
  }

  afni_setup_df <- list()

  for (l3_name in names(harvested)) {
    l3_obj <- gpa$l3_models$models[[l3_name]]
    if (!identical(l3_obj$l3_input_mode, "3dlmer")) next
    
    # harvested[[l3_name]] is a list of data.frames split by contrast
    for (con_name in names(harvested[[l3_name]])) {
      dt <- harvested[[l3_name]][[con_name]]

      analysis_ids <- dt %>%
        dplyr::select("id", "session") %>%
        dplyr::distinct()
      l3_fit <- respecify_l3_model(l3_obj, analysis_ids)

      # Drop harvested rows whose id/session were excluded due to missing L3 covariates in respecify_l3_model
      meta_keys <- l3_fit$metadata
      meta_keys$id <- as.character(meta_keys$id)
      dt$id <- as.character(dt$id)
      if (!identical(class(dt$session), class(meta_keys$session))) {
        dt$session <- as.integer(as.character(dt$session))
        meta_keys$session <- as.integer(as.character(meta_keys$session))
      }
      dt <- dplyr::semi_join(dt, meta_keys, by = c("id", "session"))
      if (nrow(dt) == 0L) {
        lg$warn("No complete covariate rows remain for AFNI 3dLMEr model '%s', contrast '%s'. Skipping.", l3_name, con_name)
        next
      }

      lmer_vars <- tryCatch(
        setdiff(all.vars(stats::as.formula(paste("~", sub("^~\\s*", "", l3_obj$lmer_formula)))), c("Subj", "id")),
        error = function(e) character(0)
      )
      model_vars <- unique(c(
        setdiff(l3_fit$model_variables, c("id", "session")),
        setdiff(lmer_vars, "session"),
        setdiff(l3_obj$random_effects_vars %||% character(0), c("id", "session", "Subj"))
      ))

      fit_model_data <- l3_fit$model_data
      transformed_session <- if ("session" %in% names(fit_model_data)) fit_model_data$session else NULL
      overlap_cols <- intersect(names(l3_fit$metadata), names(fit_model_data))
      if (length(overlap_cols) > 0L) {
        fit_model_data <- fit_model_data %>%
          dplyr::select(-dplyr::all_of(overlap_cols))
      }
      fit_model_data <- fit_model_data %>%
        dplyr::select(-dplyr::any_of("dummy"))
      subject_data <- cbind(meta_keys, fit_model_data)
      if (!is.null(transformed_session)) {
        subject_data$.transformed_session_internal <- transformed_session
      }

      datatable <- build_3dlmer_datatable(
        subject_data = subject_data,
        input_files = dt %>% dplyr::select("id", "session", "InputFile"),
        model_variables = if (".transformed_session_internal" %in% names(subject_data)) c(model_vars, ".transformed_session_internal") else model_vars
      )
      if (".transformed_session_internal" %in% names(datatable)) {
        if ("session" %in% lmer_vars) {
          datatable$session <- datatable$.transformed_session_internal
        }
        datatable$.transformed_session_internal <- NULL
      }
      validate_3dlmer_formula_datatable(
        model_formula = l3_obj$lmer_formula,
        datatable = datatable,
        context = sprintf("L3 model '%s', contrast '%s'", l3_name, con_name)
      )

      # Auto-detect qVars from the refit model data used to generate the AFNI table (including numeric session if used in lmer_formula).
      datatable_df <- as.data.frame(datatable, stringsAsFactors = FALSE)
      cand_qvars <- if ("session" %in% lmer_vars) unique(c("session", model_vars)) else model_vars
      qVars <- cand_qvars[vapply(datatable_df[, cand_qvars, drop = FALSE], function(col) is.numeric(col) && !is.factor(col), logical(1))]

      # Preserve user covariate_transform centering ("min", "max", "mean", "zscore") instead of letting 3dLMEr re-mean-center around sample mean
      qVarCenters <- NULL
      if (length(qVars) > 0L && !is.null(l3_obj$covariate_transform) && any(qVars %in% names(l3_obj$covariate_transform))) {
        qVarCenters <- vapply(qVars, function(qv) {
          if (qv %in% names(l3_obj$covariate_transform)) 0 else mean(as.numeric(datatable_df[[qv]]), na.rm = TRUE)
        }, numeric(1))
      }

      glt_codes <- emmeans_to_3dlmer_glt(
        l3_fit,
        datatable,
        qVars = qVars,
        raw_glt_codes = l3_obj$lmer_glt_codes,
        context = sprintf("L3 model '%s', contrast '%s'", l3_name, con_name)
      )

      l1_model <- dt$l1_model[1]
      l2_model <- if ("l2_model" %in% names(dt) && !is.null(dt$l2_model[1]) && !is.na(dt$l2_model[1])) dt$l2_model[1] else "none"
      l3_model <- l3_name
      l1_cope_name <- dt$l1_cope_name[1]
      l2_cope_name <- if ("l2_cope_name" %in% names(dt) && !is.null(dt$l2_cope_name[1]) && !is.na(dt$l2_cope_name[1])) dt$l2_cope_name[1] else "none"

      target_dir <- as.character(glue::glue(gpa$output_locations$afni_3dlmer_directory))
      
      prefix <- file.path(target_dir, paste0(l3_name, "_LMEr"))
      mask_file <- resolve_3dlmer_mask(
        target_dir = target_dir,
        harvested_inputs = dt,
        explicit_mask = l3_obj$lmer_mask %||% gpa$mask_file,
        requires_l2 = requires_l2,
        lg = lg
      )
      
      cmd <- build_3dlmer_command(
        prefix = prefix,
        model_formula = l3_obj$lmer_formula,
        qVars = qVars,
        qVarCenters = qVarCenters,
        glt_codes = glt_codes,
        data_table_file = "dataTable.txt", # Relative to script
        mask = mask_file,
        njobs = l3_obj$lmer_njobs %||% gpa$parallel$afni$l3_lmer_njobs %||% 1,
        ss_type = 3
      )
      
      files <- write_3dlmer_files(target_dir, datatable, cmd)
      
      afni_setup_df[[length(afni_setup_df) + 1]] <- data.frame(
        l1_model = l1_model,
        l2_model = l2_model,
        l3_model = l3_model,
        l1_cope_name = l1_cope_name,
        l2_cope_name = l2_cope_name,
        afni_script = files$script,
        mask_file = mask_file %||% NA_character_,
        output_file = prefix,
        afni_complete = FALSE,
        afni_failed = NA,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(afni_setup_df) == 0L) {
    return(list(metadata = NULL, data = NULL, id_cols = NULL))
  }

  afni_df <- do.call(rbind, afni_setup_df)

  list(
    data = afni_df,
    id_cols = c("l1_model", "l2_model", "l3_model", "l1_cope_name", "l2_cope_name")
  )
}


############
# helper function to get a cope data.frame for all level 1 models
get_l1_cope_df <- function(gpa, model_set, subj_df=NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(.data$exclude_subject == FALSE) %>%
      dplyr::select("id", "session")
  }
  if (is.null(gpa$l1_cope_names) || !checkmate::test_list(gpa$l1_cope_names)) {
    gpa <- refresh_l1_cope_names(gpa)
  }
  checkmate::assert_data_frame(model_set)
  dt <- dplyr::bind_rows(
    lapply(unique(model_set$l1_model), function(mm) {
      cope_names <- gpa$l1_cope_names[[mm]]
      if (is.null(cope_names)) {
        this_model <- gpa$l1_models$models[[mm]]
        cope_names <- extract_contrast_names(
          this_model$contrasts,
          model_name = mm,
          level_label = "L1",
          allow_empty = TRUE
        )
      }
      if (length(cope_names) == 0L) {
        stop(
          "No L1 contrasts available for model '", mm,
          "'. Ensure contrasts are defined before setup_l3_models().",
          call. = FALSE
        )
      }
      data.frame(
        l1_model = mm,
        l1_cope_number = seq_along(cope_names),
        l1_cope_name = cope_names
      )
    })
  )

  #return subject-specific rows for each cope
  dt <- dt %>% tidyr::crossing(subj_df)
  return(dt)
}

# helper function to get a cope data.frame for all level 3 models
# handles the per-subject cope numbering problem
get_l2_cope_df <- function(gpa, model_set, subj_df=NULL) {
  checkmate::assert_data_frame(model_set)

  if (is.null(gpa$l2_model_setup) ||
      !inherits(gpa$l2_model_setup, "l2_setup") ||
      is.null(gpa$l2_model_setup$fsl) ||
      !is.data.frame(gpa$l2_model_setup$fsl)) {
    stop("FSL L3 setup requires per-cope gpa$l2_model_setup$fsl rows from setup_l2_models().", call. = FALSE)
  }

  required_cols <- c(
    "id", "session", "l1_model", "l1_cope_number", "l1_cope_name",
    "l2_model", "l2_input_mode", "cope_list"
  )
  missing_cols <- setdiff(required_cols, names(gpa$l2_model_setup$fsl))
  if (length(missing_cols) > 0L) {
    stop(
      "FSL L3 setup requires per-cope L2 setup columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  l2_setup <- gpa$l2_model_setup$fsl %>%
    dplyr::filter(.data$l2_model %in% unique(model_set$l2_model))
  bad_modes <- setdiff(unique(l2_setup$l2_input_mode), c("cope_files", "l1_cope_file_passthrough"))
  if (length(bad_modes) > 0L) {
    stop("Unsupported L2 input mode in per-cope L2 setup: ", paste(bad_modes, collapse = ", "), call. = FALSE)
  }

  l2_rows <- lapply(seq_len(nrow(l2_setup)), function(ii) {
    cope_df <- l2_setup$cope_list[[ii]]
    if (is.null(cope_df) || nrow(cope_df) == 0L) return(NULL)

    if (!all(c("l2_cope_number", "l2_cope_name") %in% names(cope_df))) {
      stop("Each per-cope L2 setup row must contain cope_list entries with l2_cope_number and l2_cope_name.", call. = FALSE)
    }

    cope_df$l2_model <- l2_setup$l2_model[ii]
    cope_df$l1_model <- l2_setup$l1_model[ii]
    cope_df$l1_cope_name <- l2_setup$l1_cope_name[ii]
    cope_df$l1_cope_number <- l2_setup$l1_cope_number[ii]
    cope_df$l2_input_mode <- l2_setup$l2_input_mode[ii]
    cope_df
  })

  out <- dplyr::bind_rows(l2_rows)
  if (nrow(out) == 0L) {
    return(data.frame(
      id = character(0), session = integer(0),
      l2_cope_number = integer(0), l2_cope_name = character(0),
      l2_model = character(0), l1_model = character(0),
      l1_cope_name = character(0), l1_cope_number = integer(0),
      l2_input_mode = character(0),
      stringsAsFactors = FALSE
    ))
  }

  out
}

# helper function to get a cope data.frame for all level 3 models
get_l3_cope_df <- function(gpa, model_set, subj_df=NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(.data$exclude_subject == FALSE) %>%
      dplyr::select("id", "session")
  }
  checkmate::assert_data_frame(model_set)
  dt <- dplyr::bind_rows(
    lapply(unique(model_set$l3_model), function(mm) {
      m3 <- gpa$l3_models$models[[mm]]
      cope_names <- extract_contrast_names(
        m3$contrasts,
        model_name = mm,
        level_label = "L3",
        allow_empty = TRUE
      )
      if (length(cope_names) == 0L && identical(m3$l3_input_mode, "3dlmer")) {
        glt_nms <- names(m3$lmer_glt_codes)
        cope_names <- if (!is.null(glt_nms) && length(glt_nms) > 0L) glt_nms else "omnibus"
      }
      if (length(cope_names) == 0L) {
        stop(
          "No L3 contrasts available for model '", mm,
          "'. Ensure contrasts are defined before setup_l3_models().",
          call. = FALSE
        )
      }
      data.frame(
        l3_model = mm,
        l3_cope_number = seq_along(cope_names),
        l3_cope_name = cope_names
      )
    })
  )

  # return subject-specific rows for each cope
  dt <- dt %>% tidyr::crossing(subj_df)
  return(dt)

}


# data.table cross-join (tidyr::crossing is a bit faster and already exists)
# https://stackoverflow.com/questions/10600060/how-to-do-cross-join-in-r
# CJ.table <- function(X, Y) {
#   setkey(X[, c(k = 1, .SD)], k)[Y[, c(k = 1, .SD)], allow.cartesian = TRUE][, k := NULL]
# }

get_fsl_l3_model_df <- function(gpa, model_df, subj_df=NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(.data$exclude_subject == FALSE) %>%
      dplyr::select("id", "session")
  }
  if (!"l3_input_mode" %in% names(model_df)) {
    model_df$l3_input_mode <- vapply(
      model_df$l3_model,
      function(mm) {
        normalize_l3_input_mode(gpa$l3_models$models[[mm]]$l3_input_mode)
      },
      character(1)
    )
  }
  model_df$model_id <- seq_len(nrow(model_df))

  l1_df <- get_l1_cope_df(gpa, model_df, subj_df)

  l3_df <- get_l3_cope_df(gpa, model_df, subj_df)

  if (isTRUE(gpa$multi_run)) {
    # model_df has l1_model, l2_model, l3_model
    l2_df <- get_l2_cope_df(gpa, model_df)
    l2_join_keys <- c("l1_model", "l1_cope_number", "l1_cope_name", "l2_model")

    # Identify L2 models with l2_scope="id" (session=0L sentinel) vs real sessions
    l2_id_scope_models <- character(0)
    for (mm in unique(model_df$l2_model)) {
      if (identical(gpa$l2_models$models[[mm]]$l2_scope, "id")) {
        l2_id_scope_models <- c(l2_id_scope_models, mm)
      }
    }

    base <- model_df %>%
      tidyr::crossing(subj_df) %>%
      left_join(l1_df, by = c("id", "session", "l1_model"))

    if (length(l2_id_scope_models) > 0L) {
      # For l2_scope="id", l2_df uses session=0L; join without session key
      l2_df_id <- l2_df %>%
        dplyr::filter(.data$l2_model %in% l2_id_scope_models) %>%
        dplyr::select(-"session")
      l2_df_session <- l2_df %>%
        dplyr::filter(!.data$l2_model %in% l2_id_scope_models)

      base_id <- base %>%
        dplyr::filter(.data$l2_model %in% l2_id_scope_models) %>%
        left_join(
          l2_df_id,
          by = c("id", l2_join_keys)
        )
      base_session <- base %>%
        dplyr::filter(!.data$l2_model %in% l2_id_scope_models)

      if (nrow(base_session) > 0L) {
        base_session <- base_session %>%
          left_join(
            l2_df_session,
            by = c("id", "session", l2_join_keys)
          )
      }
      combined <- dplyr::bind_rows(base_id, base_session) %>%
        left_join(l3_df, by = c("id", "session", "l3_model"))
    } else {
      combined <- base %>%
        left_join(
          l2_df,
          by = c("id", "session", l2_join_keys)
        ) %>%
        left_join(l3_df, by = c("id", "session", "l3_model"))
    }
  } else {
    combined <- model_df %>%
      tidyr::crossing(subj_df) %>%
      left_join(l1_df, by = c("id", "session", "l1_model")) %>%
      left_join(l3_df, by = c("id", "session", "l3_model"))
  }

  return(combined)
}

get_feat_l3_inputs <- function(gpa, l3_cope_config, lg=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(l3_cope_config)

  if (is.null(lg)) { lg <- lgr::get_logger() }
  if (isTRUE(gpa$multi_run)) {
    # feat directories in $l2_model_setup
    feat_inputs <- gpa$l2_model_setup$fsl %>%
      dplyr::filter(.data$feat_complete == TRUE)
    required_cols <- c(
      "id", "session", "l1_model", "l1_cope_number", "l1_cope_name",
      "l2_model", "l2_scope", "l2_input_mode", "feat_dir",
      "passthrough_cope_file"
    )
    missing_cols <- setdiff(required_cols, names(feat_inputs))
    if (length(missing_cols) > 0L) {
      stop(
        "FSL L3 inputs require per-cope L2 setup columns: ",
        paste(required_cols, collapse = ", "),
        ". Missing: ", paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
    bad_modes <- setdiff(unique(feat_inputs$l2_input_mode), c("cope_files", "l1_cope_file_passthrough"))
    if (length(bad_modes) > 0L) {
      stop("Unsupported L2 input mode in per-cope L2 setup: ", paste(bad_modes, collapse = ", "), call. = FALSE)
    }
    l3_cope_config <- l3_cope_config %>% dplyr::select(-tidyselect::any_of("l2_input_mode"))

    # For l2_scope="id", l2_model_setup uses session=0L as sentinel because L2
    # stacks all sessions per subject. Join without session for those models,
    # then deduplicate to one row per subject (all sessions share the same FEAT).
    id_scope_models <- character(0)
    if ("l2_scope" %in% names(feat_inputs)) {
      id_scope_models <- unique(feat_inputs$l2_model[feat_inputs$l2_scope == "id"])
    }

    if (length(id_scope_models) > 0L) {
      feat_id <- feat_inputs %>%
        dplyr::filter(.data$l2_model %in% id_scope_models) %>%
        dplyr::select(-"session") %>%
        dplyr::inner_join(
          l3_cope_config %>% dplyr::filter(.data$l2_model %in% id_scope_models),
          by = c("id", "l1_model", "l1_cope_number", "l1_cope_name", "l2_model")
        ) %>%
        dplyr::distinct(.data$id, .data$l1_model, .data$l2_model, .data$l3_model, .data$l3_input_mode,
          .data$l1_cope_name, .data$l2_cope_name, .data$feat_dir, .keep_all = TRUE)

      feat_non_id <- feat_inputs %>%
        dplyr::filter(!.data$l2_model %in% id_scope_models)

      if (nrow(feat_non_id) > 0L) {
        feat_non_id <- feat_non_id %>%
          dplyr::inner_join(
            l3_cope_config %>% dplyr::filter(!.data$l2_model %in% id_scope_models),
            by = c("id", "session", "l1_model", "l1_cope_number", "l1_cope_name", "l2_model")
          )
      }

      feat_inputs <- dplyr::bind_rows(feat_id, feat_non_id)
    } else {
      #join up combination of all models with cope directories
      feat_inputs <- feat_inputs %>%
        dplyr::inner_join(
          l3_cope_config,
          by = c("id", "session", "l1_model", "l1_cope_number", "l1_cope_name", "l2_model")
        )
    }

    # Ensure any 1-run passthrough inputs have synthesized tdof_t<k>.nii.gz before L3 FLAME runs
    pt_rows <- which(feat_inputs$l2_input_mode == "l1_cope_file_passthrough")
    if (length(pt_rows) > 0L && "l1_cope_number" %in% names(feat_inputs)) {
      for (pri in pt_rows) {
        synthesize_passthrough_tdof(
          feat_dir = feat_inputs$feat_dir[pri],
          l1_cope_number = feat_inputs$l1_cope_number[pri],
          lg = lg
        )
      }
    }

    # sort out expected cope files for each model combination
    feat_inputs <- feat_inputs %>%
      dplyr::mutate(
        l2_cope_dir = ifelse(.data$l2_input_mode == "cope_files", "cope1.feat", NA_character_)
      ) %>%
      dplyr::mutate(
        cope_file = dplyr::if_else(
          .data$l2_input_mode == "l1_cope_file_passthrough",
          .data$passthrough_cope_file,
          file.path(
            .data$feat_dir,
            .data$l2_cope_dir,
            "stats",
            paste0("cope", .data$l2_cope_number, ".nii.gz")
          )
        )
      ) %>%
      dplyr::select(
        "id", "session", "l1_model", "l2_model", "l3_model", "l3_input_mode", "l1_cope_name", "l2_cope_name", "feat_dir", "cope_file"
      )

    split_on <- c("l1_cope_name", "l2_cope_name", "l1_model", "l2_model", "l3_model")
  } else {
    # feat directories in $l1_model_setup (single-run study: filter out empty/non-estimable L1 COPEs before L3!)
    l1_complete <- gpa$l1_model_setup$fsl %>%
      dplyr::filter(.data$feat_complete == TRUE)

    if (nrow(l1_complete) > 0L && all(c("feat_fsf", "feat_dir") %in% names(l1_complete))) {
      l1_for_val <- l1_complete
      if (!"run_number" %in% names(l1_for_val)) l1_for_val$run_number <- 1L
      l1_val <- build_l1_cope_validity(l1_for_val, gpa, lg = lg)
      if (nrow(l1_val) > 0L) {
        dropped_l1 <- l1_val %>% dplyr::filter(.data$valid_for_higher_level == FALSE)
        if (nrow(dropped_l1) > 0L) {
          for (di in seq_len(nrow(dropped_l1))) {
            lg$warn(
              "Dropping single-run L1 input for subject %s session %s, L1 model %s, L1 cope %s before L3 setup because cope is %s (%s).",
              dropped_l1$id[di], dropped_l1$session[di], dropped_l1$l1_model[di],
              dropped_l1$l1_cope_name[di], dropped_l1$validity_state[di],
              ifelse(is.na(dropped_l1$validity_reason[di]), "unknown", dropped_l1$validity_reason[di])
            )
          }
        }
        valid_l1_keys <- l1_val %>%
          dplyr::filter(.data$valid_for_higher_level == TRUE) %>%
          dplyr::select("id", "session", "l1_model", "l1_cope_number", "l1_cope_name") %>%
          dplyr::distinct()
        feat_inputs <- l1_complete %>%
          dplyr::inner_join(l3_cope_config, by = c("id", "session", "l1_model")) %>%
          dplyr::inner_join(valid_l1_keys, by = c("id", "session", "l1_model", "l1_cope_number", "l1_cope_name"))
      } else {
        feat_inputs <- l1_complete %>%
          dplyr::inner_join(l3_cope_config, by = c("id", "session", "l1_model"))
      }
    } else {
      feat_inputs <- l1_complete %>%
        dplyr::inner_join(l3_cope_config, by = c("id", "session", "l1_model"))
    }

    feat_inputs <- feat_inputs %>%
      dplyr::mutate(cope_file = file.path(
        .data$feat_dir,
        "stats",
        paste0("cope", .data$l1_cope_number, ".nii.gz")
      )) %>%
      dplyr::select(
        "id", "session", "l1_model", "l3_model", "l3_input_mode", "l1_cope_name", "feat_dir", "cope_file"
      )

    split_on <- c("l1_cope_name", "l1_model", "l3_model")

  }

  if (!"l3_input_mode" %in% names(feat_inputs)) {
    feat_inputs$l3_input_mode <- "per_session"
  }
  feat_inputs$l3_input_mode <- normalize_l3_input_mode(feat_inputs$l3_input_mode)
  feat_inputs$l3_session_partition <- ifelse(
    feat_inputs$l3_input_mode == "per_session",
    as.character(feat_inputs$session),
    "all"
  )
  split_on <- c(split_on, "l3_session_partition")

  if (!is.data.table(feat_inputs)) data.table::setDT(feat_inputs)
  feat_inputs <- split(feat_inputs, by = split_on)

  return(feat_inputs)
}

get_spm_contrast_file <- function(spm_dir, cope_number) {
  base <- sprintf("con_%04d.nii", cope_number)
  cand <- file.path(spm_dir, base)
  if (file.exists(cand)) return(cand)
  cand_gz <- paste0(cand, ".gz")
  if (file.exists(cand_gz)) return(cand_gz)
  cand_gz2 <- sub("\\.nii$", ".nii.gz", cand)
  if (file.exists(cand_gz2)) return(cand_gz2)
  return(cand)
}

get_spm_l3_model_df <- function(gpa, model_df, subj_df = NULL) {
  if (is.null(subj_df)) {
    subj_df <- gpa$subject_data %>%
      dplyr::filter(.data$exclude_subject == FALSE) %>%
      dplyr::select("id", "session")
  }
  if (!"l3_input_mode" %in% names(model_df)) {
    model_df$l3_input_mode <- vapply(
      model_df$l3_model,
      function(mm) {
        normalize_l3_input_mode(gpa$l3_models$models[[mm]]$l3_input_mode)
      },
      character(1)
    )
  }
  model_df$model_id <- seq_len(nrow(model_df))

  l1_df <- get_l1_cope_df(gpa, model_df, subj_df)
  l3_df <- get_l3_cope_df(gpa, model_df, subj_df)

  combined <- model_df %>%
    tidyr::crossing(subj_df) %>%
    left_join(l1_df, by = c("id", "session", "l1_model")) %>%
    left_join(l3_df, by = c("id", "session", "l3_model"))

  return(combined)
}

get_spm_l3_inputs <- function(gpa, l3_cope_config, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(l3_cope_config)

  if (is.null(lg)) { lg <- lgr::get_logger() }

  spm_inputs <- gpa$l1_model_setup$spm
  if (is.null(spm_inputs) || (is.data.frame(spm_inputs) && nrow(spm_inputs) == 0L)) {
    lg$warn("No SPM L1 model setup data found. Cannot setup SPM L3 models.")
    return(list())
  }

  if (!"spm_dir" %in% names(spm_inputs)) {
    lg$warn("No spm_dir column in gpa$l1_model_setup$spm. Cannot setup SPM L3 models.")
    return(list())
  }

  if ("spm_complete" %in% names(spm_inputs)) {
    spm_inputs <- spm_inputs %>%
      dplyr::filter(.data$spm_complete == TRUE)
  } else {
    lg$warn("No spm_complete column in gpa$l1_model_setup$spm. Proceeding without completion filter for L3 inputs.")
  }

  if ("spm_contrast_exists" %in% names(spm_inputs)) {
    spm_inputs <- spm_inputs %>%
      dplyr::filter(.data$spm_contrast_exists == TRUE)
  }

  spm_inputs <- spm_inputs %>%
    dplyr::inner_join(l3_cope_config, by = c("id", "session", "l1_model"))

  if (nrow(spm_inputs) == 0L) {
    return(list())
  }

  spm_inputs <- spm_inputs %>%
    dplyr::mutate(con_file = mapply(get_spm_contrast_file, .data$spm_dir, .data$l1_cope_number, USE.NAMES = FALSE)) %>%
    dplyr::select("id", "session", "l1_model", "l3_model", "l3_input_mode", "l1_cope_name", "spm_dir", "con_file")

  if (!is.character(spm_inputs$con_file)) {
    spm_inputs$con_file <- as.character(spm_inputs$con_file)
  }

  spm_inputs <- spm_inputs %>%
    dplyr::filter(!is.na(.data$con_file) & nzchar(.data$con_file))

  if (nrow(spm_inputs) == 0L) {
    return(list())
  }

  spm_inputs <- spm_inputs %>%
    dplyr::filter(file.exists(.data$con_file))

  if (!"l3_input_mode" %in% names(spm_inputs)) {
    spm_inputs$l3_input_mode <- "per_session"
  }
  spm_inputs$l3_input_mode <- normalize_l3_input_mode(spm_inputs$l3_input_mode)
  spm_inputs$l3_session_partition <- ifelse(
    spm_inputs$l3_input_mode == "per_session",
    as.character(spm_inputs$session),
    "all"
  )

  if (!is.data.table(spm_inputs)) data.table::setDT(spm_inputs)
  spm_inputs <- split(spm_inputs, by = c("l1_cope_name", "l1_model", "l3_model", "l3_session_partition"))

  return(spm_inputs)
}
