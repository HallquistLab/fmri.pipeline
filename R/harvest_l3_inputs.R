#' Harvest standardized L3 input tables from upstream producer backends
#'
#' @param gpa a glm_pipeline_arguments object
#' @param l3_backend target L3 backend specification or backend name
#' @param l1_model_names optional character vector of L1 models to include
#' @param l2_model_names optional character vector of L2 models to include
#' @param l3_model_names optional character vector of L3 models to process
#' @param lg optional logger
#'
#' @return a list of data.frames, one per L1/L2 contrast combination, for each L3 model
#' @keywords internal
harvest_l3_inputs <- function(gpa, l3_backend, l1_model_names = NULL, l2_model_names = NULL,
                              l3_model_names = NULL, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (is.null(lg)) lg <- lgr::get_logger("glm_pipeline/harvest_l3_inputs")

  if (is.character(l3_backend)) {
    checkmate::assert_string(l3_backend)
    l3_backend <- list(name = tolower(l3_backend))
  } else {
    checkmate::assert_list(l3_backend)
  }

  producer_backends <- backend_l3_input_provider_backends(l3_backend)
  if (length(producer_backends) == 0L) {
    producer_backends <- resolve_l3_input_producers(get_glm_backends(gpa, must_exist = FALSE), l3_backend)
  }

  if (length(producer_backends) == 0L) {
    lg$warn(
      "Could not resolve an upstream L3 input producer for backend '%s'.",
      backend_name(l3_backend)
    )
    return(NULL)
  }

  all_harvested <- list()
  for (producer_backend in producer_backends) {
    harvested <- switch(
      producer_backend,
      fsl = harvest_l3_inputs_fsl(
        gpa = gpa,
        l1_model_names = l1_model_names,
        l2_model_names = l2_model_names,
        l3_model_names = l3_model_names,
        lg = lg
      ),
      stop("harvest_l3_inputs does not support producer backend: '", producer_backend, "'")
    )

    if (is.null(harvested)) {
      lg$warn(
        "No harvested L3 inputs available from producer backend '%s' for target backend '%s'.",
        producer_backend, backend_name(l3_backend)
      )
      next
    }

    for (l3_name in names(harvested)) {
      if (is.null(all_harvested[[l3_name]])) {
        all_harvested[[l3_name]] <- harvested[[l3_name]]
      } else {
        all_harvested[[l3_name]] <- c(all_harvested[[l3_name]], harvested[[l3_name]])
      }
    }
  }

  if (length(all_harvested) == 0L) return(NULL)
  all_harvested
}

#' Harvest FSL-produced inputs for downstream L3 consumers
#'
#' @param gpa a glm_pipeline_arguments object
#' @param l1_model_names optional character vector of L1 models to include
#' @param l2_model_names optional character vector of L2 models to include
#' @param l3_model_names optional character vector of L3 models to process
#' @param lg optional logger
#'
#' @return a list of data.frames, one per L1/L2 contrast combination, for each L3 model
#' @keywords internal
harvest_l3_inputs_fsl <- function(gpa, l1_model_names = NULL, l2_model_names = NULL,
                                  l3_model_names = NULL, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_character(l1_model_names, null.ok = TRUE)
  checkmate::assert_character(l2_model_names, null.ok = TRUE)
  checkmate::assert_character(l3_model_names, null.ok = TRUE)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (is.null(lg)) lg <- lgr::get_logger("glm_pipeline/harvest_l3_inputs_fsl")

  if (is.null(l3_model_names)) {
    l3_model_names <- names(gpa$l3_models$models)
  }

  subj_df <- gpa$subject_data %>%
    dplyr::filter(.data$exclude_subject == FALSE) %>%
    dplyr::select("id", "session")

  # Support single-run multi-session longitudinal designs (multi_run == FALSE or no L2 models)
  has_l2_models <- !is.null(gpa$l2_models) && length(gpa$l2_models$models) > 0L
  if (isFALSE(gpa$multi_run) || !isTRUE(has_l2_models)) {
    this_l1_models <- l1_model_names %||% names(gpa$l1_models$models)
    if (length(this_l1_models) == 0L || length(l3_model_names) == 0L) {
      lg$warn("No L1 or L3 models available for direct L1 -> L3 harvesting.")
      return(NULL)
    }

    all_harvested <- list()
    for (this_l3 in l3_model_names) {
      model_df <- expand.grid(
        l1_model = this_l1_models,
        l3_model = this_l3,
        stringsAsFactors = FALSE
      )
      l3_cope_config <- get_fsl_l3_model_df(gpa, model_df, subj_df)
      if ("l3_cope_number" %in% names(l3_cope_config)) {
        l3_cope_config <- l3_cope_config %>%
          dplyr::filter(.data$l3_cope_number == 1L) %>%
          dplyr::select(-dplyr::any_of(c("l3_cope_number", "l3_cope_name")))
      }
      feat_inputs <- get_feat_l3_inputs(gpa, l3_cope_config, lg = lg)
      for (nm in names(feat_inputs)) {
        df_in <- feat_inputs[[nm]]
        if ("feat_dir" %in% names(df_in) && "cope_file" %in% names(df_in)) {
          reg_std <- file.path(df_in$feat_dir, "reg_standard", "stats", basename(df_in$cope_file))
          use_std <- file.exists(reg_std)
          df_in$cope_file[use_std] <- reg_std[use_std]
        }
        feat_inputs[[nm]] <- df_in %>%
          dplyr::rename(InputFile = "cope_file") %>%
          dplyr::mutate(source_backend = "fsl")
      }
      all_harvested[[this_l3]] <- c(all_harvested[[this_l3]], feat_inputs)
    }

    all_harvested <- Filter(function(x) length(x) > 0L, all_harvested)
    if (length(all_harvested) == 0L) return(NULL)
    return(all_harvested)
  }

  pairs <- resolve_l2_l3_compatible_pairs(
    gpa,
    l2_model_names = l2_model_names,
    l3_model_names = l3_model_names
  )

  if (nrow(pairs) == 0L) {
    lg$warn("No compatible L2/L3 model pairs found for harvesting.")
    return(NULL)
  }

  all_harvested <- list()
  for (i in seq_len(nrow(pairs))) {
    this_l2 <- pairs$l2_model[i]
    this_l3 <- pairs$l3_model[i]

    this_l1_models <- gpa$l2_models$models[[this_l2]]$l1_model_names
    if (is.null(this_l1_models)) {
      this_l1_models <- names(gpa$l1_models$models)
    }
    if (!is.null(l1_model_names)) {
      this_l1_models <- intersect(this_l1_models, l1_model_names)
    }
    if (length(this_l1_models) == 0L) next

    model_df <- expand.grid(
      l1_model = this_l1_models,
      l2_model = this_l2,
      l3_model = this_l3,
      stringsAsFactors = FALSE
    )

    subj_df <- gpa$subject_data %>%
      dplyr::filter(.data$exclude_subject == FALSE) %>%
      dplyr::select("id", "session")

    l3_cope_config <- get_fsl_l3_model_df(gpa, model_df, subj_df)
    if ("l3_cope_number" %in% names(l3_cope_config)) {
      l3_cope_config <- l3_cope_config %>%
        dplyr::filter(.data$l3_cope_number == 1L) %>%
        dplyr::select(-dplyr::any_of(c("l3_cope_number", "l3_cope_name")))
    }
    feat_inputs <- get_feat_l3_inputs(gpa, l3_cope_config, lg = lg)

    for (nm in names(feat_inputs)) {
      feat_inputs[[nm]] <- feat_inputs[[nm]] %>%
        dplyr::rename(InputFile = "cope_file") %>%
        dplyr::mutate(source_backend = "fsl")
    }

    all_harvested[[this_l3]] <- c(all_harvested[[this_l3]], feat_inputs)
  }

  if (length(all_harvested) == 0L) return(NULL)
  all_harvested
}

#' Compatibility wrapper for legacy internal callers
#'
#' @param gpa a glm_pipeline_arguments object
#' @param l1_model_names optional character vector of L1 models to include
#' @param l2_model_names optional character vector of L2 models to include
#' @param l3_model_names optional character vector of L3 models to process
#'
#' @return a list of harvested FSL inputs for downstream L3 consumers
#' @keywords internal
harvest_fsl_copes <- function(gpa, l1_model_names = NULL, l2_model_names = NULL, l3_model_names = NULL) {
  harvest_l3_inputs_fsl(
    gpa = gpa,
    l1_model_names = l1_model_names,
    l2_model_names = l2_model_names,
    l3_model_names = l3_model_names
  )
}
