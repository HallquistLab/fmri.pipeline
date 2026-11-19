# Harvest named contrast maps (never hard-coded AFNI sub-brick indices).
harvest_l3_inputs_afni <- function(gpa, l1_model_names = NULL, l3_model_names = NULL, lg = NULL) {
  setup <- gpa$l1_model_setup$afni
  if (is.null(setup) || nrow(setup) == 0L) return(NULL)
  if (is.null(l1_model_names)) l1_model_names <- names(gpa$l1_models$models)
  if (is.null(l3_model_names)) l3_model_names <- names(gpa$l3_models$models)
  subjects <- gpa$subject_data
  if ("exclude_subject" %in% names(subjects)) subjects <- subjects[!subjects$exclude_subject, , drop = FALSE]
  setup <- dplyr::semi_join(setup, subjects, by = c("id", "session"))
  setup <- setup[setup$l1_model %in% l1_model_names, , drop = FALSE]
  maps <- lapply(seq_len(nrow(setup)), function(i) {
    r <- setup[i, , drop = FALSE]
    if (!get_afni_glm_status(r$afni_dir)$afni_complete) return(NULL)
    mf <- utils::read.delim(file.path(r$afni_dir, "contrasts.tsv"), stringsAsFactors = FALSE)
    wanted <- rownames(gpa$l1_models$models[[r$l1_model]]$contrasts)
    mf <- mf[mf$contrast_name %in% wanted & file.exists(mf$coef_file), , drop = FALSE]
    if (nrow(mf) == 0L) return(NULL)
    data.frame(id = r$id, session = r$session, l1_model = r$l1_model,
      l1_cope_name = mf$contrast_name, l1_cope_number = mf$contrast_number,
      InputFile = mf$coef_file, mask_file = file.path(r$afni_dir, "mask.nii.gz"), source_backend = "afni")
  })
  maps <- dplyr::bind_rows(maps)
  if (!nrow(maps)) return(NULL)
  out <- list()
  for (model in l3_model_names) {
    dt <- maps
    dt$l3_model <- model
    dt$l3_input_mode <- normalize_l3_input_mode(gpa$l3_models$models[[model]]$l3_input_mode)
    dt$l3_session_partition <- ifelse(dt$l3_input_mode == "per_session", as.character(dt$session), "all")
    out[[model]] <- split(data.table::as.data.table(dt),
      by = c("l1_model", "l1_cope_name", "l3_model", "l3_session_partition"))
  }
  out
}

# Dispatch mixed models to the established 3dLMEr implementation, and ordinary
# group GLMs to 3dDeconvolve OLS with the exact pipeline design and contrasts.
setup_l3_backend_afni <- function(gpa, backend, lg, l1_model_names, l2_model_names,
                                  l3_model_names, l2_l3_pairs, subj_df, requires_l2) {
  mixed <- vapply(gpa$l3_models$models[l3_model_names], function(x)
    identical(normalize_l3_input_mode(x$l3_input_mode), "3dlmer"), logical(1))
  results <- list()
  if (any(mixed)) {
    results[[1L]] <- afni_3dlmer_setup(gpa, backend, lg, l1_model_names, l2_model_names,
      l3_model_names[mixed], l2_l3_pairs, subj_df, requires_l2)$data
  }
  providers <- get_effective_model_backends(gpa, 3L, l3_model_names, type = "producer")
  for (model in l3_model_names[!mixed]) {
    this_backend <- backend
    this_backend$l3_input_provider_backends <- providers[[model]] %||% "afni"
    if (!length(this_backend$l3_input_provider_backends)) this_backend$l3_input_provider_backends <- "afni"
    inputs <- harvest_l3_inputs(gpa, this_backend, l1_model_names, l2_model_names, model, lg)
    for (df in inputs[[model]]) {
      results[[length(results) + 1L]] <- afni_l3_model(as.data.frame(df), gpa)
    }
  }
  data <- dplyr::bind_rows(results)
  list(data = data, metadata = NULL, id_cols = "afni_script")
}

afni_l3_model <- function(input, gpa) {
  model <- unique(input$l3_model)
  checkmate::assert_string(model)
  mode <- normalize_l3_input_mode(gpa$l3_models$models[[model]]$l3_input_mode)
  # Repeated observations require the explicit subject-effects design or 3dLMEr.
  if (anyDuplicated(input[, c("id", "session")])) stop("Duplicate AFNI L3 subject/session inputs.")
  if (mode == "per_session" && length(unique(input$session)) != 1L) stop("AFNI per_session inputs must share a session.")
  if (mode == "subject_rows" && anyDuplicated(input$id)) {
    stop("AFNI subject_rows requires one input per subject; use per_session or 3dlmer for multiple sessions.")
  }
  fit <- respecify_l3_model(gpa$l3_models$models[[model]], new_data = input)
  input <- dplyr::left_join(fit$metadata, input, by = c("id", "session"))
  checkmate::assert_file_exists(input$InputFile)
  design <- as.matrix(fit$model_matrix)
  if (nrow(design) != nrow(input) || qr(design)$rank < ncol(design)) stop("AFNI L3 design is misaligned or rank deficient.")
  settings <- afni_glm_settings(gpa)
  session <- if (mode == "per_session") unique(input$session) else 0L
  l2_model <- if ("l2_model" %in% names(input)) unique(input$l2_model) else NA_character_
  l2_cope_name <- if ("l2_cope_name" %in% names(input)) unique(input$l2_cope_name) else NA_character_
  outdir <- get_output_directory(l1_model = input$l1_model[1L], l3_model = model,
    l1_contrast = input$l1_cope_name[1L], gpa = gpa, what = "l3", glm_software = "afni")
  if (!is.na(l2_model)) outdir <- file.path(outdir, paste0("L2m-", l2_model), paste0("l2c-", l2_cope_name))
  if (mode == "per_session") outdir <- file.path(outdir, paste0("ses-", session))
  masks <- if ("mask_file" %in% names(input)) input$mask_file else if ("feat_dir" %in% names(input)) {
    vapply(seq_len(nrow(input)), function(i) {
      standard <- file.path(input$feat_dir[i], "reg_standard")
      if (startsWith(input$InputFile[i], paste0(standard, "/"))) file.path(standard, "mask.nii.gz") else file.path(input$feat_dir[i], "mask.nii.gz")
    }, character(1))
  } else NULL
  if (is.null(settings$mask) && (is.null(masks) || anyNA(masks) || !all(file.exists(masks)))) {
    stop("AFNI L3 requires source masks or glm_settings$afni$mask; contrast images cannot be automasked.")
  }
  generated <- afni_write_glm(outdir, input$InputFile, design, fit$contrasts,
    mask = settings$mask, mask_inputs = if (is.null(settings$mask)) unique(masks) else NULL,
    write_residuals = settings$write_residuals, force = settings$force_l3_creation,
    cpus = gpa$parallel$afni$l3_glm_cpus %||% 1L)
  data.frame(l1_model = input$l1_model[1L], l2_model = l2_model, l3_model = model,
    l1_cope_name = input$l1_cope_name[1L], l2_cope_name = l2_cope_name, session = session,
    afni_dir = dirname(generated$script), afni_script = generated$script,
    output_file = generated$output_file, get_afni_glm_status(dirname(generated$script)))
}

get_afni_l3_status <- function(output_file, lg = NULL) {
  if (basename(output_file) == "stats.nii.gz" || file.exists(file.path(dirname(output_file), "expected_outputs.txt"))) {
    get_afni_glm_status(dirname(output_file), lg)
  } else {
    get_3dlmer_status(output_file, lg)
  }
}
