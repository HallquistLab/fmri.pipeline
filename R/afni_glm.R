# Shared AFNI command generation and status for first-level REML and group OLS.
# All user paths are shell quoted; AFNI labels are generated independently from
# model labels, with the original names retained in contrast/regressor manifests.
afni_write_glm <- function(outdir, inputs, design, contrasts, n_task = ncol(design),
                           censor = rep(1L, nrow(design)), starts = 0L, polort = -1L,
                           tr = 1, estimation = "OLS", mask = NULL, write_residuals = TRUE,
                           cpus = 1L, force = FALSE, input_files = inputs, mask_inputs = NULL) {
  design <- as.matrix(design)
  contrasts <- as.matrix(contrasts)
  checkmate::assert_matrix(design, mode = "numeric", min.rows = 2L, min.cols = 1L, any.missing = FALSE)
  checkmate::assert_matrix(contrasts, mode = "numeric", min.rows = 1L, any.missing = FALSE)
  checkmate::assert_integerish(cpus, lower = 1L, len = 1L)
  if (any(!is.finite(design)) || any(!is.finite(contrasts))) stop("AFNI design and contrasts must be finite.")
  if (is.null(colnames(contrasts)) || anyDuplicated(colnames(contrasts)) ||
      !all(colnames(contrasts) %in% colnames(design)[seq_len(n_task)])) {
    stop("AFNI contrast columns must match named task regressors.")
  }
  if (any(rowSums(abs(contrasts)) == 0)) stop("AFNI contrasts must have nonzero weights.")
  if (sum(censor) <= ncol(design)) stop("Insufficient uncensored observations for the AFNI design.")
  if (!is.null(mask)) {
    checkmate::assert_file_exists(mask)
    mask <- normalizePath(mask, mustWork = TRUE)
  }
  if (!is.null(mask_inputs)) checkmate::assert_file_exists(mask_inputs)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir <- normalizePath(outdir, mustWork = TRUE)
  script <- file.path(outdir, "run_afni_glm.sh")
  output <- file.path(outdir, "stats.nii.gz")
  input_files <- unique(c(input_files, mask, mask_inputs))
  input_files <- input_files[!is.na(input_files) & nzchar(input_files)]
  # Check spatial correspondence before AFNI can merely warn about grid differences.
  image_files <- unique(c(sub("\\[[^]]*\\]$", "", inputs), mask, mask_inputs))
  reference <- roi_image_geometry(image_files[1L])
  for (file in image_files[-1L]) {
    geometry <- roi_image_geometry(file)
    geometry$volumes <- 1L # time lengths may differ; only spatial axes must agree
    tryCatch(validate_roi_geometry(reference, geometry, image_files[1L], file),
      error = function(e) stop("AFNI inputs have incompatible spatial geometry: ", file, call. = FALSE))
  }
  signature <- list(format_version = 1L, inputs = inputs, design = design, contrasts = contrasts, n_task = n_task,
    censor = censor, starts = starts, polort = polort, tr = tr, estimation = estimation,
    mask = mask, mask_inputs = mask_inputs, write_residuals = write_residuals, cpus = cpus,
    files = file.info(input_files)[, c("size", "mtime"), drop = FALSE])
  signature_file <- file.path(outdir, "afni_setup.rds")
  if (!force && file.exists(signature_file) && file.exists(script) &&
      identical(tryCatch(readRDS(signature_file), error = function(e) NULL), signature)) {
    return(list(script = script, output_file = output))
  }
  # A changed specification invalidates previous completion, even before a rerun.
  unlink(file.path(outdir, c(".afni_complete", ".afni_failed", ".afni_start")))
  q <- function(x) shQuote(x, type = "sh")
  stim_labels <- sprintf("s%04d", seq_len(ncol(design)))
  glt_labels <- sprintf("c%04d", seq_len(nrow(contrasts)))
  contrast_names <- rownames(contrasts) %||% paste0("contrast", seq_len(nrow(contrasts)))
  manifest <- data.frame(contrast_number = seq_len(nrow(contrasts)), contrast_name = contrast_names,
    afni_label = glt_labels, coef_file = file.path(outdir, sprintf("con_%04d.nii.gz", seq_len(nrow(contrasts)))),
    tstat_file = file.path(outdir, sprintf("tstat_%04d.nii.gz", seq_len(nrow(contrasts)))), stringsAsFactors = FALSE)
  utils::write.table(manifest, file.path(outdir, "contrasts.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
  utils::write.table(data.frame(regressor = colnames(design), afni_label = stim_labels,
    baseline = seq_len(ncol(design)) > n_task), file.path(outdir, "regressors.tsv"),
    sep = "\t", row.names = FALSE, quote = TRUE)
  afni_write_1d(design, file.path(outdir, "regressors.1D"))
  afni_write_1d(censor, file.path(outdir, "censor.1D"))
  afni_write_1d(matrix(starts, nrow = 1L), file.path(outdir, "runs.1D"))
  args <- c("3dDeconvolve -overwrite", paste("-force_TR", tr), "-input input.nii.gz",
    paste("-polort", polort), "-concat runs.1D -censor censor.1D",
    "-mask mask.nii.gz -tout -float", paste("-jobs", as.integer(cpus)),
    paste("-num_stimts", ncol(design)))
  for (i in seq_len(ncol(design))) {
    args <- c(args, paste("-stim_file", i, q(paste0("regressors.1D[", i - 1L, "]")),
      "-stim_label", i, stim_labels[i], if (i > n_task) paste("-stim_base", i) else ""))
  }
  args <- c(args, paste("-num_glt", nrow(contrasts)))
  for (i in seq_len(nrow(contrasts))) {
    w <- contrasts[i, ]
    keep <- which(w != 0)
    sym <- paste0(sprintf("%+.17g", w[keep]), "*", stim_labels[match(colnames(contrasts)[keep], colnames(design))])
    args <- c(args, paste("-gltsym", q(paste("SYM:", paste(sym, collapse = " "))), "-glt_label", i, glt_labels[i]))
  }
  args <- c(args, "-x1D design.xmat.1D -bucket stats.nii.gz")
  if (estimation == "REML") {
    args <- c(args, "-x1D_stop")
    fit <- paste("3dREMLfit -overwrite -input input.nii.gz -matrix design.xmat.1D",
      "-mask mask.nii.gz -tout -Rbuck stats.nii.gz",
      if (write_residuals) "-Rerrts residuals.nii.gz" else "")
  } else {
    if (write_residuals) args <- c(args, "-errts residuals.nii.gz")
    fit <- character(0)
  }
  mask_sources <- mask_inputs %||% mask
  mask_links <- character(0)
  masks <- character(0)
  mask_commands <- if (!is.null(mask_sources)) {
    # Some AFNI readers split paths on spaces internally, even when shell quoted.
    # Local aliases avoid that parser behavior for masks supplied by the user.
    mask_links <- paste0(sprintf("source_mask_%04d", seq_along(mask_sources)),
      ifelse(grepl("\\.nii\\.gz$", mask_sources), ".nii.gz", ".nii"))
    c(paste("ln -sf", q(normalizePath(mask_sources, mustWork = TRUE)), q(mask_links)),
      if (!is.null(mask_inputs)) {
        paste("3dmask_tool -overwrite -frac 1 -prefix mask.nii.gz -input", paste(q(mask_links), collapse = " "))
      } else paste("3dcopy -overwrite", q(mask_links), "mask.nii.gz"))
  } else {
    # Automask the retained run segments of the local concatenated dataset.
    # This also avoids AFNI's whitespace splitting in external input paths.
    masks <- sprintf("run_mask_%04d.nii.gz", seq_along(starts))
    ends <- c(starts[-1L], nrow(design)) - 1L
    segments <- paste0("input.nii.gz[", starts, "..", ends, "]")
    c(paste("3dAutomask -overwrite -prefix", q(masks), q(segments)),
      paste("3dmask_tool -overwrite -frac 1 -prefix mask.nii.gz -input", paste(q(masks), collapse = " ")))
  }
  expected <- c(output, file.path(outdir, c("design.xmat.1D", "mask.nii.gz")), manifest$coef_file, manifest$tstat_file,
    if (write_residuals) file.path(outdir, "residuals.nii.gz"))
  writeLines(expected, file.path(outdir, "expected_outputs.txt"))
  lines <- c("#!/usr/bin/env bash", "set -euo pipefail", paste("cd", q(outdir)),
    "exec > afni_glm.log 2>&1", paste0("export OMP_NUM_THREADS=", as.integer(cpus)),
    "rm -f .afni_complete .afni_failed", "date -u +%FT%TZ > .afni_start",
    "trap 'code=$?; if [ \"$code\" -ne 0 ]; then rm -f .afni_complete; date -u +%FT%TZ > .afni_failed; fi' EXIT",
    paste("3dTcat -overwrite -prefix input.nii.gz", paste(q(inputs), collapse = " ")),
    mask_commands, paste(args, collapse = " \\\n  "), fit)
  for (i in seq_len(nrow(manifest))) {
    lines <- c(lines,
      paste("3dbucket -overwrite -prefix", q(basename(manifest$coef_file[i])), q(paste0("stats.nii.gz[", glt_labels[i], if (estimation == "OLS") "_GLT" else "", "#0_Coef]"))),
      paste("3dbucket -overwrite -prefix", q(basename(manifest$tstat_file[i])), q(paste0("stats.nii.gz[", glt_labels[i], if (estimation == "OLS") "_GLT" else "", "#0_Tstat]"))))
  }
  lines <- c(lines, paste("test -s", q(expected)),
    "rm -f input.nii.gz", paste("rm -f", paste(q(c(masks, mask_links)), collapse = " ")),
    "date -u +%FT%TZ > .afni_complete")
  writeLines(lines, script)
  saveRDS(signature, signature_file)
  list(script = script, output_file = output)
}

#' Inspect AFNI GLM completion and output files
#' @param afni_dir Directory generated by AFNI GLM setup.
#' @param lg Optional logger.
#' @return A one-row status data.frame.
#' @export
get_afni_glm_status <- function(afni_dir, lg = NULL) {
  checkmate::assert_string(afni_dir)
  expected_file <- file.path(afni_dir, "expected_outputs.txt")
  expected <- if (file.exists(expected_file)) readLines(expected_file, warn = FALSE) else character(0)
  outputs <- length(expected) > 0L && all(file.exists(expected)) && all(file.info(expected)$size > 0)
  failed <- file.exists(file.path(afni_dir, ".afni_failed"))
  complete <- outputs && file.exists(file.path(afni_dir, ".afni_complete")) && !failed
  data.frame(afni_complete = complete, afni_failed = failed,
             afni_started = file.exists(file.path(afni_dir, ".afni_start")),
             afni_outputs_present = outputs, to_run = !complete)
}
