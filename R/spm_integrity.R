# Shared SPM input alignment, setup signatures and verified output metadata.
spm_nonempty_files <- function(files) {
  length(files) > 0L && all(!is.na(files) & nzchar(files)) &&
    all(file.exists(files)) && all(file.info(files)$size > 0)
}

spm_file_signature <- function(files) {
  files <- unique(normalizePath(files, mustWork = TRUE))
  info <- file.info(files)[, c("size", "mtime"), drop = FALSE]
  # Hash small specifications and nuisance files; large images use size/mtime.
  info$md5 <- NA_character_
  small <- !grepl("\\.nii(\\.gz)?$", files, ignore.case = TRUE)
  info$md5[small] <- unname(tools::md5sum(files[small]))
  info$mtime[small] <- NA # rewriting identical nuisance files does not change the specification
  info
}

spm_setup_signature <- function(spec, files) {
  list(version = 1L, spec = spec, files = spm_file_signature(files))
}

spm_setup_current <- function(outdir, signature = NULL) {
  if (!file.exists(file.path(outdir, "spm_setup.rds"))) return(FALSE)
  saved <- tryCatch(readRDS(file.path(outdir, "spm_setup.rds")), error = function(e) NULL)
  if (is.null(saved) || file.exists(file.path(outdir, ".spm_invalidated"))) return(FALSE)
  if (!is.null(signature)) return(identical(saved$signature, signature))
  tryCatch(identical(saved$signature$files, spm_file_signature(rownames(saved$signature$files))),
           error = function(e) FALSE)
}

spm_invalidate_setup <- function(outdir) {
  unlink(file.path(outdir, c(".spm_complete", ".spm_fail", "spm_expected_outputs.txt",
                            "spm_contrasts.tsv", "spm_validated_token.txt")))
  writeLines("Setup must finish before this model can run.", file.path(outdir, ".spm_invalidated"))
}

spm_save_setup <- function(outdir, signature, syntax = NULL, write_residuals = FALSE) {
  saveRDS(list(signature = signature, syntax = syntax), file.path(outdir, "spm_setup.rds"))
  token <- unname(tools::md5sum(file.path(outdir, "spm_setup.rds")))
  writeLines(token, file.path(outdir, "spm_setup_token.txt"))
  helper <- system.file("matlab", "spm_write_output_manifest.m", package = "fmri.pipeline")
  stopifnot(file.exists(helper))
  esc <- function(x) gsub("'", "''", x, fixed = TRUE)
  writeLines(c(
    if (!is.null(signature$spec$settings$spm_path)) paste0("addpath('", esc(signature$spec$settings$spm_path), "');"),
    paste0("addpath('", esc(dirname(helper)), "');"),
    paste0("spm_write_output_manifest('", esc(outdir), "', ", as.integer(write_residuals), ");")),
    file.path(outdir, "validate_spm_outputs.m"))
  unlink(file.path(outdir, ".spm_invalidated"))
}

spm_read_contrast_manifest <- function(outdir) {
  path <- file.path(outdir, "spm_contrasts.tsv")
  if (!spm_nonempty_files(path)) return(NULL)
  x <- tryCatch(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character"), error = function(e) NULL)
  required <- c("contrast_number", "contrast_name", "stat", "coef_file", "stat_file")
  if (is.null(x) || !all(required %in% names(x)) || !nrow(x) ||
      anyNA(x[, required]) || anyDuplicated(x$contrast_name) || anyDuplicated(x$contrast_number)) return(NULL)
  x
}

spm_outputs_valid <- function(outdir) {
  expected_file <- file.path(outdir, "spm_expected_outputs.txt")
  manifest <- spm_read_contrast_manifest(outdir)
  if (!spm_nonempty_files(expected_file) || is.null(manifest) || !spm_setup_current(outdir)) return(FALSE)
  token <- file.path(outdir, "spm_setup_token.txt")
  validated <- file.path(outdir, "spm_validated_token.txt")
  if (!spm_nonempty_files(c(token, validated)) ||
      !identical(readLines(token, warn = FALSE), readLines(validated, warn = FALSE)) ||
      !identical(readLines(token, warn = FALSE), unname(tools::md5sum(file.path(outdir, "spm_setup.rds"))))) return(FALSE)
  expected <- readLines(expected_file, warn = FALSE)
  required_names <- file.path(outdir, "spm_expected_contrasts.txt")
  if (!spm_nonempty_files(required_names) ||
      !setequal(readLines(required_names, warn = FALSE), manifest$contrast_name)) return(FALSE)
  required <- c(file.path(outdir, "SPM.mat"), manifest$coef_file, manifest$stat_file)
  if (!all(required %in% expected) || !any(grepl("/beta_[0-9]+\\.nii(\\.gz)?$", expected))) return(FALSE)
  spm_nonempty_files(expected)
}

spm_align_inputs <- function(bdm, files, ts_files, concatenate_runs, outdir) {
  actual <- vapply(files, function(f) {
    d <- RNifti::niftiHeader(f)$dim
    if (d[1L] >= 4L) as.integer(d[5L]) else 1L
  }, integer(1))
  n <- bdm$run_volumes %||% actual
  drops <- bdm$drop_volumes %||% 0L
  checkmate::assert_integerish(n, len = length(files), lower = 1, any.missing = FALSE)
  checkmate::assert_integerish(drops, lower = 0, any.missing = FALSE)
  if (!length(drops) %in% c(1L, length(files))) stop("SPM drop_volumes must be scalar or one value per run.")
  drops <- rep_len(drops, length(files))
  offsets <- ifelse(actual == n, 0L, drops)
  if (any(actual != n + offsets)) stop("SPM image lengths do not match retained volumes and leading drops.")
  if (!is.null(ts_files)) {
    checkmate::assert_character(ts_files, any.missing = FALSE, min.chars = 1L)
    checkmate::assert_file_exists(ts_files)
    expected_length <- if (concatenate_runs) 1L else length(files)
    if (length(ts_files) != expected_length) stop("SPM nuisance files do not match the run layout.")
    for (i in seq_along(ts_files)) {
      x <- as.matrix(data.table::fread(ts_files[i], header = FALSE, data.table = FALSE))
      if (!is.numeric(x) || !ncol(x) || any(!is.finite(x))) stop("SPM confounds must contain finite numeric values: ", ts_files[i])
      nn <- if (concatenate_runs) n else n[i]
      dd <- if (concatenate_runs) drops else drops[i]
      if (nrow(x) == sum(nn + dd) && any(dd > 0L)) {
        starts <- c(0L, head(cumsum(nn + dd), -1L))
        keep <- unlist(lapply(seq_along(nn), function(j) starts[j] + dd[j] + seq_len(nn[j])))
        x <- x[keep, , drop = FALSE]
      } else if (nrow(x) != sum(nn)) {
        stop("SPM confound row count disagrees with retained/original volumes: ", ts_files[i])
      }
      aligned <- file.path(outdir, sprintf("spm_confounds_%04d.txt", i))
      utils::write.table(x, aligned, row.names = FALSE, col.names = FALSE, quote = FALSE)
      ts_files[i] <- aligned
    }
  }
  list(lengths = as.integer(n), offsets = as.integer(offsets), ts_files = ts_files)
}

# Use the same stage sequencing and completion rules for immediate and cluster runs.
spm_execute_model <- function(outdir, level, settings, compute_env = character(0)) {
  stage_names <- if (level == 1L) c("setup_glm_design.m", "run_glm.m", "estimate_glm_contrasts.m") else
    c("setup_l3_design.m", "run_l3_glm.m", "estimate_l3_contrasts.m")
  stages <- vapply(c("setup", "glm", "contrasts"), function(stage) {
    isTRUE(settings[[paste0("spm_execute_", stage)]]) ||
      (level == 3L && isTRUE(settings[[paste0("spm_execute_l3_", stage)]]))
  }, logical(1))
  if (!any(stages)) return(invisible(NULL))
  script <- file.path(outdir, "run_spm_now.sh")
  variables <- c(matlab_cmd = settings$matlab_cmd %||% "matlab", matlab_args = settings$matlab_args %||% "-batch",
    spm_setup_script = stage_names[1], spm_glm_script = stage_names[2], spm_contrast_script = stage_names[3],
    run_setup = as.integer(stages[1]), run_glm = as.integer(stages[2]), run_contrasts = as.integer(stages[3]),
    spm_level = level, cleanup_tmp = as.integer(isTRUE(settings$cleanup_tmp)))
  writeLines(c("#!/bin/bash", compute_env, paste0(names(variables), "=", shQuote(unname(variables))),
    "write_job_manifest() { :; }", readLines(system.file("scripts", "spm_runner.sh", package = "fmri.pipeline")),
    paste("spm_runner", shQuote(outdir))), script)
  if (system2("bash", shQuote(script)) != 0L) stop("SPM execution failed in ", outdir)
  invisible(NULL)
}
