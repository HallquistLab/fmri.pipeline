# FEAT specifications live beside the FSF, outside directories replaced on rerun.
feat_setup_file <- function(fsf) if (length(fsf)) paste0(fsf, ".setup.rds") else character()
feat_lock_dir <- function(fsf) paste0(fsf, ".lock")

feat_file_signature <- function(files) {
  files <- unique(normalizePath(files, winslash = "/", mustWork = FALSE))
  info <- file.info(files)[, c("size", "mtime"), drop = FALSE]
  info$md5 <- rep(NA_character_, length(files))
  small <- !grepl("\\.nii(\\.gz)?$", files, ignore.case = TRUE)
  info$md5[small] <- unname(tools::md5sum(files[small]))
  info$mtime[small] <- NA
  info
}

feat_read_setup <- function(fsf) {
  if (!file.exists(feat_setup_file(fsf))) return(NULL)
  tryCatch(suppressWarnings(readRDS(feat_setup_file(fsf))), error = function(e) NULL)
}

feat_fsf_value <- function(syntax, key) {
  prefix <- paste0("set ", key, " ")
  lines <- trimws(syntax)
  value <- trimws(substring(lines[startsWith(lines, prefix)], nchar(prefix) + 1L))
  if (!length(value)) return(NULL)
  sub('"$', "", sub('^"', "", tail(value, 1L)))
}

feat_output_directory <- function(syntax) {
  path <- feat_fsf_value(syntax, "fmri(outputdir)")
  if (is.null(path) || !nzchar(path)) stop("FEAT FSF is missing fmri(outputdir).")
  suffix <- if (identical(feat_fsf_value(syntax, "fmri(level)"), "1")) ".feat" else ".gfeat"
  if (!endsWith(path, suffix)) path <- paste0(path, suffix)
  file.path(normalizePath(dirname(path), winslash = "/", mustWork = FALSE), basename(path))
}

feat_upstream_setups <- function(inputs) {
  dirs <- sub("^(.*\\.g?feat)/.*$", "\\1", inputs)
  # A higher-level child (cope1.feat) belongs to the outer .gfeat fit.
  dirs <- sub("^(.*\\.gfeat)/.*$", "\\1", dirs)
  fsfs <- sub("\\.g?feat$", ".fsf", dirs)
  unique(fsfs[file.exists(feat_setup_file(fsfs))])
}

feat_higher_inputs <- function(copes) {
  # FEAT reads paired variance and degrees-of-freedom inputs implicitly.
  variances <- sub("/cope([0-9]+)(\\.nii(?:\\.gz)?)$", "/varcope\\1\\2", copes, perl = TRUE)
  dof <- unique(unlist(lapply(dirname(copes), function(path) {
    list.files(path, pattern = "^(dof|tdof_t[0-9]+(\\.nii(\\.gz)?)?)$", full.names = TRUE)
  })))
  unique(c(copes, variances, dof))
}

feat_setup_current <- function(fsf, seen = character()) {
  fsf <- normalizePath(fsf, winslash = "/", mustWork = FALSE)
  if (fsf %in% seen) return(FALSE)
  setup <- feat_read_setup(fsf)
  if (is.null(setup) || !identical(setup$version, 1L) || !file.exists(fsf)) return(FALSE)
  if (!identical(readLines(fsf, warn = FALSE), setup$syntax)) return(FALSE)
  if (any(!file.exists(rownames(setup$files))) ||
      !identical(setup$files, feat_file_signature(rownames(setup$files)))) return(FALSE)
  all(vapply(setup$upstream, function(up) {
    saved <- feat_read_setup(up)
    !is.null(saved) && file.exists(file.path(saved$output_dir, ".feat_complete")) &&
      feat_outputs_valid(saved$output_dir, up, c(seen, fsf))
  }, logical(1)))
}

feat_save_setup <- function(fsf, syntax, inputs, level, spec = NULL, force = FALSE) {
  # Resolve the parent even before the FSF exists (e.g. symlinked scratch).
  fsf <- file.path(normalizePath(dirname(fsf), winslash = "/"), basename(fsf))
  token <- basename(tempfile("feat_setup_"))
  feat_reserve(fsf, token)
  on.exit(feat_release(fsf, token), add = TRUE)
  output_dir <- feat_output_directory(syntax)
  # The FSF names the fit throughout setup, submission, and downstream provenance.
  suffix <- if (level == 1L) ".feat" else ".gfeat"
  if (!identical(output_dir, sub("\\.fsf$", suffix, fsf))) {
    stop("FEAT outputdir must match the configured model path. Set the FEAT output location in output_locations: ", fsf)
  }
  upstream <- feat_upstream_setups(inputs)
  files <- unique(c(inputs, feat_setup_file(upstream)))
  setup <- list(version = 1L, syntax = syntax, files = feat_file_signature(files),
    upstream = upstream, output_dir = output_dir, level = level, spec = spec)
  old <- feat_read_setup(fsf)
  changed <- !identical(old, setup) || !file.exists(fsf) ||
    !identical(readLines(fsf, warn = FALSE), syntax)
  if (changed || isTRUE(force)) {
    if (dir.exists(output_dir)) writeLines("Specification changed; previous fit retained until execution.", file.path(output_dir, ".feat_stale"))
    unlink(file.path(output_dir, c(".feat_complete", ".feat_validated.rds")))
    writeLines(syntax, fsf)
    saveRDS(setup, feat_setup_file(fsf))
  }
  invisible(changed || isTRUE(force))
}

feat_outputs_valid <- function(output_dir, fsf, seen = character()) {
  tryCatch({
    if (file.exists(file.path(output_dir, ".feat_fail")) || !feat_setup_current(fsf, seen)) return(FALSE)
    record <- file.path(output_dir, ".feat_validated.rds")
    if (!file.exists(record)) return(FALSE)
    saved <- readRDS(record)
    expected <- feat_read_setup(fsf)
    isTRUE(identical(normalizePath(output_dir, winslash = "/", mustWork = FALSE), expected$output_dir) &&
      identical(saved$setup_md5, unname(tools::md5sum(feat_setup_file(fsf)))) &&
      nrow(saved$files) > 0L && all(file.exists(rownames(saved$files))) &&
      identical(saved$files, feat_file_signature(rownames(saved$files))))
  }, error = function(e) FALSE)
}

feat_read_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- which(trimws(lines) == "/Matrix")
  if (length(start) != 1L || start >= length(lines)) stop("Invalid FEAT matrix: ", path)
  x <- as.matrix(utils::read.table(text = paste(lines[(start + 1L):length(lines)], collapse = "\n")))
  if (!is.numeric(x) || any(!is.finite(x))) stop("Nonfinite FEAT matrix: ", path)
  header <- function(key) as.integer(sub(paste0("^/", key, "[[:space:]]+"), "",
    grep(paste0("^/", key, "[[:space:]]"), lines, value = TRUE)))
  if (!identical(header("NumWaves"), ncol(x)) ||
      !identical(c(header("NumPoints"), header("NumContrasts")), nrow(x))) {
    stop("FEAT matrix dimensions do not match its header: ", path)
  }
  x
}

# Called once after estimation; subsequent status checks verify the recorded files.
feat_validate_outputs <- function(fsf) {
  if (!feat_setup_current(fsf)) stop("FEAT setup is missing or stale: ", fsf)
  setup <- feat_read_setup(fsf)
  out <- setup$output_dir
  fit <- if (setup$level == 1L) out else file.path(out, "cope1.feat")
  design_files <- file.path(fit, c("design.mat", "design.con"))
  design <- feat_read_matrix(design_files[1L])
  contrasts <- feat_read_matrix(design_files[2L])
  n <- as.integer(feat_fsf_value(setup$syntax, "fmri(ncon_real)"))
  expected_names <- vapply(seq_len(n), function(i) feat_fsf_value(setup$syntax, paste0("fmri(conname_real.", i, ")")), character(1))
  con_lines <- readLines(design_files[2L])
  actual_names <- vapply(seq_len(n), function(i) {
    line <- grep(paste0("^/ContrastName", i, "[[:space:]]"), con_lines, value = TRUE)
    if (length(line) != 1L) return(NA_character_)
    trimws(sub(paste0("^/ContrastName", i, "[[:space:]]+"), "", line))
  }, character(1))
  if (nrow(design) != as.integer(feat_fsf_value(setup$syntax, "fmri(npts)")) ||
      nrow(contrasts) != n || ncol(contrasts) != ncol(design) ||
      !identical(actual_names, expected_names)) stop("FEAT design/contrast outputs do not match setup: ", out)
  weights <- setup$spec$contrasts
  if (!is.null(weights)) {
    # FEAT appends nuisance EV columns to the contrast matrix with zero weight.
    expected <- matrix(0, nrow(contrasts), ncol(contrasts))
    if (nrow(weights) != nrow(expected) || ncol(weights) > ncol(expected)) {
      stop("FEAT contrast dimensions differ from the specification: ", out)
    }
    expected[, seq_len(ncol(weights))] <- as.matrix(weights)
    if (!isTRUE(all.equal(as.numeric(contrasts), as.numeric(expected), tolerance = 1e-6))) {
      stop("FEAT contrast weights differ from the specification: ", out)
    }
  }
  if (setup$level > 1L && !is.null(setup$spec$design)) {
    # The higher-level FSF writer rounds EV values to five decimal places.
    expected <- round(as.matrix(setup$spec$design), 5)
    if (!identical(dim(expected), dim(design)) ||
        !isTRUE(all.equal(as.numeric(design), as.numeric(expected), tolerance = 1e-6))) {
      stop("FEAT higher-level design differs from the specification: ", out)
    }
  }
  image <- function(stem) {
    files <- paste0(stem, c(".nii.gz", ".nii"))
    found <- files[file.exists(files)]
    if (length(found) != 1L) stop("Missing or ambiguous FEAT output: ", stem)
    found
  }
  mask_file <- image(file.path(fit, "mask"))
  read_image <- function(path) withCallingHandlers(RNifti::readNifti(path),
    warning = function(w) stop("Unreadable FEAT image: ", path, ": ", conditionMessage(w)))
  mask <- read_image(mask_file)
  if (any(!is.finite(mask)) || !any(mask > 0)) stop("Invalid FEAT mask: ", mask_file)
  statistics <- c("cope", "varcope", "tstat", "zstat")
  if (setup$level > 1L) statistics <- c(statistics, "tdof_t")
  stems <- unlist(lapply(statistics, function(stat) paste0(stat, seq_len(n))))
  files <- vapply(file.path(fit, "stats", stems), image, character(1))
  for (file in files) {
    x <- read_image(file)
    if (!identical(dim(x), dim(mask)) || any(!is.finite(x[mask > 0]))) stop("Invalid FEAT statistic: ", file)
    if (!isTRUE(all.equal(unname(RNifti::xform(mask)), unname(RNifti::xform(x)), tolerance = 1e-4))) {
      stop("FEAT statistic has incompatible geometry: ", file)
    }
    if (startsWith(basename(file), "varcope") && any(x[mask > 0] < 0)) stop("Negative FEAT variance: ", file)
  }
  if (setup$level == 1L) {
    dof_file <- file.path(fit, "stats", "dof")
    dof <- scan(dof_file, quiet = TRUE)
    if (!length(dof) || any(!is.finite(dof) | dof <= 0)) stop("Invalid FEAT degrees of freedom.")
    files <- c(files, dof_file)
  }
  saveRDS(list(setup_md5 = unname(tools::md5sum(feat_setup_file(fsf))),
    files = feat_file_signature(c(design_files, mask_file, files))), file.path(out, ".feat_validated.rds"))
  invisible(out)
}

feat_validate_l1_inputs <- function(files, volumes, confounds, confounds_requested) {
  checkmate::assert_integerish(volumes, lower = 1, len = length(files), any.missing = FALSE)
  for (i in seq_along(files)) {
    d <- RNifti::niftiHeader(files[i])$dim
    actual <- if (d[1L] >= 4L) d[5L] else 1L
    if (actual != volumes[i]) stop("FEAT image/design volume mismatch for ", files[i],
      ": image has ", actual, ", design has ", volumes[i], ". Apply leading-volume drops to the image before FEAT setup.")
    conf <- confounds[i]
    supplied <- !is.na(conf) && nzchar(conf)
    if (isTRUE(confounds_requested) && (!supplied || !file.exists(conf))) {
      stop("Missing requested FEAT confounds; refusing to omit nuisance adjustment: ", conf)
    }
    if (supplied) {
      if (!file.exists(conf)) stop("Missing FEAT confounds: ", conf)
      x <- as.matrix(data.table::fread(conf, header = FALSE, data.table = FALSE))
      if (!is.numeric(x) || !ncol(x) || nrow(x) != volumes[i] || any(!is.finite(x))) {
        stop("FEAT confounds must have finite numeric values and one row per retained volume: ", conf)
      }
    }
  }
  invisible(TRUE)
}

# Cleanup happens only while the model's atomic directory lock is held.
feat_prepare_run <- function(fsf, action) {
  if (!feat_setup_current(fsf)) stop("FEAT setup is missing or stale; run model setup again: ", fsf)
  setup <- feat_read_setup(fsf)
  out <- setup$output_dir
  if (dir.exists(out)) {
    # A queued cleanup must never delete an output completed since status was read.
    if (identical(action, "delete") && (file.exists(file.path(out, ".feat_complete")) ||
        file.exists(file.path(out, ".feat_stale")))) action <- "archive"
    if (identical(action, "keep")) stop("Keeping existing FEAT output; choose archive or delete to rerun: ", out)
    if (identical(action, "archive")) {
      target <- tempfile(paste0(basename(out), ".archive_"), tmpdir = dirname(out))
      if (!file.rename(out, target)) stop("Could not archive FEAT output: ", out)
    } else if (identical(action, "delete")) {
      unlink(out, recursive = TRUE)
    } else stop("Unknown FEAT folder action: ", action)
    if (dir.exists(out)) stop("FEAT output directory still exists: ", out)
  }
  invisible(out)
}

feat_reserve <- function(fsf, token) {
  lock <- feat_lock_dir(fsf)
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("FEAT model is already queued/running (or has an interrupted-job lock): ", lock,
      ". Confirm the owning job has stopped before removing the lock.")
  }
  tryCatch(writeLines(token, file.path(lock, "owner")), error = function(e) {
    unlink(lock, recursive = TRUE)
    stop(e)
  })
  invisible(lock)
}

feat_release <- function(fsf, token) {
  owner <- file.path(feat_lock_dir(fsf), "owner")
  if (file.exists(owner) && identical(readLines(owner, warn = FALSE), token)) unlink(feat_lock_dir(fsf), recursive = TRUE)
}

feat_run_action <- function(status, failed = "delete", incomplete = "delete") {
  action <- if (isTRUE(status$feat_failed)) failed else incomplete
  checkmate::assert_choice(action, c("keep", "archive", "delete"))
  if (dir.exists(status$feat_dir) && action != "keep" &&
      (isTRUE(status$feat_complete) || file.exists(file.path(status$feat_dir, ".feat_stale")))) action <- "archive"
  action
}

feat_execute_local <- function(fsf, gpa) {
  token <- basename(tempfile("feat_local_"))
  feat_reserve(fsf, token)
  # Once launched, the worker owns the lock and handles interruptions.
  launched <- FALSE
  on.exit(if (!launched) feat_release(fsf, token), add = TRUE)
  setup <- feat_read_setup(fsf)
  status <- get_feat_status(setup$output_dir, fsf)
  if (isTRUE(status$feat_complete)) return(invisible(0L))
  action <- feat_run_action(status,
    gpa$glm_settings$fsl$failed_l1_folder_action %||% "delete",
    gpa$glm_settings$fsl$incomplete_l1_folder_action %||% "delete")
  if (status$feat_dir_exists && action == "keep") return(invisible(0L))
  script <- tempfile(fileext = ".sh")
  on.exit(unlink(script), add = TRUE)
  writeLines(c("#!/bin/bash", get_compute_environment(gpa, c("fsl", "r")),
    "feat_binary=feat", "feat_use_retry=0", "feat_level=1", "write_job_manifest() { :; }",
    paste0("feat_integrity_script=", shQuote(system.file("bin", "feat_integrity.R", package = "fmri.pipeline"))),
    paste0("feat_reg_script=", shQuote(system.file("bash", "gen_feat_reg_dir", package = "fmri.pipeline"))),
    paste("source", shQuote(system.file("bash", "feat_runner.sh", package = "fmri.pipeline"))),
    paste("feat_runner", paste(shQuote(c(fsf, setup$output_dir, token, action)), collapse = " "))), script)
  launched <- TRUE
  result <- system2("bash", shQuote(script))
  if (result != 0L) stop("FEAT execution or output validation failed: ", fsf)
  invisible(result)
}
