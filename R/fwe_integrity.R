# Version the numerical/worker contract independently of the serialized plan.
fwe_implementation_version <- function() 2L

fwe_hash <- function(x) {
  path <- tempfile("fwe_hash_")
  on.exit(unlink(path))
  saveRDS(x, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

fwe_task_fingerprints <- function(tasks) {
  fields <- intersect(c("zstat_file", "correction_mask_file", "model_mask_file", "smoothness_file"), names(tasks))
  files <- unique(unlist(tasks[, fields, drop = FALSE], use.names = FALSE))
  hashes <- stats::setNames(rep(NA_character_, length(files)), files)
  present <- !is.na(files) & file.exists(files)
  hashes[present] <- unname(tools::md5sum(files[present]))
  vapply(seq_len(nrow(tasks)), function(i) {
    paths <- unlist(tasks[i, fields, drop = FALSE], use.names = FALSE)
    if (anyNA(hashes[paths])) return(NA_character_)
    fwe_hash(list(
      implementation = fwe_implementation_version(),
      inputs = hashes[paths],
      two_sided = tasks$two_sided[i],
      fwe_alpha = tasks$fwe_alpha[[i]],
      write_threshold_images = tasks$write_threshold_images[i]
    ))
  }, character(1L))
}

fwe_completion_file <- function(task) file.path(task$output_directory, "completion.rds")

fwe_task_is_complete <- function(task, artifacts) {
  tryCatch({
    if (is.na(task$input_fingerprint) || !file.exists(fwe_completion_file(task))) return(FALSE)
    receipt <- readRDS(fwe_completion_file(task))
    files <- artifacts$artifact_file[artifacts$required]
    isTRUE(receipt$success) &&
      identical(receipt$input_fingerprint, task$input_fingerprint) &&
      length(files) > 0L && all(file.exists(files)) &&
      all(file.info(files)$size > 0) &&
      identical(receipt$output_hashes, tools::md5sum(files))
  }, error = function(e) FALSE)
}

validate_ptfce_outputs <- function(task, artifacts) {
  reference <- roi_image_geometry(task$zstat_file)
  inside <- ptfce_read_mask(task$correction_mask_file, task$zstat_file)
  for (i in seq_len(nrow(artifacts))) {
    file <- artifacts$artifact_file[i]
    if (!file.exists(file) || file.info(file)$size == 0) stop("Missing or empty FWE output: ", file)
    if (artifacts$artifact_role[i] == "threshold_table") {
      tab <- utils::read.csv(file)
      checkmate::assert_names(names(tab), must.include = c("z_ptfce", "p_value", "two_sided", "number_of_resels"))
      if (!isTRUE(all.equal(sort(tab$p_value), sort(task$fwe_alpha[[1L]]), tolerance = 1e-12)) ||
          anyNA(tab$two_sided) || !all(tab$two_sided == task$two_sided) ||
          !all(is.finite(tab$z_ptfce) & tab$z_ptfce >= 0) ||
          !all(is.finite(tab$number_of_resels) & tab$number_of_resels > 0)) {
        stop("Invalid pTFCE threshold table: ", file)
      }
    } else {
      validate_roi_geometry(reference, roi_image_geometry(file), task$zstat_file, file)
      values <- RNifti::readNifti(file)
      if (any(!is.finite(values))) stop("Non-finite values in FWE output: ", file)
      if (any(values[!inside] != 0)) stop("FWE output extends outside the correction mask: ", file)
    }
  }
  invisible(TRUE)
}

# Called in a fresh process for both local and scheduled runs. Only this wrapper
# commits a receipt, after the worker exits successfully and outputs validate.
execute_fwe_task <- function(request_file, worker_script, arguments) {
  request <- readRDS(request_file)
  task <- request$task
  receipt_file <- fwe_completion_file(task)
  unlink(receipt_file)
  if (!identical(fwe_task_fingerprints(task), task$input_fingerprint)) {
    stop("FWE inputs changed after submission; replan the correction.")
  }
  unlink(request$artifacts$artifact_file)
  status <- system2(rscript_executable(), vapply(c(worker_script, arguments), shQuote, character(1L)))
  if (status != 0L) stop(structure(
    list(message = paste("FWE worker failed with exit status", status), call = NULL, status = status),
    class = c("fwe_worker_error", "error", "condition")))
  validate_ptfce_outputs(task, request$artifacts)
  if (!identical(fwe_task_fingerprints(task), task$input_fingerprint)) {
    stop("FWE inputs changed during execution; replan the correction.")
  }
  files <- request$artifacts$artifact_file[request$artifacts$required]
  receipt <- list(
    success = TRUE, input_fingerprint = task$input_fingerprint,
    output_hashes = tools::md5sum(files), completed_at = Sys.time(),
    implementation = fwe_implementation_version(),
    ptfce_version = if (requireNamespace("pTFCE", quietly = TRUE)) as.character(utils::packageVersion("pTFCE")) else NA_character_,
    worker_hash = unname(tools::md5sum(worker_script))
  )
  fwe_atomic_write(receipt_file, function(path) saveRDS(receipt, path))
  invisible(receipt)
}
