# Can also be run on a completed analysis without submitting any new jobs.
check_images <- function(paths) {
  stopifnot(length(paths) > 0L, !anyNA(paths), all(file.exists(paths)))
  for (path in paths) {
    img <- RNifti::readNifti(path)
    if (length(dim(img)) != 3L || any(is.infinite(img)) || !any(is.finite(img) & img != 0)) {
      stop("Invalid or empty 3D statistical map: ", path)
    }
  }
}

check_rows <- function(actual, expected, keys) {
  stopifnot(is.data.frame(actual), all(keys %in% names(actual)), nrow(expected) > 0L,
            !anyDuplicated(actual[keys]), nrow(actual) == nrow(expected),
            !nrow(dplyr::anti_join(expected, actual, by = keys)))
}

verify_clock_run <- function(dir, wait = FALSE) {
  prepared <- readRDS(file.path(dir, "prepared.rds"))
  engine <- prepared$glm_software
  stopifnot(length(engine) == 1L)
  location <- prepared$output_locations
  history <- jsonlite::read_json(location$project_config_json)$batch_history
  stopifnot(length(history) == 1L)
  hours <- as.numeric(Sys.getenv("GLM_TEST_TIMEOUT_HOURS", "24"))
  stopifnot(is.finite(hours), hours > 0)
  deadline <- Sys.time() + as.difftime(hours, units = "hours")
  repeat {
    jobs <- fmri.pipeline::get_tracked_job_status(sequence_id = names(history), sqlite_db = location$sqlite_db)
    failed <- jobs$status %in% c("FAILED", "FAILED_BY_EXT")
    if (any(failed)) stop(engine, " failed: ", paste(jobs$job_name[failed], jobs$job_id[failed], collapse = ", "),
                          ". Logs and scripts: ", history[[1L]]$batch_directory)
    if (any(jobs$job_name == "cleanup_glm" & jobs$status == "COMPLETED")) break
    if (!wait || Sys.time() > deadline) stop("Cleanup incomplete for ", engine, " in ", dir)
    # A process killed by Slurm may never update its SQLite row.
    active <- jobs$job_id[jobs$status %in% c("QUEUED", "STARTED")]
    if (length(active)) {
      states <- fmri.pipeline:::slurm_job_status(active)
      terminal <- sub("[+ ].*$", "", states$State) %in%
        c("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "BOOT_FAIL", "DEADLINE")
      if (any(terminal)) stop(engine, " scheduler failure: ",
                              paste(states$JobID[terminal], states$State[terminal], collapse = ", "))
    }
    Sys.sleep(30)
  }
  fit <- readRDS(location$object_cache)
  stopifnot(identical(fit$batch_run$sequence_id, names(history)))
  runs <- prepared$run_data
  subjects <- prepared$subject_data
  subject_keys <- c("id", "session")
  l1 <- fit$l1_model_setup[[engine]]
  group <- fit$l3_model_setup[[engine]]
  stopifnot(is.data.frame(group), nrow(group) == 1L, identical(group$l1_cope_name, "clock_entropy"),
            identical(group$l3_model, "l3age"))
  group_names <- rownames(prepared$l3_models$models$l3age$contrasts)
  stopifnot(length(group_names) > 0L)
  if (engine == "fsl") {
    check_rows(l1, runs, c(subject_keys, "run_number"))
    l2 <- fit$l2_model_setup$fsl
    check_rows(l2, subjects, subject_keys)
    stopifnot("feat_complete" %in% names(l1), "feat_complete" %in% names(l2),
              all(l1$feat_complete), all(l2$feat_complete), all(group$feat_complete),
              all(l2$l1_cope_name == "clock_entropy"), all(!l2$l2_passthrough),
              all(l2$n_input_files == prepared$n_expected_runs),
              all(l2$feat_actual_estimator == "fe"), identical(group$feat_requested_estimator, "flame12"),
              group$feat_actual_estimator %in% c("flame12", "flame12+flame1", "flame1"))
    maps <- fmri.pipeline::lookup_feat_outputs(fit, level = 1:3, what = c("cope", "tstat"),
      source = "setup", format = "long", include_missing = TRUE)
    expected <- c(nrow(runs), nrow(subjects), length(group_names))
    counts <- table(factor(maps$output_level, levels = 1:3), maps$statistic)
    stopifnot(all(counts[, "cope"] == expected), all(counts[, "tstat"] == expected),
              !anyDuplicated(maps$image_file), setequal(maps$l3_cope_name[maps$output_level == 3L], group_names))
    check_images(maps$image_file)
    group_inputs <- grep("^set feat_files", readLines(group$feat_fsf), value = TRUE)
    group_inputs <- sub('^[^"]*"([^"]+)".*$', '\\1', group_inputs)
    stopifnot(setequal(normalizePath(group_inputs, mustWork = TRUE),
                       normalizePath(maps$image_file[maps$output_level == 2L & maps$statistic == "cope"], mustWork = TRUE)))
    message("FEAT group estimator requested ", group$feat_requested_estimator, ", actual ", group$feat_actual_estimator)
  } else {
    check_rows(l1, subjects, subject_keys)
    spm <- engine == "spm"
    dir_col <- if (spm) "spm_dir" else "afni_dir"
    status <- paste0(engine, "_complete")
    stopifnot(status %in% names(l1), all(l1[[status]]), isTRUE(group[[status]]))
    manifest <- if (spm) "spm_contrasts.tsv" else "contrasts.tsv"
    statistic <- if (spm) "stat_file" else "tstat_file"
    verify_manifest <- function(path, expected) {
      m <- read.delim(file.path(path, manifest), stringsAsFactors = FALSE)
      stopifnot(all(c("contrast_name", "coef_file", statistic) %in% names(m)),
                !anyDuplicated(m$contrast_name), setequal(m$contrast_name, expected))
      check_images(c(m$coef_file, m[[statistic]]))
      m$coef_file
    }
    l1_maps <- vapply(l1[[dir_col]], verify_manifest, character(1), expected = "clock_entropy")
    verify_manifest(group[[dir_col]], group_names)
    if (spm) {
      stopifnot(all(l1$n_runs == prepared$n_expected_runs))
      for (i in seq_len(nrow(l1))) {
        r <- subset(runs, id == l1$id[i] & session == l1$session[i])
        r <- r[order(r$run_number), ]
        stopifnot(all(l1$run_volumes[[i]] == r$run_volumes),
                  identical(as.character(l1$l1_confound_files[[i]]), as.character(r$l1_confound_file)))
      }
      setup <- readRDS(file.path(group$spm_dir, "spm_setup.rds"))
      inputs <- rownames(setup$signature$files)
    } else {
      for (i in seq_len(nrow(l1))) {
        r <- subset(runs, id == l1$id[i] & session == l1$session[i])
        lengths <- r$run_volumes[order(r$run_number)]
        setup <- readRDS(file.path(l1$afni_dir[i], "afni_setup.rds"))
        stopifnot(setup$estimation == "REML", nrow(setup$design) == sum(lengths),
                  identical(as.integer(setup$starts), as.integer(c(0, head(cumsum(lengths), -1L)))))
      }
      setup <- readRDS(file.path(group$afni_dir, "afni_setup.rds"))
      stopifnot(setup$estimation == "OLS")
      inputs <- setup$inputs
    }
    stopifnot(length(inputs) == nrow(subjects),
              setequal(normalizePath(inputs, mustWork = TRUE), normalizePath(l1_maps, mustWork = TRUE)))
  }
  message("PASS ", engine, " (", nrow(subjects), " subjects, ", nrow(runs), " runs): ", dir)
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Supply one or more completed backend directories")
  for (path in args) verify_clock_run(path)
}
