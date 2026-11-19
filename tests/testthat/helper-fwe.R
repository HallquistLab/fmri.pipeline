write_fwe_test_image <- function(file, value = 0, dimensions = c(5L, 5L, 5L)) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  RNifti::writeNifti(array(value, dimensions), file)
  invisible(file)
}

# A valid synthetic worker result for orchestration tests. Numerical validation
# of pTFCE itself is covered separately using the real implementation.
write_fwe_test_result <- function(plan) {
  for (i in seq_len(nrow(plan$tasks))) {
    task <- plan$tasks[i, , drop = FALSE]
    artifacts <- plan$artifacts[plan$artifacts$task_id == task$task_id, , drop = FALSE]
    for (j in seq_len(nrow(artifacts))) {
      file <- artifacts$artifact_file[j]
      dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
      if (artifacts$artifact_role[j] == "threshold_table") {
        utils::write.csv(data.frame(z_ptfce = 4, p_value = task$fwe_alpha[[1L]],
          two_sided = task$two_sided, number_of_resels = 10), file, row.names = FALSE)
      } else write_fwe_test_image(file)
    }
    fmri.pipeline:::validate_ptfce_outputs(task, artifacts)
    saveRDS(list(success = TRUE, input_fingerprint = fmri.pipeline:::fwe_task_fingerprints(task),
      output_hashes = tools::md5sum(artifacts$artifact_file[artifacts$required])),
      fmri.pipeline:::fwe_completion_file(task))
  }
  invisible(plan)
}

make_fwe_execution_test_gpa <- function(root, create_smoothness = TRUE) {
  l3_feat <- file.path(
    root, "feat_l3", "L3m-l3_group", "L1m-facehouse",
    "L2m-l2_session", "l2c-overall", "FEAT_l1c-EV_face.gfeat"
  )
  image_feat_dir <- file.path(l3_feat, "cope1.feat")
  stats_dir <- file.path(image_feat_dir, "stats")
  dir.create(stats_dir, recursive = TRUE)
  for (file in file.path(stats_dir, c("zstat1.nii.gz", "zstat2.nii.gz"))) write_fwe_test_image(file)
  write_fwe_test_image(file.path(image_feat_dir, "mask.nii.gz"), 1)
  if (isTRUE(create_smoothness)) {
    writeLines(
      c("DLH 0.1", "VOLUME 1000", "RESELS 100"),
      file.path(stats_dir, "smoothness")
    )
  }

  l3_setup <- data.frame(
    l1_model = "facehouse",
    l1_cope_name = "EV_face",
    l2_model = "l2_session",
    l2_cope_name = "overall",
    l3_model = "l3_group",
    session = 0L,
    l3_input_mode = "all_sessions",
    feat_fsf = sub("\\.gfeat$", ".fsf", l3_feat),
    feat_dir = l3_feat,
    feat_complete = TRUE,
    feat_failed = FALSE,
    feat_dir_exists = TRUE,
    feat_fsf_exists = FALSE,
    stringsAsFactors = FALSE
  )
  l3_metadata <- data.frame(
    id = c("01", "01"),
    session = c(1L, 1L),
    l1_model = c("facehouse", "facehouse"),
    l1_cope_number = c(1L, 1L),
    l1_cope_name = c("EV_face", "EV_face"),
    l2_model = c("l2_session", "l2_session"),
    l2_cope_number = c(1L, 1L),
    l2_cope_name = c("overall", "overall"),
    l3_model = c("l3_group", "l3_group"),
    l3_cope_number = c(1L, 2L),
    l3_cope_name = c("group", "age"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      output_directory = root,
      scheduler = "slurm",
      lgr_threshold = "warn",
      multi_run = TRUE,
      l3_model_setup = structure(
        list(metadata = l3_metadata, fsl = l3_setup),
        class = c("l3_setup", "list")
      )
    ),
    class = c("glm_pipeline_arguments", "list")
  )
}

make_fake_ptfce_worker <- function(root, fail = FALSE) {
  worker <- file.path(root, if (fail) "fail_ptfce.R" else "fake_ptfce.R")
  if (isTRUE(fail)) {
    lines <- c("#!/usr/bin/env Rscript", "quit(save = 'no', status = 3L)")
  } else {
    lines <- c(
      "#!/usr/bin/env Rscript",
      "args <- commandArgs(trailingOnly = TRUE)",
      "one <- function(flag) args[match(flag, args) + 1L]",
      "zstat <- one('--zstat')",
      "out_dir <- one('--out_dir')",
      "dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)",
      "filename <- basename(zstat)",
      "ext <- if (grepl('\\\\.nii\\\\.gz$', filename)) '.nii.gz' else '.nii'",
      "stem <- substr(filename, 1L, nchar(filename) - nchar(ext))",
      "base <- file.path(out_dir, paste0(stem, '_ptfce'))",
      "RNifti::writeNifti(RNifti::readNifti(zstat), paste0(base, ext))",
      "start <- match('--fwep', args) + 1L",
      "end <- which(seq_along(args) > start & grepl('^--', args))[1L] - 1L",
      "alpha <- as.numeric(args[start:end])",
      "write.csv(data.frame(z_ptfce=4, p_value=alpha, two_sided='--twosided' %in% args, number_of_resels=10), paste0(base, '_zthresh.csv'), row.names=FALSE)",
      "if ('--write_thresh_imgs' %in% args) {",
      "  start <- match('--fwep', args) + 1L",
      "  stop_at <- which(seq_along(args) > start & grepl('^--', args))[1L]",
      "  end <- if (is.na(stop_at)) length(args) else stop_at - 1L",
      "  alpha <- args[start:end]",
      "  for (value in alpha) RNifti::writeNifti(RNifti::readNifti(zstat), paste0(base, '_fwep_', value, ext))",
      "}"
    )
  }
  writeLines(lines, worker)
  Sys.chmod(worker, mode = "0755")
  worker
}
