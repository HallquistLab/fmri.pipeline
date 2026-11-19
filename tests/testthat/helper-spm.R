# Synthetic output metadata for status and harvesting tests. Image validity itself
# is checked by the MATLAB integration test, not by these text fixtures.
make_spm_completed_fixture <- function(outdir, contrast_names = "cope1", numbers = seq_along(contrast_names)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  source <- file.path(dirname(outdir), paste0(basename(outdir), "_input.txt"))
  writeLines("input", source)
  signature <- fmri.pipeline:::spm_setup_signature(list(settings = list()), source)
  fmri.pipeline:::spm_save_setup(outdir, signature)
  manifest <- data.frame(contrast_number = numbers, contrast_name = contrast_names, stat = "T",
    coef_file = file.path(outdir, sprintf("con_%04d.nii", numbers)),
    stat_file = file.path(outdir, sprintf("spmT_%04d.nii", numbers)))
  outputs <- c(file.path(outdir, c("SPM.mat", "beta_0001.nii", "beta_0002.nii")),
    manifest$coef_file, manifest$stat_file)
  for (f in outputs) writeLines("synthetic output", f)
  utils::write.table(manifest, file.path(outdir, "spm_contrasts.tsv"), sep = "\t", row.names = FALSE)
  writeLines(contrast_names, file.path(outdir, "spm_expected_contrasts.txt"))
  writeLines(outputs, file.path(outdir, "spm_expected_outputs.txt"))
  file.copy(file.path(outdir, "spm_setup_token.txt"), file.path(outdir, "spm_validated_token.txt"), overwrite = TRUE)
  writeLines(rep(as.character(Sys.time()), 2L), file.path(outdir, ".spm_complete"))
  invisible(manifest)
}

make_spm_input_fixture <- function(root, volumes = c(5L, 5L), drops = c(2L, 0L), compressed = FALSE) {
  files <- vapply(seq_along(volumes), function(i) {
    run_dir <- file.path(root, paste0("run", i))
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(run_dir, if (compressed) "func.nii.gz" else "func.nii")
    RNifti::writeNifti(array(i, dim = c(2L, 2L, 2L, volumes[i])), path)
    path
  }, character(1))
  reg <- cbind(trial = 1:2, onset = c(0, 2), duration = c(1, 1), value = c(1, 1))
  attr(reg, "event") <- "A"
  structure(list(design = array(rep(list(reg), length(files)), dim = c(length(files), 1L),
    dimnames = list(paste0("run", seq_along(files)), "A")), design_concat = list(A = reg),
    run_niftis = files, run_volumes = volumes - drops, drop_volumes = drops, tr = 1), class = c("bdm", "list"))
}
