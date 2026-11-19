# Synthetic, readable FEAT outputs for status/cleanup tests; no FSL required.
make_valid_feat <- function(root, name = "model", level = 1L, inputs = character()) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  fsf <- file.path(root, paste0(name, ".fsf"))
  out <- file.path(root, paste0(name, if (level == 1L) ".feat" else ".gfeat"))
  syntax <- c(paste('set fmri(outputdir)', shQuote(out, type = "cmd")),
    paste('set fmri(level)', if (level == 1L) 1 else 2),
    'set fmri(npts) 4', 'set fmri(ncon_real) 1', 'set fmri(conname_real.1) "effect"')
  fmri.pipeline:::feat_save_setup(fsf, syntax, inputs, level)
  fit <- if (level == 1L) out else file.path(out, "cope1.feat")
  dir.create(file.path(fit, "stats"), recursive = TRUE)
  writeLines(c('/NumWaves 1', '/NumPoints 4', '/Matrix', '1', '1', '1', '1'), file.path(fit, 'design.mat'))
  writeLines(c('/ContrastName1 effect', '/NumWaves 1', '/NumContrasts 1', '/Matrix', '1'), file.path(fit, 'design.con'))
  image <- array(1, c(3L, 3L, 3L))
  RNifti::writeNifti(image, file.path(fit, 'mask.nii.gz'))
  for (stat in c('cope', 'varcope', 'tstat', 'zstat', if (level > 1L) 'tdof_t')) {
    RNifti::writeNifti(image, file.path(fit, 'stats', paste0(stat, '1.nii.gz')))
  }
  writeLines('3', file.path(fit, 'stats', 'dof'))
  fmri.pipeline:::feat_validate_outputs(fsf)
  writeLines(rep('2026-01-01 10:00:00', 2), file.path(out, '.feat_complete'))
  list(fsf = fsf, out = out, fit = fit, syntax = syntax)
}

make_feat_gpa <- function(fit) {
  gpa <- create_mock_gpa(include_l1_models = TRUE, output_directory = dirname(fit$out))
  gpa$lgr_threshold <- 'error'
  gpa$scheduler <- 'slurm'
  gpa$output_locations$sqlite_db <- NULL
  gpa$output_locations$scheduler_scripts <- file.path(dirname(fit$out), 'scripts')
  gpa$output_locations$log_directory <- file.path(dirname(fit$out), 'logs')
  gpa$parallel$fsl <- list(l1_feat_time = '0:10:00', l1_feat_memgb = 1,
    l1_feat_cpus_per_job = 1L, l1_feat_runs_per_cpu = 1L)
  gpa$glm_settings <- list(fsl = list(failed_l1_folder_action = 'delete', incomplete_l1_folder_action = 'delete'))
  gpa$l1_model_setup <- structure(list(fsl = data.frame(l1_model = 'model1',
    feat_fsf = fit$fsf, feat_dir = fit$out, feat_complete = FALSE, feat_failed = FALSE, to_run = TRUE)), class = 'l1_setup')
  gpa
}
