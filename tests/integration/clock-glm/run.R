# Rscript --vanilla tests/integration/clock-glm/run.R <new-shared-directory> [full|smoke] [fsl,spm,afni] [run|prepare]
suppressPackageStartupMessages(library(fmri.pipeline))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 4L) stop("Supply a new output directory, optional profile, backends and action.")
profile <- match.arg(if (length(args) > 1L) args[2L] else "full", c("full", "smoke"))
engines <- if (length(args) > 2L) strsplit(args[3L], ",", fixed = TRUE)[[1L]] else c("fsl", "spm", "afni")
stopifnot(length(engines) > 0L, !anyDuplicated(engines), all(engines %in% c("fsl", "spm", "afni")))
action <- match.arg(if (length(args) == 4L) args[4L] else "run", c("run", "prepare"))
data_root <- normalizePath(Sys.getenv("FMRI_TEST_DATA", "/proj/mnhallqlab/projects/fmri.pipeline_test_data"), mustWork = TRUE)
base <- readRDS(file.path(data_root, "gpa_base.rds"))
spec <- yaml::read_yaml(file.path(data_root, "sample_2_L1_spec.yaml"))
spm_path <- Sys.getenv("SPM_PATH", "/proj/mnhallqlab/lab_resources/spm12")
if ("spm" %in% engines) stopifnot(file.exists(file.path(spm_path, "spm.m")), nzchar(Sys.which("matlab")))
stopifnot(nzchar(Sys.which("sbatch")))
if (profile == "smoke") stopifnot(nzchar(Sys.which("fslroi")))
parent <- normalizePath(dirname(args[1L]), mustWork = TRUE)
if (parent == data_root || startsWith(parent, paste0(data_root, "/"))) stop("Choose an output directory outside the fixture.")
if (!dir.create(args[1L], recursive = TRUE, mode = "0700")) stop("Output must be a new directory: ", args[1L])
root <- normalizePath(args[1L])
source_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L]), mustWork = TRUE))
setwd(root)

# Use the example tables, never its historical analysis directories, caches or fitted objects.
runs <- as.data.frame(base$run_data)
stopifnot(!anyDuplicated(runs[c("id", "session", "run_number")]), nrow(runs) == 80L,
          length(unique(runs$id)) == 10L, all(table(runs$id) == 8L))
if (profile == "smoke") {
  runs <- subset(runs, id %in% head(sort(unique(id)), 4L) & run_number %in% 1:2)
  stopifnot(nrow(runs) == 8L)
}
subjects <- dplyr::semi_join(base$subject_data, runs, by = c("id", "session"))
trials <- dplyr::semi_join(base$trial_data, runs, by = c("id", "session", "run_number"))
keys <- c("id", "session", "run_number")
stopifnot(!nrow(dplyr::anti_join(runs, trials, by = keys)), !anyDuplicated(subjects[c("id","session")]))
for (i in seq_len(nrow(runs))) {
  # The saved absolute paths may refer to another mount. Preserve the fixture's
  # participant/date, preprocessing and run directory components.
  run_dir <- runs$mr_dir[i]
  src <- file.path(data_root, basename(dirname(dirname(run_dir))), basename(dirname(run_dir)), basename(run_dir))
  dst <- file.path(root, "inputs", as.character(runs$id[i]), paste0("ses-", runs$session[i]), basename(src))
  stopifnot(dir.create(dst, recursive = TRUE))
  bold <- basename(runs$run_nifti[i])
  nuisance <- basename(runs$confound_input_file[i])
  files <- file.path(src, c(bold, nuisance, base$confound_settings$motion_params_file))
  stopifnot(all(file.exists(files)))
  if (profile == "smoke") {
    # Keep acquisition timing, affine and the entire time series. Crop a small
    # central region from the already aligned MNI images for all three engines.
    shape <- RNifti::niftiHeader(files[1L])$dim[2:4]
    stopifnot(all(shape >= 17L))
    roi <- as.vector(rbind(floor((shape - 17L) / 2L), 17L))
    stopifnot(system2("fslroi", c(shQuote(files[1L]), shQuote(file.path(dst, bold)), roi),
                      env = "FSLOUTPUTTYPE=NIFTI_GZ") == 0L)
    files <- files[-1L]
  }
  stopifnot(all(file.symlink(files, file.path(dst, basename(files)))))
  runs$mr_dir[i] <- dst
  runs$run_nifti[i] <- file.path(dst, bold)
  runs$confound_input_file[i] <- file.path(dst, nuisance)
}
# One parametric contrast keeps the number of real group fits manageable.
spec$l1_models <- spec$l1_models["entropy"]
spec$l1_models$entropy$contrasts <- list(diagonal = FALSE, custom = list(clock_entropy = list(clock = 0, feedback = 0, entropy_clock = 1)))
spec$signals <- spec$signals[c("clock", "feedback", "entropy_clock")]
spec$l2_models <- list(run_mean = list(level = 2L, model_formula = "~ 1", diagonal = TRUE))
spec_file <- file.path(root, "clock.yaml")
yaml::write_yaml(spec, spec_file)
write.csv(runs[c(keys, "run_nifti")], "selected_runs.csv", row.names = FALSE)
conf <- base$confound_settings
# Exercise leading-volume removal and the four measured nuisance regressors.
# Use a fixed cohort so dropped runs/participants are detectable by the tests.
conf[c("exclude_run", "exclude_subject", "truncate_run")] <- list(NULL)
shell_config <- file.path(root, "fsl-sub.yml")
writeLines("method: shell", shell_config)
# The workers must use the same installed package and tools as the driver.
parallel <- list(pipeline_cores = 1L, finalize_cores = 1L, l1_setup_cores = 1L, l2_setup_cores = 1L,
  fsl = list(l1_feat_cpus_per_job = 2L, l1_feat_runs_per_cpu = 1L, l3_feat_cpusperjob = 2L),
  spm = list(l1_spm_cpus_per_job = 1L, l3_spm_cpus_per_job = 1L),
  afni = list(l1_glm_cpus = 2L),
  compute_environment = list(global = c(
    paste("export R_LIBS=", shQuote(paste(.libPaths(), collapse = ":")), sep = ""),
    paste("export PATH=", shQuote(Sys.getenv("PATH")), sep = ""),
    paste("export FSLSUB_CONF=", shQuote(shell_config), sep = "")),
    r = character(), spm = character(), fsl = character(), afni = character()))

prepare_and_submit <- function(engine) {
  gpa <- setup_glm_pipeline(analysis_name = paste0("clock_", engine), output_directory = file.path(root, engine),
    scheduler = "slurm", subject_data = subjects, run_data = runs, trial_data = trials,
    drop_volumes = base$drop_volumes, n_expected_runs = as.integer(nrow(runs) / nrow(subjects)),
    confound_settings = conf, use_preconvolve = TRUE, glm_software = engine, parallel = parallel,
    l1_models = NULL, l2_models = NULL, l3_models = NULL,
    glm_settings = list(spm = list(spm_path = spm_path, require_matlab = TRUE,
      l3_model_type = "mreg", concatenate_runs = FALSE, fmri_t = 16L, fmri_t0 = 8L),
      afni = list(estimation = "REML", write_residuals = FALSE), fsl = list()),
    additional = list(feat_l1_args = list(thresh = 0),
      feat_l2_args = list(thresh = 0, featwatcher_yn = 0),
      feat_l3_args = list(thresh = 0, featwatcher_yn = 0)))
  gpa <- build_l1_models(gpa, from_spec_file = spec_file)
  if (engine == "fsl") gpa <- build_l2_models(gpa, from_spec_file = spec_file)
  gpa <- build_l3_models(gpa, from_spec_file = spec_file)
  stopifnot(identical(rownames(gpa$l1_models$models$entropy$contrasts), "clock_entropy"))
  # Finalize synchronously so the exact cohort, confounds and design inputs are
  # checked before queuing expensive fits. Other stages run via normal Slurm jobs.
  gpa <- finalize_pipeline_configuration(gpa)
  stopifnot(isTRUE(gpa$finalize_complete), nrow(gpa$run_data) == nrow(runs),
            !any(gpa$run_data$exclude_run), !any(gpa$subject_data$exclude_subject),
            all(file.exists(gpa$run_data$l1_confound_file)))
  for (i in seq_len(nrow(runs))) {
    nuisance <- as.matrix(read.table(gpa$run_data$l1_confound_file[i]))
    stopifnot(ncol(nuisance) == length(conf$l1_confound_regressors),
              nrow(nuisance) == gpa$run_data$run_volumes[i],
              gpa$run_data$orig_volumes[i] - gpa$run_data$run_volumes[i] == base$drop_volumes,
              all(is.finite(nuisance)))
  }
  saveRDS(gpa, file.path(root, engine, "prepared.rds"))
  if (action == "prepare") return(invisible(TRUE))
  run_glm_pipeline(gpa, l1_model_names = "entropy", l2_model_names = if (engine == "fsl") "run_mean" else "none",
                   l3_model_names = "l3age")
}
source(file.path(source_dir, "validate.R"))
run_each <- function(selected, fn) {
  vapply(selected, function(engine) tryCatch({
    fn(engine)
    "ok"
  }, error = function(e) {
    message(engine, ": ", conditionMessage(e))
    conditionMessage(e)
  }), character(1))
}
status <- run_each(engines, prepare_and_submit)
if (action == "run") {
  submitted <- names(status)[status == "ok"]
  status[submitted] <- run_each(submitted, function(engine) verify_clock_run(file.path(root, engine), wait = TRUE))
}
write.csv(data.frame(backend = names(status), action, status = unname(status)), "results.csv", row.names = FALSE)
if (any(status != "ok")) stop("Clock GLM suite failed. See ", file.path(root, "results.csv"))
