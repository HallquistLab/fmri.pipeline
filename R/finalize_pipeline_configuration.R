#' This is a small helper function to validate the glm_model_arguments list structure.
#' It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
#' N.B. gpa is a shorthand abbreviation for glm_model_arguments, to save typing
#'
#' @param gpa A \code{glm_pipeline_arguments} object setup by \code{setup_glm_pipeline}
#' @param refinalize A logical indicating whether to force checks and finalize steps on an object that
#'   was previously finalized.
#' @importFrom stringr str_count fixed
#' @importFrom magrittr %>%
#' @importFrom lgr get_logger
#' @importFrom DBI dbConnect dbIsValid dbWriteTable dbDisconnect
#' @importFrom RSQLite SQLite
#' @export
finalize_pipeline_configuration <- function(gpa, refinalize = FALSE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_logical(refinalize, len = 1L)
  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  if (isTRUE(gpa$finalize_complete) && isFALSE(refinalize)) {
    lg$debug("In finalize_pipeline_configuration, finalization of gpa already complete. Returning object unchanged.")
    return(gpa)
  }

  lg$debug("In finalize_pipeline_configuration, opening SQLite connection to: %s", gpa$output_locations$sqlite_db)
  sqlite_con <- DBI::dbConnect(RSQLite::SQLite(), gpa$output_locations$sqlite_db)

  # final checks on compute environment now that we're running inside the compute environment
  test_compute_environment(gpa, what = unique(c("r", gpa$glm_software)), stop_on_fail = TRUE)
  if ("afni" %in% gpa$glm_software) test_afni_glm_compute_environment(gpa)
  if ("spm" %in% gpa$glm_software) {
    test_spm_compute_environment(gpa, stop_on_fail = TRUE)
  }

  # l1 models must be specified to get started (hard enforcement)
  if (!checkmate::test_class(gpa$l1_models, "l1_model_set")) {
    msg <- "Could not find valid $l1_models specification in gpa. Be sure to run build_l1_models() before proceeding!"
    lg$error(msg)
    stop(msg)
  }

  # l2 models are optional and depend on whether it's a multi-run setup
  if (isTRUE(gpa$multi_run) && "fsl" %in% get_level_backend_names(gpa, 2L) &&
      !checkmate::test_class(gpa$l2_models, "hi_model_set")) {
    msg <- "Could not find valid $l2_models specification in gpa. We will proceed, but weird things may happen! I suggest running build_l2_models()."
    lg$warn(msg)
    warning(msg)
  }

  # l3 models are optional but almost always should be in place
  if (!checkmate::test_class(gpa$l3_models, "hi_model_set")) {
    msg <- "Could not find valid $l3_models specification in gpa. We will proceed, but weird things may happen! I suggest running build_l3_models()."
    lg$warn(msg)
    warning(msg)
  }

  # new approach: use internal model names for creating output directories at subject level
  # default to <analysis_name>/<l1_model>
  # add suffix if using preconvolution approach
  gpa$l1_models$models <- lapply(gpa$l1_models$models, function(mm) {
    mm$outdir <- file.path(gpa$analysis_name, paste0(mm$name, ifelse(gpa$use_preconvolve, "_preconvolve", "")))
    return(mm)
  })

  if (!is.null(gpa$run_number_regex)) {
    if (stringr::str_count(gpa$run_number_regex, stringr::fixed("(")) != 1L) {
      stop(
        "run_number_regex: ", gpa$run_number_regex,
        " must have exactly one opening parenthesis, denoting start of run number capture"
      )
    }

    if (stringr::str_count(gpa$run_number_regex, stringr::fixed(")")) != 1L) {
      stop(
        "run_number_regex: ", gpa$run_number_regex,
        " must have exactly one closing parenthesis, denoting end of run number capture"
      )
    }
  }

  # setup l1 cope names cache from current L1 model contrasts
  gpa <- refresh_l1_cope_names(gpa, lg = lg)


  if (is.null(gpa$center_l3_predictors)) gpa$center_l3_predictors <- TRUE
  if (is.null(gpa$bad_ids)) gpa$bad_ids <- c()
  if (is.null(gpa$scheduler)) gpa$scheduler <- "slurm" # HPC batch system
  gpa$scheduler <- tolower(gpa$scheduler)
  checkmate::assert_subset(gpa$scheduler, c("slurm", "sbatch", "torque", "qsub", "local", "sh"))

  if (is.null(gpa$zthresh)) gpa$zthresh <- 3.09 # 1-tailed p=.001 for z stat
  if (is.null(gpa$clustsize)) gpa$clustsize <- 50 # arbitrary reasonable lower bound on cluster size
  if (is.null(gpa$glm_software)) gpa$glm_software <- "fsl" # default to FSL FEAT

  if (is.null(gpa$log_json)) gpa$log_json <- TRUE # whether to write JSON log files
  if (is.null(gpa$log_txt)) gpa$log_txt <- TRUE # whether to write text log files

  if (is.null(gpa$l1_setup_log)) {
    l1_setup_log <- paste0(names(gpa$l1_models$models), "_l1setup") %>% setNames(names(gpa$l1_models$models))
    lg$debug("l1_setup_log is %s", l1_setup_log)
    gpa$l1_setup_log <- l1_setup_log
  }
  if (is.null(gpa$l1_execution_log)) {
    l1_execution_log <- paste0(names(gpa$l1_models$models), "_l1execution") %>% setNames(names(gpa$l1_models$models))
    lg$debug("l1_execution_log is %s", l1_execution_log)
    gpa$l1_execution_log <- l1_execution_log
  }

  if (is.null(gpa$n_expected_runs)) gpa$n_expected_runs <- 1 # assume single run case

  # remove bad ids before running anything further
  if (!is.null(gpa$bad_ids) && length(gpa$bad_ids) > 0L) {
    lg$info("Removing the following IDs from data structure before beginning analysis: %s", paste(gpa$bad_ids, collapse = ", "))
    gpa$subject_data <- gpa$subject_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$run_data <- gpa$run_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
    gpa$trial_data <- gpa$trial_data %>% dplyr::filter(!id %in% gpa$bad_ids) # remove bad ids
  }

  # build design matrix default arguments
  if (is.null(gpa$additional$bdm_args)) {
    gpa$additional$bdm_args <- list(
      baseline_coef_order = 2, center_values = TRUE,
      plot = FALSE, convolve_wi_run = TRUE, output_directory = "run_timing"
    )
  }

  # default settings for feat l1
  if (is.null(gpa$additional$feat_l1_args)) {
    gpa$additional$feat_l1_args <- list(z_thresh = 1.96, prob_thresh = .05, paradigm_hp = 120)
  }

  # identify and validate niftis for each run
  gpa <- lookup_nifti_inputs(gpa)

  #####
  # handle GLM settings and defaults
  if (is.null(gpa$glm_settings) || gpa$glm_settings[1L] == "default") {
    lg$info("Using default settings for GLM implementation")
    gpa$glm_settings <- list(
      fsl = list(),
      afni = list(),
      spm = list()
    )
  }

  fsl_defaults <- list(
    force_l1_creation = FALSE, # whether to overwrite existing level 1 setup files (e.g., .fsf)
    failed_l1_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l1_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    force_l2_creation = FALSE, # whether to overwrite existing level 2 setup files (e.g., .fsf)
    failed_l2_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l2_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    force_l3_creation = FALSE, # whether to overwrite existing level 3 setup files (e.g., .fsf)
    failed_l3_folder_action = "delete", # whether to 'delete', 'archive', or 'skip' failed folders before running feat l1 fsf jobs
    incomplete_l3_folder_action = "archive", # whether to 'delete', 'archive', or 'skip' incomplete folders before running feat l1 fsf jobs
    auto_retry_l3_excessive_outliers = TRUE, # whether to retry L3 FLAME failures with robust_yn=0 when excessive outliers are detected
    replace_l1_nifti_symlink = TRUE # whether to replace filtered_func_data with symlink to the (same) input data
  )

  gpa$glm_settings$fsl <- populate_defaults(gpa$glm_settings$fsl, fsl_defaults)

  spm_defaults <- list(
    hpf = 100,
    hrf_derivs = "none",
    cvi = "AR(1)",
    estimation_method = "Classical",
    write_residuals = FALSE,
    fmri_t = NULL,
    fmri_t0 = NULL,
    condition_contrasts = TRUE,
    unit_contrasts = TRUE,
    effects_of_interest_F = TRUE,
    l2_projection_model = NULL,
    l2_projection_interactions = NULL,
    l2_projection_interaction_unit_contrasts = TRUE,
    l2_projection_interaction_contrast_modes = NULL,
    l1_session_mode = "separate",
    l3_model_type = "flexible_factorial",
    l3_use_one_sample_when_intercept_only = TRUE,
    spm_execute_setup = FALSE,
    spm_execute_glm = FALSE,
    spm_execute_contrasts = FALSE,
    spm_execute_l3_setup = FALSE,
    spm_execute_l3_glm = FALSE,
    spm_execute_l3_contrasts = FALSE,
    run_l1_setup = TRUE,
    run_l1_glm = TRUE,
    run_l1_contrasts = TRUE,
    run_l3_setup = TRUE,
    run_l3_glm = TRUE,
    run_l3_contrasts = TRUE,
    print_spm_run_instructions = FALSE,
    require_matlab = FALSE,
    matlab_cmd = "matlab",
    matlab_args = "-batch",
    matlab_timeout = 120,
    concatenate_runs = FALSE,
    generate_qsub = TRUE,
    execute_qsub = FALSE,
    cleanup_tmp = TRUE,
    nifti_tmpdir = NULL,
    spm_path = "/proj/mnhallqlab/lab_resources/spm12",
    force_l1_creation = FALSE,
    force_l3_creation = FALSE
  )

  gpa$glm_settings$spm <- populate_defaults(gpa$glm_settings$spm, spm_defaults)

  gpa$glm_settings$afni <- afni_glm_settings(gpa)

  # initialize GLM backend specs and resolved functions
  gpa <- initialize_glm_backends(gpa)

  # normalize and validate L2/L3 model signature fields now that model objects are finalized
  gpa <- normalize_longitudinal_model_signatures(gpa, lg)

  # process confound settings
  gpa <- finalize_confound_settings(gpa, lg)

  # populate subject exclusions
  gpa <- calculate_subject_exclusions(gpa)

  # cache gpa object to file
  res <- tryCatch(saveRDS(gpa, file = gpa$output_locations$object_cache), error = function(e) {
    lg$error("Could not save gpa object to file: %s", gpa$output_locations$object_cache)
    return(NULL)
  })

  # save subject, run, and trial data to the database, too
  lg$info("Writing run_data to sqlite db: %s", gpa$output_locations$sqlite_db)
  DBI::dbWriteTable(conn = sqlite_con, name = "run_data", value = gpa$run_data, overwrite = TRUE)

  lg$info("Writing subject_data to sqlite db: %s", gpa$output_locations$sqlite_db)
  DBI::dbWriteTable(conn = sqlite_con, name = "subject_data", value = gpa$subject_data, overwrite = TRUE)

  lg$info("Writing trial_data to sqlite db: %s", gpa$output_locations$sqlite_db)
  DBI::dbWriteTable(conn = sqlite_con, name = "trial_data", value = gpa$trial_data, overwrite = TRUE)

  lg$debug("Setting finalize_complete to TRUE")
  gpa$finalize_complete <- TRUE

  try(DBI::dbDisconnect(sqlite_con)) # close SQLite connection

  return(gpa)
}

disable_unavailable_spike_volumes <- function(gpa, lg) {
  if (is.null(gpa$confound_settings$spike_volumes)) {
    return(gpa)
  }

  has_motion <- any(isTRUE(gpa$run_data$motion_params_present), na.rm = TRUE)
  if (!isTRUE(has_motion)) {
    lg$warn(
      "spike_volumes was requested, but no motion parameter files were found for any run; spike regressors will not be generated"
    )
    gpa$confound_settings$spike_volumes <- NULL
  }

  gpa
}

setup_parallel_settings <- function(gpa, lg = NULL) {
  checkmate::assert_class(lg, "Logger")

  # ---- PARALLELISM SETUP
  parse_mem_to_gb <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    x <- trimws(tolower(as.character(x[1L])))
    if (!nzchar(x)) return(NA_real_)
    rr <- regexec("^([0-9]+(?:\\.[0-9]+)?)\\s*([a-z]*)$", x)
    mm <- regmatches(x, rr)[[1]]
    if (length(mm) != 3L) return(NA_real_)
    val <- suppressWarnings(as.numeric(mm[2L]))
    if (is.na(val)) return(NA_real_)
    unit <- mm[3L]
    if (unit %in% c("", "g", "gb")) {
      mult <- 1
    } else if (unit %in% c("m", "mb")) {
      mult <- 1 / 1024
    } else if (unit %in% c("k", "kb")) {
      mult <- 1 / (1024 * 1024)
    } else if (unit %in% c("t", "tb")) {
      mult <- 1024
    } else {
      mult <- NA_real_
    }
    if (is.na(mult)) return(NA_real_)
    val * mult
  }

  enforce_min_mem <- function(value, min_gb, style = c("plain_gb", "suffix_g"), label = "memory") {
    style <- match.arg(style)
    gb <- parse_mem_to_gb(value)
    if (is.na(gb)) {
      lg$warn("Could not parse %s value '%s'. Using %d GB.", label, as.character(value), as.integer(min_gb))
      gb <- min_gb
    }
    if (gb < min_gb) {
      lg$info("Increasing %s from %.2f GB to %d GB minimum.", label, gb, as.integer(min_gb))
      gb <- min_gb
    }
    rounded <- as.integer(ceiling(gb))
    if (style == "suffix_g") paste0(rounded, "G") else as.character(rounded)
  }

  specify_cores <- function(gpa, field_name, default=1L) {
    # l1_setup_cores defines how many cores to use when looping over subjects within setup_l1_models
    if (is.null(gpa$parallel[[field_name]]) || gpa$parallel[[field_name]] == "default") {
      gpa$parallel[[field_name]] <- default
    }

    checkmate::assert_integerish(gpa$parallel[[field_name]], lower = 1)
    return(gpa)
  }

  # pipeline_cores: number of cores used in push_pipeline when looping over l1 model variants
  gpa <- specify_cores(gpa, "pipeline_cores", ifelse(is.null(gpa$l1_models$models), 4, length(gpa$l1_models$models)))

  # l1_setup_cores defines how many cores to use when looping over models in setup_l1_models
  # default to 4 cores in setup_lvl1_models
  gpa <- specify_cores(gpa, "l1_setup_cores", 4)

  # default to 4GB per core for l1 setup
  if (is.null(gpa$parallel$l1_setup_memgb)) gpa$parallel$l1_setup_memgb <- paste0(4 * gpa$parallel$l1_setup_cores, "G")

  # l2_setup_cores defines how many cores to use when looping over models in setup_l2_models
  gpa <- specify_cores(gpa, "l2_setup_cores", 4)

  # finalize_cores defines how many cores to use when finalizing the pipeline before setting up and running models
  gpa <- specify_cores(gpa, "finalize_cores", 4)

  if (is.null(gpa$parallel$slurm)) gpa$parallel$slurm <- list()
  if (is.null(gpa$parallel$torque)) gpa$parallel$torque <- list()

  if (gpa$scheduler == "slurm") {
    if (is.null(gpa$parallel$sched_args)) {
      # Jan 2022: turns out Longleaf doesn't want us to use a partition by default!
      # gpa$parallel$sched_args <- c("-p general")
      # lg$info("Using default SLURM scheduler arguments: ")
      # lg$info("Argument: %s", gpa$parallel$sched_args)
    } else {
      checkmate::assert_character(gpa$parallel$sched_args)
    }
  } else if (gpa$scheduler == "torque") {
    # gpa$parallel$sched_args <- c("-A mnh5174_c_g_sc_default", "-W group_list=mnh5174_collab")
    if (is.null(gpa$parallel$sched_args)) {
      gpa$parallel$sched_args <- c("-j oe", "-m n")
      lg$info("Using default PBS scheduler arguments: ")
      lg$info("Argument: %s", gpa$parallel$sched_args)
    }
  }

  # time for finalize_pipeline_configuration in run_glm_pipeline
  if (is.null(gpa$parallel$finalize_time)) gpa$parallel$finalize_time <- "6:00:00" # 6.0 hours (includes run truncation, which is slow)
  if (is.null(gpa$parallel$l1_setup_time)) gpa$parallel$l1_setup_time <- "4:00:00" # 4.0 hours

  # 14 hours for l2 setup and execution to clear scheduler (all jobs)
  if (is.null(gpa$parallel$l2_setup_run_time)) gpa$parallel$l2_setup_run_time <- "14:00:00"

  # 80 hours for all L3 analyses to clear scheduler (all jobs)
  if (is.null(gpa$parallel$l3_setup_run_time)) gpa$parallel$l3_setup_run_time <- "80:00:00"
  if (is.null(gpa$parallel$compute_environment$global) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default global compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$global <- "module use /proj/mnhallqlab/sw/modules"
  }

  if (is.null(gpa$parallel$compute_environment$r) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default R compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$r <- c(
      "module unload r",
      "module load r/4.5.1"
    )
  }

  # fsl_parallel_defaults <- list(
  #   l1_feat_time = "8:00:00", # 8 hours
  # )

  # number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(gpa$parallel$fsl$l1_feat_alljobs_time)) gpa$parallel$fsl$l1_feat_alljobs_time <- "72:00:00" # 3 days for all jobs
  if (is.null(gpa$parallel$fsl$l1_feat_time)) gpa$parallel$fsl$l1_feat_time <- "10:00:00" # 10 hours
  if (is.null(gpa$parallel$fsl$l1_feat_memgb)) gpa$parallel$fsl$l1_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l1_feat_cpus_per_job)) gpa$parallel$fsl$l1_feat_cpus_per_job <- 8
  if (is.null(gpa$parallel$fsl$l1_feat_runs_per_cpu)) gpa$parallel$fsl$l1_feat_runs_per_cpu <- 2 # 2 runs per cpu in a job
  if (is.null(gpa$parallel$fsl$l2_feat_time)) gpa$parallel$fsl$l2_feat_time <- "2:00:00" # 2 hours
  if (is.null(gpa$parallel$fsl$l2_feat_memgb)) gpa$parallel$fsl$l2_feat_memgb <- "12" # 12 GB by default
  if (is.null(gpa$parallel$fsl$l3_feat_time)) gpa$parallel$fsl$l3_feat_time <- "24:00:00" # 24 hours
  if (is.null(gpa$parallel$fsl$l3_feat_memgb)) gpa$parallel$fsl$l3_feat_memgb <- "32" # 32 GB by default
  if (is.null(gpa$parallel$fsl$l3_feat_cpusperjob)) gpa$parallel$fsl$l3_feat_cpusperjob <- 16 # cpus used to process all slices

  if (is.null(gpa$parallel$spm)) gpa$parallel$spm <- list()
  if (is.null(gpa$parallel$spm$l1_spm_alljobs_time)) gpa$parallel$spm$l1_spm_alljobs_time <- "24:00:00"
  if (is.null(gpa$parallel$spm$l1_spm_time)) gpa$parallel$spm$l1_spm_time <- "8:00:00"
  if (is.null(gpa$parallel$spm$l1_spm_memgb)) gpa$parallel$spm$l1_spm_memgb <- "12"
  if (is.null(gpa$parallel$spm$l1_spm_cpus_per_job)) gpa$parallel$spm$l1_spm_cpus_per_job <- 4
  if (is.null(gpa$parallel$spm$l1_spm_runs_per_cpu)) gpa$parallel$spm$l1_spm_runs_per_cpu <- 1
  if (is.null(gpa$parallel$spm$l3_spm_alljobs_time)) gpa$parallel$spm$l3_spm_alljobs_time <- "24:00:00"
  if (is.null(gpa$parallel$spm$l3_spm_time)) gpa$parallel$spm$l3_spm_time <- "8:00:00"
  if (is.null(gpa$parallel$spm$l3_spm_memgb)) gpa$parallel$spm$l3_spm_memgb <- "16"

  if (is.null(gpa$parallel$afni)) gpa$parallel$afni <- list()
  gpa$parallel$afni <- populate_defaults(gpa$parallel$afni, list(
    l1_glm_alljobs_time = "24:00:00", l1_glm_time = "8:00:00", l1_glm_memgb = "12", l1_glm_cpus = 4L,
    l3_glm_time = "8:00:00", l3_glm_memgb = "16", l3_glm_cpus = 1L
  ))
  if (is.null(gpa$parallel$afni$l3_lmer_time)) gpa$parallel$afni$l3_lmer_time <- "48:00:00"
  if (is.null(gpa$parallel$afni$l3_lmer_memgb)) gpa$parallel$afni$l3_lmer_memgb <- "32"
  if (is.null(gpa$parallel$afni$l3_lmer_njobs)) gpa$parallel$afni$l3_lmer_njobs <- 8

  gpa$parallel$l1_setup_memgb <- enforce_min_mem(
    gpa$parallel$l1_setup_memgb,
    min_gb = 16,
    style = "suffix_g",
    label = "parallel$l1_setup_memgb"
  )

  gpa$parallel$fsl$l1_feat_memgb <- enforce_min_mem(
    gpa$parallel$fsl$l1_feat_memgb,
    min_gb = 16,
    style = "plain_gb",
    label = "parallel$fsl$l1_feat_memgb"
  )
  gpa$parallel$fsl$l2_feat_memgb <- enforce_min_mem(
    gpa$parallel$fsl$l2_feat_memgb,
    min_gb = 16,
    style = "plain_gb",
    label = "parallel$fsl$l2_feat_memgb"
  )
  gpa$parallel$fsl$l3_feat_memgb <- enforce_min_mem(
    gpa$parallel$fsl$l3_feat_memgb,
    min_gb = 16,
    style = "plain_gb",
    label = "parallel$fsl$l3_feat_memgb"
  )

  gpa$parallel$spm$l1_spm_memgb <- enforce_min_mem(
    gpa$parallel$spm$l1_spm_memgb,
    min_gb = 32,
    style = "plain_gb",
    label = "parallel$spm$l1_spm_memgb"
  )
  gpa$parallel$spm$l3_spm_memgb <- enforce_min_mem(
    gpa$parallel$spm$l3_spm_memgb,
    min_gb = 32,
    style = "plain_gb",
    label = "parallel$spm$l3_spm_memgb"
  )
  if (is.null(gpa$parallel$spm$l3_spm_cpus_per_job)) gpa$parallel$spm$l3_spm_cpus_per_job <- 4
  if (is.null(gpa$parallel$spm$l3_spm_runs_per_cpu)) gpa$parallel$spm$l3_spm_runs_per_cpu <- 1

  if (is.null(gpa$parallel$compute_environment$fsl) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default FSL compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$fsl <- c(
      "module unload fsl", # remove any current fsl module
      "module load fsl/6.0.7.18" # load latest version (2021)
    )
  }

  if (is.null(gpa$parallel$compute_environment$afni) && isTRUE(grepl("(longleaf|ll\\.unc\\.edu)", gpa$nodename))) {
    lg$info("Using default AFNI compute environment for UNC Longleaf")
    gpa$parallel$compute_environment$afni <- c(
      "module unload afni", # remove any current afni module
      "module load afni/26.0.00" # load latest version (2026)
    )
  }

  if (is.null(gpa$parallel$fsl$slurm_l1_array)) {

  }

  return(gpa)
}

# helper function for settings up $output_directory and $output_locations
setup_output_locations <- function(gpa, lg = NULL) {
  checkmate::assert_class(lg, "Logger")

  # sort out file locations
  if (is.null(gpa$output_directory) || gpa$output_directory == "default") {
    gpa$output_directory <- file.path(getwd(), gpa$analysis_name)
  }

  lg$info("Output directory for this analysis will be: %s", gpa$output_directory)

  # see quickstart.Rmd > output settings

  # build out ability to consolidate outputs in one folder, to use specific paths for some outputs, etc.
  # if user specifies gpa$output_directory that matches gpa$analysis_name, don't at this as subfolder

  if (length(unique(gpa$run_data$session)) == 1L) {
    feat_sub_directory <- file.path("{gpa$output_directory}", "feat_l1", "sub-{id}")
    feat_l2_sub_directory <- file.path("{gpa$output_directory}", "feat_l2", "sub-{id}")
    spm_sub_directory <- file.path("{gpa$output_directory}", "spm_l1", "sub-{id}")
    afni_sub_directory <- file.path("{gpa$output_directory}", "afni_l1", "sub-{id}")
  } else {
    feat_sub_directory <- file.path("{gpa$output_directory}", "feat_l1", "sub-{id}", "ses-{session}")
    feat_l2_sub_directory <- file.path("{gpa$output_directory}", "feat_l2", "sub-{id}", "ses-{session}")
    spm_sub_directory <- file.path("{gpa$output_directory}", "spm_l1", "sub-{id}", "ses-{session}")
    afni_sub_directory <- file.path("{gpa$output_directory}", "afni_l1", "sub-{id}", "ses-{session}")
  }

  output_defaults <- list(
    # default to BIDS-style consolidated output
    consolidated = TRUE,
    feat_sub_directory = feat_sub_directory,
    feat_ses_directory = feat_sub_directory, # no difference in defaults
    feat_l1_directory = file.path(feat_sub_directory, "{l1_model}"),
    # organize L2 outputs by the lower-level cope because FEAT receives one L1 cope at a time
    feat_l2_directory = file.path(
      feat_l2_sub_directory,
      "L1m-{l1_model}", "l1c-{l1_cope_label}", "L2m-{l2_model}"
    ),
    # for pooled L2 models (l2_scope='id'), avoid synthetic ses-0 folders
    feat_l2_id_scope_directory = file.path(
      "{gpa$output_directory}", "feat_l2", "sub-{id}",
      "L1m-{l1_model}", "l1c-{l1_cope_label}", "L2m-{l2_model}"
    ),
    spm_sub_directory = spm_sub_directory,
    spm_l1_directory = file.path(spm_sub_directory, "{l1_model}"),
    # for pooled SPM L1 models (l1_session_mode='pooled'), avoid synthetic ses-0 folders
    spm_l1_id_scope_directory = file.path("{gpa$output_directory}", "spm_l1", "sub-{id}", "{l1_model}"),
    # old default:
    # spm_l3_directory = file.path("{gpa$output_directory}", "spm_l3", "L1m-{l1_model}", "l1c-{l1_cope_name}", "L3m-{l3_model}"),
    spm_l3_directory = file.path("{gpa$output_directory}", "spm_l3", "L3m-{l3_model}", "L1m-{l1_model}", "l1c-{l1_cope_name}"),
    afni_sub_directory = afni_sub_directory,
    afni_l1_directory = file.path(afni_sub_directory, "{l1_model}"),
    afni_l3_directory = file.path("{gpa$output_directory}", "afni_l3", "L3m-{l3_model}",
                                  "L1m-{l1_model}", "l1c-{l1_cope_name}"),
    # old default:
    # feat_l3_directory = ifelse(isTRUE(gpa$multi_run),
    #   file.path("{gpa$output_directory}", "feat_l3", "L1m-{l1_model}", "L2m-{l2_model}_l2c-{l2_contrast}", "L3m-{l3_model}"),
    #   file.path("{gpa$output_directory}", "feat_l3", "L1m_{l1_model}", "L3m-{l3_model}")
    # )
    feat_l3_directory = ifelse(isTRUE(gpa$multi_run),
      file.path(
        "{gpa$output_directory}", "feat_l3", "L3m-{l3_model}",
        "L1m-{l1_model}", "L2m-{l2_model}", "l2c-{l2_contrast}"
      ),
      file.path(
        "{gpa$output_directory}", "feat_l3", "L3m-{l3_model}",
        "L1m-{l1_model}"
      )
    ),
    feat_l3_fsf = "FEAT_l1c-{l1_contrast}.fsf",
    # old default:
    # feat_l3_combined_filename = ifelse(isTRUE(gpa$multi_run),
    #   file.path("{gpa$output_directory}", "feat_l3_combined", "L1m-{l1_model}", "l1c-{l1_cope_name}", "L2m-{l2_model}_L3m-{l3_model}_stats"),
    #   file.path("{gpa$output_directory}", "feat_l3_combined", "L1m-{l1_model}", "l1c-{l1_cope_name}", "L3m-{l3_model}_stats")
    # )
    feat_l3_combined_filename = ifelse(isTRUE(gpa$multi_run), # settings for combining FEAT L3 models into a smaller set of AFNI files for visualization
      file.path(
        "{gpa$output_directory}", "feat_l3_combined", "L3m-{l3_model}",
        "L1m-{l1_model}", "L2m-{l2_model}", "l1c-{l1_cope_name}",
        "stats"
      ),
      file.path(
        "{gpa$output_directory}", "feat_l3_combined", "L3m-{l3_model}",
        "L1m-{l1_model}", "l1c-{l1_cope_name}", "stats"
      )
    ),
    feat_l3_combined_briknames = ifelse(isTRUE(gpa$multi_run),
      "l2c-{l2_cope_name}_l3c-{l3_cope_name}",
      "l3c-{l3_cope_name}"
    ),
    # old default:
    # spm_l3_combined_filename = file.path("{gpa$output_directory}", "spm_l3_combined", "L1m-{l1_model}", "l1c-{l1_cope_name}", "L3m-{l3_model}_stats"),
    spm_l3_combined_filename = file.path(
      "{gpa$output_directory}", "spm_l3_combined", "L3m-{l3_model}",
      "L1m-{l1_model}", "l1c-{l1_cope_name}", "stats"
    ),
    spm_l3_combined_briknames = "l3c-{l3_cope_name}",
    # old default:
    # afni_3dlmer_directory = file.path(
    #   "{gpa$output_directory}", "afni_3dlmer",
    #   "L1m-{l1_model}", "l1c-{l1_cope_name}",
    #   "L2m-{l2_model}_l2c-{l2_cope_name}",
    #   "L3m-{l3_model}"
    # ),
    afni_3dlmer_directory = file.path(
      "{gpa$output_directory}", "afni_3dlmer", "L3m-{l3_model}",
      "L1m-{l1_model}", "L2m-{l2_model}",
      "l1c-{l1_cope_name}", "l2c-{l2_cope_name}"
    ),
    scheduler_scripts = file.path(gpa$output_directory, "scheduler_scripts"),
    sqlite_db = file.path(gpa$output_directory, paste0(gpa$analysis_name, ".sqlite")),
    project_config_json = file.path(gpa$output_directory, "project_config.json"),
    object_cache = file.path(gpa$output_directory, paste0(gpa$analysis_name, ".rds")),
    log_directory = file.path(gpa$output_directory, "logs"),
    setup_l1_log_txt = file.path(gpa$output_directory, "logs", "setup_l1_models.txt"),
    setup_l1_log_json = file.path(gpa$output_directory, "logs", "setup_l1_models.json"),
    setup_l2_log_txt = file.path(gpa$output_directory, "logs", "setup_l2_models.txt"),
    setup_l2_log_json = file.path(gpa$output_directory, "logs", "setup_l2_models.json"),
    setup_l2_l1_cope_validity_tsv = file.path(gpa$output_directory, "logs", "setup_l2_l1_cope_validity.tsv"),
    setup_l3_log_txt = file.path(gpa$output_directory, "logs", "setup_l3_models.txt"),
    setup_l3_log_json = file.path(gpa$output_directory, "logs", "setup_l3_models.json")
  )

  if (checkmate::test_string(gpa$output_locations) && gpa$output_locations[1L] == "default") {
    gpa$output_locations <- output_defaults
  }

  miss_fields <- setdiff(names(output_defaults), names(gpa$output_locations))
  if (length(miss_fields) > 0L) {
    for (mm in miss_fields) {
      lg$info("Populating missing $output_locations field: %s with default: %s", mm, output_defaults[[mm]])
      gpa$output_locations[[mm]] <- output_defaults[[mm]]
    }
  }

  # ensure log directory exists if configured
  log_dir <- gpa$output_locations$log_directory
  if (checkmate::test_string(log_dir) && !dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
  }

  return(gpa)
}

#' Helper function to populate confound information for pipeline files
#'
#' @param gpa a glm_pipeline_arguments object for population
#' @param lg a Logger object for logging results of confound processing
#' @keywords internal
#' @importFrom parallel mclapply
finalize_confound_settings <- function(gpa, lg) {
  checkmate::assert_class(lg, "Logger")
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  # validate confound settings
  default_exclude_run <- "mean(framewise_displacement) > 0.5 | max(framewise_displacement) > 6"
  confound_defaults <- list(
    motion_params_file = NULL,
    motion_params_colnames = NULL,
    confound_input_file = NULL,
    l1_confound_regressors = NULL, # column names in motion_params_file and/or confound_input_file ("all" to include all columns)
    exclude_run = NULL,
    truncate_run = NULL, # example: framewise_displacement > 1 & time > last_onset
    spike_volumes = NULL,
    na_strings = getOption("datatable.na.strings", "NA") # default na.strings argument for data.table::fread calls
  )

  user_specified_exclude_run <- !is.null(gpa$confound_settings) && "exclude_run" %in% names(gpa$confound_settings)
  if (is.null(gpa$confound_settings)) {
    lg$info("Using default settings for confounds and exclusions")
    lg$info("No motion parameters, confounds, spike regressors, or run exclusions are assumed by default")
  } else {
    checkmate::assert_string(gpa$confound_settings$exclude_run, null.ok = TRUE)
    checkmate::assert_string(gpa$confound_settings$exclude_subject, null.ok = TRUE)
  }

  gpa$confound_settings <- populate_defaults(gpa$confound_settings, confound_defaults)

  #####
  # populate confounds
  if ("motion_params_file" %in% names(gpa$run_data)) {
    lg$info("motion_params_file column already in run_data. Not using motion_params_file confounds specification.")
  } else if (!is.null(gpa$confound_settings$motion_params_file)) {
    checkmate::assert_string(gpa$confound_settings$motion_params_file)
    gpa$run_data$motion_params_file <- gpa$confound_settings$motion_params_file # this gets expanded by get_mr_abspath in get_l1_confounds
  }

  if (!"motion_params_file" %in% names(gpa$run_data)) {
    gpa$run_data$motion_params_file <- gpa$run_data$motion_params_present <- NA_character_
  } else {
    gpa$run_data$motion_params_present <- file.exists(get_mr_abspath(gpa$run_data, "motion_params_file"))
  }

  if ("confound_input_file" %in% names(gpa$run_data)) {
    lg$info("confound_input_file column already in run_data. Not using confound_input_file confounds specification.")
  } else if (!is.null(gpa$confound_settings$confound_input_file)) {
    checkmate::assert_string(gpa$confound_settings$confound_input_file)
    gpa$run_data$confound_input_file <- gpa$confound_settings$confound_input_file # this gets expanded by get_mr_abspath in get_l1_confounds
  }

  # determine whether confound input files are present
  if (!"confound_input_file" %in% names(gpa$run_data)) {
    gpa$run_data$confound_input_file <- gpa$run_data$confound_input_file_present <- NA_character_
  } else {
    gpa$run_data$confound_input_file_present <- file.exists(get_mr_abspath(gpa$run_data, "confound_input_file"))
  }

  gpa <- disable_unavailable_spike_volumes(gpa, lg)

  if (isTRUE(user_specified_exclude_run) &&
      is.character(gpa$confound_settings$exclude_run) &&
      length(gpa$confound_settings$exclude_run) == 1L &&
      tolower(gpa$confound_settings$exclude_run) == "default") {
    has_fd <- FALSE

    if ("motion_params_file" %in% names(gpa$run_data)) {
      has_fd <- any(isTRUE(gpa$run_data$motion_params_present), na.rm = TRUE)
    }

    if (!has_fd && !is.null(gpa$confound_settings$confound_input_colnames)) {
      has_fd <- "framewise_displacement" %in% gpa$confound_settings$confound_input_colnames
    }

    if (!has_fd && !is.null(gpa$confound_settings$motion_params_colnames)) {
      has_fd <- "framewise_displacement" %in% gpa$confound_settings$motion_params_colnames
    }

    if (!has_fd && "confound_input_file" %in% names(gpa$run_data)) {
      confound_files <- get_mr_abspath(gpa$run_data, "confound_input_file")
      confound_files <- confound_files[!is.na(confound_files) & nzchar(confound_files)]
      confound_files <- confound_files[file.exists(confound_files)]
      if (length(confound_files) > 0L) {
        confound_header <- tryCatch(
          data.table::fread(confound_files[1L], nrows = 0, data.table = FALSE),
          error = function(e) NULL
        )
        if (!is.null(confound_header)) {
          has_fd <- "framewise_displacement" %in% names(confound_header)
        }
      }
    }

    if (isTRUE(has_fd)) {
      gpa$confound_settings$exclude_run <- default_exclude_run
      lg$info("Using default exclude_run: %s", default_exclude_run)
    } else {
      gpa$confound_settings$exclude_run <- NULL
      lg$warn("exclude_run was set to 'default' but framewise_displacement was not found; no runs will be excluded by default")
    }
  } else if (is.null(gpa$confound_settings$exclude_run)) {
    lg$info("exclude_run is unset; no runs will be excluded by default")
  }

  rhs_to_vars <- function(str) {
    if (is.null(str)) {
      NULL # return NULL if input is NULL
    } else if (checkmate::test_string(str)) {
      if (!grepl("^\\s*~", str)) str <- paste("~", str)
      all.vars(as.formula(paste("~", str)))
    } else if (checkmate::test_formula(str)) {
      all.vars(str)
    } else {
      stop("rhs_to_vars input is not a string or formula")
    }
  }

  # figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_exclusion_columns <- rhs_to_vars(gpa$confound_settings$exclude_run)

  # figure out all confound columns that will be used in the pipeline
  gpa$confound_settings$run_truncation_columns <- rhs_to_vars(gpa$confound_settings$truncate_run)

  # TODO: Should this become 'id_exclusion_columns' and should we support session versus subject exclusion
  # (E.g., in longitudinal analysis)
  gpa$confound_settings$subject_exclusion_columns <- rhs_to_vars(gpa$confound_settings$exclude_subject)

  gpa$confound_settings$all_confound_columns <- unique(c(
    gpa$confound_settings$l1_confound_regressors,
    gpa$confound_settings$run_exclusion_columns,
    gpa$confound_settings$run_truncation_columns,
    gpa$confound_settings$subject_exclusion_columns
  ))

  # handle lookup and creation of all confound files

  # numeric row number of each input to aid in tracking
  gpa$run_data$input_number <- seq_len(nrow(gpa$run_data))

  # populate confounds in SQLite database and calculate run exclusions
  # TODO: allow external $exclude_run from user, add internal calculated $calc_exclude run

  # confound_info <- parallel::mclapply(seq_len(nrow(gpa$run_data)), function(ii) {
  #   # this should add rows to the SQLite data for a subject if not yet present, or just return those rows if they exist
  #   l1_info <- get_l1_confounds(
  #     id = gpa$run_data$id[ii], session = gpa$run_data$session[ii], run_number = gpa$run_data$run_number[ii],
  #     gpa = gpa, drop_volumes = gpa$drop_volumes
  #   )[c("l1_confound_file", "exclude_run")]
  #   return(l1_info)
  # }, mc.cores = gpa$parallel$finalize_cores)

  # add onset + offset + isi data to run_data so that get_l1_confounds can handle run truncation appropriately
  gpa <- populate_last_events(gpa, lg)

  # for each run, calculate confounds, exclusions, and truncation
  run_list <- lapply(seq_len(nrow(gpa$run_data)), function(ii) {
    get_l1_confounds(run_df = gpa$run_data[ii, , drop = FALSE], gpa = gpa)
  })

  if (length(unique(sapply(run_list, length))) == 1L) {
    gpa$run_data <- data.table::rbindlist(run_list) |> data.table::setDF()
  } else {
    msg <- "Lengths of confound run_df elements have different lengths. Cannot recombine"
    lg$error(msg)
    stop(msg)
  }

  return(gpa)
}
