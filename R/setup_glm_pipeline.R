# register bg_ custom classes so that DBI writeTable works with this as a formal (S4?) class
setOldClass(c("bg_subject_data", "data.frame"))
setOldClass(c("bg_run_data", "data.frame"))
setOldClass(c("bg_session_data", "data.frame"))
setOldClass(c("bg_block_data", "data.frame"))
setOldClass(c("bg_trial_data", "data.frame"))

#' Main worker function for setting up an analysis pipeline
#'
#' @param analysis_name A character string providing a useful name for identifying this analysis. Practically, this
#'   influences the top-level folder name of the group analysis outputs, as well as the name of .RData objects
#'   saved by this function to the output directory for the analysis.
#' @param scheduler Which HPC scheduler system should be used for queueing jobs. Options are 'slurm', 'torque', or 'local'.
#' @param subject_data A data.frame containing Level 3 covariates and subject/session metadata. For cross-sectional
#'   analyses this is typically one row per subject. For longitudinal analyses, this must contain one row per
#'   \code{id}/\code{session}; person-level covariates such as sex, baseline age, or treatment arm should be repeated
#'   across that subject's sessions. Columns from \code{subject_data} can be used as covariates in group
#'   (aka 'level 3') analyses. If \code{NULL}, then this will be distilled from \code{trial_data} by looking for
#'   variables that vary at the same level as \code{vm["id"]} and \code{vm["session"]}. Do not provide
#'   \code{session_data}; \code{setup_glm_pipeline()} derives \code{gpa$session_data} internally from
#'   \code{run_data} and \code{subject_data}.
#' @param run_data A data.frame containing all run-level data such as run condition or run number. Columns from
#'   \code{run_data} can be used as covariates in subject (aka 'level 2') analyses. If \code{NULL}, this will be
#'   distilled from \code{trial_data} by looking for variables that vary at the same level as \code{vm["run_number"]}.
#'   Columns that are constant within each \code{id}/\code{session} are copied into the internally derived
#'   \code{gpa$session_data}.
#' @param trial_data A data.frame containing trial-level statistics for all subjects. Data should be stacked in long
#'   form such that each row represents a single trial for a given subject and the number of total rows is subjects x trials.
#'   If you wish, you can pass a single trial-level data frame that also contains all run-level and subject-level covariates
#'   (i.e., a combined long format, where variables at different levels are all included as columns). In this case,
#'   \code{setup_glm_pipeline} will detect which variables occur at each level and parse these accordingly into
#'   \code{subject_data} and \code{run_data}; \code{session_data} is still constructed internally.
#' @param ppi_data An optional `data.frame` containing physiological signals (e.g., deconvolved time series in seed ROIs) that
#'   can be used in setting up the level 1 models for a PPI analysis. Specifically, if this `data.frame` is provided, you will be
#'   able to create interaction regressors for any convolved event-related signal and the columns of `ppi_data`. This data.frame
#'   must include subject id, session id (if relevant), and run_number as columns. All additional columns should be eligible
#'   physiological signals for PPI analysis.
#' @param vm A named character vector containing key identifying columns in \code{subject_data} and
#'   \code{trial_data}. Minimally, this vector should contain the elements 'id'
#' @param bad_ids An optional vector of ids in \code{subject_data} and \code{trial_data} that should be excluded from analysis.
#' @param fmri_file_regex A character string containing a Perl-compatible regular expression for the subfolder and filename
#'   within the \code{mr_dir} field in \code{subject_data}.
#' @param group_output_directory Group-level output directory policy.
#' @param output_locations Output location template settings.
#' @param fmri_path_regex Optional regular expression for extracting fMRI paths.
#' @param run_number_regex Optional regular expression for extracting run numbers.
#' @param tr The repetition time of the scanning sequence in seconds. Used for setting up design matrices. If this is NULL, the
#'   function will look for a \code{tr} field in the \code{run_data} object, which specifies TR at run level (e.g., if varies).
#' @param output_directory The output directory for all configuration and job submission files for this analysis.
#'   Default is a subfolder called "glm_out" in the current working directory.
#' @param drop_volumes The number of volumes to drop from the fMRI data and convolved regressors prior to analysis.
#'   Default is 0.
#' @param use_preconvolve For FSL, pass pipeline-convolved EVs (\code{TRUE}, the default)
#'   or let FEAT convolve onset/duration/value files (\code{FALSE}). AFNI uses pipeline-convolved
#'   task regressors. SPM builds conditions and parametric modulators from event timing.
#' @param l1_models An \code{l1_model_set} object containing all level 1 (run) models to be included in GLM pipeline.
#'   If "prompt" is passed, \code{setup_glm_pipeline} will call \code{build_l1_models} to setup l1 models
#'   interactively. Optionally, this argument can be \code{NULL} if you want to setup l1 models later, though
#'   the resulting object will not be functional within the pipeline until l1 models are provided.
#' @param l2_models An \code{hi_model_set} object containing all level 2 (subject) models to be included in GLM pipeline.
#'   If "prompt" is passed, \code{setup_glm_pipeline} will call \code{build_l2_models} to setup l2 models
#'   interactively. Use \code{NULL} to defer setup or omit L2 for a native AFNI or SPM fit.
#'   SPM can optionally use an L2 specification to project run-level predictors into L1.
#' @param l3_models An \code{hi_model_set} object containing all level 3 (sample) models to be included in GLM pipeline.
#'   If "prompt" is passed, \code{setup_glm_pipeline} will call \code{build_l3_models} to setup l3 models
#'   interactively. Optionally, this argument can be \code{NULL} if you want to setup l3 models later, though
#'   the resulting object will not be functional within the pipeline until l3 models are provided.
#' @param glm_software Default GLM engine(s), case insensitive. Use \code{"fsl"} (default)
#'   for run-wise FEAT and an explicit multi-run L2, \code{"afni"} for shared task effects
#'   across runs within a visit, or \code{"spm"} for an SPM first-level model and group analysis.
#'   A vector enables multiple engines. Set \code{gpa$level_backends} to assign them to levels.
#' @param n_expected_runs Number of expected runs per subject. Used to determine 2- versus 3-level analysis
#'   (for FSL), and for providing feedback about subjects who have unexpected numbers of runs.
#' @param additional A list of additional metadata that will be added to the \code{glm.pipeline} object returned
#'   by the function. This can be useful if there are other identifiers that you want for long-term storage or
#'   off-shoot functions.
#' @param glm_settings A named list of engine settings, e.g.,
#'   \code{list(afni = list(estimation = "OLS"))} or \code{list(spm = list(spm_path = "/path/to/spm12"))}.
#'   The default AFNI estimator is REML. Configure SPM with your local installation path.
#' @param confound_settings Confound-file and censoring settings.
#' @param parallel Parallelization and scheduler settings.
#' @param lgr_threshold The logging threshold used to determine whether to output messages of different severity to
#'   the screen and to log files. Default is "info", which produces all messages, warnings, and errors, but not debug
#'   or trace statements. To output only concerning errors, change to "error". See: 
#'   \url{https://s-fleck.github.io/lgr/articles/lgr.html} for details.
#'
#' @seealso \code{\link{run_glm_pipeline}}, \code{vignette("glm-engines", package = "fmri.pipeline")}
#' @importFrom checkmate assert_subset assert_data_frame assert_number assert_integerish 
#' @importFrom checkmate assert_list assert_logical test_string test_class
#' @importFrom dplyr mutate_at group_by select vars inner_join filter count is_grouped_df ungroup n_distinct
#' @importFrom stats aggregate
#' @importFrom tidyselect everything
#' @importFrom utils capture.output menu
#' @export
setup_glm_pipeline <- function(analysis_name = "glm_analysis", scheduler = "slurm",
                               output_directory = file.path(getwd(), analysis_name),
                               subject_data = NULL, run_data = NULL, trial_data = NULL,
                               ppi_data = NULL,
                               group_output_directory = "default",
                               output_locations = "default",
                               vm = c(
                                 id = "id", session = "session", run_number = "run_number", block_number = "block_number", trial = "trial",
                                 run_trial = "run_trial", mr_dir = "mr_dir", run_nifti = "run_nifti",
                                 exclude_subject = "exclude_subject"
                               ),
                               bad_ids = NULL, tr = NULL,
                               fmri_file_regex = ".*\\.nii(\\.gz)?", fmri_path_regex = NULL,
                               run_number_regex = ".*run-*([0-9]+).*", drop_volumes = 0L,
                               l1_models = "prompt", l2_models = "prompt", l3_models = "prompt",
                               glm_software = "fsl", n_expected_runs = 1L,
                               use_preconvolve = TRUE,
                               glm_settings = "default",
                               confound_settings = list(
                                 motion_params_file = NULL,
                                 motion_params_colnames = NULL,
                                 confound_input_file = NULL,
                                 confound_input_colnames = NULL, # names of confound columns -- if null, we will attempt to find a header row
                                 l1_confound_regressors = NULL, # column names in motion_params_file and/or confound_input_file ("all" to include all columns)
                                 exclude_run = NULL,
                                 truncate_run = NULL, # "framewise_displacement > 1 & volume > last_onset"
                                 exclude_subject = NULL,
                                 spike_volumes = NULL
                               ),
                               parallel = list(
                                 # number of cores used when looping over l1 setup of design matrices and syntax for each subject
                                 l1_setup_cores = 4L,

                                 # number of cores used  when looping over l1 model variants in push_pipeline
                                 pipeline_cores = "default"
                               ), additional = list(
                                 # additional feat level 1 settings (uses internal FSL nomenclature)
                                 feat_l1_args = list(z_thresh = 1.96, prob_thresh = .05)
                               ),
                               lgr_threshold = "info") {

  checkmate::assert_string(analysis_name) # must be scalar string
  checkmate::assert_subset(scheduler, c("slurm", "sbatch", "torque", "qsub", "local", "sh"), empty.ok = FALSE)
  checkmate::assert_data_frame(subject_data, null.ok = TRUE)
  checkmate::assert_data_frame(run_data, null.ok = TRUE)
  checkmate::assert_data_frame(trial_data)
  checkmate::assert_character(vm, unique = TRUE) # all values of vm must refer to distinct columns

  # would be insane to have super-long TR (suggests milliseconds, not seconds passed in)
  checkmate::assert_number(tr, lower = 0.01, upper = 20, null.ok = TRUE)
  checkmate::assert_string(fmri_file_regex, null.ok = TRUE)
  checkmate::assert_string(run_number_regex, null.ok = TRUE)
  checkmate::assert_integerish(drop_volumes)
  checkmate::assert_character(glm_software)
  checkmate::assert_logical(use_preconvolve, null.ok = FALSE)

  glm_software <- tolower(glm_software)
  checkmate::assert_subset(glm_software, c("fsl", "spm", "afni"))
  checkmate::assert_integerish(n_expected_runs, lower = 1L, null.ok = TRUE)
  if (checkmate::test_string(lgr_threshold)) {
    checkmate::assert_subset(lgr_threshold, c("off", "fatal", "error", "warn", "info", "debug", "trace", "all"))
  } else if (checkmate::test_number(lgr_threshold)) {
    checkmate::assert_integerish(lgr_threshold, lower = 0, len = 1L, all.missing = TRUE)
    lgr_threshold <- as.integer(lgr_threshold)
  } else {
    lgr_threshold <- "info" # default to info in case of weird input
  }

  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")
  lg$set_threshold(lgr_threshold)

  if (!basename(output_directory) == analysis_name) {
    lg$info("Appending analysis_name %s to output_directory %s", analysis_name, output_directory)
    output_directory <- file.path(output_directory, analysis_name)
  }

  # setup output directory, if needed
  if (!dir.exists(output_directory)) {
    lg$info("Setting up output directory for pipeline: %s", output_directory)
    dir.create(output_directory, recursive = TRUE)
  }

  # always use the full path internally
  output_directory <- normalizePath(output_directory)

  # validate and fill in variable mapping vector (if user only passes some fields)
  default_vm <- c(
    id = "id", session = "session",  run_number = "run_number", block_number = "block_number", trial = "trial", 
    run_trial = "run_trial", mr_dir = "mr_dir", run_nifti = "run_nifti"
  )

  default_vm[names(vm)] <- vm # override defaults with user inputs
  vm <- default_vm # reassign full vm

  # validate completeness of trial data
  trial_data <- validate_input_data(trial_data, vm, lg, level="trial")

  # create run data, if needed
  if (is.null(run_data)) {
    lg$info("Distilling run_data object from trial_data by finding variables that vary at run level")
    idcols <- vm[c("id", "session", "run_number")]

    variation_df <- trial_data %>%
      dplyr::select(-dplyr::all_of(idcols)) %>%
      aggregate(by = trial_data[, idcols], FUN = n_distinct)

    # should include the id and run columns
    one_cols <- names(which(sapply(variation_df, function(col) {
      all(col == 1)
    }) == TRUE))

    keepcols <- union(idcols, one_cols)
    lg$info("Retaining columns: %s", paste(keepcols, collapse = ", "))

    # Retain first row for each combination of idcols, keep all columns that vary at run level
    run_data <- trial_data[!duplicated(trial_data[, keepcols]), keepcols, drop=FALSE]
  }

  # check completeness of run data and correct any data expectation problems
  run_data <- validate_input_data(run_data, vm, lg, level="run")

  # handle setup of TR as global versus run-level
  if (is.null(tr) && is.null(run_data[["tr"]])) {
    stop("A TR was not provided in the setup_glm_pipeline call or in run_data. Please provide one or the other.")
  } else if (!is.null(tr) && is.null(run_data[["tr"]])) {
    # single TR applies to all runs. Copy it to run_data as a variable
    run_data$tr <- tr
  } else if (!is.null(tr) && !is.null(run_data[["tr"]])) {
    warning(
      "TR provided in both the setup_glm_pipeline call and in run_data. ",
      "We will use the values in run_data and ignore the tr argument to setup_glm_pipeline."
    )
  }

  # final test on TR
  if (!checkmate::test_numeric(run_data$tr, lower=0.01, upper=20)) {
    stop("Please provide a TR in seconds.")
  }

  # create subject data
  if (is.null(subject_data)) {
    lg$info("Distilling subject_data object from trial_data by finding variables that vary at subject level")
    idcols <- vm[c("id", "session")]

    variation_df <- trial_data %>%
      dplyr::select(-dplyr::all_of(idcols)) %>%
      aggregate(by = trial_data[, idcols], FUN = n_distinct)

    # should include the id column itself
    one_cols <- names(which(sapply(variation_df, function(col) {
      all(col == 1)
    }) == TRUE))

    keepcols <- union(idcols, one_cols)
    lg$info("Retaining columns: ", paste(one_cols, collapse = ", "))

    # Retain first row for each combination of idcols, keep all columns that vary at subject level
    subject_data <- trial_data[!duplicated(trial_data[, keepcols]), keepcols, drop = FALSE]
  }

  # check completeness of subject data and correct any data expectation problems
  subject_data <- validate_input_data(subject_data, vm, lg, level="subject")

  # handle optional PPI data
  ppi_data <- validate_input_data(ppi_data, vm, lg, level="ppi")

  # convert all names to internal conventions from this point forward
  names(trial_data) <- names_to_internal(trial_data, vm)
  names(run_data) <- names_to_internal(run_data, vm)
  names(subject_data) <- names_to_internal(subject_data, vm)
  names(ppi_data) <- names_to_internal(ppi_data, vm)

  # whether to run a 2-level or 3-level analysis
  multi_run <- ifelse(length(unique(trial_data$run_number)) > 1L, TRUE, FALSE)

  # TODO: should probably look at names in subject, run, and trial data to make sure they all line up

  # compare ids in subject_data, run_data, and trial_data
  subject_data_ids <- unique(subject_data$id)
  run_data_ids <- unique(run_data$id)
  trial_data_ids <- unique(trial_data$id)

  all_ids <- named_list(subject_data_ids, run_data_ids, trial_data_ids)
  match_ids <- Reduce(intersect, all_ids) # only keep ids present at all three levels
  union_ids <- Reduce(union, all_ids) # only keep ids present at all three levels

  # print all pairwise differences in ids
  setdiff_list_combn(all_ids)

  if (length(match_ids) == 0L) {
    msg <- "No ids are in common across subject_data, run_data, and trial_data! We cannot proceed with setup."
    lg$error(msg)
    stop(msg)
  } else if (length(match_ids) < length(union_ids)) {
    lg$warn("The ids in subject_data, run_data, and trial_data are not identical. Only the ids in common will be analyzed!")
    lg$warn(glue("Number of non-matching ids: {length(union_ids) - length(match_ids)}"))
    lg$warn(glue("Dropped ids: {paste(setdiff(union_ids, match_ids), collapse=', ')}"))
  }

  subject_data <- subject_data %>% dplyr::filter(.data$id %in% !!match_ids)
  run_data <- run_data %>% dplyr::filter(.data$id %in% !!match_ids)
  trial_data <- trial_data %>% dplyr::filter(.data$id %in% !!match_ids)

  # enforce that subject id + session must be unique (only one row per combination)
  subj_counts <- subject_data %>% dplyr::count(.data$id, .data$session)
  if (any(subj_counts$n > 1)) {
    subj_dupes <- subj_counts %>% dplyr::filter(.data$n > 1)
    msg <- "At least one id + session combination in subject_data is duplicated. All rows in subject_data must represent unique observations!"
    lg$error(msg)
    lg$error("Problematic entries: ")
    lg$error("%s", capture.output(print(subject_data %>% dplyr::inner_join(subj_dupes, by=c("id", "session")))))
    stop(msg)
  }

  # force as.character so that class attributes always match between missing runs and valid runs that are truncated (affects rbindlist)
  # https://github.com/Rdatatable/data.table/issues/3911
  if ("mr_dir" %in% names(subject_data)) subject_data$mr_dir <- as.character(subject_data$mr_dir)
  if ("mr_dir" %in% names(run_data)) run_data$mr_dir <- as.character(run_data$mr_dir)
  if ("run_nifti" %in% names(run_data)) run_data$run_nifti <- as.character(run_data$run_nifti)
  if ("confound_input_file" %in% names(run_data)) run_data$confound_input_file <- as.character(run_data$confound_input_file)

  # distill a session-level dataset so that session predictors can be exposed explicitly in pooled L2 models
  session_data <- derive_session_data(run_data = run_data, subject_data = subject_data, lg = lg)

  if (!is.null(l1_models)) {
    if (checkmate::test_string(l1_models)) {
      if (l1_models == "prompt") l1_models <- build_l1_models(trial_data = trial_data)
      else if (checkmate::test_file_exists(l1_models)) l1_models <- build_l1_models(from_spec_file = l1_models) # treat input as YAML spec file
      else stop("Don't know how to interpret l1_models argument: ", l1_models)
    } else if (!checkmate::test_class(l1_models, "l1_model_set")) {
      stop("l1_models argument is not of class l1_model_set. Use build_l1_models to create this.")
    }
  } # else allow nulls in case user wants to specify things later

  if (!is.null(l2_models)) {
    if (checkmate::test_string(l2_models)) {
      if (l2_models != "prompt" && !checkmate::test_file_exists(l2_models)) {
        stop("Don't know how to interpret l2_models argument: ", l2_models)
      }
    } else if (!checkmate::test_class(l2_models, "hi_model_set")) {
      stop("l2_models argument is not of class hi_model_set. Use build_l2_models to create this.")
    }
  } # else allow nulls in case user wants to specify things later
  
  if (!is.null(l3_models)) {
    if (checkmate::test_string(l3_models)) {
      if (l3_models != "prompt" && !checkmate::test_file_exists(l3_models)) {
        stop("Don't know how to interpret l3_models argument: ", l3_models)
      }
    } else if (!checkmate::test_class(l3_models, "hi_model_set")) {
      stop("l3_models argument is not of class hi_model_set. Use build_l3_models to create this.")
    }
  } # else allow nulls in case user wants to specify things later

  gpa <- list(
    # metadata
    analysis_name = analysis_name,
    scheduler = scheduler,
    output_directory = output_directory,
    subject_data = subject_data,
    run_data = run_data,
    session_data = session_data,
    trial_data = trial_data,
    ppi_data = ppi_data,
    vm = vm,
    bad_ids = bad_ids,
    tr = tr,
    multi_run = multi_run, # 2- or 3-level analysis
    glm_settings = glm_settings,
    confound_settings = confound_settings,
    n_expected_runs = n_expected_runs,
    output_locations = output_locations,

    # l1 analysis details
    fmri_file_regex = fmri_file_regex,
    fmri_path_regex = fmri_path_regex,
    run_number_regex = run_number_regex,
    drop_volumes = drop_volumes,
    use_preconvolve = use_preconvolve,
    l1_models = l1_models,

    # l2 analysis details
    l2_models = l2_models,

    # l3 analysis details
    l3_models = l3_models,
    glm_software = glm_software,
    finalize_complete = FALSE,
    lgr_threshold = lgr_threshold
  )

  # validate and populate any other pipeline details before execution
  # gpa <- finalize_pipeline_configuration(gpa)

  class(gpa) <- c("list", "glm_pipeline_arguments")

  # populate $output_locations
  gpa <- setup_output_locations(gpa, lg)

  # copy in settings passed by user
  gpa$parallel <- parallel

  # add name of node/host on which this is run (useful for logic about different compute environments)
  info <- Sys.info()
  gpa$nodename <- info["nodename"]
  gpa$sys_info <- info # populate full system information to object

  # populate $parallel
  gpa <- setup_parallel_settings(gpa, lg)

  if (checkmate::test_string(gpa$l2_models)) {
    if (gpa$l2_models == "prompt") {
      gpa <- build_l2_models(gpa)
    } else if (checkmate::test_file_exists(gpa$l2_models)) {
      gpa <- build_l2_models(gpa, from_spec_file = gpa$l2_models)
    }
  }

  if (checkmate::test_string(gpa$l3_models)) {
    if (gpa$l3_models == "prompt") {
      gpa <- build_l3_models(gpa)
    } else if (checkmate::test_file_exists(gpa$l3_models)) {
      gpa <- build_l3_models(gpa, from_spec_file = gpa$l3_models)
    }
  }

  # ensure glm_settings is a list before any compute environment checks
  if (is.null(gpa$glm_settings) || identical(gpa$glm_settings, "default")) {
    gpa$glm_settings <- list(fsl = list(), afni = list(), spm = list())
  } else if (is.list(gpa$glm_settings)) {
    if (is.null(gpa$glm_settings$fsl)) gpa$glm_settings$fsl <- list()
    if (is.null(gpa$glm_settings$afni)) gpa$glm_settings$afni <- list()
    if (is.null(gpa$glm_settings$spm)) gpa$glm_settings$spm <- list()
  }

  # normalize longitudinal model signatures early for downstream compatibility checks
  gpa <- normalize_longitudinal_model_signatures(gpa, lg)

  # initial checks on compute environment
  test_compute_environment(gpa, what = unique(c("r", gpa$glm_software)), stop_on_fail = FALSE)
  if ("spm" %in% gpa$glm_software) {
    test_spm_compute_environment(gpa, stop_on_fail = FALSE)
  }

  return(gpa)
}

#' Derive the internal session-level data.frame from run- and subject-level inputs
#'
#' \code{session_data} is not a user-supplied setup object. \code{setup_glm_pipeline()}
#' derives it after validating and normalizing \code{run_data} and \code{subject_data}.
#'
#' @noRd
derive_session_data <- function(run_data, subject_data, lg) {
  checkmate::assert_data_frame(run_data)
  checkmate::assert_data_frame(subject_data)
  checkmate::assert_subset(c("id", "session", "run_number"), names(run_data))
  checkmate::assert_subset(c("id", "session"), names(subject_data))

  run_cols <- setdiff(names(run_data), c("id", "session", "run_number"))
  session_key <- interaction(run_data$id, run_data$session, drop = TRUE, lex.order = TRUE)
  constant_cols <- run_cols[vapply(run_cols, function(cc) {
    split_vals <- split(run_data[[cc]], session_key)
    all(vapply(split_vals, function(xx) {
      uu <- unique(xx[!is.na(xx)])
      length(uu) <= 1L
    }, logical(1)))
  }, logical(1))]

  if (length(constant_cols) > 0L) {
    run_session_data <- run_data %>%
      dplyr::group_by(.data$id, .data$session) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(constant_cols), function(xx) {
          non_missing <- xx[!is.na(xx)]
          if (length(non_missing) > 0L) non_missing[[1L]] else xx[[1L]]
        }),
        .groups = "drop"
      )
  } else {
    run_session_data <- run_data %>% dplyr::distinct(.data$id, .data$session)
  }

  subject_session_data <- subject_data %>% dplyr::distinct(.data$id, .data$session, .keep_all = TRUE)
  shared_cols <- intersect(
    setdiff(names(run_session_data), c("id", "session")),
    setdiff(names(subject_session_data), c("id", "session"))
  )

  if (length(shared_cols) > 0L) {
    cmp <- dplyr::left_join(
      run_session_data %>% dplyr::select("id", "session", dplyr::all_of(shared_cols)),
      subject_session_data %>% dplyr::select("id", "session", dplyr::all_of(shared_cols)),
      by = c("id", "session"), suffix = c(".run", ".subject")
    )

    conflicting <- character(0)
    for (cc in shared_cols) {
      run_col <- cmp[[paste0(cc, ".run")]]
      subj_col <- cmp[[paste0(cc, ".subject")]]
      mismatch <- !is.na(run_col) & !is.na(subj_col) & run_col != subj_col
      if (any(mismatch)) conflicting <- c(conflicting, cc)
    }

    if (length(conflicting) > 0L) {
      lg$warn(
        "Session-level variables conflict between run_data and subject_data for %s. Preferring subject_data values.",
        paste(conflicting, collapse = ", ")
      )
    }
  }

  session_data <- dplyr::full_join(
    run_session_data, subject_session_data,
    by = c("id", "session"), suffix = c(".run", ".subject")
  )

  if (length(shared_cols) > 0L) {
    for (cc in shared_cols) {
      run_col <- paste0(cc, ".run")
      subj_col <- paste0(cc, ".subject")
      session_data[[cc]] <- dplyr::coalesce(session_data[[subj_col]], session_data[[run_col]])
    }

    drop_cols <- c(paste0(shared_cols, ".run"), paste0(shared_cols, ".subject"))
    session_data <- session_data %>% dplyr::select(-dplyr::all_of(drop_cols))
  }

  session_data <- session_data %>% dplyr::arrange(.data$id, .data$session)
  class(session_data) <- c("bg_session_data", class(session_data))

  if (length(constant_cols) > 0L) {
    lg$info(
      "Derived session_data from run_data using %d columns constant within id/session.",
      length(constant_cols)
    )
  } else {
    lg$info("No run-level columns were constant within id/session while deriving session_data.")
  }

  lg$info("Session data contains %d rows and %d columns.", nrow(session_data), ncol(session_data))
  return(session_data)
}

#' Configure the compute environment for pipeline jobs
#'
#' @details Set shell initialization commands in \code{gpa$parallel$compute_environment}.
#'   Use \code{global} for commands shared by all jobs and the \code{fsl},
#'   \code{afni}, or \code{spm} entry for commands specific to an engine.
#'   For SPM, make the configured MATLAB-compatible executable available here
#'   and set the SPM installation path in \code{gpa$glm_settings$spm$spm_path}.
#'
#' @param gpa a \code{glm_pipeline_arguments} object
#' @param preselect_action optional menu action to apply without prompting.
#' @return a modified gpa object containing settings for the compute environment
#' @export
setup_compute_environment <- function(gpa, preselect_action = NULL) {

  if(!is.null(preselect_action)) {
    if (preselect_action %in% c(0L, 5L)){
      cat(c("The environment is configured.", "\n"))
    }
  } else {
    cat(c(
      "\n",
      "We will now handle setup of the compute environment on your system.",
      "The commands you provide here are included at the beginning of scripts produced",
      "by the pipeline in order to configure the environment for your system.",
      "For example, load the FSL or AFNI commands used by your jobs, or the",
      "MATLAB-compatible executable used for SPM. On a cluster, use the module",
      "commands for your site.", "\n"
    ), sep = "\n")
  }
  sanitize_compute_environment <- function(ll = NULL) {
    if (!is.list(ll)) ll <- list()
    if (!is.null(ll$global)) checkmate::assert_character(ll$global)
    if (!is.null(ll$afni)) checkmate::assert_character(ll$afni)
    if (!is.null(ll$fsl)) checkmate::assert_character(ll$fsl)
    if (!is.null(ll$spm)) checkmate::assert_character(ll$spm)
    return(ll)
  }

  gpa$parallel$compute_environment <- sanitize_compute_environment(gpa$parallel$compute_environment)

  get_lines <- function(extant = NULL, level = NULL) {
    lines <- extant # return starts with current input
    add_more <- 1
    while (add_more != 3L) {
      cat(c(
        "\n-----------",
        glue("Current commands in the {level} compute environment: "),
        paste("  ", lines, collapse = "\n"),
        "\n-----------"
      ), sep = "\n")


      add_more <- menu(c("Add command", "Delete command", "Done with command setup"),
        title = paste(toupper(level), "compute environment menu")
      )

      if (add_more == 3L) {
        break
      } else if (add_more == 2L) {
        if (length(lines) == 0L) {
          cat("No existing commands to delete\n")
        } else {
          which_del <- menu(lines, title = "Which command would you like to delete?")
          if (which_del == 0L) {
            cat("Canceling out of deletion\n")
          } else {
            cat("Deleting '", lines[which_del], "'\n", sep = "")
            lines <- lines[-which_del]
          }
        }
      } else if (add_more == 1L) {
        cat(c(
          paste("Please enter commands to be included in the", level, "compute environment."),
          "These will be included in shell scripts and should be executable commands such as",
          "export PATH=/local/dir:$PATH or module load c3d/1.1.0.\n",
          "When you are finished entering commands, just press return on a blank line."
        ), sep = "\n")

        cmds <- scan(what = "character", sep = "\n")
        cat(
          glue("Adding the following commands to the {level} compute environment:"),
          paste("  ", cmds, collapse = "\n"),
          sep = "\n"
        )
        lines <- c(lines, cmds)
      }
    }

    return(lines)
  }

  gpa$parallel$compute_environment <- sanitize_compute_environment(gpa$parallel$compute_environment)

  setup_global_compute <- function(gpa) {
    # global
    cat(c(
      "Let's define any commands you want run at the beginning of *all* scripts (i.e., global).",
      "For example, on a cluster, you might have a command like module use /proj/mnhallqlab/sw/modules",
      "to point the system toward a local directory containing modules for custom programs.",
      "\n",
      "On a local machine, this might include commands like export PATH=/my/local/bin:$PATH",
      "to include certain folders on the system path."
    ), sep = "\n")

    gpa$parallel$compute_environment$global <- get_lines(gpa$parallel$compute_environment$global, "global")
    return(gpa)
  }

  setup_afni_compute <- function(gpa) {
    cat(c(
      "Please provide any commands needed to make AFNI available for jobs that require this.",
      "For example, on a cluster, you might have a command like module load afni/23.0.07."
    ), sep = "\n")

    gpa$parallel$compute_environment$afni <- get_lines(gpa$parallel$compute_environment$afni, "afni")
    return(gpa)
  }

  setup_fsl_compute <- function(gpa) {
    cat(c(
      "Please provide any commands needed to make FSL available for jobs that require this.",
      "For example, on a cluster, you might have a command like module load fsl/6.0.6."
    ), sep = "\n")

    gpa$parallel$compute_environment$fsl <- get_lines(gpa$parallel$compute_environment$fsl, "fsl")
    return(gpa)
  }

  setup_spm_compute <- function(gpa) {
    cat(c(
      "Please provide any commands needed to make MATLAB or Octave available for SPM jobs.",
      "For example, on a cluster, you might have a command like module load matlab/R2017b,",
      "or module load octave/7.3.0."
    ), sep = "\n")

    gpa$parallel$compute_environment$spm <- get_lines(gpa$parallel$compute_environment$spm, "spm")
    return(gpa)
  }

  setup_r_compute <- function(gpa) {
    cat(c(
      "Please provide any commands needed to make R available for jobs that require this.",
      "For example, on a cluster, you might have a command like module load r/4.2.1."
    ), sep = "\n")

    gpa$parallel$compute_environment$r <- get_lines(gpa$parallel$compute_environment$r, "R")
    return(gpa)
  }

  action <- 1
  if(!is.null(preselect_action)) {
    if (preselect_action %in% c(0L, 5L)){
      action <- preselect_action
    } else {
      print("Invalid preselected action. Please select from the menu.")
      action <- 1
    }
  }
  while (!action %in% c(0L, 5L)) {
    cat(c(
      "Here are the current compute environment settings:",
      "\nGlobal (included in all job scripts)",
      paste("  ", gpa$parallel$compute_environment$global, collapse = "\n"),
      "\nAFNI",
      paste("  ", gpa$parallel$compute_environment$afni, collapse = "\n"),
      "\nFSL",
      paste("  ", gpa$parallel$compute_environment$fsl, collapse = "\n"),
      "\nSPM (MATLAB/Octave)",
      paste("  ", gpa$parallel$compute_environment$spm, collapse = "\n"),
      "\nR",
      paste("  ", gpa$parallel$compute_environment$r, collapse = "\n"), "\n"
    ), sep = "\n")

    action <- menu(c("Modify global", "Modify AFNI", "Modify FSL", "Modify SPM (MATLAB/Octave)", "Modify R", "Finish setup"), title = "What would you like to do?")
    if (action == 1L) {
      gpa <- setup_global_compute(gpa)
    } else if (action == 2L) {
      gpa <- setup_afni_compute(gpa)
    } else if (action == 3L) {
      gpa <- setup_fsl_compute(gpa)
    } else if (action == 4L) {
      gpa <- setup_spm_compute(gpa)
    } else if (action == 5L) {
      gpa <- setup_r_compute(gpa)
    } else if (action == 6L) {
      break
    }
  }
  return(gpa)
}



#' helper function to verify contents of subject, run, block, and trial data
#' @importFrom dplyr is_grouped_df ungroup
#' @keywords internal
validate_input_data <- function(df, vm, lg, level = "trial") {
  if (is.null(df)) return(NULL) # allow quick NULL return on NULL input
  checkmate::assert_data_frame(df)

  # Pipeline metadata uses ordinary data.frame subsetting semantics. Copy only
  # incoming data.tables so validation cannot mutate the caller by reference.
  if (data.table::is.data.table(df)) {
    df <- as.data.frame(df)
  }

  # having a grouped data.frame can cause problems (e.g., having the grouping variable unexpectedly come back as a column in select)
  if (is_grouped_df(df)) df <- df %>% ungroup()

  # enforce id column in all input data
  if (!vm["id"] %in% names(df)) {
    msg <- sprintf("Cannot find expected id column %s in %s data. Is your variable mapping (vm) correct?", vm["id"], level)
    lg$error(msg)
    stop(msg)
  }

  # enforce that id is stored as a character vector so that joins between data.frames operate on a consistent character type
  # regardless of how the id was coded in the input data
  lg$debug("Converting id to character in %s data", level)
  df[[vm["id"]]] <- as.character(df[[vm["id"]]])

  if (!vm["session"] %in% names(df)) {
    lg$debug("Adding session = 1 to %s data", level)
    df[[vm["session"]]] <- 1L
  }

  # for trial and run data, code default run_number of 1, if missing
  if (level != "subject")  {
    if (!vm["run_number"] %in% names(df)) {
      lg$debug("Adding run_number = 1 to %s data", level)
      df[[vm["run_number"]]] <- 1L
    }
    if (!checkmate::test_integerish(df[[vm["run_number"]]], lower=1L, upper=1000L)) {
      msg <- sprintf("%s must contain positive integers", vm["run_number"])
      lg$error(msg)
      stop(msg)
    }
  }

  # verify that block data contain block-level variation
  if (level == "block") {
    if (!vm["block_number"] %in% names(df)) {
      msg <- sprintf("Cannot find expected block column %s in %s data. Is your variable mapping (vm) correct?", vm["block_number", level])
      lg$error(msg)
      stop(msg)
    }
  }

  # add relevant class information so that data summary functions work as expected
  # this is a bit poorly designed because we haven't renamed the columns to the internal standard yet, so for a second, the columns
  # could be wrong. But we need to defer the renaming until after the validation to line up the data.frames with each other above
  if (level == "subject" & !inherits(df, "bg_subject_data")) {
    class(df) <- c("bg_subject_data", class(df))
  } else if (level == "run" & !inherits(df, "bg_run_data")) {
    class(df) <- c("bg_run_data", class(df))
  } else if (level == "block" & !inherits(df, "bg_block_data")) {
    class(df) <- c("bg_block_data", class(df))
  } else if (level == "trial" & !inherits(df, "bg_trial_data")) {
    class(df) <- c("bg_trial_data", class(df)) 
  } else if (level == "ppi" & !inherits(df, "bg_ppi_data")) {
    class(df) <- c("bg_ppi_data", class(df))
  }

  return(df)
}

#' internal function to validate the contents of a block data.frame
#' @importFrom dplyr group_by summarise n n_distinct filter
#' @keywords internal
validate_block_data <- function(df) {
  checkmate::assert_class(df, "bg_block_data")

  # check that a) there are no duplicate block_numbers in each run counts of 
  counts <- df %>%
    dplyr::group_by(.data$id, .data$session, .data$run_number) %>%
    dplyr::summarise(n_blocks = dplyr::n_distinct(.data$block_number), n_rows = dplyr::n(), duplicates = any(duplicated(.data$block_number)), .groups = "drop")

  if (all(counts$n_rows) == 1L) {
    msg <- "All runs in block_data contain only a single block. It is unclear why block_data are provided."
    lg$error(msg)
    stop(msg)
  }
  
  if (any(counts$n_blocks) == 1L) {
    single_block <- counts %>%
      dplyr::filter(.data$n_blocks == 1L)

    msg <- sprintf("%d runs contained only a single block. Please verify that this matches your expectations!", nrow(single_block))
    lg$warn(msg)
    lg$warn("%s", capture.output(print(single_block)))
    warning(msg)
  }

  if (any(counts$duplicates)) {
    dupes <- counts %>%
      dplyr::filter(.data$duplicates == TRUE)
    
    msg <- "Duplicate blocks identified for some runs. Unable to proceed until this is resolved."
    lg$error(msg)
    lg$error("%s", capture.output(print(dupes)))
    stop(msg)
  }

}

#' Summarize the contents of the subject data
#' @param object A \code{bg_subject_data} object containing subject/session-level data.
#' @param ... ignored.
#' @details This is typically run using \code{summary(gpa$trial_data)}
#' @export
summary.bg_subject_data <- function(object, ...) {


}

#' Summarize the contents of the subject data
#' @param object A \code{bg_block_data} object containing block-level data.
#' @param ... ignored.
#' @details This is typically run using \code{summary(gpa$block_data)}
#' @export
summary.bg_block_data <- function(object, ...) {
  cat(
    "Summary of block data:",
    "Number of ids: ", length(unique(object$id)), "\n"
  )
}
