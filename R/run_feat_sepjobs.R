#' function to submit a set of jobs on a cluster to estimate many Feat level 1 models
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#' @param level GLM level to run.
#' @param model_names an optional level 1 model name used to subset the runs to submit to feat.
#'   If not provided, all level 1 models in \code{gpa} will be submitted for feat estimation.
#' @param rerun a logical indicating whether to re-run an existing directory. Default: FALSE
#' @param wait_for a parent job that should complete before these jobs commence
#' @return a vector of job ids for all files that were submitted to the cluster
#' @importFrom glue glue
#' @export
run_feat_sepjobs <- function(gpa, level=1L, model_names=NULL, rerun=FALSE, wait_for="") {
  # this version of the FSL LVL1 feat estimation creates multiple qsub scripts in a temporary directory
  # where each script has a number of feat calls that are forked, then the script waits for completion
  ## This circumvents the ICS limit on multiple nodes in R using a SNOW cluster.
  ## The primary calculation is how many files there are to run relative to processors and files per processor (chunking)

  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower=1L, upper=3L, len=1L)
  checkmate::assert_character(model_names, null.ok=TRUE)
  checkmate::assert_logical(rerun)

  lg <- lgr::get_logger(paste0("glm_pipeline/l", level, "_estimation"))
  lg$set_threshold(gpa$lgr_threshold)

  # retreive sqlite_db
  tracking_sqlite_db <- gpa$output_locations$sqlite_db
  # retrieve cmd line update function
  upd_job_status_path <- system.file("bin/upd_job_status.R", package = "fmri.pipeline")
  log_flame_fallbacks_path <- system.file("bin/log_flame_runner_fallbacks.R", package = "fmri.pipeline")

  # TODO: support named subsetting in model_names, like list(l1_model="pe_only", l2_model=c("l2_1", "l2_2"))
  # level-specific copy-paste is a bit clunky here, but at least it's clear...
  if (level == 1) {
    if (!checkmate::test_class(gpa$l1_model_setup, "l1_setup")) {
      lg$error("In run_feat_sepjobs, did not find an l1_setup object in gpa$l1_model_setup.")
      lg$error("Make sure to run setup_l1_models before running run_feat_sepjobs")
      return(NULL)
    }

    feat_queue <- gpa$l1_model_setup$fsl

    feat_time   <- gpa$parallel$fsl$l1_feat_time
    feat_memgb  <- gpa$parallel$fsl$l1_feat_memgb
    feat_cpus   <- gpa$parallel$fsl$l1_feat_cpus_per_job # number of cpus per scheduled job
    runsperproc <- gpa$parallel$fsl$l1_feat_runs_per_cpu # number of feat calls per processor
  } else if (level == 2) {
    if (!checkmate::test_class(gpa$l2_model_setup, "l2_setup")) {
      lg$error("In run_feat_sepjobs, did not find an l2_setup object in gpa$l2_model_setup.")
      lg$error("Make sure to run setup_l2_models before running run_feat_sepjobs")
      return(NULL)
    }

    feat_queue <- gpa$l2_model_setup$fsl
    if ("l2_passthrough" %in% names(feat_queue)) {
      feat_queue <- feat_queue %>% dplyr::filter(!(.data$l2_passthrough %in% TRUE))
    }

    feat_time <- gpa$parallel$fsl$l2_feat_time
    feat_memgb <- gpa$parallel$fsl$l2_feat_memgb
    feat_cpus <- 8 # TODO: populate parallel settings for l2
    runsperproc <- 2 # number of feat calls per processor
  } else if (level == 3) {
    if (!checkmate::test_class(gpa$l3_model_setup, "l3_setup")) {
      lg$error("In run_feat_sepjobs, did not find an gpa$l3_model_setup object.")
      lg$error("Make sure to run setup_l3_models before running run_feat_sepjobs")
      return(NULL)
    }

    feat_queue <- gpa$l3_model_setup$fsl

    feat_time <- gpa$parallel$fsl$l3_feat_time
    # memory is for total job, not per cpu at l3
    feat_memgb <- ceiling(as.numeric(gpa$parallel$fsl$l3_feat_memgb) / as.numeric(gpa$parallel$fsl$l3_feat_cpusperjob))
    feat_cpus <- gpa$parallel$fsl$l3_feat_cpusperjob
    runsperproc <- 1 # number of feat calls per processor
  }

  model_column <- paste0("l", level, "_model")
  required_queue_columns <- c(
    model_column, "feat_fsf", "feat_dir", "feat_complete", "feat_failed", "to_run"
  )
  if (!is.data.frame(feat_queue)) {
    stop(
      sprintf(
        "gpa$l%d_model_setup$fsl is not a data.frame. Make sure setup_l%d_models() produced FSL inputs before calling run_feat_sepjobs().",
        level, level
      ),
      call. = FALSE
    )
  }
  missing_queue_columns <- setdiff(required_queue_columns, names(feat_queue))
  if (length(missing_queue_columns) > 0L) {
    stop(
      sprintf(
        paste0(
          "gpa$l%d_model_setup$fsl is missing required column(s): %s. ",
          "The setup step did not produce a valid FSL job table; inspect the setup_l%d_models log for the upstream error."
        ),
        level, paste(missing_queue_columns, collapse = ", "), level
      ),
      call. = FALSE
    )
  }
  if (!is.null(model_names[1L])) {
    feat_queue <- feat_queue %>%
      dplyr::filter(.data[[model_column]] %in% model_names)
  }

  fail_action <- gpa$glm_settings$fsl[[glue("failed_l{level}_folder_action")]]
  incomplete_action <- gpa$glm_settings$fsl[[glue("incomplete_l{level}_folder_action")]]
  if (is.null(fail_action)) fail_action <- "delete" #should be populated in finalize step, but just in case
  if (is.null(incomplete_action)) incomplete_action <- "delete"

  # Refresh from disk before considering any cleanup or submission. A lock covers
  # the entire queued/running lifetime; no output is removed by this function.
  feat_output_directory <- file.path(gpa$output_locations$scheduler_scripts, paste0("feat_l", level))
  statuses <- lapply(seq_len(nrow(feat_queue)), function(i)
    get_feat_status(feat_queue$feat_dir[i], feat_queue$feat_fsf[i], lg = lg))
  feat_job_df <- dplyr::bind_rows(statuses)
  if (!nrow(feat_job_df)) return(invisible(character(0)))
  feat_job_df$to_run <- !feat_job_df$feat_complete | isTRUE(rerun)
  feat_job_df$action <- vapply(seq_len(nrow(feat_job_df)), function(i)
    feat_run_action(feat_job_df[i, ], fail_action, incomplete_action), character(1))
  kept <- feat_job_df$to_run & feat_job_df$feat_dir_exists & feat_job_df$action == "keep"
  if (any(kept)) lg$warn("Keeping existing FEAT output without resubmission: %s", feat_job_df$feat_dir[kept])
  feat_job_df$to_run[kept] <- FALSE
  for (i in which(feat_job_df$to_run)) {
    fsf <- feat_job_df$feat_fsf[i]
    if (dir.exists(feat_lock_dir(fsf))) stop("FEAT model is already queued/running or locked after interruption: ", fsf)
    if (!feat_setup_current(fsf)) stop("FEAT setup is missing or stale; run model setup again: ", fsf)
    if (!identical(feat_read_setup(fsf)$output_dir,
        normalizePath(feat_job_df$feat_dir[i], winslash = "/", mustWork = FALSE))) {
      stop("FEAT output directory differs from the submitted setup: ", fsf)
    }
  }
  to_run <- feat_job_df$feat_fsf[feat_job_df$to_run]

  if (length(to_run) == 0L) {
    lg$warn("No Feat level %d .fsf files to execute.", level)
    return(invisible(character(0)))
  }

  lg$info("About to run the following fsf files in parallel:")
  lg$info("  File: %s", to_run)

  child_log_directory <- Sys.getenv("SLURM_SUBMIT_DIR", unset = "")
  if (!nzchar(child_log_directory)) child_log_directory <- Sys.getenv("PBS_O_WORKDIR", unset = "")
  if (!nzchar(child_log_directory)) child_log_directory <- getwd()

  if (gpa$scheduler %in% c("slurm", "sbatch")) {
    file_suffix <- ".sbatch"
    preamble <- function(job_name, extra, batch_file) {
      c(
        "#!/bin/bash",
        "#SBATCH -N 1", #always single node for now
        paste0("#SBATCH -n ", feat_cpus),
        paste0("#SBATCH --time=", feat_time),
        paste0("#SBATCH --mem-per-cpu=", feat_memgb, "G"),
        ifelse(any(nzchar(wait_for)), paste0("#SBATCH --dependency=afterok:", paste(wait_for, collapse=":")), ""), # allow job dependency on upstream setup
        sched_args_to_header(gpa), # analysis-level SBATCH directives
        paste("#SBATCH -J", scheduler_safe_token(job_name, max_chars = 80L)),
        scheduler_output_directives(gpa$scheduler, child_log_directory, job_name = job_name, extra = extra),
        "",
        "job_id=$SLURM_JOB_ID",
        paste0("job_name=", shQuote(job_name)),
        "scheduler_name=slurm",
        paste0("batch_directory=", shQuote(child_log_directory)),
        paste0("batch_file=", shQuote(batch_file)),
        scheduler_runtime_log_assignment(gpa$scheduler, child_log_directory, job_name = job_name, extra = extra),
        "",
        get_compute_environment(gpa, c("fsl", "r")), # r needed for job tracking sqlite
        "",
        "cd \"$SLURM_SUBMIT_DIR\""
      )
    }

  } else if (gpa$scheduler %in% c("torque", "qsub")) {
    file_suffix <- ".pbs"
    preamble <- function(job_name, extra, batch_file) {
      c(
        "#!/bin/bash",
        paste0("#PBS -l nodes=1:ppn=", feat_cpus),
        paste0("#PBS -l pmem=", feat_memgb, "gb"),
        ifelse(any(nzchar(wait_for)), paste0("#PBS -W depend=afterok:", paste(wait_for, collapse=":")), ""), # allow job dependency on upstream setup
        paste0("#PBS -l walltime=", feat_time),
        sched_args_to_header(gpa), # analysis-level PBS directives
        paste("#PBS -N", substr(scheduler_safe_token(job_name, max_chars = 80L), 1L, 15L)),
        scheduler_output_directives(gpa$scheduler, child_log_directory, job_name = job_name, extra = extra),
        "",
        "job_id=$PBS_JOBID",
        paste0("job_name=", shQuote(job_name)),
        "scheduler_name=pbs",
        paste0("batch_directory=", shQuote(child_log_directory)),
        paste0("batch_file=", shQuote(batch_file)),
        scheduler_runtime_log_assignment(gpa$scheduler, child_log_directory, job_name = job_name, extra = extra),
        "",
        get_compute_environment(gpa, c("fsl", "r")),
        "",
        "cd \"$PBS_O_WORKDIR\""
      )
    }
  }

  if (!file.exists(feat_output_directory)) {
    lg$debug("Creating l%d working directory: %s", level, feat_output_directory)
    dir.create(feat_output_directory, recursive = TRUE)
  }

  estimation_log_dir <- gpa$output_locations$log_directory
  if (is.null(estimation_log_dir) || !nzchar(estimation_log_dir)) {
    estimation_log_dir <- file.path(gpa$output_directory, "logs")
  }
  dir.create(estimation_log_dir, recursive = TRUE, showWarnings = FALSE)
  estimation_log_txt <- if (isTRUE(gpa$log_txt)) {
    file.path(estimation_log_dir, sprintf("l%d_estimation.txt", level))
  } else {
    "NULL"
  }
  estimation_log_json <- if (isTRUE(gpa$log_json)) {
    file.path(estimation_log_dir, sprintf("l%d_estimation.json", level))
  } else {
    "NULL"
  }

  auto_retry_l3_excessive_outliers <- isTRUE(gpa$glm_settings$fsl$auto_retry_l3_excessive_outliers) ||
    is.null(gpa$glm_settings$fsl$auto_retry_l3_excessive_outliers)

  if (level == 3) {
    njobs <- length(to_run) # parallel across slices, one job per model
    feat_binary <- system.file("bin/feat_parallel", package = "fmri.pipeline")
    stopifnot(file.exists(feat_binary))
    if (isTRUE(auto_retry_l3_excessive_outliers)) {
      feat_retry_binary <- system.file("bash/run_feat_with_outlier_retry.sh", package = "fmri.pipeline")
      stopifnot(file.exists(feat_retry_binary))
      lg$info("L3 auto-retry for excessive outlier failures is enabled.")
    } else {
      feat_retry_binary <- NA_character_
      lg$info("L3 auto-retry for excessive outlier failures is disabled.")
    }
    df <- data.frame(fsf = to_run, job = 1:njobs, stringsAsFactors = FALSE)
  } else {
    njobs <- ceiling(length(to_run) / (feat_cpus * runsperproc))
    feat_binary <- "feat"
    feat_retry_binary <- NA_character_
    # use length.out on rep to ensure that the vectors align even if chunks are uneven wrt files to run
    df <- data.frame(
      fsf = to_run,
      job = rep(1:njobs, each = feat_cpus * runsperproc, length.out = length(to_run)), stringsAsFactors = FALSE
    )
  }

  df <- df[order(df$job), ]

  submission_id <- basename(tempfile(pattern = "job"))
  joblist <- rep(NA_character_, njobs)
  for (j in seq_len(njobs)) {
    outfile <- paste0(feat_output_directory, "/featsep_l", level, "_", j, "_", submission_id, file_suffix)
    job_name <- paste0("featsep_l", level, "_", j)
    cat(preamble(job_name, submission_id, outfile), file=outfile, sep="\n")
    thisrun <- with(df, fsf[job==j])
    rows <- match(thisrun, feat_job_df$feat_fsf)
    variables <- c(feat_binary = feat_binary, feat_level = level,
      feat_use_retry = as.integer(level == 3L && isTRUE(auto_retry_l3_excessive_outliers)),
      feat_retry_binary = feat_retry_binary, feat_fallback_logger = log_flame_fallbacks_path,
      feat_log_txt = estimation_log_txt, feat_log_json = estimation_log_json,
      feat_log_threshold = as.character(gpa$lgr_threshold),
      feat_integrity_script = system.file("bin", "feat_integrity.R", package = "fmri.pipeline"),
      feat_reg_script = system.file("bash", "gen_feat_reg_dir", package = "fmri.pipeline"))
    tracking_cmd <- if (!is.null(tracking_sqlite_db)) paste("Rscript", shQuote(upd_job_status_path),
      "--job_id \"$job_id\" --sqlite_db", shQuote(tracking_sqlite_db)) else "true"
    cat(c(paste0(names(variables), "=", shQuote(unname(variables))), job_manifest_shell_function(),
      readLines(system.file("bash", "feat_runner.sh", package = "fmri.pipeline")),
      "pids=()", "job_failed=0",
      "feat_killed() { for pid in \"${pids[@]}\"; do kill -TERM \"$pid\" 2>/dev/null; done; wait; exit 1; }",
      "trap feat_killed TERM INT HUP", paste(tracking_cmd, "--status STARTED")),
      file = outfile, sep = "\n", append = TRUE)
    for (i in rows) {
      args <- c(feat_job_df$feat_fsf[i], feat_job_df$feat_dir[i], submission_id, feat_job_df$action[i])
      if (level == 3L) args <- c(args, as.character(feat_cpus))
      cat(c(paste("feat_runner", paste(shQuote(args, type = "sh"), collapse = " "), "&"), "pids+=(\"$!\")"),
        file = outfile, sep = "\n", append = TRUE)
    }
    cat(c("for pid in \"${pids[@]}\"; do", "  if ! wait \"$pid\"; then job_failed=1; fi", "done",
      "if [ $job_failed -ne 0 ]; then job_status=FAILED; else job_status=COMPLETED; fi",
      paste(tracking_cmd, "--status \"${job_status}\""), "exit $job_failed"),
      file = outfile, sep = "\n", append = TRUE)

    # gather tracking arguments
    if (!is.null(tracking_sqlite_db)) {
      tracking_args <- list(job_name = job_name,
                            batch_directory = feat_output_directory,
                            n_nodes = 1,
                            n_cpus = feat_cpus,
                            wall_time = feat_time,
                            scheduler_options = gpa$parallel$sched_args)
      if (level == 3) {
        tracking_args$mem_total <- feat_memgb
      } else {
        tracking_args$mem_per_cpu <- feat_memgb
      }
      if (any(nzchar(wait_for))) tracking_args$parent_job_id <- wait_for
    } else {
      tracking_args = NULL
    }

    reserved <- character()
    joblist[j] <- tryCatch({
      for (i in rows) {
        fsf <- feat_job_df$feat_fsf[i]
        feat_reserve(fsf, submission_id)
        reserved <- c(reserved, fsf)
        current <- get_feat_status(feat_job_df$feat_dir[i], fsf, lg = lg)
        if ((!isTRUE(rerun) && isTRUE(current$feat_complete)) ||
            !identical(feat_run_action(current, fail_action, incomplete_action), feat_job_df$action[i])) {
          stop("FEAT status changed while reserving the batch. Retry submission: ", fsf)
        }
      }
      job <- cluster_job_submit(outfile, scheduler = gpa$scheduler, tracking_sqlite_db = tracking_sqlite_db, tracking_args = tracking_args)
      if (length(job) != 1L || is.na(job) || !nzchar(job)) stop("FEAT scheduler submission did not return a job ID.")
      job
    }, error = function(e) {
      for (fsf in reserved) feat_release(fsf, submission_id)
      stop(e)
    })
    #joblist[j] <- "dummy"
  }

  # write the list of separate feat qsub jobs that are now queued (so that LVL2 can wait on these)
  # should also return this to the caller as a function?
  writeLines(joblist, con = file.path(feat_output_directory, paste0("sep_l", level, "_jobs.txt")))

  return(joblist)
}
