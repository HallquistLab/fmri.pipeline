#' This function submits a single script to a high-performance cluster using a scheduler (Slurm or TORQUE).
#' It accepts a vector of arguments to be passed to the scheduler and
#' a vector of environment variables that should be passed to the compute node at job execution.
#'
#' The function returns the jobid of the scheduled job.
#'
#' @param script A script that should be executed by the scheduler. This can contain scheduler directives, but in the
#'      case of conflicts, the directives passed with \code{sched_args} will take precedence.
#' @param scheduler Which scheduler to use for job submission. Options are 'qsub', 'torque', 'sbatch', 'slurm', or 'sh'.
#'      The terms 'qsub' and 'torque' are aliases (where 'torque' submits via the qsub command). Likewise for 'sbatch'
#'      and 'slurm'. The scheduler 'sh' does not submit to any scheduler at all, but instead executes the command
#'      immediately via sh.
#' @param sched_args A character vector of arguments to be included in the scheduling command. On TORQUE, these
#'      will typically begin with '-l' such as '-l wall_time=10:00:00'.
#' @param env_variables A named character vector containing environment variables and their values to be passed
#'      to the \code{script} at execution time. This is handled by the '-v' directive on TORQUE clusters and
#'      by '--export' on Slurm clusters. The names of this vector are the environment variable names and
#'      the values of the vector are the environment variable values to be passed in.
#'      If you want to propagate the current value of an environment variable to the compute node at runtime,
#'      use NA as the value of the element in \code{env_variables}. See examples.
#' @param export_all Whether to export all environment variables to the compute node at runtime. Default: FALSE
#' @param echo Whether to echo the job submission command to the terminal at the time it is scheduled. Default: TRUE.
#' @param fail_on_error Whether to stop execution of the script (TRUE), or issue a warning (FALSE) if the job
#'      submission fails. Defaults to FALSE (i.e., issue a warning).
#' @param wait_jobs a character string of jobs or process ids that should complete before this job is executed
#' @param wait_signal on torque or slurm clusters, the signal that should indicate that parent jobs have finished.
#' @param repolling_interval The number of seconds to wait before rechecking job status (used only for local scheduler)
#' @param tracking_sqlite_db Optional SQLite database used for job tracking.
#' @param tracking_args Optional job-tracking metadata passed through to the tracking database.
#'
#' @return A character string containing the jobid of the scheduled job, or NULL
#'   after a failed submission when fail_on_error is FALSE.
#'
#' @examples
#' \dontrun{
#'   #simple PBS submission
#'   cluster_job_submit(
#'     "myscript.bash",
#'     scheduler = "torque",
#'     sched_args = c("-l walltime=10:00:00", "-l nodes=1:ppn=20"),
#'     env_variables = c(RUN_INDEX = 2, MODEL_NAME = "FSE21")
#'   )
#'
#'   #To forward environment variables without explicitly providing values. Note that these must
#'   #  be in R's system environment (cf. Sys.getenv) at execution time to forward as expected.
#'   cluster_job_submit('myscript.sbatch', scheduler="slurm",
#'      sched_args=c('-p general', '-N 1', '-n 12', '--mem=10g', '-t 02-00:00:00'),
#'      env_variables=c(RUN_INDEX=2, R_HOME=NA, JAVA_HOME=NA))
#' }
#'
#' @author Michael Hallquist
#' @importFrom tools file_path_sans_ext
#' @importFrom checkmate assert_character assert_subset
#' @importFrom utils write.csv
#' @export
cluster_job_submit <- function(script, scheduler="slurm", sched_args=NULL,
                           env_variables=NULL, export_all=FALSE, echo=TRUE, fail_on_error=FALSE, 
                           wait_jobs=NULL, wait_signal="afterok", repolling_interval=60, 
                           tracking_sqlite_db=NULL, tracking_args=list()) {

  checkmate::assert_string(script)
  checkmate::assert_file_exists(script)
  checkmate::assert_string(scheduler)
  checkmate::assert_subset(scheduler, c("qsub", "torque", "sbatch", "slurm", "sh", "local"))
  checkmate::assert_logical(export_all, max.len = 1L)
  checkmate::assert_logical(echo, max.len = 1L)
  checkmate::assert_logical(fail_on_error, max.len=1L)
  if (!is.null(wait_jobs)) checkmate::assert_character(wait_jobs)
  if (is.character(tracking_args) || is.numeric(tracking_args)) tracking_args <- as.list(tracking_args) # coerce tracking args to list
  if (length(tracking_args) > 0 && is.null(tracking_sqlite_db)) {
    warning("Tracking arguments provided to `cluster_job_submit` but `tracking_sqlite_db` is NULL")
  }

  if (scheduler %in% c("torque", "qsub")) {
    scheduler <- "qsub" # simpler internal tracking
    if (isTRUE(export_all)) {
      sched_args <- c("-V", sched_args)
    } # directive to export all environment variables to script
  } else if (scheduler %in% c("slurm", "sbatch")) {
    scheduler <- "sbatch"
    if (isTRUE(export_all)) {
      env_variables <- c(ALL = NA, env_variables)
    } # directive to export all environment variables to script
  } else if (scheduler %in% c("sh", "local")) {
    scheduler <- "sh"
    if (!is.null(sched_args)) {
      message("Omitting scheduler arguments for sh/local execution")
    }
    sched_args <- NULL # not relevant
  }

  # Scheduler arguments are just pasted together with spaces.
  # Thus, arguments like '--mem=5g' and '-n 12' are not handled differently
  if (!is.null(sched_args)) { sched_args <- paste(sched_args, collapse=" ") }

  #subfunction to handle variable=value and variable combinations
  paste_args <- function(str_vec) {
    nms <- names(str_vec)
    stopifnot(!is.null(nms))
    sapply(seq_along(str_vec), function(x) {
      if (is.na(str_vec[x])) {
        return(nms[x]) #just the name of the env variable (forwards from environment)
      } else {
        # force argument to be quoted to avoid problems with spaces
        val <- shQuote(str_vec[x], type = "sh")
        return(paste0(nms[x], "=", val))
      }
    })
  }

  #pass through environment variables
  if (!is.null(env_variables)) {
    env_variables <- paste_args(env_variables) #convert to appropriate name-value pairs
    if (scheduler == "qsub") {
      env_variables <- paste("-v", paste(env_variables, collapse=","))
    } else if (scheduler == "sbatch") {
      env_variables <- paste0("--export=", paste(env_variables, collapse=","))
    } else if (scheduler == "sh") {
      env_variables <- paste(sapply(env_variables, function(x) {
        ifelse(grepl("=", x, fixed = TRUE), x, paste0(x, "=\"$", x, "\""))
      }), collapse = " ")
    }

    sched_args <- paste(sched_args, env_variables)
  }

  if (!is.null(wait_jobs)) {
    jcomb <- paste(wait_jobs, collapse = ":") # multiple jobs are separated by colons
    if (scheduler == "qsub") {
      sched_args <- paste0(sched_args, " -W depend=", wait_signal, ":", jcomb)
    } else if (scheduler == "sbatch") {
      sched_args <- paste0(sched_args, " --dependency=", wait_signal, ":", jcomb)
    }
  }

  # use unique temp files to avoid parallel collisions in job tracking
  sub_stdout <- paste0(tempfile(), "_", tools::file_path_sans_ext(basename(script)), "_stdout") 
  sub_stderr <- paste0(tempfile(), "_", tools::file_path_sans_ext(basename(script)), "_stderr")
  sub_pid <- paste0(tempfile(), "_", tools::file_path_sans_ext(basename(script)), "_pid")

  if (scheduler == "sh") {
    # if an R script is provided, execute with Rscript --vanilla as command
    if (grepl(".+\\.R$", script, ignore.case = TRUE)) {
      bin <- "Rscript --vanilla"
    } else {
      bin <- "sh"
    }

    # for local scheduler, we need to hold jobs manually by waiting for relevant parents to complete
    if (!is.null(wait_jobs)) {
      message("Waiting for the following jobs to finish: ", paste(wait_jobs, collapse=","))
      wait_for_job(wait_jobs, repolling_interval = repolling_interval, scheduler = scheduler)
      # ZV: grab return code here for parent job failing?
    }

    # for direct execution, need to pass environment variables by prepending
    cmd <- paste(env_variables, bin, script)
    if (isTRUE(echo)) cat(cmd, "\n") # echo command to terminal
    # submit the job script and return the jobid by forking to background and returning PID
    jobres <- system(paste(cmd, ">", sub_stdout, "2>", sub_stderr, "& echo $! >", sub_pid), wait = FALSE)
    Sys.sleep(.05) #sometimes the pid file is not in place when file.exists executes -- add a bit of time to ensure that it reads
    jobid <- if (file.exists(sub_pid)) scan(file = sub_pid, what = "char", sep = "\n", quiet = TRUE) else ""
  } else {
    cmd <- paste(scheduler, sched_args, script)
    if (isTRUE(echo)) cat(cmd, "\n")
    # submit the job script and return the jobid
    jobres <- system2(scheduler, args = paste(sched_args, script), stdout = sub_stdout, stderr = sub_stderr)
    jobid <- if (file.exists(sub_stdout)) scan(file = sub_stdout, what = "char", sep = "\n", quiet = TRUE) else ""
  }

  joberr <- if (file.exists(sub_stderr)) {
    paste(scan(file = sub_stderr, what = "char", sep = "\n", quiet = TRUE), collapse = ". ")
  } else {
    ""
  }

  if (jobres != 0) {
    jobid <- NULL
    if (isTRUE(fail_on_error)) {
      stop("Job submission failed: ", script, ", error: ", joberr, ", errcode: ", jobres)
    } else {
      warning("Job submission failed: ", script, ", error: ", joberr, ", errcode: ", jobres)
      return(NULL) # no job exists to track or attach command metadata to
    }
  } else {
    jobid <- sub("Submitted batch job ", "", jobid, fixed = TRUE) # replace irrelevant details if needed
    if (!is.null(tracking_args)) tracking_args$status <- "QUEUED" # on successful submission, tracking status defaults to "QUEUED"
  }

  if (!is.null(tracking_args)) {
    # populate tracking arguments that are undefined and can be defined at this point
    tracking_args <- populate_list_arg(tracking_args, "batch_file", script)
    tracking_args <- populate_list_arg(tracking_args, "scheduler", scheduler)
    tracking_args <- populate_list_arg(tracking_args, "scheduler_options", sched_args, append = TRUE)
  }

  # once a job_id has been generated, we add it to the tracking db
  # if a job has a NULL id or the tracking_sqlite_db arg is NULL, the function will return invisible NULL
  insert_tracked_job(sqlite_db = tracking_sqlite_db, job_id = jobid, tracking_args = tracking_args)

  if (!is.null(wait_jobs)) {
    # add any parent jobs using the wait_jobs argument (defaults to last parent id in list)
    add_tracked_job_parent(sqlite_db = tracking_sqlite_db, job_id = jobid, parent_job_id = wait_jobs[length(wait_jobs)], child_level = 1) 
  } else if (!is.null(tracking_args$parent_job_id)) {
    # in the case that a parent job id is passed in through the tracking_args list
    add_tracked_job_parent(sqlite_db = tracking_sqlite_db, job_id = jobid, parent_job_id = tracking_args$parent_job_id, child_level = 1) 
  }
  
  attr(jobid, "cmd") <- cmd # add the command executed as an attribute

  return(jobid)
}



#' helper function to submit a set of shell jobs that are independent of one another
#'
#' @param job_list a list or character vector where each element represents an independent job to execute in a shell environment
#' @param commands_per_cpu how many elements from \code{job_list} are executed by each core within a single job
#' @param cpus_per_job how many cpus/cores are requested for each job
#' @param memgb_per_command amount of memory (RAM) requested for each command (in GB)
#' @param time_per_job amount of time requested for each job
#' @param time_per_command amount of time requested for each individual command
#' @param fork_jobs if TRUE, all jobs within a single batch will run simultaneously using the fork (&) approach.
#'   Each task runs in its own subshell with `set -e`. The batch waits for every task and returns
#'   the first nonzero task status in submission order, or the cleanup status if all tasks succeed.
#' @param pre user-specified code to include in the job script prior to job_list (e.g., module load commands)
#' @param post user-specified cleanup code run after all tasks finish; its success cannot mask a task failure.
#' @param sched_args scheduler directives passed to \code{cluster_job_submit()}.
#' @param env_variables named environment variables passed to \code{cluster_job_submit()}.
#' @param wait_jobs jobs that must complete before submitted jobs run.
#' @param scheduler scheduler backend used for submission.
#' @param job_out_dir the directory where job scripts should be written
#' @param job_script_prefix the filename prefix for each job script
#' @param log_file a csv log file containing the job ids and commands that were executed
#' @param debug a logical indicating whether to actually submit the jobs (TRUE) or just create the scripts for inspection (FALSE)
#'
#' @importFrom tidyr unnest
#' @importFrom checkmate assert_multi_class
#'
#' @export
cluster_submit_shell_jobs <- function(job_list, commands_per_cpu = 1L, cpus_per_job = 8L, memgb_per_command = 8, time_per_job="1:00:00",
  time_per_command = NULL, fork_jobs = TRUE, pre=NULL, post=NULL, sched_args = NULL, env_variables = NULL, wait_jobs = NULL, scheduler="slurm",
  job_out_dir=getwd(), job_script_prefix="job", log_file="cluster_submit_jobs.csv", debug = FALSE)
{
  checkmate::assert_multi_class(job_list, c("list", "character"))
  checkmate::assert_integerish(commands_per_cpu, len = 1L, lower = 1, upper = 1e4)
  checkmate::assert_integerish(cpus_per_job, len=1L, lower=1, upper=1e3)
  checkmate::assert_number(memgb_per_command, lower = 0.01, upper = 1e3) # max 1TB
  checkmate::assert_logical(fork_jobs, len = 1L)
  if (!is.null(time_per_command)) {
    if (!is.null(time_per_job)) {
      message("Using time_per_command specification, even though time_per_job was also provided")
    }

    # convert to lubridate object for math
    time_per_command <- dhms(time_per_command)

    if (isTRUE(fork_jobs)) {
      # if we fork, then the total job duration is just the number of commands per cpu (default 1) multiplied
      # by the per-command time. Note that at present, the forking does not respect a job limit equal to cpus_per_job.
      # So, if you say commands_per_cpu=3 and cpus_per_job=8, then 24 jobs will launch at once. Need a jobs | wc -l waiting approach
      time_per_job <- time_per_command * commands_per_cpu
    } else {
      # if we do not fork, the total time increases by a factor of cpus_per_job
      # this is useful if each command is a multithreaded job (e.g., MCMC sampling multiple chains).
      time_per_job <- time_per_command * commands_per_cpu * cpus_per_job
    }

    # convert back to time string from period
    time_per_job <- as_dhms(time_per_job)
  } else {
    time_per_job <- validate_dhms(time_per_job)
  }

  n_jobs <- ceiling(length(job_list) / (cpus_per_job * commands_per_cpu))

  # how every job submission script begins
  preamble <- c(
    "#!/bin/bash",
    ""
  )

  if (scheduler == "slurm") {
    file_suffix <- ".sbatch"
    preamble <- c(preamble,
      "#SBATCH -N 1", #always single node for now
      paste0("#SBATCH -n ", cpus_per_job),
      paste0("#SBATCH --time=", time_per_job),
      paste0("#SBATCH --mem-per-cpu=", memgb_per_command, "G"),
      'cd "$SLURM_SUBMIT_DIR" || exit $?',
      ""
    )
  } else if (scheduler == "torque") {
    file_suffix <- ".pbs"
    preamble <- c(
      preamble,
      paste0("#PBS -l nodes=1:ppn=", cpus_per_job),
      paste0("#PBS -l pmem=", memgb_per_command, "gb"),
      paste0("#PBS -l walltime=", time_per_job),
      'cd "$PBS_O_WORKDIR" || exit $?',
      ""
    )
  } else if (scheduler == "local") {
    file_suffix <- ".bash"
  }

  # add any user-specified code that should precede each job execution
  # Fail setup immediately, but collect task failures without abandoning siblings.
  preamble <- c(preamble, "set -e", pre, "set +e", "batch_status=0", "batch_pids=", "")

  job_df <- data.frame(
    job_number = rep(seq_len(n_jobs), each = cpus_per_job * commands_per_cpu, length.out = length(job_list)),
    job_id = NA_character_,
    stringsAsFactors = FALSE
  )

  job_df$cmd <- job_list # post-assignment keeps as list, if relevant

  submitted_jobs <- rep(NA_character_, n_jobs)

  for (j in seq_len(n_jobs)) {
    this_run <- job_df$cmd[job_df$job_number == j]

    outfile <- file.path(job_out_dir, paste0(job_script_prefix, "_", j, file_suffix))
    cat(preamble, file = outfile, sep = "\n")

    # Each logical task gets a subshell with errexit, so a later successful
    # command in a multi-command task cannot mask an earlier failure.
    for (task in as.list(this_run)) {
      task_syntax <- c("(", "set -e", task, if (isTRUE(fork_jobs)) ") &" else ")")
      if (isTRUE(fork_jobs)) {
        task_syntax <- c(task_syntax, 'batch_pids="$batch_pids $!"')
      } else {
        task_syntax <- c(task_syntax, "batch_child_status=$?",
          'if [ "$batch_status" -eq 0 ]; then batch_status=$batch_child_status; fi')
      }
      cat(task_syntax, file = outfile, sep = "\n", append = TRUE)
    }

    if (isTRUE(fork_jobs)) {
      cat(c(
        'for batch_pid in $batch_pids; do',
        '  if wait "$batch_pid"; then',
        '    :',
        '  else',
        '    batch_child_status=$?',
        '    if [ "$batch_status" -eq 0 ]; then batch_status=$batch_child_status; fi',
        '  fi',
        'done'
      ), file = outfile, sep = "\n", append = TRUE)
    }

    # Run cleanup after every task has finished, retaining any earlier failure.
    if (!is.null(post)) {
      cat(c("(", "set -e", post, ")", "batch_post_status=$?",
        'if [ "$batch_status" -eq 0 ]; then batch_status=$batch_post_status; fi'),
        file = outfile, sep = "\n", append = TRUE)
    }
    cat('exit "$batch_status"\n', file = outfile, append = TRUE)

    # in debug mode, we only create the scripts, but we do not submit them
    if (isTRUE(debug)) {
      job_id <- paste0("dummy", j)
    } else {
      job_id <- cluster_job_submit(outfile, scheduler = scheduler, sched_args = sched_args, env_variables = env_variables, wait_jobs = wait_jobs)
    }

    job_df$job_id[job_df$job_number == j] <- job_id
    submitted_jobs[j] <- job_id

  }

  # write job log to file
  if (is.list(job_list)) {
    write.csv(job_df %>% unnest(cols="cmd"), file = log_file, row.names = FALSE)
  } else {
    write.csv(job_df, file = log_file, row.names = FALSE)
  }

  return(submitted_jobs)
}
