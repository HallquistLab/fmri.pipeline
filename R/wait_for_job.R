#' This function pauses execution of an R script while a scheduled qsub job is not yet complete.
#'
#' It is intended to give you control over job dependencies within R when the formal PBS
#' depend approach is insufficient, especially in the case of a script that spawns child jobs that
#' need to be scheduled or complete before the parent script should continue.
#'
#' @param job_ids One or more job ids of existing PBS or slurm jobs, or process ids of a local process for
#'   \code{scheduler="sh"}.
#' @param repolling_interval How often to recheck the job status, in seconds. Default: 30
#' @param max_wait How long to wait on the job before giving up, in seconds. Default: 24 hours (86,400 seconds)
#' @param scheduler What scheduler is used for job execution.
#'   Options: c("torque", "qsub", "slurm", "sbatch", "sh", "local")
#' @param quiet If \code{TRUE}, \code{wait_for_job} will not print out any status updates on jobs. If \code{FALSE},
#'   the function prints out status updates for each tracked job so that the user knows what's holding up progress.
#' @param stop_on_timeout Whether to stop if \code{max_wait} is reached.
#'
#' @return Invisibly returns TRUE when all jobs complete successfully, or FALSE
#'   when any job fails or is cancelled (or a non-stopping timeout occurs).
#'
#' @details Note that for the \code{scheduler} argument, "torque" and "qsub" are the same;
#'   "slurm" and "sbatch" are the same, and "sh" and "local" are the same.
#' @examples
#' \dontrun{
#' # example on qsub/torque cluster
#' wait_for_job("7968857.torque01.util.production.int.aci.ics.psu.edu", scheduler = "torque")
#'
#' # example of waiting for two jobs on slurm cluster
#' wait_for_job(c("24147864", "24147876"), scheduler = "slurm")
#'
#' # example of waiting for two jobs on local machine
#' wait_for_job(c("9843", "9844"), scheduler = "local")
#' }
#'
#' @author Michael Hallquist
#' @importFrom dplyr full_join if_else bind_rows
#' @importFrom utils type.convert
#' @export
wait_for_job <- function(job_ids, repolling_interval = 60, max_wait = 60 * 60 * 24,
                         scheduler = "local", quiet = TRUE, stop_on_timeout = TRUE) {
  checkmate::assert_number(repolling_interval, lower = 0.1, upper = 2e5)
  checkmate::assert_number(max_wait, lower = 1, upper = 1814400) # 21 days
  scheduler <- tolower(scheduler) # ignore case
  checkmate::assert_subset(scheduler, c("torque", "qsub", "slurm", "sbatch", "sh", "local"))

  job_complete <- FALSE
  wait_start <- Sys.time()

  get_job_status <- function() { # use variables in parent environment
    if (scheduler %in% c("slurm", "sbatch")) {
      status <- slurm_job_status(job_ids)
      state <- vapply(status$State, function(x) {
        # sacct can append a cancelling user ("CANCELLED by 123") or a
        # truncation marker. Compare the base state, not that display text.
        x <- sub("[+[:space:]].*$", "", toupper(trimws(x)))
        switch(x,
          "BOOT_FAIL" = "failed",
          "CANCELLED" = "failed",
          "COMPLETED" = "complete",
          "DEADLINE" = "failed",
          "FAILED" = "failed",
          "NODE_FAIL" = "failed",
          "OUT_OF_MEMORY" = "failed",
          "PENDING" = "queued",
          "PREEMPTED" = "failed",
          "RUNNING" = "running",
          "REQUEUED" = "queued",
          "REVOKED" = "failed",
          "SUSPENDED" = "suspended",
          "TIMEOUT" = "failed",
          "MISSING" = "missing",
          "unknown"
        )
      }, character(1))    
    } else if (scheduler %in% c("sh", "local")) {
      status <- local_job_status(job_ids)
      state <- vapply(status$STAT, function(x) {
        switch(x,
          "C" = "complete",
          "I" = "running",
          "R" = "running",
          "S" = "running",
          "T" = "suspended",
          "U" = "running",
          "Z" = "failed",
          "unknown"
        )
      }, character(1))     
    } else if (scheduler %in% c("torque", "qsub")) {
      # QSUB
      status <- torque_job_status(job_ids)
      state <- status$State

      # no need for additional mapping in simple torque results
      # state <- sapply(status$State, function(x) {
      #   switch(x,
      #     "C" = "complete",
      #     "R" = "running",
      #     "Q" = "queued",
      #     "H" = "suspended",
      #     "W" = "suspended", # waiting
      #     stop("Unable to understand job state: ", x)
      #   )
      # })
    } else {
      stop("unknown scheduler: ", scheduler)
    }
    return(state)
  }

  ret_code <- NULL # should be set to TRUE if all jobs complete and FALSE if any job fails 

  while (job_complete == FALSE) {
    status <- get_job_status()

    # update wait time
    wait_total <- as.numeric(difftime(Sys.time(), wait_start, units = "secs"))

    # Debugging
    # cat("Wait so far: ", wait_total, "\n")

    if (any(status == "running")) {
      if (isFALSE(quiet)) {
        cat("Job(s) still running:", paste(job_ids[status == "running"], collapse = ", "), "\n")
      }
    }

    if (any(status == "queued")) {
      if (isFALSE(quiet)) {
        cat("Job(s) still queued:", paste(job_ids[status == "queued"], collapse = ", "), "\n")
      }
    }

    if (any(status == "suspended")) {
      if (isFALSE(quiet)) {
        cat("Job(s) suspended:", paste(job_ids[status == "suspended"], collapse = ", "), "\n")
      }
    }

    if (any(status == "missing")) {
      if (isFALSE(quiet)) {
        cat("Job(s) missing from scheduler response:", paste(job_ids[status == "missing"], collapse = ", "), "\n")
      }
    }

    if (wait_total > max_wait) {
      if (isTRUE(stop_on_timeout)) {
        stop("Maximum wait time: ", max_wait, " exceeded. Stopping execution of parent script because something is wrong.")
      } else {
        return(FALSE)
      }
    } else if (all(status %in% c("failed", "complete"))) {
      job_complete <- TRUE # drop out of this loop
      if (isFALSE(quiet)) {
        cat("All jobs have finished.\n")
      }
      if (any(status == "failed")) {
        cat("The following jobs(s) failed:", paste(job_ids[status == "failed"], collapse = ", "), "\n")
        ret_code <- FALSE
      } else {
        ret_code <- TRUE
      }
    } else {
      Sys.sleep(repolling_interval) # wait and repoll jobs
    }
  }

  return(invisible(ret_code))
}

# calls sacct with a job list
slurm_job_status <- function(job_ids = NULL, user = NULL, sacct_format = "jobid,submit,timelimit,start,end,state") {
  if (!is.null(job_ids)) {
    jstring <- paste("-j", paste(job_ids, collapse = ","))
  } else {
    jstring <- ""
  }

  if (!is.null(user)) {
    ustring <- paste("-u", paste(user, collapse = ","))
  } else {
    ustring <- ""
  }

  # -P specifies a parsable output separated by pipes
  # -X avoids printing subsidiary jobs within each job id
  #cmd <- paste("sacct", jstring, ustring, "-X -P -o", sacct_format)
  cmd <- paste(jstring, ustring, "-X -P -o", sacct_format)
  # cat(cmd, "\n")
  res <- system2("sacct", args = cmd, stdout = TRUE)

  # Create base job ID data frame
  df_base <- data.frame(JobID = job_ids, stringsAsFactors = FALSE)
  df_empty <- df_base %>%
    mutate(
      Submit = NA_character_,
      Timelimit = NA_character_,
      Start = NA_character_,
      End = NA_character_,
      State = "MISSING",
      stringsAsFactors = FALSE
    )

  # handle non-zero exit status -- return empty data
  if (!is.null(attr(res, "status"))) {
    warning("sacct call generated non-zero exit status")
    print(cmd)
    return(df_empty)
  }

  if (length(res) == 1L) {
    # data.table::fread will break down (see Github issue )
    out <- readr::read_delim(I(res), delim = "|", show_col_types = FALSE)
  } else {
    out <- data.table::fread(text = res, data.table=FALSE)
  }

  # Check required fields
  if (!all(c("JobID", "State") %in% names(out))) {
    warning("Missing required columns in sacct output")
    return(df_empty)
  }

  out$JobID <- as.character(out$JobID)

  merged <- merge(df_base, out, by = "JobID", all.x = TRUE)

  if ("State" %in% names(merged)) {
    merged$State[is.na(merged$State)] <- "MISSING"
  } else {
    merged$State <- "MISSING"
  }

  return(merged)
}

# torque does not keep information about completed jobs available in qstat or qselect
# thus, need to log when a job is listed as queued, so that it 'going missing' is evidence of it being completed
torque_job_status <- function(job_ids, user = NULL) {
  #res <- system2("qstat", args = paste("-f", paste(job_ids, collapse=" "), "| grep -i 'job_state'"), stdout = TRUE)

  q_jobs <- system2("qselect", args = "-u $USER -s QW", stdout = TRUE) # queued jobs
  r_jobs <- system2("qselect", args = "-u $USER -s EHRT", stdout = TRUE) # running jobs
  c_jobs <- system2("qselect", args = "-u $USER -s C", stdout = TRUE) # complete jobs
  m_jobs <- setdiff(job_ids, c(q_jobs, r_jobs, c_jobs)) # missing jobs
  #state <- c("queued", "running", "complete", "missing")
  state <- c("queued", "running", "complete", "complete")

  # TORQUE clusters only keep jobs with status C (complete) for a limited period of time. After that, the job comes back as missing.
  # Because of this, if one job finishes at time X and another finishes at time Y, job X will be 'missing' if job Y takes a very long time.
  # Thus, we return any missing jobs as complete, which could be problematic if they are truly missing immediately after submission (as happened with slurm).
  # Ideally, we would track a job within wait_for_job such that it can be missing initially, then move into running, then move into complete.

  j_list <- list(q_jobs, r_jobs, c_jobs, m_jobs)
  state_list <- list()
  for (ii in seq_along(j_list)) {
    if (length(j_list[[ii]]) > 0L) {
      state_list[[state[ii]]] <- data.frame(JobID = j_list[[ii]], State = state[ii])
    }
  }

  state_df <- bind_rows(state_list)

  if (!is.null(attr(q_jobs, "status"))) {
    warning("qselect call generated non-zero exit status")
    return(data.frame(JobID = job_ids, State = "missing"))
  }

  #job_state <- sub(".*job_state = ([A-z]).*", "\\1", res, perl = TRUE)

  return(state_df)
}

local_job_status <- function(job_ids = NULL, user = NULL,
                             ps_format = "user,pid,state,time,etime,%cpu,%mem,comm,xstat") {
  job_ids <- type.convert(job_ids, as.is = T) # convert to integers
  checkmate::assert_integerish(job_ids)

  if (!is.null(job_ids)) {
    jstring <- paste("-p", paste(job_ids, collapse = ","))
  } else {
    jstring <- ""
  }

  if (!is.null(user)) {
    ustring <- paste("-u", paste(user, collapse = ","))
  } else {
    ustring <- ""
  }

  # cat(paste("ps", jstring, ustring, "-o", ps_format), sep = "\n")
  res <- suppressWarnings(system2("ps", args = paste(jstring, ustring, "-o", ps_format), stdout = TRUE)) # intern=TRUE)

  # need to trap res of length 1 (just header row) to avoid data.table bug.
  if (!is.null(attr(res, "status")) && attr(res, "status") != 0) {
    hrow <- strsplit(res, "\\s+")[[1]]
    dt <- data.frame(matrix(NA, nrow = length(job_ids), ncol = length(hrow)))
    names(dt) <- hrow
    dt$PID <- as.integer(job_ids)
  } else {
    stopifnot(length(res) > 1)
    # fread and any other parsing can break down with consecutive spaces in body of output.
    # This happens with lstart and start, avoid these for now.
    # header <- gregexpr("\\b", res[1], perl = T)
    # l2 <- gregexpr("\\b", res[2], perl=T)
    dt <- data.table::fread(text = res)
  }

  # fix difference in column naming between FreeBSD and *nux (make all like FreeBSD)
  data.table::setnames(dt, c("S", "COMMAND"), c("STAT", "COMM"), skip_absent = TRUE)

  # build df that fills in missing jobs (completed/killed)
  all_dt <- data.frame(PID = as.integer(job_ids)) %>%
    dplyr::full_join(dt, by = "PID") %>%
    mutate(STAT = substr(STAT, 1, 1)) # only care about first character of state

  all_dt$STAT[is.na(all_dt$STAT)] <- "C" # complete

  return(all_dt)
}
