# Execution of compiled familywise-error correction plans

fwe_run_schema_version <- function() 1L

refresh_fwe_method <- function(method, plan) {
  method <- as_fwe_method_adapter(method)
  UseMethod("refresh_fwe_method", method)
}

#' @keywords internal
#' @export
refresh_fwe_method.default <- function(method, plan) {
  stop(
    "No FWE refresh adapter is implemented for method: ", method$name,
    call. = FALSE
  )
}

ptfce_missing_inputs <- function(tasks) {
  vapply(seq_len(nrow(tasks)), function(ii) {
    missing <- character(0)
    if (!file.exists(tasks$zstat_file[ii])) missing <- c(missing, "zstat_file")
    if (!file.exists(tasks$correction_mask_file[ii])) {
      missing <- c(missing, "correction_mask_file")
    }
    if ("model_mask_file" %in% names(tasks) && !file.exists(tasks$model_mask_file[ii])) {
      missing <- c(missing, "model_mask_file")
    }
    if (!file.exists(tasks$smoothness_file[ii])) {
      missing <- c(missing, "smoothness_file")
    }
    paste(missing, collapse = ",")
  }, character(1L))
}

fwe_complete_tasks <- function(tasks, artifacts) {
  vapply(seq_len(nrow(tasks)), function(i) {
    fwe_task_is_complete(tasks[i, , drop = FALSE],
      artifacts[artifacts$task_id == tasks$task_id[i], , drop = FALSE])
  }, logical(1L))
}

#' @keywords internal
#' @export
refresh_fwe_method.fwe_method_ptfce <- function(method, plan) {
  plan$analyses$model_mask_exists <- file.exists(plan$analyses$model_mask_file)
  plan$analyses$smoothness_exists <- file.exists(plan$analyses$smoothness_file)
  plan$artifacts$exists <- file.exists(plan$artifacts$artifact_file)
  plan$tasks$missing_inputs <- ptfce_missing_inputs(plan$tasks)
  plan$tasks$input_fingerprint <- fwe_task_fingerprints(plan$tasks)
  plan$tasks$status <- ifelse(
    nzchar(plan$tasks$missing_inputs), "blocked", "ready"
  )
  plan$tasks$status[fwe_complete_tasks(plan$tasks, plan$artifacts)] <- "complete"
  plan
}

#' Refresh the state of a compiled FWE plan
#'
#' Recheck method inputs and expected artifacts without rerunning a correction.
#' This is especially useful after scheduler-submitted jobs have finished.
#'
#' @param plan an \code{fwe_plan} or path to an RDS plan written by
#'   \code{write_fwe_plan()}.
#' @return a validated \code{fwe_plan} with refreshed task and artifact state.
#' @export
refresh_fwe_plan <- function(plan) {
  if (checkmate::test_string(plan) && file.exists(plan)) {
    plan <- read_fwe_plan(plan)
  }
  validate_fwe_plan(plan)
  plan <- refresh_fwe_method(plan$spec$method, plan)
  class(plan) <- c("fwe_plan", "list")
  validate_fwe_plan(plan)
  plan
}

build_fwe_commands <- function(method, plan, task_rows, worker_script = NULL) {
  method <- as_fwe_method_adapter(method)
  UseMethod("build_fwe_commands", method)
}

#' @keywords internal
#' @export
build_fwe_commands.default <- function(method, plan, task_rows,
                                       worker_script = NULL) {
  stop(
    "No FWE execution adapter is implemented for method: ", method$name,
    call. = FALSE
  )
}

resolve_ptfce_worker_script <- function(worker_script = NULL) {
  if (is.null(worker_script)) {
    worker_script <- system.file(
      "bin", "ptfce_zstat.R", package = "fmri.pipeline"
    )
    if (!nzchar(worker_script)) {
      stop(
        "Cannot locate bin/ptfce_zstat.R in the fmri.pipeline installation.",
        call. = FALSE
      )
    }
  }
  checkmate::assert_file_exists(worker_script)
  normalizePath(worker_script, mustWork = TRUE)
}

fwe_shell_command <- function(executable, arguments) {
  tokens <- c(executable, arguments)
  paste(vapply(tokens, shQuote, character(1L), type = "sh"), collapse = " ")
}

#' @keywords internal
#' @export
build_fwe_commands.fwe_method_ptfce <- function(
    method, plan, task_rows, worker_script = NULL) {
  worker_script <- resolve_ptfce_worker_script(worker_script)
  wrapper <- system.file("bin", "fwe_task.R", package = "fmri.pipeline")
  checkmate::assert_file_exists(wrapper)
  rscript <- rscript_executable()
  arguments <- lapply(seq_len(nrow(task_rows)), function(ii) {
    alpha <- task_rows$fwe_alpha[[ii]]
    task_arguments <- c(
      "--zstat", task_rows$zstat_file[ii],
      "--mask", task_rows$correction_mask_file[ii],
      if ("model_mask_file" %in% names(task_rows)) c("--model_mask", task_rows$model_mask_file[ii]),
      "--fsl_smoothest", task_rows$smoothness_file[ii],
      "--out_dir", task_rows$output_directory[ii],
      "--fwep", vapply(alpha, format_fwe_probability, character(1L)),
      if (isTRUE(task_rows$two_sided[ii])) "--twosided" else "--onesided",
      if (isTRUE(task_rows$write_threshold_images[ii])) {
        "--write_thresh_imgs"
      }
    )
    task_arguments <- task_arguments[
      !is.na(task_arguments) & nzchar(task_arguments)
    ]
    c(wrapper, "--request", file.path(task_rows$output_directory[ii], "request.rds"),
      "--worker", worker_script, task_arguments)
  })
  commands <- vapply(
    arguments,
    function(task_arguments) fwe_shell_command(rscript, task_arguments),
    character(1L)
  )

  result <- data.frame(
    task_id = task_rows$task_id,
    executable = rep(rscript, nrow(task_rows)),
    command = unname(commands),
    stringsAsFactors = FALSE
  )
  result$arguments <- I(arguments)
  result
}

normalize_fwe_run_scheduler <- function(scheduler) {
  checkmate::assert_string(scheduler)
  scheduler <- tolower(scheduler)
  scheduler <- switch(scheduler,
    sbatch = "slurm",
    qsub = "torque",
    sh = "local",
    scheduler
  )
  checkmate::assert_subset(
    scheduler, c("local", "slurm", "torque"), empty.ok = FALSE
  )
  scheduler
}

normalize_fwe_run_environment <- function(env_variables) {
  if (is.null(env_variables)) return(character(0))
  checkmate::assert_character(env_variables, names = "unique")
  if (is.null(names(env_variables)) || any(!nzchar(names(env_variables)))) {
    stop("env_variables must be a named character vector.", call. = FALSE)
  }
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*$", names(env_variables)))) {
    stop("env_variables contains an invalid environment variable name.", call. = FALSE)
  }
  specified <- !is.na(env_variables)
  values <- vapply(
    env_variables[specified], shQuote, character(1L), type = "sh"
  )
  paste0(names(env_variables)[specified], "=", values)
}

with_fwe_run_environment <- function(env_variables, code) {
  if (length(env_variables) == 0L) return(force(code))

  variable_names <- names(env_variables)
  previous_values <- Sys.getenv(variable_names, unset = NA_character_)
  names(previous_values) <- variable_names
  on.exit({
    previously_unset <- is.na(previous_values)
    if (any(previously_unset)) {
      Sys.unsetenv(variable_names[previously_unset])
    }
    if (any(!previously_unset)) {
      do.call(
        Sys.setenv,
        as.list(previous_values[!previously_unset])
      )
    }
  }, add = TRUE)

  do.call(Sys.setenv, as.list(env_variables))
  force(code)
}

new_fwe_run <- function(plan, scheduler, dry_run, force, execution,
                        job_ids, started_at, finished_at) {
  job_ids <- as.character(job_ids)
  job_ids <- unique(job_ids[!is.na(job_ids) & nzchar(job_ids)])
  value <- list(
    schema_version = fwe_run_schema_version(),
    plan_name = plan$name,
    method = plan$method,
    scheduler = scheduler,
    dry_run = dry_run,
    force = force,
    execution = execution,
    job_ids = unname(job_ids),
    result = NULL,
    plan = plan,
    started_at = started_at,
    finished_at = finished_at
  )
  class(value) <- c("fwe_run", "list")
  value
}

finalize_fwe_run <- function(run, persist_results) {
  if (isTRUE(persist_results)) {
    result <- collect_fwe_results(run, require_complete = FALSE, persist = TRUE)
    run$plan <- result$plan
    run$result <- result
  }
  run
}

stop_fwe_execution <- function(message, run) {
  condition <- structure(
    list(message = message, call = NULL, run = run),
    class = c("fwe_execution_error", "error", "condition")
  )
  stop(condition)
}

#' Execute a compiled FWE plan
#'
#' Execute ready tasks locally or submit them through Slurm or TORQUE. Completed
#' tasks are skipped unless \code{force = TRUE}. Local execution is synchronous
#' and returns a refreshed plan; scheduler execution returns submitted job IDs
#' and can later be refreshed with \code{refresh_fwe_plan()}.
#'
#' @param plan an \code{fwe_plan} or path to an RDS plan.
#' @param scheduler execution backend. Defaults to the scheduler compiled into
#'   the plan.
#' @param force rerun tasks whose required artifacts already exist.
#' @param dry_run compile commands without creating directories, executing, or
#'   submitting jobs.
#' @param allow_blocked skip blocked tasks instead of rejecting the run.
#' @param worker_script optional method worker override. Primarily useful for
#'   development and testing.
#' @param commands_per_cpu,cpus_per_job,memgb_per_command,time_per_command
#'   resource settings passed to \code{cluster_submit_shell_jobs()} for Slurm
#'   and TORQUE execution.
#' @param sched_args,env_variables,wait_jobs optional scheduler settings passed
#'   to \code{cluster_submit_shell_jobs()}. Named \code{env_variables} are also
#'   applied during local execution.
#' @param tracking_sqlite_db,tracking_args optional SQLite database and metadata
#'   for registering submitted correction jobs in pipeline job tracking.
#' @param stop_on_error raise an \code{fwe_execution_error} if a local command
#'   fails, does not produce all required artifacts, or a scheduler submission
#'   fails. The condition contains the run report in its \code{run} field.
#' @param persist_results if \code{TRUE}, real executions and submissions update
#'   the persistent artifact manifest and RDS result snapshot. Dry runs never
#'   write result files.
#'
#' @return an \code{fwe_run} containing the execution report and refreshed
#'   \code{fwe_plan}.
#' @export
run_fwe_plan <- function(
    plan, scheduler = NULL, force = FALSE, dry_run = FALSE,
    allow_blocked = FALSE, worker_script = NULL,
    commands_per_cpu = NULL, cpus_per_job = NULL,
    memgb_per_command = NULL, time_per_command = NULL,
    sched_args = NULL, env_variables = NULL, wait_jobs = NULL,
    stop_on_error = TRUE, persist_results = TRUE,
    tracking_sqlite_db = NULL, tracking_args = list()) {
  if (checkmate::test_string(plan) && file.exists(plan)) {
    plan <- read_fwe_plan(plan)
  }
  validate_fwe_plan(plan)
  checkmate::assert_logical(force, len = 1L, any.missing = FALSE)
  checkmate::assert_logical(dry_run, len = 1L, any.missing = FALSE)
  checkmate::assert_logical(allow_blocked, len = 1L, any.missing = FALSE)
  checkmate::assert_logical(stop_on_error, len = 1L, any.missing = FALSE)
  checkmate::assert_logical(persist_results, len = 1L, any.missing = FALSE)
  compute <- normalize_fwe_compute(plan$spec$compute)
  commands_per_cpu <- commands_per_cpu %||% compute$commands_per_cpu
  cpus_per_job <- cpus_per_job %||% compute$cpus_per_job
  memgb_per_command <- memgb_per_command %||% compute$memgb_per_command
  time_per_command <- time_per_command %||% compute$time_per_command
  if (is.null(scheduler)) scheduler <- plan$scheduler
  scheduler <- normalize_fwe_run_scheduler(scheduler)
  started_at <- Sys.time()
  plan <- refresh_fwe_plan(plan)

  missing_inputs <- nzchar(plan$tasks$missing_inputs)
  blocked_for_run <- missing_inputs &
    (plan$tasks$status != "complete" | isTRUE(force))
  if (any(blocked_for_run) && !isTRUE(allow_blocked)) {
    blocked <- plan$tasks$task_id[blocked_for_run]
    stop(
      "Cannot execute an FWE plan with blocked tasks: ",
      paste(blocked, collapse = ", "),
      ". Use allow_blocked = TRUE to skip them.",
      call. = FALSE
    )
  }

  execution <- data.frame(
    task_id = plan$tasks$task_id,
    plan_status = plan$tasks$status,
    execution_status = ifelse(
      plan$tasks$status == "complete" & !force, "skipped_complete",
      ifelse(missing_inputs, "skipped_blocked", "pending")
    ),
    command = NA_character_,
    script_file = NA_character_,
    log_file = file.path(
      plan$output_directory, "logs", paste0(plan$tasks$task_id, ".log")
    ),
    exit_status = NA_integer_,
    job_id = NA_character_,
    stringsAsFactors = FALSE
  )
  run_rows <- which(execution$execution_status == "pending")

  if (length(run_rows) > 0L) {
    command_table <- build_fwe_commands(
      plan$spec$method,
      plan,
      plan$tasks[run_rows, , drop = FALSE],
      worker_script = worker_script
    )
    execution$command[run_rows] <- command_table$command[
      match(execution$task_id[run_rows], command_table$task_id)
    ]
  }

  if (isTRUE(dry_run)) {
    execution$execution_status[run_rows] <- "planned"
    return(new_fwe_run(
      plan, scheduler, dry_run, force, execution, character(0),
      started_at, Sys.time()
    ))
  }

  if (length(run_rows) == 0L) {
    run <- new_fwe_run(
      plan, scheduler, dry_run, force, execution, character(0),
      started_at, Sys.time()
    )
    return(finalize_fwe_run(run, persist_results))
  }

  output_directories <- plan$tasks$output_directory[
    match(execution$task_id[run_rows], plan$tasks$task_id)
  ]
  invisible(vapply(output_directories, dir.create, logical(1L),
                   recursive = TRUE, showWarnings = FALSE))
  jobs_directory <- file.path(plan$output_directory, "jobs")
  logs_directory <- file.path(plan$output_directory, "logs")
  dir.create(jobs_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_directory, recursive = TRUE, showWarnings = FALSE)

  for (ii in run_rows) {
    task <- plan$tasks[ii, , drop = FALSE]
    unlink(fwe_completion_file(task))
    request <- list(task = task, artifacts = plan$artifacts[
      plan$artifacts$task_id == task$task_id, , drop = FALSE])
    fwe_atomic_write(file.path(task$output_directory, "request.rds"),
      function(path) saveRDS(request, path))
  }

  if (identical(scheduler, "local")) {
    # Validate names and values before changing the current process.
    normalize_fwe_run_environment(env_variables)
    local_env <- env_variables[!is.na(env_variables)]
    # A child R process must not run the parent R CMD check harness.
    local_env <- c(
      local_env[names(local_env) != "R_TESTS"],
      R_TESTS = ""
    )
    with_fwe_run_environment(local_env, {
      for (ii in run_rows) {
        script_file <- file.path(
          jobs_directory, paste0(execution$task_id[ii], ".bash")
        )
        writeLines(
          c("#!/bin/bash", "set -euo pipefail", execution$command[ii]),
          script_file
        )
        Sys.chmod(script_file, mode = "0755")
        execution$script_file[ii] <- script_file
        command_row <- match(
          execution$task_id[ii], command_table$task_id
        )
        process_arguments <- vapply(
          command_table$arguments[[command_row]],
          shQuote,
          character(1L)
        )
        execution$exit_status[ii] <- suppressWarnings(system2(
          command_table$executable[command_row], process_arguments,
          stdout = execution$log_file[ii],
          stderr = execution$log_file[ii]
        ))
      }
    })
    plan <- refresh_fwe_plan(plan)
    refreshed_status <- plan$tasks$status[
      match(execution$task_id[run_rows], plan$tasks$task_id)
    ]
    execution$execution_status[run_rows] <- ifelse(
      execution$exit_status[run_rows] != 0L, "failed",
      ifelse(refreshed_status == "complete", "completed", "incomplete")
    )
  } else {
    submission_log <- file.path(jobs_directory, "submission.csv")
    job_ids <- cluster_submit_shell_jobs(
      job_list = execution$command[run_rows],
      commands_per_cpu = commands_per_cpu,
      cpus_per_job = cpus_per_job,
      memgb_per_command = memgb_per_command,
      time_per_job = NULL,
      time_per_command = time_per_command,
      fork_jobs = TRUE,
      sched_args = sched_args,
      env_variables = env_variables,
      wait_jobs = wait_jobs,
      scheduler = scheduler,
      job_out_dir = jobs_directory,
      job_script_prefix = paste0("job_", fs::path_sanitize(plan$name)),
      log_file = submission_log,
      debug = FALSE
    )
    submitted <- utils::read.csv(submission_log, stringsAsFactors = FALSE)
    for (job_id in unique(stats::na.omit(submitted$job_id))) {
      if (!nzchar(job_id)) next
      insert_tracked_job(tracking_sqlite_db, as.character(job_id), utils::modifyList(list(
        job_name = paste0("fwe_", plan$name), scheduler = scheduler,
        batch_directory = jobs_directory, n_cpus = cpus_per_job
      ), tracking_args))
    }
    submitted_rows <- match(execution$command[run_rows], submitted$cmd)
    execution$job_id[run_rows] <- submitted$job_id[submitted_rows]
    execution$execution_status[run_rows] <- ifelse(
      is.na(execution$job_id[run_rows]) |
        !nzchar(execution$job_id[run_rows]),
      "submission_failed", "submitted"
    )
  }

  if (!exists("job_ids", inherits = FALSE)) {
    job_ids <- execution$job_id
  }
  result <- new_fwe_run(
    plan, scheduler, dry_run, force, execution, job_ids,
    started_at, Sys.time()
  )
  result <- finalize_fwe_run(result, persist_results)
  failures <- execution$execution_status %in%
    c("failed", "incomplete", "submission_failed")
  if (any(failures) && isTRUE(stop_on_error)) {
    stop_fwe_execution(
      paste0(
        "FWE execution failed or produced incomplete artifacts for task(s): ",
        paste(execution$task_id[failures], collapse = ", ")
      ),
      result
    )
  }
  result
}

#' @keywords internal
#' @export
print.fwe_run <- function(x, ...) {
  statuses <- table(x$execution$execution_status)
  status_text <- paste(
    paste0(names(statuses), "=", as.integer(statuses)), collapse = ", "
  )
  cat(sprintf("<fwe_run> %s\n", x$plan_name))
  cat(sprintf("  method: %s\n", x$method))
  cat(sprintf("  scheduler: %s\n", x$scheduler))
  cat(sprintf("  dry run: %s\n", x$dry_run))
  cat(sprintf("  tasks: %s\n", status_text))
  invisible(x)
}
