#' Submit AFNI GLMs to the configured scheduler
#' @param gpa A glm_pipeline_arguments object.
#' @param level Analysis level (1 or 3). Runs are combined within L1.
#' @param model_names Optional model subset.
#' @param rerun Whether to rerun complete models.
#' @param wait_for Optional parent job IDs.
#' @return Submitted job IDs.
#' @export
run_afni_sepjobs <- function(gpa, level = 1L, model_names = NULL, rerun = FALSE, wait_for = "") {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_choice(level, c(1L, 3L))
  checkmate::assert_character(model_names, null.ok = TRUE)
  checkmate::assert_flag(rerun)
  checkmate::assert_character(wait_for, any.missing = FALSE)
  queue <- gpa[[paste0("l", level, "_model_setup")]]$afni
  if (is.null(queue) || !nrow(queue)) return(NULL)
  model_col <- paste0("l", level, "_model")
  if (!is.null(model_names)) queue <- queue[queue[[model_col]] %in% model_names, , drop = FALSE]
  if (!nrow(queue)) return(NULL)
  jobs <- character(0)
  is_glm <- basename(queue$afni_script) == "run_afni_glm.sh"
  if (any(is_glm & !file.exists(file.path(dirname(queue$afni_script), "expected_outputs.txt")))) {
    stop("Missing AFNI GLM output manifest; rerun model setup.")
  }
  if (any(!is_glm)) {
    if (level != 3L) stop("Missing AFNI GLM output manifest; run setup_l1_models again.")
    mixed_gpa <- gpa
    mixed_gpa$l3_model_setup$afni <- queue[!is_glm, , drop = FALSE]
    jobs <- run_3dlmer_sepjobs(mixed_gpa, rerun = rerun, wait_for = wait_for)
  }
  queue <- queue[is_glm, , drop = FALSE]
  resources <- gpa$parallel$afni
  wall_time <- resources[[paste0("l", level, "_glm_time")]] %||% "8:00:00"
  mem <- resources[[paste0("l", level, "_glm_memgb")]] %||% "12"
  cpus <- resources[[paste0("l", level, "_glm_cpus")]] %||% 1L
  for (i in seq_len(nrow(queue))) {
    row <- queue[i, , drop = FALSE]
    if (!rerun && get_afni_glm_status(dirname(row$afni_script))$afni_complete) next
    jobs <- c(jobs, submit_afni_glm_job(gpa, row$afni_script,
      paste0("afni_l", level, "_", row[[model_col]]), wall_time, mem, cpus, wait_for))
  }
  jobs
}

submit_afni_glm_job <- function(gpa, script, job_name, wall_time, mem, cpus, wait_for) {
  checkmate::assert_file_exists(script)
  checkmate::assert_integerish(cpus, lower = 1L, len = 1L)
  checkmate::assert_choice(gpa$scheduler, c("slurm", "sbatch", "torque", "qsub", "local", "sh"))
  q <- function(x) shQuote(x, type = "sh")
  log_dir <- scheduler_child_log_directory(gpa)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  slurm <- gpa$scheduler %in% c("slurm", "sbatch")
  pbs <- gpa$scheduler %in% c("torque", "qsub")
  token <- basename(tempfile("afni_glm_"))
  batch <- file.path(dirname(script), paste0(token, if (slurm) ".sbatch" else if (pbs) ".pbs" else ".sh"))
  deps <- wait_for[nzchar(wait_for)]
  header <- if (slurm) c("#SBATCH -N 1", "#SBATCH --ntasks=1", paste0("#SBATCH --cpus-per-task=", cpus),
    paste0("#SBATCH --time=", wall_time), paste0("#SBATCH --mem=", mem, "G"),
    paste("#SBATCH -J", scheduler_safe_token(job_name, max_chars = 80L)),
    if (length(deps)) paste0("#SBATCH --dependency=afterok:", paste(deps, collapse = ":"))) else if (pbs) c(
    paste0("#PBS -l nodes=1:ppn=", cpus), paste0("#PBS -l mem=", mem, "gb"),
    paste0("#PBS -l walltime=", wall_time), paste("#PBS -N", substr(scheduler_safe_token(job_name), 1L, 15L)),
    if (length(deps)) paste0("#PBS -W depend=afterok:", paste(deps, collapse = ":"))) else character(0)
  tracking <- gpa$output_locations$sqlite_db
  updater <- system.file("bin/upd_job_status.R", package = "fmri.pipeline")
  track_cmd <- if (!is.null(tracking)) paste("Rscript", q(updater), '--job_id "$job_id" --sqlite_db', q(tracking), '--status "$job_status"') else ":"
  lines <- c("#!/usr/bin/env bash", header,
    if (slurm || pbs) sched_args_to_header(gpa),
    scheduler_output_directives(gpa$scheduler, log_dir, job_name = job_name, extra = token),
    if (slurm) "job_id=$SLURM_JOB_ID" else if (pbs) "job_id=$PBS_JOBID" else "job_id=$$",
    paste0("job_name=", q(job_name)), paste0("batch_directory=", q(log_dir)), paste0("batch_file=", q(batch)),
    paste0("scheduler_name=", if (slurm) "slurm" else if (pbs) "pbs" else "local"),
    scheduler_runtime_log_assignment(gpa$scheduler, log_dir, job_name = job_name, extra = token),
    job_manifest_shell_function(),
    # The trap preserves estimation failure even when bookkeeping succeeds.
    "finish_afni_job() {", "  exit_code=$?", "  trap - EXIT",
    '  if [ "$exit_code" -eq 0 ]; then job_status=COMPLETED; else job_status=FAILED; fi',
    paste("  write_job_manifest", q(dirname(script)), "afni_glm", '"$job_status"', q(script)),
    paste(" ", track_cmd), '  exit "$exit_code"', "}", "trap finish_afni_job EXIT",
    "set -e", get_compute_environment(gpa, c("afni", "r")), "job_status=STARTED", track_cmd,
    paste("bash", q(script)))
  writeLines(lines, batch)
  cluster_job_submit(batch, scheduler = gpa$scheduler, fail_on_error = TRUE,
    wait_jobs = if (length(deps) && !slurm && !pbs) deps else NULL,
    tracking_sqlite_db = tracking, tracking_args = if (!is.null(tracking)) list(
      job_name = job_name, batch_directory = dirname(script), n_nodes = 1L, n_cpus = cpus,
      wall_time = wall_time, mem_total = mem, scheduler_options = gpa$parallel$sched_args,
      parent_job_id = deps) else list())
}
