#' Submit SPM jobs to the cluster for unattended MATLAB/Octave execution
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing model specification
#' @param level level of analysis (1 or 3)
#' @param model_names optional model names used to subset SPM runs
#' @param rerun logical indicating whether to rerun existing directories
#' @param wait_for optional parent job ids that should complete before these jobs commence
#' @return a vector of job ids for all scripts that were submitted
#' @importFrom glue glue
#' @export
run_spm_sepjobs <- function(gpa, level = 1L, model_names = NULL, rerun = FALSE, wait_for = "") {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1L, upper = 3L, len = 1L)
  checkmate::assert_character(model_names, null.ok = TRUE)
  checkmate::assert_logical(rerun)

  lg <- lgr::get_logger(paste0("glm_pipeline/spm_l", level, "_estimation"))
  lg$set_threshold(gpa$lgr_threshold)
  upd_job_status_path <- system.file("bin/upd_job_status.R", package = "fmri.pipeline")

  if (level == 1L) {
    if (!checkmate::test_class(gpa$l1_model_setup, "l1_setup")) {
      lg$error("In run_spm_sepjobs, did not find an l1_setup object in gpa$l1_model_setup.")
      lg$error("Make sure to run setup_l1_models before running run_spm_sepjobs")
      return(NULL)
    }

    spm_queue <- gpa$l1_model_setup$spm
    if (!is.null(model_names[1L])) {
      spm_queue <- spm_queue %>% dplyr::filter(l1_model %in% !!model_names)
    }

    spm_time <- gpa$parallel$spm$l1_spm_time
    spm_memgb <- gpa$parallel$spm$l1_spm_memgb
    spm_cpus <- gpa$parallel$spm$l1_spm_cpus_per_job
    runsperproc <- gpa$parallel$spm$l1_spm_runs_per_cpu

    run_setup <- isTRUE(gpa$glm_settings$spm$run_l1_setup)
    run_glm <- isTRUE(gpa$glm_settings$spm$run_l1_glm)
    run_contrasts <- isTRUE(gpa$glm_settings$spm$run_l1_contrasts)
    cleanup_tmp <- isTRUE(gpa$glm_settings$spm$cleanup_tmp)

    script_map <- c(
      setup = "setup_glm_design.m",
      glm = "run_glm.m",
      contrasts = "estimate_glm_contrasts.m"
    )

    spm_output_directory <- file.path(gpa$output_locations$scheduler_scripts, "spm_l1")
  } else if (level == 3L) {
    if (!checkmate::test_class(gpa$l3_model_setup, "l3_setup")) {
      lg$error("In run_spm_sepjobs, did not find a l3_setup object in gpa$l3_model_setup.")
      lg$error("Make sure to run setup_l3_models before running run_spm_sepjobs")
      return(NULL)
    }

    spm_queue <- gpa$l3_model_setup$spm
    if (!is.null(model_names[1L])) {
      spm_queue <- spm_queue %>% dplyr::filter(l3_model %in% !!model_names)
    }

    spm_time <- gpa$parallel$spm$l3_spm_time
    spm_memgb <- gpa$parallel$spm$l3_spm_memgb
    spm_cpus <- gpa$parallel$spm$l3_spm_cpus_per_job
    runsperproc <- gpa$parallel$spm$l3_spm_runs_per_cpu

    run_setup <- isTRUE(gpa$glm_settings$spm$run_l3_setup)
    run_glm <- isTRUE(gpa$glm_settings$spm$run_l3_glm)
    run_contrasts <- isTRUE(gpa$glm_settings$spm$run_l3_contrasts)
    cleanup_tmp <- isTRUE(gpa$glm_settings$spm$cleanup_tmp)

    script_map <- c(
      setup = "setup_l3_design.m",
      glm = "run_l3_glm.m",
      contrasts = "estimate_l3_contrasts.m"
    )

    spm_output_directory <- file.path(gpa$output_locations$scheduler_scripts, "spm_l3")
  }

  if (is.null(spm_queue) || nrow(spm_queue) == 0L) {
    lg$warn("No SPM jobs found for level %d.", level)
    return(NULL)
  }

  if (!"spm_dir" %in% names(spm_queue)) {
    stop("SPM queue is missing spm_dir column. Ensure setup_l1_models/setup_l3_models ran successfully.")
  }

  # Re-read disk state: the queue may predate a failure or a removed output.
  spm_job_df <- dplyr::bind_rows(lapply(unique(spm_queue$spm_dir), get_spm_status, lg = lg))
  spm_job_df$spm_dir <- unique(spm_queue$spm_dir)
  spm_job_df$to_run <- !spm_job_df$spm_complete

  if (isTRUE(rerun)) {
    lg$info("rerun = TRUE in run_spm_sepjobs. All SPM directories will be marked for job execution.")
    spm_job_df$to_run <- TRUE
  }

  to_run_dirs <- spm_job_df %>%
    dplyr::filter(to_run == TRUE) %>%
    dplyr::pull("spm_dir")

  if (length(to_run_dirs) == 0L) {
    lg$warn("No SPM level %d directories to execute.", level)
    return(NULL)
  }

  current <- vapply(to_run_dirs, spm_setup_current, logical(1))
  if (any(!current)) stop("SPM setup is missing or stale; rerun setup_l", level,
                         "_models before submitting: ", paste(to_run_dirs[!current], collapse = ", "))

  if (!file.exists(spm_output_directory)) {
    lg$debug("Creating spm_l%d working directory: %s", level, spm_output_directory)
    dir.create(spm_output_directory, recursive = TRUE)
  }

  if (is.null(gpa$parallel$compute_environment$spm) || length(gpa$parallel$compute_environment$spm) == 0L) {
    lg$warn("No SPM compute environment configured. MATLAB/Octave may not be available on compute nodes.")
  }

  matlab_cmd <- gpa$glm_settings$spm$matlab_cmd
  if (is.null(matlab_cmd) || !nzchar(matlab_cmd)) matlab_cmd <- "matlab"
  matlab_args <- gpa$glm_settings$spm$matlab_args
  if (is.null(matlab_args) || !nzchar(matlab_args)) matlab_args <- "-batch"

  child_log_directory <- Sys.getenv("SLURM_SUBMIT_DIR", unset = "")
  if (!nzchar(child_log_directory)) child_log_directory <- Sys.getenv("PBS_O_WORKDIR", unset = "")
  if (!nzchar(child_log_directory)) child_log_directory <- getwd()

  # Scheduler header
  if (gpa$scheduler %in% c("slurm", "sbatch")) {
    file_suffix <- ".sbatch"
    preamble <- function(job_name, extra, batch_file) {
      c(
        "#!/bin/bash",
        "#SBATCH -N 1",
        paste0("#SBATCH -n ", spm_cpus),
        paste0("#SBATCH --time=", spm_time),
        paste0("#SBATCH --mem-per-cpu=", spm_memgb, "G"),
        ifelse(any(nzchar(wait_for)), paste0("#SBATCH --dependency=afterok:", paste(wait_for, collapse=":")), ""),
        sched_args_to_header(gpa),
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
        get_compute_environment(gpa, c("spm", "r")),
        "",
        "cd $SLURM_SUBMIT_DIR"
      )
    }
  } else if (gpa$scheduler %in% c("torque", "qsub")) {
    file_suffix <- ".pbs"
    preamble <- function(job_name, extra, batch_file) {
      c(
        "#!/bin/bash",
        paste0("#PBS -l nodes=1:ppn=", spm_cpus),
        paste0("#PBS -l pmem=", spm_memgb, "gb"),
        ifelse(any(nzchar(wait_for)), paste0("#PBS -W depend=afterok:", paste(wait_for, collapse=":")), ""),
        paste0("#PBS -l walltime=", spm_time),
        sched_args_to_header(gpa),
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
        get_compute_environment(gpa, c("spm", "r")),
        "",
        "cd $PBS_O_WORKDIR"
      )
    }
  } else {
    file_suffix <- ".sh"
    preamble <- function(job_name, extra, batch_file) {
      c(
        "#!/bin/bash",
        get_compute_environment(gpa, c("spm", "r")),
        "job_id=$$",
        paste0("job_name=", shQuote(job_name)),
        "scheduler_name=local",
        paste0("batch_directory=", shQuote(child_log_directory)),
        paste0("batch_file=", shQuote(batch_file)),
        scheduler_runtime_log_assignment(gpa$scheduler, child_log_directory, job_name = job_name, extra = extra)
      )
    }
  }

  njobs <- ceiling(length(to_run_dirs) / (spm_cpus * runsperproc))
  df <- data.frame(
    spm_dir = to_run_dirs,
    job = rep(1:njobs, each = spm_cpus * runsperproc, length.out = length(to_run_dirs)),
    stringsAsFactors = FALSE
  )

  submission_id <- basename(tempfile(pattern = "job"))
  joblist <- rep(NA_character_, njobs)

  tracking_sqlite_db <- gpa$output_locations$sqlite_db
  for (j in seq_len(njobs)) {
    outfile <- file.path(spm_output_directory, paste0("spmsep_l", level, "_", j, "_", submission_id, file_suffix))
    job_name <- paste0("spmsep_l", level, "_", j)
    cat(preamble(job_name, submission_id, outfile), file = outfile, sep = "\n")

    thisrun <- with(df, spm_dir[job == j])
    cat(
      "",
      job_manifest_shell_function(),
      "",
      "job_failed=0",
      paste0("matlab_cmd=", shQuote(matlab_cmd)),
      paste0("matlab_args=", shQuote(matlab_args)),
      paste0("spm_setup_script=", shQuote(script_map[["setup"]])),
      paste0("spm_glm_script=", shQuote(script_map[["glm"]])),
      paste0("spm_contrast_script=", shQuote(script_map[["contrasts"]])),
      "",
      paste0("cleanup_tmp=", ifelse(cleanup_tmp, "1", "0")),
      paste0("run_setup=", as.integer(run_setup)),
      paste0("run_glm=", as.integer(run_glm)),
      paste0("run_contrasts=", as.integer(run_contrasts)),
      paste0("spm_level=", level),
      readLines(system.file("scripts", "spm_runner.sh", package = "fmri.pipeline")),
      "",
      "function spm_killed() {",
      "  kill_time=$( date )",
      "  if [ -n \"${odir}\" ]; then",
      "    echo $kill_time > \"${odir}/.spm_fail\"",
      paste0("    write_job_manifest \"${odir}\" \"spm_l", level, "\" \"FAILED\" \"${odir}\""),
      "  fi",
      paste("  Rscript", upd_job_status_path, "--job_id", "\"$job_id\"", "--sqlite_db", tracking_sqlite_db, "--status", "FAILED"),
      "  exit 1",
      "}",
      "trap spm_killed SIGTERM",
      paste("Rscript", upd_job_status_path, "--job_id", "\"$job_id\"", "--sqlite_db", tracking_sqlite_db, "--status", "STARTED"),
      "",
      sep = "\n",
      file = outfile,
      append = TRUE
    )

    for (spm_dir in thisrun) {
      if (!dir.exists(spm_dir)) {
        lg$warn("Skipping missing spm_dir: %s", spm_dir)
        next
      }
      cat(
        paste0("spm_runner ", shQuote(spm_dir)),
        "if [ $? -ne 0 ]; then job_failed=1; fi",
        "",
        sep = "\n",
        file = outfile,
        append = TRUE
      )
    }

    cat(
      "if [ $job_failed -ne 0 ]; then",
      "  job_status=\"FAILED\"",
      "else",
      "  job_status=\"COMPLETED\"",
      "fi",
      paste("Rscript", upd_job_status_path, "--job_id", "\"$job_id\"", "--sqlite_db", tracking_sqlite_db, "--status", "\"${job_status}\""),
      "if [ \"${job_status}\" = \"FAILED\" ]; then",
      "  exit 1",
      "fi",
      "",
      sep = "\n",
      file = outfile,
      append = TRUE
    )

    tracking_args <- list(
      job_name = job_name,
      batch_directory = spm_output_directory,
      n_nodes = 1,
      n_cpus = spm_cpus,
      wall_time = spm_time,
      mem_per_cpu = spm_memgb,
      scheduler_options = gpa$parallel$sched_args
    )

    if (any(nzchar(wait_for))) tracking_args$parent_job_id <- wait_for

    joblist[j] <- cluster_job_submit(
      outfile,
      scheduler = gpa$scheduler,
      tracking_sqlite_db = tracking_sqlite_db,
      tracking_args = tracking_args,
      fail_on_error = TRUE
    )
  }

  writeLines(joblist, con = file.path(spm_output_directory, paste0("spmsep_l", level, "_jobs.txt")))

  return(joblist)
}

#' Combine SPM L3 outputs into AFNI BRIK+HEAD files for visualization
#'
#' @param gpa a glm_pipeline_arguments object having a populated l3_model_setup field.
#' @param spm_l3_combined_filename a glue expression for the path and filename prefix.
#' @param spm_l3_combined_briknames a glue expression for naming the subbriks in the AFNI output.
#' @param template_brain an optional filename for the MNI template that should be used as an underlay in AFNI.
#' @details This function mirrors the logic of combine_feat_l3_to_afni, but for SPM L3 outputs. It detects
#'   contrast maps (con_*.nii) and statistic maps (spmT_*.nii / spmF_*.nii) and concatenates them into AFNI
#'   BRIK/HEAD files with informative labels.
#' @return A data.frame describing the images that were combined.
#' @export
#' @importFrom glue glue_data
#' @importFrom tidyr unnest
combine_spm_l3_to_afni <- function(gpa, spm_l3_combined_filename=NULL, spm_l3_combined_briknames=NULL, template_brain=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (!is.null(template_brain)) {
    checkmate::assert_file_exists(template_brain)
    template_ext <- sub(".*(\\.nii(\\.gz)*)", "\\1", template_brain, perl=TRUE)
  }

  if (is.null(spm_l3_combined_filename) && !is.null(gpa$output_locations$spm_l3_combined_filename)) {
    spm_l3_combined_filename <- gpa$output_locations$spm_l3_combined_filename
  }

  if (is.null(spm_l3_combined_briknames) && !is.null(gpa$output_locations$spm_l3_combined_briknames)) {
    spm_l3_combined_briknames <- gpa$output_locations$spm_l3_combined_briknames
  }

  if (is.null(spm_l3_combined_filename) || is.null(spm_l3_combined_briknames)) {
    stop("spm_l3_combined_filename and spm_l3_combined_briknames must be provided (or present in gpa$output_locations).")
  }

  checkmate::assert_string(spm_l3_combined_filename)
  checkmate::assert_string(spm_l3_combined_briknames)

  if (is.null(gpa$l3_model_setup$spm) || nrow(gpa$l3_model_setup$spm) == 0L) {
    stop("No SPM L3 model setup found in gpa$l3_model_setup$spm.")
  }

  lg <- lgr::get_logger("glm_pipeline/combine_spm_l3_to_afni")
  lg$set_threshold(gpa$lgr_threshold)

  meta_df <- gpa$l3_model_setup$spm %>%
    dplyr::select("l1_model", "l1_cope_name", "l3_model", "spm_dir")

  meta_df$contrast_df <- lapply(meta_df$spm_dir, spm_l3_collect_contrasts, lg = lg)
  meta_df <- meta_df %>% unnest(cols="contrast_df", keep_empty = FALSE)

  if (nrow(meta_df) == 0L) {
    lg$warn("No SPM L3 contrasts found for AFNI combination.")
    return(meta_df)
  }

  meta_df <- meta_df %>%
    dplyr::mutate(
      afni_out = glue_data(., !!spm_l3_combined_filename),
      afni_briks = glue_data(., !!spm_l3_combined_briknames)
    )

  meta_split <- split(meta_df, meta_df$afni_out)

  lapply(meta_split, function(ss) {
    # Mirror FEAT behavior by forcing +tlrc outputs for group-level maps.
    afni_out <- sub("(\\+tlrc|\\+orig)*$", "+tlrc", ss$afni_out[1])
    afni_dir <- dirname(afni_out)
    if (!dir.exists(afni_dir)) {
      dir.create(afni_dir, recursive = TRUE)
      if (!is.null(template_brain)) file.symlink(template_brain, file.path(afni_dir, paste0("template_brain", template_ext)))
    }

    # build long-form list of images (contrast then statistic)
    long_df <- rbind(
      data.frame(
        image_type = "con",
        stat = ss$stat,
        df1 = ss$df1,
        df2 = ss$df2,
        afni_briks = ss$afni_briks,
        nii_file = ss$con_file,
        stringsAsFactors = FALSE
      ),
      data.frame(
        image_type = ifelse(ss$stat == "F", "spmF", "spmT"),
        stat = ss$stat,
        df1 = ss$df1,
        df2 = ss$df2,
        afni_briks = ss$afni_briks,
        nii_file = ss$stat_file,
        stringsAsFactors = FALSE
      )
    )

    long_df <- long_df %>%
      dplyr::filter(!is.na(.data$nii_file) & file.exists(.data$nii_file))

    if (nrow(long_df) == 0L) {
      lg$warn("No images found for AFNI combination in: %s", afni_out)
      return(NULL)
    }

    tcatcall <- paste("3dTcat -overwrite -prefix", afni_out, paste(long_df$nii_file, collapse = " "))
    tcat_ret <- run_afni_command(tcatcall)
    tcat_status <- suppressWarnings(as.integer(tcat_ret[1]))
    if (is.na(tcat_status) || tcat_status != 0L) {
      msg <- sprintf("AFNI command failed with exit status %s: %s", as.character(tcat_ret[1]), tcatcall)
      lg$error(msg)
      stop(msg)
    }

    # tack on statistic type as suffix to sub-brik name (avoid ambiguity)
    long_df$afni_briks <- paste(long_df$afni_briks, long_df$image_type, sep = "_")

    refit_parts <- c("3drefit -fbuc")
    stat_idx <- which(long_df$image_type %in% c("spmT", "spmF")) - 1
    if (length(stat_idx) > 0L) {
      for (ii in seq_along(stat_idx)) {
        row <- long_df[stat_idx[ii] + 1, , drop = FALSE]
        if (row$stat == "T" && is.finite(row$df1)) {
          refit_parts <- c(refit_parts, paste("-substatpar", stat_idx[ii], "fitt", row$df1))
        } else if (row$stat == "F" && is.finite(row$df1) && is.finite(row$df2)) {
          refit_parts <- c(refit_parts, paste("-substatpar", stat_idx[ii], "fift", row$df1, row$df2))
        }
      }
    }

    refitcall <- paste0(
      paste(refit_parts, collapse = " "),
      " -relabel_all_str '", paste(long_df$afni_briks, collapse = " "), "' ", afni_out
    )
    refit_ret <- run_afni_command(refitcall)
    refit_status <- suppressWarnings(as.integer(refit_ret[1]))
    if (is.na(refit_status) || refit_status != 0L) {
      msg <- sprintf("AFNI command failed with exit status %s: %s", as.character(refit_ret[1]), refitcall)
      lg$error(msg)
      stop(msg)
    }
  })

  return(meta_df)
}

spm_l3_collect_contrasts <- function(spm_dir, lg = NULL) {
  checkmate::assert_string(spm_dir)
  if (is.null(lg)) lg <- lgr::get_logger()

  xcon_df <- spm_l3_read_xcon(spm_dir, lg = lg)
  if (is.null(xcon_df) || nrow(xcon_df) == 0L) {
    xcon_df <- spm_l3_parse_contrast_script(spm_dir, lg = lg)
  }
  if (is.null(xcon_df) || nrow(xcon_df) == 0L) {
    xcon_df <- spm_l3_contrasts_from_files(spm_dir, lg = lg)
  }

  if (is.null(xcon_df) || nrow(xcon_df) == 0L) {
    return(data.frame())
  }

  # ensure files are set
  xcon_df$con_file <- ifelse(
    !is.na(xcon_df$con_file), xcon_df$con_file,
    sapply(xcon_df$l3_cope_number, function(ii) spm_find_file(spm_dir, sprintf("con_%04d.nii", ii)))
  )

  xcon_df$stat_file <- ifelse(
    !is.na(xcon_df$stat_file), xcon_df$stat_file,
    sapply(seq_len(nrow(xcon_df)), function(i) {
      stat <- xcon_df$stat[i]
      base <- ifelse(stat == "F", sprintf("spmF_%04d.nii", xcon_df$l3_cope_number[i]),
                     sprintf("spmT_%04d.nii", xcon_df$l3_cope_number[i]))
      spm_find_file(spm_dir, base)
    })
  )

  return(xcon_df)
}

spm_l3_read_xcon <- function(spm_dir, lg = NULL) {
  spm_mat <- file.path(spm_dir, "SPM.mat")
  if (!file.exists(spm_mat)) return(NULL)
  if (is.null(lg)) lg <- lgr::get_logger()

  res <- tryCatch({
    mat <- R.matlab::readMat(spm_mat, fixNames = FALSE)
    spm <- mat$SPM
    if (is.null(spm) || is.null(spm$xCon)) return(NULL)
    spm_xcon_to_df(spm$xCon, spm_dir)
  }, error = function(e) {
    lg$debug("Failed to read SPM.mat contrasts in %s: %s", spm_dir, as.character(e))
    return(NULL)
  })

  return(res)
}

spm_xcon_to_df <- function(xcon, spm_dir) {
  if (is.null(xcon)) return(NULL)

  # case 1: struct fields are named lists
  if (is.list(xcon) && !is.null(names(xcon)) && all(c("name", "STAT") %in% names(xcon))) {
    n <- length(xcon$name)
    out <- lapply(seq_len(n), function(i) {
      name <- spm_matlab_to_string(spm_get_index(xcon$name, i))
      stat <- spm_matlab_to_string(spm_get_index(xcon$STAT, i))
      df <- spm_get_index(xcon$df, i)
      vcon <- spm_get_struct_field(xcon$Vcon, i, "fname")
      vspm <- spm_get_struct_field(xcon$Vspm, i, "fname")
      data.frame(
        l3_cope_number = i,
        l3_cope_name = ifelse(is.na(name), paste0("contrast_", i), name),
        stat = ifelse(is.na(stat) || stat == "", "T", stat),
        df1 = ifelse(length(df) >= 1, as.numeric(df[1]), NA_real_),
        df2 = ifelse(length(df) >= 2, as.numeric(df[2]), NA_real_),
        con_file = ifelse(is.na(vcon), NA_character_, spm_find_file(spm_dir, vcon)),
        stat_file = ifelse(is.na(vspm), NA_character_, spm_find_file(spm_dir, vspm)),
        stringsAsFactors = FALSE
      )
    })
    return(dplyr::bind_rows(out))
  }

  # case 2: list of structs
  if (is.list(xcon) && length(xcon) > 0L && is.list(xcon[[1]])) {
    out <- lapply(seq_along(xcon), function(i) {
      xi <- xcon[[i]]
      name <- spm_matlab_to_string(xi$name)
      stat <- spm_matlab_to_string(xi$STAT)
      df <- xi$df
      vcon <- spm_matlab_to_string(spm_get_struct_field(xi$Vcon, 1, "fname"))
      vspm <- spm_matlab_to_string(spm_get_struct_field(xi$Vspm, 1, "fname"))
      data.frame(
        l3_cope_number = i,
        l3_cope_name = ifelse(is.na(name), paste0("contrast_", i), name),
        stat = ifelse(is.na(stat) || stat == "", "T", stat),
        df1 = ifelse(length(df) >= 1, as.numeric(df[1]), NA_real_),
        df2 = ifelse(length(df) >= 2, as.numeric(df[2]), NA_real_),
        con_file = ifelse(is.na(vcon), NA_character_, spm_find_file(spm_dir, vcon)),
        stat_file = ifelse(is.na(vspm), NA_character_, spm_find_file(spm_dir, vspm)),
        stringsAsFactors = FALSE
      )
    })
    return(dplyr::bind_rows(out))
  }

  return(NULL)
}

spm_l3_parse_contrast_script <- function(spm_dir, lg = NULL) {
  if (is.null(lg)) lg <- lgr::get_logger()
  script <- file.path(spm_dir, "estimate_l3_contrasts.m")
  if (!file.exists(script)) return(NULL)

  lines <- readLines(script, warn = FALSE)
  rx <- "consess\\{(\\d+)\\}\\.(tcon|fcon)\\.name\\s*=\\s*'([^']*)'"
  m <- regmatches(lines, gregexpr(rx, lines, perl = TRUE))
  hits <- unlist(m)
  if (length(hits) == 0L) return(NULL)

  parse_hit <- function(x) {
    parts <- sub(rx, "\\1|\\2|\\3", x, perl = TRUE)
    parts <- strsplit(parts, "\\|", fixed = FALSE)[[1]]
    data.frame(
      l3_cope_number = as.integer(parts[1]),
      l3_cope_name = parts[3],
      stat = ifelse(parts[2] == "fcon", "F", "T"),
      df1 = NA_real_,
      df2 = NA_real_,
      con_file = NA_character_,
      stat_file = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  out <- dplyr::bind_rows(lapply(hits, parse_hit))
  out <- out[order(out$l3_cope_number), , drop = FALSE]
  return(out)
}

spm_l3_contrasts_from_files <- function(spm_dir, lg = NULL) {
  if (is.null(lg)) lg <- lgr::get_logger()
  con_files <- list.files(spm_dir, pattern = "^con_\\d+\\.nii(\\.gz)?$", full.names = TRUE)
  spmT_files <- list.files(spm_dir, pattern = "^spmT_\\d+\\.nii(\\.gz)?$", full.names = TRUE)
  spmF_files <- list.files(spm_dir, pattern = "^spmF_\\d+\\.nii(\\.gz)?$", full.names = TRUE)

  get_num <- function(x, prefix) as.integer(sub(paste0("^", prefix, "_(\\d+).*"), "\\1", basename(x)))
  nums <- unique(c(get_num(con_files, "con"), get_num(spmT_files, "spmT"), get_num(spmF_files, "spmF")))
  nums <- nums[!is.na(nums)]
  if (length(nums) == 0L) return(NULL)
  nums <- sort(nums)

  out <- lapply(nums, function(ii) {
    stat <- if (any(grepl(sprintf("^spmF_%04d", ii), basename(spmF_files)))) "F" else "T"
    data.frame(
      l3_cope_number = ii,
      l3_cope_name = paste0("contrast_", ii),
      stat = stat,
      df1 = NA_real_,
      df2 = NA_real_,
      con_file = NA_character_,
      stat_file = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(out)
}

spm_find_file <- function(spm_dir, fname) {
  if (is.null(fname) || !nzchar(fname)) return(NA_character_)
  cand <- file.path(spm_dir, fname)
  if (file.exists(cand)) return(cand)
  cand_gz <- paste0(cand, ".gz")
  if (file.exists(cand_gz)) return(cand_gz)
  cand_gz2 <- sub("\\.nii$", ".nii.gz", cand)
  if (file.exists(cand_gz2)) return(cand_gz2)
  return(cand)
}

spm_get_index <- function(x, i) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) {
    if (length(x) >= i) return(x[[i]])
  }
  if (is.matrix(x)) {
    if (nrow(x) >= i) return(x[i, , drop = TRUE])
  }
  if (length(x) >= i) return(x[i])
  return(NULL)
}

spm_get_struct_field <- function(x, i, field) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.null(names(x)) && field %in% names(x)) {
    return(spm_get_index(x[[field]], i))
  }
  if (is.list(x) && length(x) >= i) {
    xi <- x[[i]]
    if (is.list(xi) && field %in% names(xi)) return(xi[[field]])
  }
  return(NULL)
}

spm_matlab_to_string <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.character(x)) {
    if (is.matrix(x)) return(paste(x, collapse = ""))
    if (length(x) == 1L) return(x)
    if (all(nchar(x) == 1)) return(paste(x, collapse = ""))
    return(x[1])
  }
  if (is.numeric(x)) {
    x <- x[!is.na(x) & x != 0]
    if (length(x) == 0L) return(NA_character_)
    return(paste(intToUtf8(x), collapse = ""))
  }
  if (is.list(x) && length(x) == 1L) return(spm_matlab_to_string(x[[1]]))
  return(as.character(x[1]))
}

#' This is a wrapper around the spm_extract_anatomical_rois.m script in the inst directory
#'
#' @param l1spmdirs character vector of level 1 SPM directories containing SPM.mat files
#' @param masks character vector of NIfTI mask images for each anatomical ROI of interest
#' @param threshold p-value threshold applied to contrast within mask before extraction
#' @param threshdesc multiple comparisons correction on p-value. 'none' or 'FWE'
#' @param session which session (run) to use for extracting time series
#' @param extent exclude clusters having fewer than voxels than extent
#' @param adjust_F_index index of F-test in SPM.mat to adjust for all effects of interest
#' @param contrast_index index of t-test contrast in SPM.mat that is of interest
#' @param ncores number of cores to use in a parallel approach; function parallelizes over \code{l1spmdirs}
#' @param spm_path path to spm12 installation; added to MATLAB path at runtime
#' @param matlab_path location of MATLAB binary; used with matlabr for run_matlab_code()
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach %dopar%
#' @importFrom iterators iter
#' @importFrom parallel makeCluster stopCluster
#' @importFrom matlabr have_matlab run_matlab_code
#'
#' @author Michael Hallquist
#' 
#' @export
#' 
spm_extract_anatomical_rois <- function(l1spmdirs, masks, threshold=0.2, threshdesc='none', session=1, extent=0,
                                        adjust_F_index=1, contrast_index=NULL, ncores=1,
                                        spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12",
                                        matlab_path="/opt/aci/sw/matlab/R2017b/bin") {
  
  if (missing(l1spmdirs)) { stop("Need to pass in a character vector of all level 1 directories containing SPM.mat files.") }
  stopifnot(all(dir.exists(l1spmdirs)))
  stopifnot(all(file.exists(masks)))
  if (is.null(names(masks))) {
    names(masks) <- sub("(\\.hdr|\\.nii|\\.img)(\\.gz)*$", "", basename(masks), perl=TRUE) #develop basic naming scheme
  }

  stopifnot(dir.exists(spm_path))
  stopifnot(dir.exists(matlab_path))

  #set the matlab path for matlabr
  options(matlab.path=matlab_path)

  if (!have_matlab()) { stop("Unable to find MATLAB installation") }
  
  if (ncores > 1) {
    cl <- makeCluster(ncores)
    on.exit(try(stopCluster(cl)))
    registerDoParallel(cl)
  } else {
    registerDoSEQ()
  }

  nifti_tmpdir <- tempdir()
  gzipped <- grepl(".nii.gz$", masks)
  dir.create(nifti_tmpdir, showWarnings=FALSE)
  masks[gzipped] <- sapply(masks[gzipped], function(x) {
    tmpout <- tempfile(fileext=".nii", tmpdir=nifti_tmpdir)
    system(paste0("gunzip -c ", x, " > ", tmpout))
    return(tmpout)
  })

  #add single quotes around mask strings if not present
  mask_names <- sub("^'?([^']+)'?", "'\\1'", names(masks), perl=TRUE)
  masks <- sub("^'?([^']+)'?", "'\\1'", masks, perl=TRUE)

  #TODO: support multi-session data

  spm_preamble <- c(
    ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
    "spm('defaults', 'fmri');",
    "spm_jobman('initcfg');",
    ""    
  )

  #contains matlab scripts for extraction
  matlab_scripts <- system.file("matlab", package = "dependlab")
  stopifnot(dir.exists(matlab_scripts))

  res <- foreach(dd=iter(l1spmdirs), .packages="matlabr") %dopar% {
    #set the matlab path for matlabr in each worker
    options(matlab.path=matlab_path)

    m_string <- c(spm_preamble,
      paste0("addpath('", matlab_scripts, "');"),
      "cfg = struct();",
      "masks = {",
      masks,
      "};",
      "names = {",
      mask_names,
      "};",
      "",
      paste0("threshold = ", threshold, ";"),
      paste0("threshdesc = ", sub("^'?([^']+)'?", "'\\1'", threshdesc, perl=TRUE), ";"),
      paste0("session = ", session, ";"),
      paste0("extent = ", extent, ";"),
      paste0("adjust_F_index = ", ifelse(is.null(adjust_F_index), "[]", adjust_F_index), "; %adjust time series for all effects of interest"),
      paste0("contrast_index = ", ifelse(is.null(contrast_index), "[]", contrast_index), "; %the contrast of interest"),
      "",
      "for jj = 1 : numel(masks)",
      paste0("  cfg(jj).target_dir = fullfile('", dd, "');"),
      "  cfg(jj).mask = masks{jj};",
      "  cfg(jj).adjust = adjust_F_index;",
      "  cfg(jj).session = session;",
      "  cfg(jj).name = names{jj};"
    )

    if (!is.null(contrast_index)) {
      m_string <- c(m_string,
      "  cfg(jj).contrast = contrast_index;",
      "  cfg(jj).threshold = threshold;",
      "  cfg(jj).threshdesc = threshdesc;",
      "  cfg(jj).extent = extent;"
      )
    }

    m_string <- c(m_string,
      "end",
      "spm_extract_anatomical_rois(cfg);"
    )      

    run_matlab_code(m_string, endlines = FALSE, verbose = TRUE, add_clear_all = FALSE)

    #this argument to run_matlab_code adds the paths at the end, not the beginning, which doesn't help
    #, paths_to_add = matlab_scripts)    
  }
}
