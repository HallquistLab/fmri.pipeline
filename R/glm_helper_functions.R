#' Wrapper for running an AFNI command safely within R
#'
#' @param args AFNI command string to be run
#' @param afnidir Location of AFNI installation. If NULL, the function will search the environment for AFNIDIR or afni on PATH
#' @param stdout File target for redirecting stdout. If NULL, stdout will not be written to file.
#' @param stderr File target for redirecting stderr. If NULL, stderr will not be written to file.
#' @param echo Whether to print AFNI command to the screen. Default: TRUE
#' @param omp_num_threads sets the number of OpenMP threads used for this AFNI command (if supported)
#' @param ... Arguments passed through to the `system` command.
#'
#' @return The exit status of the executed AFNI command. 0 for success, non-zero for failure
#'
#' @details
#'
#' This command ensures that AFNI commands are run in an environment that is setup correctly
#'
#' @author Michael Hallquist
#' @export
#'
#' @examples
#'
#' \dontrun{
#' run_afni_command("3dcopy test_data copy_data")
#' }
run_afni_command <- function(args, afnidir=NULL, stdout=NULL, stderr=NULL, echo = TRUE, omp_num_threads=1L, ...) {
  checkmate::assert_string(args, null.ok = FALSE)
  checkmate::assert_string(afnidir, null.ok = TRUE)
  if (!is.null(afnidir)) checkmate::assert_directory_exists(afnidir)

  checkmate::assert_string(stdout, null.ok = TRUE)
  checkmate::assert_string(stderr, null.ok = TRUE)
  checkmate::assert_logical(echo, len = 1L)
  checkmate::assert_integerish(omp_num_threads, lower = 1, len = 1)

  if (is.null(afnidir)) {
    # look for AFNIDIR in system environment if not passed in
    afnidir <- Sys.getenv("AFNIDIR")
    if (isFALSE(nzchar(afnidir))) {
      # look for location of afni binary on PATH
      afni_loc <- system("command -v afni", intern = TRUE)
      exit_code <- attr(afni_loc, "status")
      if (!is.null(exit_code) && exit_code == 1) {
        warning("Could not find afni using AFNIDIR or system PATH. Defaulting to ", paste0(normalizePath("~/"), "/abin"))
        afnidir <- paste0(normalizePath("~/"), "/abin")
      } else {
        afnidir <- dirname(afni_loc)
      }
    }
  }

  checkmate::assert_file_exists(file.path(afnidir, "afni")) # make sure afni actually exists in this folder

  cur_threads <- Sys.getenv("OMP_NUM_THREADS")
  if (omp_num_threads > 1L) {
    if (isTRUE(echo)) cat("Setting OMP_NUM_THREADS to:", omp_num_threads, "\n")
    Sys.setenv(OMP_NUM_THREADS = omp_num_threads) # setup OMP threads
  }
  Sys.setenv(AFNIDIR=afnidir) #export to R environment
  afnisetup <- paste0("AFNIDIR=", afnidir, "; PATH=${AFNIDIR}:${PATH}; DYLD_FALLBACK_LIBRARY_PATH=${AFNIDIR}; ${AFNIDIR}/")
  afnicmd  <- paste0(afnisetup, args)
  if (!is.null(stdout)) afnicmd <- paste(afnicmd, ">", stdout)
  if (!is.null(stderr)) afnicmd <- paste(afnicmd, "2>", stderr)
  if (isTRUE(echo)) cat("AFNI command:", afnicmd, "\n")
  retcode <- system(afnicmd, ...)
  if (omp_num_threads > 1L && cur_threads != "") {
    Sys.setenv(OMP_NUM_THREADS = cur_threads) # reset threads
  }
  return(retcode)
}


#' Wrapper for running an FSL command safely within R
#'
#' @param args FSL command string to be run
#' @param fsldir Location of FSL installation. If NULL, the function will search the environment for FSLDIR or FSL commands in the PATH
#' @param stdout File target for redirecting stdout. If NULL, stdout will not be captured
#' @param stderr File target for redirecting stderr. If NULL, stderr will not be captured
#' @param echo Whether to print FSL command to the screen. Default: TRUE
#' @param ... Arguments passed through to the `system` command.
#'
#' @return The exit status of the executed FSL command. 0 for success, non-zero for failure
#'
#' @details
#'
#' This command ensures that FSL command are run in an environment with FSL setup correctly.
#'
#' @author Michael Hallquist
#' @export
#'
#' @examples
#' \dontrun{
#' run_fsl_command("fslmaths test_data copy_data")
#' }
run_fsl_command <- function(args, fsldir = NULL, stdout = NULL, stderr = NULL, echo = TRUE, ...) {
  checkmate::assert_string(args, null.ok = FALSE)
  checkmate::assert_string(fsldir, null.ok = TRUE)
  if (!is.null(fsldir)) checkmate::assert_directory_exists(fsldir)

  checkmate::assert_string(stdout, null.ok = TRUE)
  checkmate::assert_string(stderr, null.ok = TRUE)
  checkmate::assert_logical(echo, len = 1L)

  if (is.null(fsldir)) {
    # look for FSLDIR in system environment if not passed in
    fsldir <- Sys.getenv("FSLDIR")
    if (isFALSE(nzchar(fsldir))) {
      # check for FSLDIR in .bashrc or .profile
      bashrc_fsldir <- ""
      if (file.exists("~/.profile")) {
        bashrc_fsldir <- system("source ~/.profile && echo $FSLDIR", intern = TRUE)
      }

      if (nzchar(bashrc_fsldir) && file.exists("~/.bashrc")) {
        bashrc_fsldir <- system("source ~/.bashrc && echo $FSLDIR", intern = TRUE)
      }

      # Fallback: look for location of fsl feat on PATH
      if (nzchar(bashrc_fsldir)) {
        feat_loc <- system("command -v feat", intern = TRUE)
        exit_code <- attr(feat_loc, "status")
        if (!is.null(exit_code) && exit_code == 1) {
          warning("Could not find FSL using FSLDIR or system PATH. Defaulting to Defaulting to /usr/local/fsl.")
          fsldir <- "/usr/local/fsl"
        } else {
          fsldir <- dirname(dirname(feat_loc))
        }
      }
    }
  }

  checkmate::assert_file_exists(file.path(fsldir, "bin", "feat")) # make sure FSL actually exists in this folder

  #Sys.setenv(LD_LIBRARY_PATH="/gpfs/group/mnh5174/default/sw/openblas/lib")
  Sys.setenv(FSLDIR=fsldir) #export to R environment
  fslsetup <- paste0("FSLDIR=", fsldir, "; PATH=${FSLDIR}/bin:${PATH}; . ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/")
  fslcmd <- paste0(fslsetup, args)
  if (!is.null(stdout)) fslcmd <- paste(fslcmd, ">", stdout)
  if (!is.null(stderr)) fslcmd <- paste(fslcmd, "2>", stderr)
  if (isTRUE(echo)) cat("FSL command: ", fslcmd, "\n")
  retcode <- system(fslcmd, ...)
  return(retcode)
}

# Internal function used to get the commands that should be included in batch submission scripts
# to setup the compute environment.
# @param gpa a glm_pipeline_arguments object
# @param what which programs to pull from the compute environment. If not specified, it pulls all
#   programs. At present, the options are "global" (only global setup), "fsl", "afni", "spm", and "r".
get_compute_environment <- function(gpa, what="all") {
  # always include global
  compute_string <- gpa$parallel$compute_environment$global

  if (any(what %in% c("all", "fsl"))) {
    compute_string <- c(compute_string, gpa$parallel$compute_environment$fsl)
  }

  if (any(what %in% c("all", "afni"))) {
    compute_string <- c(compute_string, gpa$parallel$compute_environment$afni)
  }

  if (any(what %in% c("all", "spm")) && "spm" %in% gpa$glm_software) {
    compute_string <- c(compute_string, gpa$parallel$compute_environment$spm)
  }

  if (any(what %in% c("all", "r"))) {
    compute_string <- c(compute_string, gpa$parallel$compute_environment$r)
  }

  return(compute_string)
}

gpa_uses_3dlmer <- function(gpa) {
  if (is.null(gpa$l3_models) || !is.list(gpa$l3_models) || is.null(gpa$l3_models$models)) {
    return(FALSE)
  }

  modes <- vapply(
    gpa$l3_models$models,
    function(mm) {
      mode <- mm$l3_input_mode
      if (is.null(mode) || length(mode) == 0L || is.na(mode[1L])) return("")
      tolower(trimws(as.character(mode[1L])))
    },
    character(1)
  )

  any(modes == "3dlmer")
}

# write a small test script for the programs used in the compute environment
# currently capable of testing FSL, R, and AFNI
test_compute_environment <- function(gpa, what="all", stop_on_fail=TRUE) {
  #"set -x",
  what <- unique(tolower(what))
  uses_3dlmer <- gpa_uses_3dlmer(gpa)
  if (isTRUE(uses_3dlmer) && any(what %in% c("all", "afni"))) {
    what <- unique(c(what, "r"))
  }

  prog_str <- c(
    "#!/bin/bash", "checks_passed=1", "",
    "echo 'Checking the compute environment settings'",
    "echo '-----------------------------------------'",
    get_compute_environment(gpa, what = what)
  )

  if (any(what %in% c("all", "afni"))) {
    prog_str <- c(
      prog_str,
      "afni_loc=$(command -v afni)
      [ $? -eq 0 ] && afni_exists=1 || afni_exists=0
      if [ $afni_exists -eq 0 ]; then
        echo 'Cannot find afni!'
        checks_passed=0
      else
        echo \"afni found: ${afni_loc}\"
      fi"
    )

    if (isTRUE(uses_3dlmer)) {
      prog_str <- c(
        prog_str,
        "lmer_loc=$(command -v 3dLMEr)
        [ $? -eq 0 ] && lmer_exists=1 || lmer_exists=0
        if [ $lmer_exists -eq 0 ]; then
          echo 'Cannot find 3dLMEr, which is required by l3_input_mode=3dlmer models!'
          checks_passed=0
        else
          echo \"3dLMEr found: ${lmer_loc}\"
        fi"
      )
    }
  }

  if (any(what %in% c("all", "fsl"))) {
    prog_str <- c(
      prog_str,
      "feat_loc=$(command -v feat)
      [ $? -eq 0 ] && feat_exists=1 || feat_exists=0
      if [ $feat_exists -eq 0 ]; then
        echo 'Cannot find FSL feat!'
        checks_passed=0
      else
        echo \"feat found: ${feat_loc}\"
      fi
      if [ -n \"$FSLDIR\" ]; then
        echo \"FSLDIR set: $FSLDIR\"
      else
        echo 'FSLDIR is not set. This is necessary for FSL setup.'
        checks_passed=0
      fi
      "
    )
  }

  if (any(what %in% c("all", "r"))) {
    prog_str <- c(
      prog_str,
      "r_loc=$(command -v R)
      [ $? -eq 0 ] && r_exists=1 || r_exists=0
      if [ $r_exists -eq 0 ]; then
        echo 'Cannot find R!'
        checks_passed=0
      else
        echo \"R found: ${r_loc}\"
      fi
      "
    )
  }

  if (isTRUE(uses_3dlmer) && any(what %in% c("all", "afni", "r"))) {
    lmer_r_packages <- NULL
    if (!is.null(gpa$glm_settings) && !is.null(gpa$glm_settings$afni)) {
      lmer_r_packages <- gpa$glm_settings$afni$lmer_r_packages
    }
    if (is.null(lmer_r_packages) || length(lmer_r_packages) == 0L) {
      lmer_r_packages <- c("lme4", "lmerTest", "phia", "snow")
    }
    lmer_r_packages <- unique(as.character(lmer_r_packages))
    lmer_r_packages <- lmer_r_packages[nzchar(lmer_r_packages)]
    lmer_pkg_expr <- paste(shQuote(lmer_r_packages), collapse = ", ")
    prog_str <- c(
      prog_str,
      sprintf(
        "rscript_loc=$(command -v Rscript)
        [ $? -eq 0 ] && rscript_exists=1 || rscript_exists=0
        if [ $rscript_exists -eq 0 ]; then
          echo 'Cannot find Rscript for AFNI 3dLMEr R package checks!'
          checks_passed=0
        else
          echo \"Rscript found for AFNI 3dLMEr checks: ${rscript_loc}\"
          Rscript -e \"pkgs <- c(%s); missing <- pkgs[!vapply(pkgs, function(pkg) requireNamespace(pkg, quietly = TRUE), logical(1))]; if (length(missing)) stop('Missing R package(s) needed by AFNI 3dLMEr: ', paste(missing, collapse = ', ')); cat('AFNI 3dLMEr R package check passed: ', paste(pkgs, collapse = ', '), '\\n', sep = '')\"
          if [ $? -ne 0 ]; then
            checks_passed=0
          fi
        fi",
        lmer_pkg_expr
      )
    )
  }

  if (any(what %in% c("all", "spm")) && "spm" %in% gpa$glm_software) {
    matlab_cmd <- gpa$glm_settings$spm$matlab_cmd
    if (is.null(matlab_cmd) || !nzchar(matlab_cmd)) matlab_cmd <- "matlab"
    prog_str <- c(
      prog_str,
      sprintf(
        "%s_loc=$(command -v %s)\n      [ $? -eq 0 ] && %s_exists=1 || %s_exists=0\n      if [ $%s_exists -eq 0 ]; then\n        echo 'Cannot find %s!'\n        checks_passed=0\n      else\n        echo \"%s found: ${%s_loc}\"\n      fi",
        matlab_cmd, matlab_cmd, matlab_cmd, matlab_cmd, matlab_cmd, matlab_cmd, matlab_cmd, matlab_cmd
      )
    )
  }

  prog_str <- c(
    prog_str,
    "
    if [ $checks_passed -eq 0 ]; then
      echo 'One or more checks failed.'
      exit 1
    else
      echo 'All checks passed. Please still verify that the program locations match your expectations!'
    fi
    "
  )

  check_script <- tempfile(pattern="chk")
  writeLines(prog_str, check_script)
  res <- suppressWarnings(system(paste("bash -lc", shQuote(paste("source", check_script))), intern = TRUE))
  cat(res, sep = "\n")
  exit_code <- attr(res, "status")
  if ("spm" %in% gpa$glm_software && isTRUE(gpa$glm_settings$spm$require_matlab)) {
    matlab_cmd <- gpa$glm_settings$spm$matlab_cmd
    if (is.null(matlab_cmd) || !nzchar(matlab_cmd)) matlab_cmd <- "matlab"
    if (!is.null(exit_code) && exit_code != 0L) {
      stop(sprintf("MATLAB command '%s' not available; spm.require_matlab is TRUE.", matlab_cmd))
    }
  }
  if (!is.null(exit_code) && exit_code == 1) {
    if (stop_on_fail) {
      stop("One or more problems were detected in the compute environment. Setup cannot continue.")
    } else {
      warning("One or more problems were detected in the compute environment. The pipeline may fail if these are not resolved.")
    }
  }
}

# Local SPM test using the configured MATLAB/Octave command and compute environment.
# This is intended to catch missing SPM paths before submitting jobs.
test_spm_compute_environment <- function(gpa, stop_on_fail = TRUE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_logical(stop_on_fail, len = 1L)
  lg <- lgr::get_logger("glm_pipeline/test_spm_compute_environment")

  matlab_cmd <- gpa$glm_settings$spm$matlab_cmd
  if (is.null(matlab_cmd) || !nzchar(matlab_cmd)) matlab_cmd <- "matlab"
  matlab_args <- gpa$glm_settings$spm$matlab_args
  if (is.null(matlab_args) || !nzchar(matlab_args)) matlab_args <- "-batch"
  matlab_timeout <- gpa$glm_settings$spm$matlab_timeout
  if (is.null(matlab_timeout) || !is.numeric(matlab_timeout) || length(matlab_timeout) != 1L || is.na(matlab_timeout)) {
    matlab_timeout <- 120
  }

  spm_path <- gpa$glm_settings$spm$spm_path
  if (is.null(spm_path) || !nzchar(spm_path)) {
    msg <- "spm_path is not set; cannot validate SPM availability."
    if (stop_on_fail) stop(msg) else warning(msg)
    return(invisible(FALSE))
  }

  compute_env <- get_compute_environment(gpa, c("spm"))
  spm_cmd <- paste0(
    "try; ",
    "addpath('", spm_path, "'); ",
    "spm('defaults','fmri'); ",
    "spm_get_defaults('cmdline',1); ",
    "spm_jobman('initcfg'); ",
    "catch ME; disp(getReport(ME,'extended')); exit(1); end; exit(0);"
  )
  matlab_call <- paste(matlab_cmd, matlab_args, shQuote(spm_cmd))
  cmd <- if (!is.null(compute_env) && length(compute_env) > 0L) {
    paste(c(compute_env, matlab_call), collapse = " && ")
  } else {
    matlab_call
  }
  # Use coreutils timeout if available since R.utils::withTimeout cannot reliably interrupt system().
  timeout_bin <- Sys.which("timeout")
  if (nzchar(timeout_bin)) {
    cmd <- paste(timeout_bin, matlab_timeout, "bash -lc", shQuote(cmd))
  } else {
    cmd <- paste("bash -lc", shQuote(cmd))
  }

  lg$info("Starting SPM compute environment check (timeout: %ss)", matlab_timeout)
  res <- tryCatch(
    R.utils::withTimeout(
      suppressWarnings(system(cmd, intern = TRUE)),
      timeout = matlab_timeout,
      onTimeout = "error"
    ),
    TimeoutException = function(e) {
      msg <- sprintf("SPM compute environment check timed out after %ss. MATLAB may be hanging.", matlab_timeout)
      if (stop_on_fail) stop(msg) else warning(msg)
      return(structure(character(0), status = 124L))
    }
  )
  lg$info("Finished SPM compute environment check")
  exit_code <- attr(res, "status")
  if (!is.null(exit_code) && exit_code != 0L) {
    msg <- paste("SPM compute environment check failed:", paste(res, collapse = "\n"))
    if (stop_on_fail) stop(msg) else warning(msg)
    return(invisible(FALSE))
  }

  invisible(TRUE)
}

#' internal function for getting a file extension that may include a compressed ending
#' this is adapted from neurobase::file_imgext, but extended for other file types
#' @keywords internal
file_ext <- function(file, withdot = TRUE) {
  file <- tolower(file)
  matches <- grepl("^.*\\.(csv|dat|hdr|img|brik|head|nii|txt|tsv|yaml|json|1d)(\\.gz|\\.bz2|\\.zip|\\.xz)*$", file)
  ext <- rep(NA, length = length(file)) # return NA for inputs that can't be parsed
  ext[matches] <- sub("^(.*)\\.(csv|dat|hdr|img|brik|head|nii|txt|tsv|yaml|json|1d)(\\.gz|\\.bz2|\\.zip|\\.xz)*$", "\\2\\3", file[matches])
  if (isTRUE(withdot)) ext[matches] <- paste0(".", ext[matches])
  return(ext)
}


#' internal function for returning a file (and path) without its extension
#' at present, it returns NA if no recognized extension is there
#' @keywords internal
file_sans_ext <- function(file, withdot = TRUE) {
  matches <- grepl("^.*\\.(csv|dat|hdr|img|brik|head|nii|txt|tsv|yaml|json|1d)(\\.gz|\\.bz2|\\.zip|\\.xz)*$", file, ignore.case = TRUE)
  fout <- rep(NA, length=length(file)) # return NA for inputs that can't be parsed
  fout[matches] <- sub("^(.*)\\.(csv|dat|hdr|img|brik|head|nii|txt|tsv|yaml|json|1d)(\\.gz|\\.bz2|\\.zip|\\.xz)*$", "\\1", file[matches], ignore.case = TRUE)
  return(fout)
}

#' helper function for generating motion regressors from raw 6-parameter motion coregistration
#'
#' @param motion_params_file file containing 6 motion parameters
#' @param col.names names of columns in \code{motion_params_file}, in order from left to right
#' @param drop_volumes number of volumes to drop from beginning of motion params
#' @param last_volume final volume to include from motion params (if truncated at end). If \code{NULL},
#'   then the end of the time series is not truncated.
#' @param rot_units The units of the rotation parameters. Default is "rad" for radians
#' @param tra_units The units of the translation parameter. Only support "mm" right now millimeters
#'
#' @importFrom data.table fread setnames
#' @importFrom checkmate assert_file_exists assert_subset assert_integerish
#' @keywords internal
generate_motion_regressors <- function(motion_params_file = "motion.par",
                                       col.names = c("rx", "ry", "rz", "tx", "ty", "tz"),
                                       demean = FALSE, drop_volumes = 0L, last_volume=NULL,
                                       rot_units="rad", tra_units="mm",
                                       na.strings="NA", lg=NULL) {

  checkmate::assert_class(lg, "Logger")
  checkmate::assert_file_exists(motion_params_file)
  if (is.null(col.names)) {
    # explicit defaults in case of null
    lg$debug("In generate_motion_regressors, defaulting column order rx, ry, rz, tx, ty, tz")
    col.names <- c("rx", "ry", "rz", "tx", "ty", "tz")
  }
  checkmate::assert_subset(col.names, c("rx", "ry", "rz", "tx", "ty", "tz"))
  checkmate::assert_logical(demean, len = 1L)
  checkmate::assert_integerish(drop_volumes, lower = 0)
  checkmate::assert_subset(rot_units, c("rad", "deg"))
  checkmate::assert_subset(tra_units, c("mm"))
  checkmate::assert_character(na.strings)

  # read raw motion parameters file
  mot <- data.table::fread(motion_params_file, na.strings = na.strings)

  if (ncol(mot) != 6L) {
    lg$error("In generate_motion_regressors, motion file %s has %d columns instead of 6!", motion_params_file, ncol(mot))
    return(invisible(NULL))
  } else {
    data.table::setnames(mot, col.names)
  }

  if (is.null(last_volume)) { last_volume <- nrow(mot) }
  checkmate::assert_integerish(last_volume, upper=nrow(mot))

  # subset motion regressors to match drops and truncations in fmri data
  mot <- mot[(1 + drop_volumes):last_volume, ]

  mot_deriv <- mot[, lapply(.SD, function(x) { c(0, diff(x)) })]
  data.table::setnames(mot_deriv, paste0("d", names(mot))) #add delta to names
  mot <- cbind(mot, mot_deriv)

  mot_quad <- mot[, lapply(.SD, function(x) { x^2 })]
  data.table::setnames(mot_quad, paste0("q", names(mot))) #add delta to names
  mot <- cbind(mot, mot_quad)

  # https://wiki.cam.ac.uk/bmuwiki/FMRI
  if (rot_units=="rad") {
    framewise_displacement <- apply(mot_deriv[, c("drx", "dry", "drz")], 1, function(x) 50 * sum(abs(x))) +
      apply(mot_deriv[, c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
  } else if (rot_units=="deg") {
    framewise_displacement <- apply(mot_deriv[, c("drx", "dry", "drz")], 1, function(x) 50 * (pi / 180) * sum(abs(x))) +
      apply(mot_deriv[, c("dtx", "dty", "dtz")], 1, function(x) sum(abs(x)))
  } else {
    stop("unknown")
  }
  mot <- cbind(mot, framewise_displacement = framewise_displacement)

  # Demean all columns, if requested
  # Note that this is not used internally at present to avoid surprising results when testing motion thresholds (e.g., framewise_displacement)
  if (isTRUE(demean)) {
    mot <- mot[, lapply(.SD, function(x) { x - mean(x, na.rm=TRUE) }) ]
  }

  ## just PCA motion on the current run
  ## mregressors <- pca_motion(run_nifti[r], runlengths[r], motion_parfile="motion.par", numpcs=3, drop_volumes=drop_volumes)$motion_pcs_concat

  return(as.data.frame(mot))
}

#' compute spike regressors from a data.frame or matrix of motion parameters and a set of
#'   expressions that are evaluated against this parameter matrix.
#' 
#' @param mot a data.frame or matrix with volumes on rows and named motion parameters on columns
#' @param spike_volume a character vector of expressions to evaluate against \code{mot}. Resulting
#'   columns will be prefixed with the names of each expression. If the expressions are unnamed,
#'   prefixes will be expr1_, expr2_, etc.
#'
#' @return a matrix of spike regressors (volumes on rows, spike regressors on columns)
#' @importFrom dplyr lead lag
#' @keywords internal
compute_spike_regressors <- function(mot = NULL, spike_volume = NULL, lg=NULL) {
  if (is.null(mot)) return(NULL)
  if (is.null(spike_volume)) return(NULL)

  checkmate::assert_character(spike_volume, null.ok = TRUE)
  checkmate::assert_class(lg, "Logger")

  spikes <- do.call(cbind, lapply(seq_along(spike_volume), function(ii) {
    has_bounds <- grepl(";", spike_volume[ii], fixed = TRUE)
    if (isTRUE(has_bounds)) {
      esplit <- strsplit(spike_volume[ii], "\\s*;\\s*", perl = TRUE)[[1L]]
      stopifnot(length(esplit) == 2L)
      spike_bounds <- as.integer(eval(parse(text = paste0("c(", esplit[1], ")"))))
      expr <- esplit[2L]
    } else {
      spike_bounds <- 0L # only volume where expr is TRUE
      expr <- spike_volume[ii]
    }

    spike_vec <- tryCatch(with(mot, eval(parse(text = expr))),
      error = function(e) {
        lg$error(
          "Problem evaluating spike regressors for file %s, expr: %s",
          motion_params_file, expr
        )
        return(NULL)
      }
    )

    which_spike <- which(spike_vec)
    if (length(which_spike) == 0L) return(NULL) #don't attempt to generate spikes if no TRUE values

    # will evaluate to NULL if which_spike is empty
    spike_df <- do.call(cbind, lapply(which_spike, function(xx) {
      vec <- rep(0, nrow(mot))
      vec[xx] <- 1
      return(vec)
    }))

    colnames(spike_df) <- paste0("spike_", seq_len(ncol(spike_df)))

    if (!is.null(spike_df) && !identical(spike_bounds, 0L)) {
      # apply shifts
      shifts <- spike_bounds[spike_bounds != 0L]
      res <- do.call(cbind, lapply(shifts, function(ss) {
        shift_mat <- apply(spike_df, 2, function(col) {
          if (ss < 0) {
            dplyr::lead(col, abs(ss), default = 0)
          } else {
            dplyr::lag(col, ss, default = 0)
          }
        })
        colnames(shift_mat) <- paste0("spike_", ifelse(ss < 0, "m", "p"), abs(ss), "_", 1:ncol(shift_mat))
        return(shift_mat)
      }))

      spike_df <- cbind(spike_df, res)
    }

    if (is.null(names(spike_volume)[ii]) || names(spike_volume)[ii] == "") {
      colnames(spike_df) <- paste0("expr", ii, "_", colnames(spike_df))
    } else {
      colnames(spike_df) <- paste0(names(spike_volume)[ii], "_", colnames(spike_df))
    }
    return(spike_df)

  }))

  # we don't want to generate any duplicate columns (linear dependency problems) 
  # due to overlapping expressions or ranges of evaluation
  spikes <- spikes[, !duplicated(spikes, MARGIN=2), drop = FALSE]

  return(spikes)
}


#' small helper function for compressing motion parameters using PCA
#'
#' @param motion_df volumes x motion parameters data frame for PCA compression
#' @param num_pcs number of principal components to extract
#' @param zscore whether to standardize motion parameters prior to PCA (this is a good idea)
#' @param verbose whether to print the variance explained by PCs and the multiple correlations of
#'   among motion parameters
#'
#' @return A matrix of PCA-compressed motion regressors
#' @importFrom psych smc
#' @importFrom stats prcomp
#' @keywords internal
pca_motion <- function(motion_df, num_pcs=3L, zscore=TRUE, verbose=FALSE) {
  checkmate::assert_data_frame(motion_df)
  checkmate::assert_integerish(num_pcs, lower=1, upper=50)
  checkmate::assert_logical(zscore)
  checkmate::assert_logical(verbose)

  #compute the PCA decomposition of motion parameters and their derivatives
  pc <- prcomp(motion_df, retx=TRUE, scale.=zscore)
  cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))

  if (isTRUE(verbose)) message("first", num_pcs, "motion principal components account for: ", round(cumvar[num_pcs], 3))
  mregressors <- pc$x[, 1:num_pcs] #cf Churchill et al. 2012 PLoS ONE
  attr(mregressors, "variance.explained") <- cumvar[num_pcs]

  if (isTRUE(verbose)) {
    cat("Multiple correlation of motion parameters:\n\n")
    print(round(sqrt(psych::smc(motion_df)), 2))
  }

  return(as.data.frame(mregressors))
}

generateRunMask <- function(mr_files, outdir=getwd(), outfile="runmask") {
  if (file.exists(file.path(outdir, paste0(outfile, ".nii.gz")))) { return(invisible(NULL)) }
  ##generate mask of mr_files where temporal min is > 0 for all runs
  for (f in seq_along(mr_files)) {
    run_fsl_command(paste0("fslmaths ", mr_files[f], " -Tmin -bin ", outdir, "/tmin", f))#, fsldir="/usr/local/ni_tools/fsl")
  }

  ##sum mins together over runs and threshold at number of runs
  run_fsl_command(paste0("fslmaths ", paste(paste0(outdir, "/tmin", seq_along(mr_files)), collapse=" -add "), " ", outdir, "/tminsum"))#, fsldir="/usr/local/ni_tools/fsl")
  run_fsl_command(paste0("fslmaths ", outdir, "/tminsum -thr ", length(mr_files), " -bin ", outdir, "/", outfile))#, fsldir="/usr/local/ni_tools/fsl")
  run_fsl_command(paste0("imrm ", outdir, "/tmin*"))#, fsldir="/usr/local/ni_tools/fsl") #cleanup 
}

#compute the mean of each cluster in roimask
#within a 4d array: x, y, z, subbrik
#
#roimask is an x,y,z array with integer mask values
#ni4d is a concatenated set of cope values, either across subjects
# or, in the case of beta series, within a subject (over trials)
get_cluster_means <- function(roimask, ni4d) {
  #ni4d should be a 4s array in which the fourth dimension is either subject (single beta) or a beta series (single subject)
  maskvals <- sort(unique(as.vector(roimask)))
  maskvals <- maskvals[!maskvals == 0]

  roimats <- sapply(maskvals, function(v) {
    mi <- which(roimask==v, arr.ind=TRUE)
    nsubbriks <- dim(ni4d)[4]
    nvox <- nrow(mi)
    mi4d <- cbind(pracma::repmat(mi, nsubbriks, 1), rep(1:nsubbriks, each=nvox))

    mat <- matrix(ni4d[mi4d], nrow=nvox, ncol=nsubbriks) #need to manually reshape into matrix from vector
    #for each subject, compute huber m-estimator of location/center Winsorizing at 2SD across voxels (similar to voxel mean)
    #clusavg <- apply(mat, 2, function(x) { MASS::huber(x, k=2)$mu })
    clusavg <- apply(mat, 2, mean)

    return(clusavg)
  })
}

# internal function to extract mean beta series
# NOT USED AT THE MOMENT
get_beta_series <- function(inputs, roimask, n_bs=50) {
  #inputs <- inputs[1:5] #speed up testing

  beta_res <- foreach(i=iter(seq_along(inputs)), .packages=c("oro.nifti", "dplyr", "abind"), .export=c("get_cluster_means", "mat2df")) %dopar% {
  #beta_res <- lapply(1:length(inputs), function(i) {
    run_dirs <- list.files(path=inputs[i], pattern="FEAT_LVL1_run\\d+\\.feat", recursive=FALSE, full.names=TRUE)
    run_betas <- list()

    for (r in seq_along(run_dirs)) {
      runnum <- as.numeric(sub(".*/FEAT_LVL1_run(\\d+)\\.feat", "\\1", run_dirs[r], perl=TRUE))
      copes <- file.path(run_dirs[r], "stats", paste0("cope", 1:n_bs, ".nii.gz"))

      if (!all(file.exists(copes))) {
        cat("File listing in dir:", run_dirs[r], list.files(run_dirs[r], recursive=FALSE), sep="\n", append=TRUE)
        cat("Cannot find all ", length(copes), " beta series copes in: ",
          run_dirs[r], "\n",
          file = "beta_series_errors.txt", append = TRUE
        )
        next #don't error for now, just omit
      }

      #in testing, it is faster to use readNIfTI to read each cope file individually (7.5s) than to use fslmerge + readNIfTI on the 4d (13s)
      #concat_file <- tempfile()
      #system.time(run_fsl_command(paste("fslmerge -t", concat_file, paste(copes, collapse=" "))))
      #system.time(cout <- readNIfTI(concat_file, reorient=FALSE)@.Data)

      cout <- do.call(abind, list(along=4, lapply(copes, function(x) { readNIfTI(x, reorient=FALSE)@.Data })))
      beta_series_cluster_means <- get_cluster_means(roimask, cout)
      beta_melt <- mat2df(as.matrix(beta_series_cluster_means), value_name="bs_value", varnames=c("trial", "cluster_number")) %>%
        mutate(feat_input_id=i, run=runnum)

      run_betas[[r]] <- beta_melt
    }

    return(do.call(rbind, run_betas))

  }

  return(do.call(rbind, beta_res))
}

#' Internal helper function for printing a summary of a contrast matrix
#' @param cmat a contrast matrix where row names are the contrast names and column
#'   names are the column names of the corresponding design matrix.
#' @keywords internal
summarize_contrasts <- function(cmat) {
  sapply(seq_len(nrow(cmat)), function(x) {
    con_name <- rownames(cmat)[x]
    cvec <- cmat[x, ]
    nzcols <- which(cvec != 0)
    cols <- colnames(cmat)[nzcols]
    cat("Contrast: ", con_name, "\n")
    for (ii in seq_along(cols)) {
      cat(cols[ii], "=", cvec[nzcols[ii]], "\n")
    }
    cat("-----\n")
  })
}

#' helper function to rename columns of input data.frame to internal nomenclature
#'   based on the variable mapping (vm) vector
#' 
#' @param df a data.frame containing columns to be renamed to internal standards
#' @param vm a named vector of columns in \code{df} that identify internal constructs
#'   such as id, session, and run.
#'
#' @return a modified version of \code{df} with column names modified to use internal names
#' @keywords internal
names_to_internal <- function(df, vm) {
  # look for naming collisions
  nm <- names(df)
  nv <- names(vm)
  conflict <- intersect(nm[nm %in% names(vm)], nv[nv != vm])

  if (length(conflict) > 0L) {
    new_names <- paste0("..", conflict, "..")
    cat(
      "To avoid naming conflicts, renaming these columns:", paste(conflict, collapse = ", "),
      "to:", paste0(new_names, collapse = ", "), "\n"
    )

    nm[which(nm %in% conflict)] <- new_names
  }

  # convert variable names to internal constructs via the variable mapping
  vm_present <- vm[which(vm %in% nm)]
  if (length(vm_present) > 0L) {
    vm_pos <- sapply(vm_present, function(x) which(nm == x)[1])
    nm[vm_pos] <- names(vm_present)
  }

  return(nm)
}

#' small helper function to pull absolute paths to a given column in run_data or trial_data
#'
#' @param mr_df a data.frame from a gpa object, following standard variable mapping nomenclature
#' @param col a character string denoting the column in \code{mr_df} to be used for looking up
#'   absolute paths
#'
#' @details Note that if a given value of the requested column is an absolute path, it will
#'   not be combined with $mr_dir to generate the combined path. Thus, the \code{col} in
#'   \code{mr_df} can contain a mixture of relative an absolute paths. The relative paths will
#'   be combined with 
#' @return Absolute paths to all files in the specified \code{col}
#' @importFrom R.utils getAbsolutePath
#' @keywords internal
get_mr_abspath <- function(mr_df, col="run_nifti") {
  checkmate::assert_data_frame(mr_df)
  checkmate::assert_string(col)
  checkmate::assert_subset(col, names(mr_df))

  sapply(seq_len(nrow(mr_df)), function(ii) {
    if (!"mr_dir" %in% names(mr_df)) {
      # this should be rare (or ideally, not happen at all)
      # but if there is no mr_dir in the mr_df, we can only use the column itself
      wd <- NULL # defaults to getwd()
    } else {
      this_dir <- mr_df$mr_dir[ii]
      if (is.na(this_dir)) {
        wd <- getwd() # defaults to getwd()
      } else {
        wd <- this_dir
      }
    }

    R.utils::getAbsolutePath(mr_df[[col]][ii], workDirectory = wd, expandTilde = TRUE)
  })

}



#' small helper function to populate subject exclusions in $run_data and $subject_data
#'
#' @param gpa a \code{glm_pipeline_arguments} object that already has the \code{$exclude_run}
#'   column populated in \code{$l1_model_setup}.
#' @return a modified copy of \code{gpa} where $exclude_subject has been added to $run_data
#'   and $subject_data. This function also adds $n_good_runs to $subject_data which is useful
#'   if we want to enforce a lower bound on the number of runs used to define subject exclusion.
#' @keywords internal
#' @importFrom dplyr filter count left_join
#' @importFrom checkmate assert_class assert_data_table
#' @importFrom lgr get_logger
#' @author Michael Hallquist
calculate_subject_exclusions <- function(gpa) {
  lg <- lgr::get_logger("glm_pipeline/setup_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(gpa$run_data)
  checkmate::assert_subset(c("id", "session", "run_number", "exclude_run"), names(gpa$run_data)) # verify that exclude_run is populated

  if ("exclude_subject" %in% names(gpa$subject_data) && "exclude_subject" %in% names(gpa$run_data)) {
    lg$debug("exclude_subject already calculated and populated in $subject_data and $run_data")
    return(gpa)
  }

  # if no subject exclusion string is provided, then keep all subjects
  if (is.null(gpa$confound_settings$exclude_subject)) {
    gpa$run_data$exclude_subject <- FALSE
    gpa$subject_data$exclude_subject <- FALSE
    return(gpa)
  }

  n_good_runs_df <- gpa$run_data %>%
    dplyr::group_by(.data$id, .data$session) %>%
    dplyr::summarize(n_good_runs = sum(.data$exclude_run == FALSE)) %>% # sum of non-excluded runs
    dplyr::ungroup() %>%
    dplyr::select("id", "session", "n_good_runs")

  gpa$subject_data <- gpa$subject_data %>%
    dplyr::left_join(n_good_runs_df, by = c("id", "session"))

  gpa$run_data <- gpa$run_data %>%
    dplyr::left_join(n_good_runs_df, by = c("id", "session"))

  # evaluate subject exclusion in the context of the broader $subject_data data.frame
  gpa$subject_data$exclude_subject <- sapply(seq_len(nrow(gpa$subject_data)), function(ss) {
    s_info <- gpa$subject_data[ss, , drop = FALSE]
    tryCatch(with(s_info, eval(parse(text = gpa$confound_settings$exclude_subject))),
      error = function(e) {
        lg$error(
          "Problem evaluating subject exclusion for subject: %s, session: %s, expr: %s",
          s_info$id, s_info$session,
          gpa$confound_settings$exclude_subject
        )
        lg$error("Defaulting to retaining this subject.")
        return(FALSE)
      }
    )
  })

  # propagate subject exclusion down to $run_data
  gpa$run_data <- gpa$run_data %>%
    dplyr::left_join(
      gpa$subject_data %>% dplyr::select("id", "session", "exclude_subject"),
      by = c("id", "session")
    )

  return(gpa)
}

#' Internal function to add last_onset and last_offset columns to gpa$run_data based on events in l1_models
#' @param gpa a \code{glm_pipeline_arguments} object containing valid $run_data and $l1_models objects
#' @param lg a Logger object for logging results
#' @details The last_onset and last_offset columns are calculated for each run based on the timing of all
#'   events in the $l1_models$events list. These are then used to facilitate run truncation if the user
#'   requests truncation after a final onset or offset.
#' @return a modified copy of gpa with $run_data populated with last_offset and last_onset columns (times in seconds)
#' @keywords internal
#' @importFrom dplyr summarize group_by left_join
#' @importFrom utils capture.output
populate_last_events <- function(gpa, lg) {
  # get all events as a long data.frame
  m_events <- data.table::rbindlist(lapply(gpa$l1_models$events, function(this_event) this_event$data))

  # default to 0 ISI if column is missing (shouldn't happen under normal circumstances)
  if (!"isi" %in% names(m_events)) m_events[, isi := 0]

  # convert NA values for ISI to zero with message
  na_isi <- is.na(m_events$isi)

  if (all(na_isi)) {
    lg$debug("In populate_last_events, all ISIs are NA. Converting to 0")
    m_events[, isi := 0]
  } else if (any(na_isi)) {
    lg$debug("Changing NA ISIs to 0 for the following cases: ")
    lg$debug("%s", capture.output(print(m_events[na_isi, c("id", "session", "run_number", "event", "trial")], nrows = 1e5)))
    m_events$isi[na_isi] <- 0
  }

  # The remaining processing follows ordinary data.frame/dplyr conventions.
  # m_events is a fresh rbindlist result, so conversion by reference is safe.
  data.table::setDF(m_events)

  last_events <- m_events %>%
    dplyr::group_by(.data$id, .data$session, .data$run_number) %>%
    dplyr::summarize(
      last_event_idx = which.max(.data$onset),
      last_event = .data$event[which.max(.data$onset)],
      last_onset = max(.data$onset, na.rm = TRUE),
      last_offset = max(.data$onset + .data$duration, na.rm = TRUE),
      last_isi = .data$isi[which.max(.data$onset)],
      .groups="drop")

  gpa$run_data <- gpa$run_data %>%
    dplyr::left_join(last_events, by=c("id", "session", "run_number"))

  return(gpa)
}

# slurm_job_array <- function(job_name = "slurm_array") {
# }

#' helper function to convert the $sched_args field to a vector
#'  of directives that can be included dynamically in the header of
#'  PBS or SBATCH scripts.
#'
#' @param gpa a \code{glm_pipeline_arguments} object containing the
#'   $parallel$sched_args field
#'
#' @return a character vector where the $sched_args elements are converted to
#'   corresponding # SBATCH or # PBS directives for inclusion in scripts
#' @keywords internal
sched_args_to_header <- function(gpa) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")

  #if there are no scheduler arguments, return nothing
  if (is.null(gpa$parallel$sched_args)) return(NULL)

  is_slurm <- gpa$scheduler == "slurm"
  directives <- sapply(gpa$parallel$sched_args, function(x) {
    ifelse(isTRUE(is_slurm), paste("#SBATCH", x), paste("#PBS", x))
  }, USE.NAMES = FALSE)
  return(directives)
}

parse_emm_spec_string <- function(vv) {
  checkmate::assert_string(vv)
  parts <- strsplit(vv, "\\|", perl = TRUE)[[1L]]
  focal_raw <- trimws(parts[1L])
  by_raw <- if (length(parts) > 1L) trimws(paste(parts[-1L], collapse = "+")) else ""

  focal_vars <- trimws(unlist(strsplit(focal_raw, "[:*+]", perl = TRUE)))
  focal_vars <- focal_vars[nzchar(focal_vars)]
  by_vars <- if (nzchar(by_raw)) {
    bv <- trimws(unlist(strsplit(by_raw, "[:*+,\\s]+", perl = TRUE)))
    bv[nzchar(bv)]
  } else {
    character(0)
  }

  focal_formula_rhs <- paste(focal_vars, collapse = " * ")
  list(
    raw = vv,
    focal_raw = focal_raw,
    focal_vars = focal_vars,
    focal_formula = stats::as.formula(paste("~", focal_formula_rhs)),
    by_vars = if (length(by_vars) > 0L) by_vars else NULL
  )
}

sanitize_contrast_token <- function(x) {
  x <- as.character(x)
  # Disambiguate valence signs like CS+ vs CS- before make.names turns both into 'CS.'
  has_plus_or_trailing_minus <- grepl("\\+", x) | grepl("(?<=[A-Za-z0-9])-(?=\\s|$|\\))", x, perl = TRUE)
  if (any(has_plus_or_trailing_minus)) {
    x[has_plus_or_trailing_minus] <- gsub("\\+", "_pos", x[has_plus_or_trailing_minus])
    x[has_plus_or_trailing_minus] <- gsub("(?<=[A-Za-z0-9])-(?=\\s|$|\\))", "_neg", x[has_plus_or_trailing_minus], perl = TRUE)
  }
  x
}

sanitize_contrast_diff_label <- function(contrast_chr) {
  # Split on emmeans ' - ' separator first so minus signs inside factor levels can be preserved distinctly
  pieces <- strsplit(as.character(contrast_chr), "\\s+-\\s+", perl = TRUE)
  out <- vapply(pieces, function(p) {
    p <- gsub("^\\((.*)\\)$", "\\1", trimws(p))
    p <- sanitize_contrast_token(p)
    paste(make.names(p, unique = FALSE), collapse = "_M_")
  }, character(1))
  out
}

format_by_prefix <- function(edata, by_vars) {
  if (is.null(by_vars) || length(by_vars) == 0L) return(rep("", nrow(edata)))
  present_by <- intersect(by_vars, names(edata))
  if (length(present_by) == 0L) return(rep("", nrow(edata)))
  apply(edata[, present_by, drop = FALSE], 1, function(row_vals) {
    vals <- sanitize_contrast_token(row_vals)
    paste0(paste(present_by, make.names(vals, unique = FALSE), sep = ".", collapse = "_"), "_")
  })
}

quiet_emm <- function(expr) {
  withCallingHandlers(
    suppressMessages(expr),
    warning = function(w) {
      if (grepl("essentially perfect fit", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

mask_empty_cell_marginal_slices <- function(lmfit, ee, focal_vars, by_vars = NULL, weights = "cells", label = "", ignore_blocking_vars = character(0)) {
  if (is.null(lmfit) || identical(weights, "equal")) return(ee)
  rg <- tryCatch(quiet_emm(emmeans::ref_grid(lmfit, nesting = NULL)), error = function(e) NULL)
  if (is.null(rg) || !".wgt." %in% names(rg@grid)) return(ee)

  vtypes <- extract_lm_var_types(lmfit)
  cat_in_rg <- setdiff(intersect(names(rg@levels), vtypes$cat_vars), ignore_blocking_vars)
  if (length(cat_in_rg) == 0L) return(ee)

  grid_check <- rg@grid
  if (length(ignore_blocking_vars) > 0L && any(ignore_blocking_vars %in% names(grid_check))) {
    grid_check <- stats::aggregate(
      list(.wgt. = grid_check$.wgt.),
      by = grid_check[, cat_in_rg, drop = FALSE],
      FUN = function(x) sum(x, na.rm = TRUE)
    )
  }

  if (!any(grid_check$.wgt. == 0, na.rm = TRUE)) return(ee)

  conditioning_vars <- unique(c(focal_vars, by_vars))
  marginalized_cats <- setdiff(cat_in_rg, conditioning_vars)
  if (length(marginalized_cats) == 0L) return(ee)

  cond_cols <- intersect(conditioning_vars, intersect(names(ee@grid), names(grid_check)))
  bad_rows <- rep(FALSE, nrow(ee@grid))
  if (length(cond_cols) == 0L) {
    # Overall mean marginalizing across all categorical cells when at least one cell is empty
    bad_rows[] <- TRUE
  } else {
    for (rr in seq_len(nrow(ee@grid))) {
      match_idx <- rep(TRUE, nrow(grid_check))
      for (cc in cond_cols) {
        match_idx <- match_idx & (as.character(grid_check[[cc]]) == as.character(ee@grid[[cc]][rr]))
      }
      if (any(match_idx) && any(grid_check$.wgt.[match_idx] == 0, na.rm = TRUE)) {
        bad_rows[rr] <- TRUE
      }
    }
  }

  if (any(bad_rows)) {
    warning(
      sprintf(
        "Dropping marginal contrast(s) for '%s' under weights='%s' because one or more marginalized cells (%s) have zero observations.",
        if (nzchar(label)) label else "overall",
        weights,
        paste(marginalized_cats, collapse = ", ")
      ),
      call. = FALSE
    )
    ee@linfct[bad_rows, ] <- NA_real_
  }
  ee
}

align_custom_contrasts_to_colnames <- function(raw_custom, c_colnames) {
  if (is.null(raw_custom) || length(raw_custom) == 0L) return(NULL)
  cmat <- parse_custom_contrast_spec(raw_custom, reg_names = NULL)
  if (is.null(cmat) || nrow(cmat) == 0L) return(NULL)

  if (is.null(colnames(cmat))) {
    if (ncol(cmat) != length(c_colnames)) {
      warning(
        sprintf(
          "Dropping unnamed custom contrast matrix because ncol (%d) does not match model regressors (%d).",
          ncol(cmat), length(c_colnames)
        ),
        call. = FALSE
      )
      return(NULL)
    }
    colnames(cmat) <- c_colnames
    return(cmat)
  }

  missing_cols <- setdiff(colnames(cmat), c_colnames)
  keep_rows <- rep(TRUE, nrow(cmat))
  if (length(missing_cols) > 0L) {
    uses_missing <- apply(cmat[, missing_cols, drop = FALSE], 1, function(r) any(abs(r) > 1e-8))
    if (any(uses_missing)) {
      warning(
        sprintf(
          "Dropping custom contrast(s) (%s) that reference omitted/aliased regressor(s): %s.",
          paste(rownames(cmat)[uses_missing], collapse = ", "),
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
      keep_rows <- !uses_missing
    }
  }

  cmat <- cmat[keep_rows, , drop = FALSE]
  if (nrow(cmat) == 0L) return(NULL)

  out <- matrix(
    0,
    nrow = nrow(cmat),
    ncol = length(c_colnames),
    dimnames = list(rownames(cmat), c_colnames)
  )
  shared_cols <- intersect(colnames(cmat), c_colnames)
  if (length(shared_cols) > 0L) {
    out[, shared_cols] <- cmat[, shared_cols, drop = FALSE]
  }
  out
}

#' helper function to generate a contrast matrix from an lm() object
#'   and a set of user-specified contrasts using emmeans
#'
#' @param mobj a model object created by build_l<X>_models
#' @param lmfit an optional lm() object used for emmeans calculations. If provided, this object
#'   will be used instead of mobj$lmfit (the parent lm on the overall dataset).
#' @param allow_interactive_dedupe Optional logical flag controlling whether duplicate contrasts
#'   prompt interactively (defaults to \code{FALSE} on refits, otherwise \code{rlang::is_interactive()}).
#'
#' @return a modified copy of the model object \code{mobj} with $contrast_list and $contrasts
#'   fully populated
#' @keywords internal
#' @importFrom emmeans emmeans emtrends
#' @importFrom graphics pairs
#' @importFrom stats coef na.omit
#' @importFrom utils menu
get_contrasts_from_spec <- function(mobj, lmfit = NULL, allow_interactive_dedupe = NULL) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "l1_wi_spec", "hi_model_spec")) # verify that we have an object of known structure
  is_refit <- !is.null(lmfit)
  if (is.null(allow_interactive_dedupe)) {
    allow_interactive_dedupe <- if (isTRUE(is_refit)) FALSE else rlang::is_interactive()
  }

  # Preserve custom contrasts before resetting cached contrast_list on refits
  preserved_custom <- mobj$contrast_spec$custom %||% mobj$contrast_list$custom
  if (!is.null(preserved_custom)) {
    mobj$contrast_spec$custom <- preserved_custom
  }

  if (is.null(lmfit)) {
    lmfit <- mobj$lmfit # use parent lmfit
  } else {
    # override existing contrasts to force lmfit-specific respecification
    mobj$contrast_list <- list()
    mobj$contrasts <- NULL
  }
  spec <- mobj$contrast_spec
  contrast_list <- mobj$contrast_list

  checkmate::assert_class(lmfit, "lm", null.ok=TRUE)
  checkmate::assert_list(spec)
  checkmate::assert_list(contrast_list)

  # prefix for contrast names
  prefix <- ifelse(inherits(mobj, "l1_wi_spec"), paste0(mobj$signal_name, "."), "")

  if (is.null(lmfit)) {
    c_colnames <- spec$regressors
  } else {
    c_colnames <- names(coef(lmfit)) # prefer model-specific regressors (if they vary by subject)
  }

  # for within-subject factors, always add the signal name as a prefix (cf. get_regressors_from_signal)
  if (nzchar(prefix) && !all(startsWith(c_colnames, prefix))) {
    c_colnames <- paste0(prefix, c_colnames)
  }

  c_custom <- align_custom_contrasts_to_colnames(preserved_custom, c_colnames)

  # Refresh cat_vars dynamically if lmfit is present and cat_vars was empty (e.g., ~ 0 + A:B or refit)
  if (!is.null(lmfit) && length(spec$cat_vars) == 0L) {
    spec$cat_vars <- extract_lm_var_types(lmfit)$cat_vars
    mobj$contrast_spec$cat_vars <- spec$cat_vars
  }

  # If no contrasts were specified in spec or contrast_list for L2/L3 models,
  # default to diagonal contrasts and warn the user.
  if (inherits(mobj, "hi_model_spec")) {
    has_list_contrasts <- any(vapply(contrast_list, function(x) {
      if (is.null(x)) return(FALSE)
      if (is.matrix(x)) return(nrow(x) > 0L)
      length(x) > 0L
    }, logical(1))) || (!is.null(c_custom) && nrow(c_custom) > 0L)

    has_spec_contrasts <- any(c(
      isTRUE(spec$diagonal),
      length(spec$cond_means) > 0L,
      length(spec$pairwise_diffs) > 0L,
      length(spec$interaction_diffs) > 0L,
      isTRUE(spec$cell_means),
      isTRUE(spec$overall_response),
      length(spec$simple_slopes) > 0L,
      length(spec$slope_diffs) > 0L
    ))

    if (!isTRUE(has_list_contrasts) && !isTRUE(has_spec_contrasts)) {
      warning(
        "No valid contrasts were specified; adding diagonal contrasts by default.",
        call. = FALSE
      )
      spec$diagonal <- TRUE
      mobj$contrast_spec$diagonal <- TRUE
    }
  }

  ### add diagonal contrasts
  c_diagonal <- contrast_list$diagonal
  if (isTRUE(spec$diagonal) && is.null(c_diagonal)) {
    diag_mat <- diag(length(c_colnames))
    # c_colnames already carries `prefix` for l1_wi_spec; avoid double-prefixing EV_<sig>.<sig>.<coef>
    rownames(diag_mat) <- paste0("EV_", c_colnames)
    colnames(diag_mat) <- c_colnames

    if (isTRUE(mobj$suppress_subject_ev_diagonals)) {
      keep_diag <- !grepl("^EV_id", rownames(diag_mat))
      diag_mat <- diag_mat[keep_diag, , drop = FALSE]
    }

    c_diagonal <- diag_mat
  }

  ### guard emmeans contrasts if lmfit is missing
  c_cond_means <- contrast_list$cond_means
  if (is.null(lmfit) && (length(spec$cond_means) > 0L || length(spec$pairwise_diffs) > 0L ||
    length(spec$interaction_diffs) > 0L || isTRUE(spec$cell_means) ||
    isTRUE(spec$overall_response) || length(spec$simple_slopes) > 0L || length(spec$slope_diffs) > 0L)) {
    warning(
      "emmeans-based contrasts requested, but no model is available. Dropping emmeans contrasts.",
      call. = FALSE
    )
    spec$cond_means <- character(0)
    spec$pairwise_diffs <- character(0)
    spec$interaction_diffs <- character(0)
    spec$cell_means <- FALSE
    spec$overall_response <- FALSE
    spec$simple_slopes <- list()
    spec$slope_diffs <- list()
  }

  eff_weights <- spec$weights %||% "cells"
  blocking_vars <- if (isTRUE(mobj$suppress_subject_ev_diagonals)) "id" else character(0)

  ### add condition means, if requested
  if (length(spec$cond_means) > 0L && is.null(c_cond_means)) {
    for (vv in spec$cond_means) {
      ps <- parse_emm_spec_string(vv)
      ee <- quiet_emm(emmeans(lmfit, specs = ps$focal_formula, by = ps$by_vars, weights = eff_weights, nesting = NULL))
      ee <- mask_empty_cell_marginal_slices(lmfit, ee, focal_vars = ps$focal_vars, by_vars = ps$by_vars, weights = eff_weights, label = vv, ignore_blocking_vars = blocking_vars)
      edata <- suppressWarnings(summary(ee))
      econ <- ee@linfct
      econ[abs(econ) < 1e-8] <- 0

      focal_vals <- apply(edata[, ps$focal_vars, drop = FALSE], 1, function(r) {
        paste(make.names(sanitize_contrast_token(r), unique = FALSE), collapse = ".")
      })
      by_pref <- format_by_prefix(edata, ps$by_vars)
      focal_label <- gsub("\\s+", "", ps$focal_raw)
      enames <- paste0(by_pref, paste(focal_label, focal_vals, sep = "."))

      which_na <- is.na(edata$emmean) | apply(econ, 1, anyNA)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      if (nrow(econ) > 0L) {
        rownames(econ) <- paste0(prefix, enames)
        c_cond_means <- rbind(c_cond_means, econ)
      }
    }
    if (!is.null(c_cond_means) && nrow(c_cond_means) > 0L) {
      colnames(c_cond_means) <- c_colnames
    }
  }

  ### add cell means for all factors, if requested
  c_cell_means <- contrast_list$cell_means
  if (isTRUE(spec$cell_means) && is.null(c_cell_means) && length(spec$cat_vars) > 0L) {
    ee <- quiet_emm(emmeans(
      lmfit,
      as.formula(paste("~", paste(spec$cat_vars, collapse = "*"))),
      weights = eff_weights,
      nesting = NULL
    ))
    edata <- suppressWarnings(summary(ee))

    econ <- ee@linfct
    econ[abs(econ) < 1e-8] <- 0
    enames <- apply(edata[, spec$cat_vars, drop = FALSE], 1, function(x) {
      vals <- make.names(sanitize_contrast_token(x), unique = FALSE)
      paste(names(x), vals, sep = ".", collapse = "_")
    })

    which_na <- is.na(edata$emmean) | apply(econ, 1, anyNA)
    if (any(which_na)) {
      econ <- econ[!which_na, , drop = FALSE]
      enames <- enames[!which_na]
    }

    if (nrow(econ) > 0L) {
      rownames(econ) <- paste0(prefix, enames)
      c_cell_means <- rbind(c_cell_means, econ)
      colnames(c_cell_means) <- c_colnames
    }
  }

  ### add overall response mean
  c_overall <- contrast_list$overall
  if (isTRUE(spec$overall_response) && is.null(c_overall)) {
    ee <- quiet_emm(emmeans(lmfit, ~1, weights = eff_weights, nesting = NULL))
    ee <- mask_empty_cell_marginal_slices(lmfit, ee, focal_vars = character(0), by_vars = NULL, weights = eff_weights, label = "overall_response", ignore_blocking_vars = blocking_vars)
    edata <- suppressWarnings(summary(ee))
    econ <- ee@linfct
    econ[abs(econ) < 1e-8] <- 0
    which_na <- is.na(edata$emmean) | apply(econ, 1, anyNA)
    if (!all(which_na)) {
      econ <- econ[!which_na, , drop = FALSE]
      rownames(econ) <- paste0(prefix, "overall")
      colnames(econ) <- c_colnames
      c_overall <- econ
    }
  }

  ### add pairwise differences, if requested
  c_pairwise_diffs <- contrast_list$pairwise_diffs
  if (length(spec$pairwise_diffs) > 0L && is.null(c_pairwise_diffs)) {
    for (vv in spec$pairwise_diffs) {
      ps <- parse_emm_spec_string(vv)
      ee <- quiet_emm(emmeans(lmfit, specs = ps$focal_formula, by = ps$by_vars, weights = eff_weights, nesting = NULL))
      ee <- mask_empty_cell_marginal_slices(lmfit, ee, focal_vars = ps$focal_vars, by_vars = ps$by_vars, weights = eff_weights, label = vv, ignore_blocking_vars = blocking_vars)
      pp <- pairs(ee)
      edata <- suppressWarnings(summary(pp))
      econ <- pp@linfct
      econ[abs(econ) < 1e-8] <- 0
      by_pref <- format_by_prefix(edata, ps$by_vars)
      focal_label <- gsub("\\s+", "", ps$focal_raw)
      enames <- paste0(by_pref, paste(focal_label, sanitize_contrast_diff_label(edata$contrast), sep = "."))

      which_na <- is.na(edata$estimate) | apply(econ, 1, anyNA)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      if (nrow(econ) > 0L) {
        rownames(econ) <- paste0(prefix, enames)
        c_pairwise_diffs <- rbind(c_pairwise_diffs, econ)
      }
    }
    if (!is.null(c_pairwise_diffs) && nrow(c_pairwise_diffs) > 0L) {
      colnames(c_pairwise_diffs) <- c_colnames
    }
  }

  ### add factorial interaction (difference-of-differences) contrasts, if requested
  c_interaction_diffs <- contrast_list$interaction_diffs
  if (length(spec$interaction_diffs) > 0L && is.null(c_interaction_diffs)) {
    for (vv in spec$interaction_diffs) {
      ps <- parse_emm_spec_string(vv)
      if (length(ps$focal_vars) < 2L) {
        warning(sprintf("interaction_diffs entry '%s' has fewer than 2 interacting factors; skipping.", vv), call. = FALSE)
        next
      }
      ee <- quiet_emm(emmeans(lmfit, specs = ps$focal_formula, by = ps$by_vars, weights = eff_weights, nesting = NULL))
      ee <- mask_empty_cell_marginal_slices(lmfit, ee, focal_vars = ps$focal_vars, by_vars = ps$by_vars, weights = eff_weights, label = vv, ignore_blocking_vars = blocking_vars)
      ic <- emmeans::contrast(ee, interaction = "pairwise")
      edata <- suppressWarnings(summary(ic))
      econ <- ic@linfct
      econ[abs(econ) < 1e-8] <- 0

      pw_cols <- grep("_pairwise$", names(edata), value = TRUE)
      if (length(pw_cols) == 0L) pw_cols <- intersect(names(edata), paste0(ps$focal_vars, "_pairwise"))
      int_labels <- apply(edata[, pw_cols, drop = FALSE], 1, function(r) {
        paste(sanitize_contrast_diff_label(r), collapse = "_x_")
      })
      by_pref <- format_by_prefix(edata, ps$by_vars)
      focal_tag <- paste(ps$focal_vars, collapse = "_x_")
      enames <- paste0(by_pref, focal_tag, ".", int_labels)

      which_na <- is.na(edata$estimate) | apply(econ, 1, anyNA)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      if (nrow(econ) > 0L) {
        rownames(econ) <- paste0(prefix, enames)
        c_interaction_diffs <- rbind(c_interaction_diffs, econ)
      }
    }
    if (!is.null(c_interaction_diffs) && nrow(c_interaction_diffs) > 0L) {
      colnames(c_interaction_diffs) <- c_colnames
    }
  }

  ### add per-cell model-predicted simple slopes
  c_simple_slopes <- contrast_list$simple_slopes
  if (length(spec$simple_slopes) > 0L && is.null(c_simple_slopes)) {
    for (vv in seq_along(spec$simple_slopes)) {
      trend_var <- names(spec$simple_slopes)[vv]
      comb_list <- spec$simple_slopes[[vv]]
      if (is.character(comb_list)) comb_list <- as.list(comb_list)
      for (comb in comb_list) {
        comb_vars <- trimws(unlist(strsplit(as.character(comb), "[:*+]", perl = TRUE)))
        comb_vars <- comb_vars[nzchar(comb_vars)]
        ee <- quiet_emm(emtrends(
          lmfit,
          specs = as.formula(paste("~", paste(comb_vars, collapse = "*"))),
          var = trend_var,
          weights = eff_weights,
          nesting = NULL
        ))
        ee <- mask_empty_cell_marginal_slices(lmfit, ee, focal_vars = comb_vars, by_vars = NULL, weights = eff_weights, label = paste0(trend_var, ":", paste(comb_vars, collapse = "*")), ignore_blocking_vars = blocking_vars)
        edata <- suppressWarnings(summary(ee))
        econ <- ee@linfct
        econ[abs(econ) < 1e-8] <- 0
        enames <- paste(
          trend_var,
          apply(edata[, comb_vars, drop = FALSE], 1, function(r) paste(make.names(sanitize_contrast_token(r), unique = FALSE), collapse = ".")),
          sep = "."
        )

        trend_col <- paste(trend_var, "trend", sep = ".")
        which_na <- is.na(edata[[trend_col]]) | apply(econ, 1, anyNA)
        if (any(which_na)) {
          econ <- econ[!which_na, , drop = FALSE]
          enames <- enames[!which_na]
        }

        if (nrow(econ) > 0L) {
          rownames(econ) <- paste0(prefix, enames)
          c_simple_slopes <- rbind(c_simple_slopes, econ)
        }
      }
    }
    if (!is.null(c_simple_slopes) && nrow(c_simple_slopes) > 0L) {
      colnames(c_simple_slopes) <- c_colnames
    }
  }

  ### add pairwise / interaction differences in simple slopes
  c_slope_diffs <- contrast_list$slope_diffs
  if (length(spec$slope_diffs) > 0L && is.null(c_slope_diffs)) {
    for (vv in seq_along(spec$slope_diffs)) {
      trend_var <- names(spec$slope_diffs)[vv]
      comb_list <- spec$slope_diffs[[vv]]
      if (is.character(comb_list)) comb_list <- as.list(comb_list)
      for (comb in comb_list) {
        comb_str <- paste(as.character(comb), collapse = ":")
        ps <- parse_emm_spec_string(comb_str)
        ee <- quiet_emm(emtrends(
          lmfit,
          specs = ps$focal_formula,
          by = ps$by_vars,
          var = trend_var,
          weights = eff_weights,
          nesting = NULL
        ))
        ee <- mask_empty_cell_marginal_slices(lmfit, ee, focal_vars = ps$focal_vars, by_vars = ps$by_vars, weights = eff_weights, label = paste0(trend_var, ":", comb_str), ignore_blocking_vars = blocking_vars)
        if (length(ps$focal_vars) >= 2L && is.null(ps$by_vars)) {
          pp <- emmeans::contrast(ee, interaction = "pairwise")
          edata <- suppressWarnings(summary(pp))
          pw_cols <- grep("_pairwise$", names(edata), value = TRUE)
          diff_lbl <- apply(edata[, pw_cols, drop = FALSE], 1, function(r) paste(sanitize_contrast_diff_label(r), collapse = "_x_"))
          enames <- paste(trend_var, paste(ps$focal_vars, collapse = "_x_"), diff_lbl, sep = ".")
        } else {
          pp <- pairs(ee)
          edata <- suppressWarnings(summary(pp))
          by_pref <- format_by_prefix(edata, ps$by_vars)
          enames <- paste0(by_pref, paste(trend_var, paste(ps$focal_vars, collapse = "_"), sanitize_contrast_diff_label(edata$contrast), sep = "."))
        }
        econ <- pp@linfct
        econ[abs(econ) < 1e-8] <- 0
        which_na <- is.na(edata$estimate) | apply(econ, 1, anyNA)
        if (any(which_na)) {
          econ <- econ[!which_na, , drop = FALSE]
          enames <- enames[!which_na]
        }
        if (nrow(econ) > 0L) {
          rownames(econ) <- paste0(prefix, enames)
          c_slope_diffs <- rbind(c_slope_diffs, econ)
        }
      }
    }
    if (!is.null(c_slope_diffs) && nrow(c_slope_diffs) > 0L) {
      colnames(c_slope_diffs) <- c_colnames
    }
  }

  # walk through all within-subject factor contrasts and return as a single matrix
  c_within <- get_wi_contrast_matrix(mobj, c_colnames)

  # combine each element
  cmat_full <- rbind(
    c_diagonal,
    c_cond_means,
    c_cell_means,
    c_overall,
    c_pairwise_diffs,
    c_interaction_diffs,
    c_simple_slopes,
    c_slope_diffs,
    c_within,
    c_custom
  )

  if (is.null(cmat_full)) {
    # no contrasts specified! Not good, but for now, just create an empty matrix and let upstream handle bothering the user
    cmat_full <- matrix(numeric(0), nrow=0, ncol=length(c_colnames), dimnames=list(NULL, c_colnames))
  } else if (!is.null(rownames(cmat_full)) && anyDuplicated(rownames(cmat_full))) {
    rownames(cmat_full) <- make.unique(rownames(cmat_full), sep = "_")
  }

  # handle contrasts from the full matrix that should be deleted
  # these are processed dynamically here so that the same set can be dropped if a contrast is respecified for a data subset or subject
  if (!is.null(spec$delete)) {
    if (!is.character(spec$delete)) {
      warning("contrast_spec$delete must be a character vector; ignoring invalid values.", call. = FALSE)
      spec$delete <- NULL
    } else {
      spec$delete <- spec$delete[!is.na(spec$delete) & nzchar(spec$delete)]
      if (length(spec$delete) == 0L) spec$delete <- NULL
    }
  }

  if (!is.null(spec$delete)) {
    which_del <- match(spec$delete, rownames(cmat_full))
    if (!any(is.na(which_del))) {
      cmat <- cmat_full[-1*which_del, , drop=FALSE]
    } else {
      if (!isTRUE(is_refit)) {
        warning("Some contrasts to be dropped were not found in the matrix: ", paste(spec$delete[is.na(which_del)], collapse = ", "))
      }
      which_del <- na.omit(which_del)
      cmat <- if (length(which_del) > 0L) cmat_full[-1 * na.omit(which_del), , drop = FALSE] else cmat_full
    }
  } else {
    cmat <- cmat_full
  }

  # check for duplicated contrast rows and resolve naming collisions.
  # In non-interactive contexts (or automated per-subject/L3 refits), retain the first row
  # deterministically and drop the rest without modifying global spec$delete on refits.
  dupe_list <- get_dupe_rows(cmat)
  drop_rows <- c()
  if (length(dupe_list) > 0L) {
    is_interactive <- isTRUE(allow_interactive_dedupe) && rlang::is_interactive()
    if (is.null(rownames(cmat))) {
      rownames(cmat) <- paste0("contrast_", seq_len(nrow(cmat)))
    }
    if (isTRUE(is_interactive)) {
      message("Duplicate contrasts found in your matrix. This can occur for many benign reasons, but we need you to say what you want
      to call these contrasts in the output so that it is clear to you.\n")
    } else if (!isTRUE(is_refit)) {
      warning(
        "Duplicate contrasts found in non-interactive mode; keeping the first name from each duplicate set.",
        call. = FALSE
      )
    }
    
    for (dd in dupe_list) {
      if (!isTRUE(is_refit)) {
        message("These contrasts are the same: ", paste(rownames(cmat)[dd], collapse=", "))
      }
      if (isTRUE(is_interactive)) {
        whichkeep <- menu(rownames(cmat)[dd], title = "Which name should we use?")
        if (!is.numeric(whichkeep) || length(whichkeep) != 1L || is.na(whichkeep) ||
            whichkeep < 1L || whichkeep > length(dd)) {
          whichkeep <- 1L
          warning(
            "Invalid duplicate-contrast selection; defaulting to first contrast name: ",
            rownames(cmat)[dd][whichkeep],
            call. = FALSE
          )
        }
      } else {
        whichkeep <- 1L
        if (!isTRUE(is_refit)) {
          message(
            "Non-interactive session detected while resolving duplicate contrasts; keeping '",
            rownames(cmat)[dd][whichkeep],
            "' and dropping: ",
            paste(rownames(cmat)[dd][-whichkeep], collapse = ", ")
          )
        }
      }
      drop_rows <- c(drop_rows, dd[-whichkeep])
    }

    if (!isTRUE(is_refit)) {
      # add these back into model specification so that the same decisions are made if the spec is used to repopulate a new matrix
      mobj$contrast_spec$delete <- unique(c(mobj$contrast_spec$delete, rownames(cmat)[drop_rows]))
    }
  } 

  dupes <- rep(FALSE, nrow(cmat))
  dupes[drop_rows] <- TRUE
  
  cmat <- cmat[!dupes, , drop=FALSE] # drop duplicated contrasts

  # parentheses in contrast names generate problems in running through FEAT (chokes on submission)
  if (!is.null(cmat) && !is.null(rownames(cmat))) {
    rownames(cmat) <- gsub("(Intercept)", "Intercept", rownames(cmat), fixed = TRUE)
    if (anyDuplicated(rownames(cmat))) {
      rownames(cmat) <- make.unique(rownames(cmat), sep = "_")
    }
  }

  # populate any updates to the contrast_list object based on new calculations
  contrast_list$diagonal <- c_diagonal
  contrast_list$cond_means <- c_cond_means
  contrast_list$cell_means <- c_cell_means
  contrast_list$overall <- c_overall
  contrast_list$pairwise_diffs <- c_pairwise_diffs
  contrast_list$interaction_diffs <- c_interaction_diffs
  contrast_list$simple_slopes <- c_simple_slopes
  contrast_list$slope_diffs <- c_slope_diffs
  contrast_list$custom <- c_custom

  # populate contrasts back into model object (spec should not be updated)
  mobj$contrast_list <- contrast_list
  mobj$contrasts <- cmat
  return(mobj)
}

# helper function to identify duplicate contrasts
get_dupe_rows <- function(mat, enforce_rownames=FALSE) {
  if (!inherits(mat, "matrix")) return(list()) # return empty list if we are not passed a contrast matrix
  
  # rows to check -- start with all, then winnow
  rtc <- seq_len(nrow(mat))
  
  # if enforce_rownames is FALSE, drop rownames from local copy so that identical does not fail
  if (isFALSE(enforce_rownames)) rownames(mat) <- NULL
  
  dupe_list<- list()
  while(length(rtc) > 1L) {
    last_rtc <- length(rtc)
    to_compare <- mat[rtc[last_rtc],,drop=FALSE] # always compare preceding rows to last remaining row
    rdupes <- sapply(1:(length(rtc) - 1), function(i) identical(mat[rtc[i],,drop=FALSE], to_compare))
    to_drop <- last_rtc # always drop last eligible row after checking for dupes of it
    if (any(rdupes)) {
      # rdupes will always be 1 less than length of rtc. Add FALSE as last element so that rtc subset never gets scalar rtc[TRUE], which would grab all elements
      dupe_list <- c(dupe_list, list(c(rtc[c(rdupes, FALSE)], rtc[last_rtc])))
      to_drop <- c(to_drop, which(rdupes))
    }
    rtc <- rtc[-to_drop] # drop rows that have been processed
  }
  
  return(dupe_list)
}



get_wi_contrast_matrix <- function(mobj, c_colnames) {
  # this function is only intended to crawl over within-subject contrasts of an l1 model
  if (!checkmate::test_class(mobj, "l1_model_spec") || is.null(mobj$wi_models)) {
    return(invisible(NULL)) # quietly return NULL if we don't have an l1 model
  }

  # get contrast matrix for each wi_model
  clist <- lapply(mobj$wi_models, "[[", "contrasts")
  clist <- Filter(function(cc) is.matrix(cc) && nrow(cc) > 0L, clist)
  if (length(clist) == 0L) return(NULL)

  cnames <- make.unique(do.call(c, lapply(clist, rownames)), sep = "_")
  cmat <- matrix(0, nrow = length(cnames), ncol = length(c_colnames), dimnames = list(cnames, c_colnames))

  # fill by row position (not rowname) so duplicate raw labels can never overwrite or leave all-zero rows
  row_offset <- 0L
  for (cc in clist) {
    idx <- row_offset + seq_len(nrow(cc))
    shared_cols <- intersect(colnames(cc), c_colnames)
    if (length(shared_cols) > 0L) {
      cmat[idx, shared_cols] <- cc[, shared_cols, drop = FALSE]
    }
    row_offset <- row_offset + nrow(cc)
  }

  return(cmat)
}

# get_run_lmfits <- function() {

# }

#' Build run-level L2 model data augmented with session-level predictors
#'
#' @param gpa A \code{glm_pipeline_arguments} object.
#' @param lg Optional logger.
#' @return A data.frame suitable for L2 model fitting/respecification.
#' @keywords internal
#' @importFrom stats setNames
compose_l2_model_data <- function(gpa, lg = NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(gpa$run_data)

  if (is.null(lg)) {
    lg <- lgr::get_logger("glm_pipeline/l2_setup")
  }

  l2_data <- gpa$run_data
  origin <- setNames(rep("run", ncol(l2_data)), names(l2_data))

  if (is.null(gpa$session_data)) {
    attr(l2_data, "l2_var_origin") <- origin
    return(l2_data)
  }

  checkmate::assert_data_frame(gpa$session_data)
  checkmate::assert_subset(c("id", "session"), names(gpa$session_data))

  session_dupes <- gpa$session_data %>%
    dplyr::count(.data$id, .data$session, name = "n") %>%
    dplyr::filter(.data$n > 1L)
  if (nrow(session_dupes) > 0L) {
    stop("session_data must have unique id/session rows", call. = FALSE)
  }

  sdat <- gpa$session_data
  session_cols <- setdiff(names(sdat), c("id", "session"))
  if (length(session_cols) == 0L) {
    attr(l2_data, "l2_var_origin") <- origin
    return(l2_data)
  }

  conflicts <- intersect(session_cols, names(l2_data))
  redundant_conflicts <- character(0)
  rename_map <- character(0)
  if (length(conflicts) > 0L) {
    columns_equivalent <- function(run_col, session_col) {
      na_match <- is.na(run_col) & is.na(session_col)
      keep <- !na_match
      if (!any(keep)) return(TRUE)

      run_keep <- run_col[keep]
      session_keep <- session_col[keep]

      if (is.factor(run_keep) || is.factor(session_keep) ||
          is.character(run_keep) || is.character(session_keep)) {
        return(identical(as.character(run_keep), as.character(session_keep)))
      }

      isTRUE(all.equal(run_keep, session_keep, check.attributes = FALSE))
    }

    for (cc in conflicts) {
      cmp <- dplyr::left_join(
        l2_data %>% dplyr::select("id", "session", run_value = dplyr::all_of(cc)),
        sdat %>% dplyr::select("id", "session", session_value = dplyr::all_of(cc)),
        by = c("id", "session")
      )

      if (columns_equivalent(cmp$run_value, cmp$session_value)) {
        redundant_conflicts <- c(redundant_conflicts, cc)
      }
    }

    if (length(redundant_conflicts) > 0L) {
      sdat <- sdat %>% dplyr::select(-dplyr::all_of(redundant_conflicts))
      session_cols <- setdiff(session_cols, redundant_conflicts)
      conflicts <- setdiff(conflicts, redundant_conflicts)
      lg$info(
        "Skipping session-level columns already present in run_data with identical values: %s",
        paste(redundant_conflicts, collapse = ", ")
      )
    }
  }

  if (length(session_cols) == 0L) {
    attr(l2_data, "l2_var_origin") <- origin
    return(l2_data)
  }

  if (length(conflicts) > 0L) {
    existing <- names(l2_data)
    for (cc in conflicts) {
      base_name <- paste0("session_", cc)
      candidate <- base_name
      idx <- 1L
      while (candidate %in% c(existing, session_cols, unname(rename_map))) {
        candidate <- paste0(base_name, "_", idx)
        idx <- idx + 1L
      }
      rename_map[cc] <- candidate
      existing <- c(existing, candidate)
    }

    for (cc in names(rename_map)) {
      names(sdat)[names(sdat) == cc] <- rename_map[[cc]]
    }

    lg$warn(
      "Renaming session-level columns to avoid collisions with run-level columns: %s",
      paste(sprintf("%s->%s", names(rename_map), unname(rename_map)), collapse = ", ")
    )
  }

  l2_data <- dplyr::left_join(l2_data, sdat, by = c("id", "session"))
  new_cols <- setdiff(names(l2_data), names(origin))
  if (length(new_cols) > 0L) {
    origin[new_cols] <- "session"
  }

  attr(l2_data, "l2_var_origin") <- origin
  if (length(rename_map) > 0L) attr(l2_data, "session_rename_map") <- rename_map
  l2_data
}

deparse_formula_str <- function(f) {
  form <- paste(deparse(f, width.cutoff = 500L), collapse = " ")
  gsub("\\s+", " ", form, perl = FALSE)
}

apply_factor_coding_to_col <- function(fac, coding) {
  if (!is.factor(fac) || nlevels(fac) < 2L || is.null(coding) || !nzchar(coding)) {
    return(fac)
  }
  coding <- tolower(trimws(as.character(coding)[1L]))
  cfun <- switch(coding,
    treatment = stats::contr.treatment,
    contr.treatment = stats::contr.treatment,
    sum = stats::contr.sum,
    contr.sum = stats::contr.sum,
    helmert = stats::contr.helmert,
    contr.helmert = stats::contr.helmert,
    poly = stats::contr.poly,
    contr.poly = stats::contr.poly,
    stop("Unsupported factor_coding '", coding, "'. Allowed: treatment, sum, helmert, poly.", call. = FALSE)
  )
  stats::contrasts(fac) <- cfun(levels(fac))
  fac
}

fit_lm_with_canonical_levels <- function(model_formula, model_data, canonical_levels = NULL, factor_coding = NULL) {
  rhs_vars <- setdiff(all.vars(stats::delete.response(stats::terms(model_formula))), "dummy")
  fac_vars <- rhs_vars[vapply(rhs_vars, function(v) v %in% names(model_data) && is.factor(model_data[[v]]), logical(1))]

  needs_zero_weight_pad <- FALSE
  if (length(fac_vars) > 0L) {
    for (fv in fac_vars) {
      obs_levs <- unique(as.character(stats::na.omit(model_data[[fv]])))
      all_levs <- canonical_levels[[fv]] %||% levels(model_data[[fv]])
      if (length(obs_levs) < 2L || !all(all_levs %in% obs_levs)) {
        needs_zero_weight_pad <- TRUE
        break
      }
    }
  }

  f_clean <- stats::as.formula(deparse_formula_str(model_formula), env = environment())
  if (!isTRUE(needs_zero_weight_pad)) {
    lmfit <- stats::lm(f_clean, data = model_data)
    return(list(lmfit = lmfit, obs_idx = seq_len(nrow(model_data)), padded = FALSE))
  }

  # Build zero-weight grid over canonical factor levels so stats::model.matrix never drops
  # unused levels or errors on 1-level factors, while keeping exact OLS estimates on observed rows.
  lev_list <- stats::setNames(lapply(fac_vars, function(fv) {
    canonical_levels[[fv]] %||% levels(model_data[[fv]])
  }), fac_vars)
  pad_grid <- expand.grid(lev_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  obs_n <- nrow(model_data)
  # emmeans needs rows representing every factor level, even when their weights
  # are zero. Its continuous-covariate reductions are unweighted, so copy each
  # observed covariate row equally often rather than inventing zeros or means.
  # This preserves both the empirical distribution and discrete covariate values.
  pad_n <- obs_n * ceiling(nrow(pad_grid) / obs_n)
  pad_df <- model_data[rep(seq_len(obs_n), length.out = pad_n), , drop = FALSE]
  grid_idx <- rep(seq_len(nrow(pad_grid)), length.out = pad_n)
  for (fv in fac_vars) {
    pad_df[[fv]] <- factor(pad_grid[[fv]][grid_idx], levels = lev_list[[fv]])
  }
  if ("dummy" %in% names(pad_df)) pad_df$dummy <- 0

  df_aug <- rbind(model_data, pad_df)
  for (fv in fac_vars) {
    df_aug[[fv]] <- factor(as.character(df_aug[[fv]]), levels = lev_list[[fv]])
    if (!is.null(factor_coding) && length(factor_coding) > 0L) {
      code_this <- if (!is.null(names(factor_coding)) && fv %in% names(factor_coding)) {
        factor_coding[[fv]]
      } else if (is.null(names(factor_coding)) && length(factor_coding) == 1L) {
        factor_coding[[1L]]
      } else {
        NULL
      }
      df_aug[[fv]] <- apply_factor_coding_to_col(df_aug[[fv]], code_this)
    }
  }
  df_aug$.obs_weight_internal <- c(rep(1, obs_n), rep(0, nrow(pad_df)))

  lmfit <- suppressWarnings(stats::lm(f_clean, data = df_aug, weights = .obs_weight_internal))
  list(lmfit = lmfit, obs_idx = seq_len(obs_n), padded = TRUE)
}

prune_rank_deficient_design_and_contrasts <- function(X_obs_full, coef_vec, cmat, lg = NULL, context = NULL) {
  if (nrow(X_obs_full) == 0L) {
    empty_dmat <- matrix(numeric(0), nrow = 0L, ncol = 0L)
    empty_cmat <- matrix(numeric(0), nrow = 0L, ncol = 0L)
    return(list(dmat = empty_dmat, cmat = empty_cmat, dropped_contrasts = rownames(cmat) %||% character(0)))
  }

  keep_cols <- !is.na(coef_vec)
  if (any(keep_cols)) {
    X_sub <- X_obs_full[, keep_cols, drop = FALSE]
    qr_sub <- qr(X_sub)
    if (qr_sub$rank < ncol(X_sub)) {
      piv_keep <- qr_sub$pivot[seq_len(qr_sub$rank)]
      sub_names <- colnames(X_sub)[sort(piv_keep)]
      keep_cols <- colnames(X_obs_full) %in% sub_names
    }
  }

  if (all(keep_cols) && qr(X_obs_full)$rank == ncol(X_obs_full)) {
    return(list(dmat = X_obs_full, cmat = cmat, dropped_contrasts = character(0)))
  }

  X_reduced <- X_obs_full[, keep_cols, drop = FALSE]
  attr(X_reduced, "model_runs") <- attr(X_obs_full, "model_runs")

  if (is.null(cmat) || nrow(cmat) == 0L) {
    empty_c <- matrix(numeric(0), nrow = 0L, ncol = ncol(X_reduced), dimnames = list(NULL, colnames(X_reduced)))
    return(list(dmat = X_reduced, cmat = empty_c, dropped_contrasts = character(0)))
  }

  # Align contrast columns to X_obs_full if needed
  if (!identical(colnames(cmat), colnames(X_obs_full))) {
    cmat_aligned <- matrix(0, nrow = nrow(cmat), ncol = ncol(X_obs_full), dimnames = list(rownames(cmat), colnames(X_obs_full)))
    shared <- intersect(colnames(cmat), colnames(X_obs_full))
    if (length(shared) > 0L) cmat_aligned[, shared] <- cmat[, shared, drop = FALSE]
    cmat <- cmat_aligned
  }

  nb <- estimability::nonest.basis(X_obs_full)
  estble <- vapply(seq_len(nrow(cmat)), function(rr) {
    row_vec <- as.numeric(cmat[rr, ])
    if (anyNA(row_vec) || all(abs(row_vec) <= 1e-12)) return(FALSE)
    isTRUE(estimability::is.estble(row_vec, nb, tol = 1e-5))
  }, logical(1))

  dropped <- if (!is.null(rownames(cmat))) rownames(cmat)[!estble] else paste0("contrast_", which(!estble))
  if (length(dropped) > 0L && !is.null(lg)) {
    lg$warn(
      "Pruning %d non-estimable contrast(s)%s due to missing/aliased design columns: %s",
      length(dropped),
      if (!is.null(context) && nzchar(context)) paste0(" in ", context) else "",
      paste(dropped, collapse = ", ")
    )
  }

  cmat_reduced <- cmat[estble, keep_cols, drop = FALSE]
  list(dmat = X_reduced, cmat = cmat_reduced, dropped_contrasts = dropped)
}

#' Recalculate an L2 model by grouping units derived from available run data
#'
#' @param mobj an \code{l1_model_spec} or \code{hi_model_spec} object containing the GLM model to run
#' @param data The run-level data frame containing data for all ids and sessions. This is split into chunks by \code{split_on}.
#' @param split_on Character vector defining grouping columns used to respecify the model. Must include \code{"id"} and may include
#'   \code{"session"} (default: \code{c("id", "session")}).
#' @param aggregated_session Integer session value assigned when grouping across sessions (i.e., when \code{split_on = "id"}). This keeps
#'   output schemas consistent for downstream joins and cope metadata. Default is \code{0L}.
#'
#' @return a modified copy of \code{mobj} where the \code{$by_subject} field has been added
#'
#' @details The function adds the \code{$by_subject} field, which contains the design matrices and contrasts
#'   for each grouping unit in \code{data}, based on the runs that remain available for that unit. For example, if
#'   a subject is missing runs (or they are dropped from analysis), some contrasts may change or drop out of the model.
#'
#' The \code{$by_subject} field is a keyed data.table containing list elements for \code{cope_list}
#'   (mapping cope numbers to contrast names), \code{contrasts}, and the design matrix for each unit.
#'   It also includes \code{grouping_scope} (\code{"id_session"} or \code{"id"}).
#'
#' @keywords internal
#' @importFrom checkmate assert_data_frame assert_multi_class assert_subset
#' @importFrom data.table data.table
#' @importFrom stats model.matrix terms
respecify_l2_models_by_subject <- function(mobj, data, split_on = c("id", "session"), aggregated_session = 0L) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure
  checkmate::assert_data_frame(data)
  checkmate::assert_subset(c("id", "session"), names(data))
  checkmate::assert_character(split_on, min.len = 1L, any.missing = FALSE)
  checkmate::assert_subset(split_on, c("id", "session"))
  if (!"id" %in% split_on) {
    stop("split_on must include 'id' for respecify_l2_models_by_subject", call. = FALSE)
  }
  checkmate::assert_integerish(aggregated_session, len = 1L)
  split_on <- unique(split_on)
  lg <- lgr::get_logger("glm_pipeline/l2_setup")

  # create a nested data.table object for fitting each grouping unit
  data <- data.table(data, key = split_on)
  data$dummy <- seq_len(nrow(data))
  dsplit <- data[, .(dt = list(.SD)), by = split_on]

  if (!"session" %in% names(dsplit)) {
    dsplit$session <- as.integer(aggregated_session)
  }

  # use model formula from parent object
  model_formula <- terms(mobj$lmfit)
  canonical_levels <- mobj$factor_levels
  if (is.null(canonical_levels) && !is.null(mobj$model_data)) {
    fac_cols <- names(mobj$model_data)[vapply(mobj$model_data, is.factor, logical(1))]
    if (length(fac_cols) > 0L) {
      canonical_levels <- stats::setNames(lapply(fac_cols, function(fc) levels(mobj$model_data[[fc]])), fac_cols)
    }
  }

  cope_list <- list()
  contrast_list <- list()
  model_matrix_list <- list()
  n_l2_copes <- rep(NA_integer_, nrow(dsplit))
  for (vv in seq_len(nrow(dsplit))) {
    model_data <- as.data.frame(dsplit[[vv, "dt"]])
    # Restore grouping columns removed from .SD by data.table `by = split_on`
    for (sc in split_on) {
      if (!sc %in% names(model_data)) model_data[[sc]] <- dsplit[[sc]][vv]
    }
    # Match initial specification: remove incomplete observations before computing
    # transformations, then apply settings within the retained subject/session.
    complete <- stats::complete.cases(model_data[, all.vars(model_formula), drop = FALSE])
    model_data <- model_data[complete, , drop = FALSE]

    if (nrow(model_data) == 0L) {
      cope_list[[vv]] <- data.frame(
        id = character(0), session = integer(0),
        l2_cope_number = integer(0), l2_cope_name = character(0)
      )
      contrast_list[[vv]] <- matrix(numeric(0), nrow = 0L, ncol = 0L)
      model_matrix_list[[vv]] <- matrix(numeric(0), nrow = 0L, ncol = 0L)
      n_l2_copes[vv] <- 0L
      next
    }

    run_meta_cols <- intersect(c("id", "session", "run_number"), names(model_data))
    kept_run_meta <- if (length(run_meta_cols) > 0L) model_data[, run_meta_cols, drop = FALSE] else NULL

    model_data <- apply_model_covariate_settings(model_data, mobj, lg)
    fit_res <- fit_lm_with_canonical_levels(
      model_formula = model_formula,
      model_data = model_data,
      canonical_levels = canonical_levels,
      factor_coding = mobj$factor_coding %||% mobj$contrast_spec$factor_coding
    )
    lmfit <- fit_res$lmfit

    mm <- get_contrasts_from_spec(mobj, lmfit, allow_interactive_dedupe = FALSE)
    X_full <- model.matrix(lmfit)[fit_res$obs_idx, , drop = FALSE]
    attr(X_full, "model_runs") <- kept_run_meta

    pruned <- prune_rank_deficient_design_and_contrasts(
      X_obs_full = X_full,
      coef_vec = stats::coef(lmfit),
      cmat = mm$contrasts,
      lg = lg,
      context = sprintf("L2 subject %s session %s", dsplit$id[vv], dsplit$session[vv])
    )
    dmat_unit <- pruned$dmat
    cmat_unit <- pruned$cmat

    n_contrasts <- if (is.null(cmat_unit)) 0L else nrow(cmat_unit)
    contrast_names <- if (n_contrasts > 0L) rownames(cmat_unit) else character(0)
    if (is.null(contrast_names)) {
      contrast_names <- if (n_contrasts > 0L) paste0("contrast_", seq_len(n_contrasts)) else character(0)
    }
    cope_df <- data.frame(
      id = rep(dsplit$id[vv], n_contrasts),
      session = rep(dsplit$session[vv], n_contrasts),
      l2_cope_number = if (n_contrasts > 0L) seq_len(n_contrasts) else integer(0),
      l2_cope_name = contrast_names
    )

    cope_list[[vv]] <- cope_df
    contrast_list[[vv]] <- cmat_unit
    model_matrix_list[[vv]] <- dmat_unit
    n_l2_copes[vv] <- n_contrasts
  }

  dsplit[, cope_list := list(cope_list)]
  dsplit[, contrasts := list(contrast_list)]
  dsplit[, model_matrix := list(model_matrix_list)]
  dsplit[, n_l2_copes := n_l2_copes]
  dsplit[, grouping_scope := ifelse("session" %in% split_on, "id_session", "id")]
  dsplit[, dt := NULL] # no longer need original split data

  mobj$by_subject <- dsplit
  return(mobj)
}

#' Helper function to re-calculate the l3 design matrix for a model based on available data
#'
#' @param mobj a \code{hi_model_spec} object containing the L3 GLM model to run
#' @param new_data The run-level data frame containing data for all ids and sessions. This will be split into individual chunks
#'
#' @return a modified copy of \code{mobj} where the $by_subject field has been added
#'
#' @details The function adds the $by_subject field, which contains the design matrices and contrasts
#'   for each subject and session in \code{data} based on the available data for that session. For example, if
#'   a subject is missing a few runs (or these are dropped from analysis), then some contrasts may change or drop out
#'   of the model.
#'
#' The $by_subject field is a keyed data.table object containing list elements for the cope_list (mapping cope numbers
#'   to contrast names), the contrasts, and the design matrix for each session.
#'
#' @keywords internal
#' @importFrom checkmate assert_data_frame assert_multi_class assert_subset
#' @importFrom data.table data.table

mobj_refit_lm <- function(mobj, new_data) {

}

#' Format a clear diagnostic for aliased higher-level model coefficients
#'
#' @param mobj a model object populated by \code{mobj_fit_lm}
#' @param model_formula optional formula string used for display
#' @return a character string suitable for an error or interactive prompt
#' @keywords internal
format_hi_model_alias_message <- function(mobj, model_formula = NULL) {
  checkmate::assert_class(mobj, "hi_model_spec")

  model_label <- sprintf("L%d model", mobj$level %||% NA_integer_)
  if (!is.null(mobj$name) && nzchar(mobj$name)) {
    model_label <- sprintf("%s '%s'", model_label, mobj$name)
  }

  aliased_terms <- mobj$aliased_terms %||% character(0)
  alias_complete <- mobj$alias_complete
  formula_label <- model_formula %||% mobj$model_formula
  if (!is.null(formula_label)) formula_label <- paste(as.character(formula_label), collapse = " ")

  coef_to_model_term <- function(coef_names) {
    if (is.null(mobj$model_matrix) || is.null(mobj$lmfit)) return(coef_names)
    assign_idx <- attr(mobj$model_matrix, "assign")
    term_labels <- c("(Intercept)", attr(stats::terms(mobj$lmfit), "term.labels"))
    coef_map <- stats::setNames(term_labels[assign_idx + 1L], colnames(mobj$model_matrix))
    unique(stats::na.omit(coef_map[intersect(coef_names, names(coef_map))]))
  }

  relation_lines <- character(0)
  if (!is.null(alias_complete) && length(aliased_terms) > 0L) {
    for (coef_name in aliased_terms) {
      alias_row <- alias_complete[coef_name, , drop = FALSE]
      dep_cols <- colnames(alias_complete)[abs(as.numeric(alias_row[1L, ])) > 1e-10]
      lhs_terms <- coef_to_model_term(coef_name)
      rhs_terms <- setdiff(coef_to_model_term(dep_cols), lhs_terms)
      if (length(rhs_terms) == 0L) rhs_terms <- dep_cols
      relation_lines <- c(
        relation_lines,
        sprintf(
          "  - %s (%s) is perfectly explained by: %s",
          coef_name,
          paste(lhs_terms, collapse = ", "),
          paste(rhs_terms, collapse = ", ")
        )
      )
    }
  }

  data_hint_lines <- character(0)
  model_data <- mobj$model_data
  model_vars <- setdiff(mobj$model_variables, c("dummy"))
  if (!is.null(model_data) && length(model_vars) > 1L) {
    present_vars <- intersect(model_vars, names(model_data))
    is_num <- vapply(model_data[, present_vars, drop = FALSE], function(x) is.numeric(x) || is.integer(x), logical(1))
    is_cat <- vapply(model_data[, present_vars, drop = FALSE], function(x) is.factor(x) || is.character(x), logical(1))
    num_vars <- present_vars[is_num]
    cat_vars <- present_vars[is_cat]

    for (num_var in num_vars) {
      for (cat_var in cat_vars) {
        dd <- model_data[, c(num_var, cat_var), drop = FALSE]
        dd <- dd[stats::complete.cases(dd), , drop = FALSE]
        if (nrow(dd) == 0L) next
        vals_by_level <- split(dd[[num_var]], as.character(dd[[cat_var]]))
        one_value_per_level <- all(vapply(vals_by_level, function(x) length(unique(x)) <= 1L, logical(1)))
        if (isTRUE(one_value_per_level) && length(unique(dd[[num_var]])) > 1L && length(vals_by_level) > 1L) {
          data_hint_lines <- c(
            data_hint_lines,
            sprintf(
              "  - %s has exactly one value within each %s level, so it is redundant with %s.",
              num_var, cat_var, cat_var
            )
          )
        }
      }
    }
  }
  data_hint_lines <- unique(data_hint_lines)

  lines <- c(
    sprintf("%s contains aliased coefficients; the fixed-effects design is rank deficient.", model_label),
    "",
    "Aliased coefficient(s):",
    paste0("  - ", aliased_terms)
  )

  if (!is.null(formula_label) && nzchar(formula_label)) {
    lines <- c(lines, "", "Model formula:", paste0("  ", formula_label))
  }

  if (length(relation_lines) > 0L) {
    lines <- c(lines, "", "Perfect redundancy detected in the design matrix:", relation_lines)
  }

  if (length(data_hint_lines) > 0L) {
    lines <- c(lines, "", "Data pattern suggesting the cause:", data_hint_lines)
  }

  paste(c(
    lines,
    "",
    "Please remove or recode redundant predictors and respecify the model."
  ), collapse = "\n")
}

# Apply the same predictor settings during initial specification and L2/L3 refits.
# Callers remove incomplete model rows before computing transformation parameters.
apply_model_covariate_settings <- function(data, mobj, lg) {
  if (is.null(lg)) lg <- lgr::get_logger()

  # 1. Re-apply num2fac conversions first so numeric/integer columns flagged as factors stay categorical on refits
  if (!is.null(mobj$num2fac) && length(mobj$num2fac) > 0L) {
    for (vv in intersect(mobj$num2fac, names(data))) {
      levs <- mobj$factor_levels[[vv]]
      if (!is.null(levs) && length(levs) > 0L) {
        data[[vv]] <- factor(as.character(data[[vv]]), levels = as.character(levs))
      } else {
        data[[vv]] <- factor(data[[vv]])
      }
    }
  }

  # 2. Enforce canonical factor levels where recorded on mobj
  if (!is.null(mobj$factor_levels) && length(mobj$factor_levels) > 0L) {
    for (vv in intersect(names(mobj$factor_levels), names(data))) {
      if (is.factor(data[[vv]]) || is.character(data[[vv]]) || vv %in% mobj$num2fac) {
        levs <- as.character(mobj$factor_levels[[vv]])
        obs_vals <- unique(as.character(stats::na.omit(data[[vv]])))
        all_levs <- unique(c(levs, obs_vals))
        data[[vv]] <- factor(as.character(data[[vv]]), levels = all_levs)
      }
    }
  }

  # 3. Continuous covariate transformations (guarding zero-SD / 1-row zscore against NaN)
  if (!is.null(mobj$covariate_transform)) {
    for (cc in seq_along(mobj$covariate_transform)) {
      cname <- names(mobj$covariate_transform)[cc]
      if (!cname %in% names(data) || is.factor(data[[cname]])) next
      ctype <- mobj$covariate_transform[cc]
      x <- as.numeric(data[[cname]])
      if (all(is.na(x))) next
      if (ctype == "mean") {
        data[[cname]] <- x - mean(x, na.rm = TRUE)
      } else if (ctype == "zscore") {
        s <- stats::sd(x, na.rm = TRUE)
        if (is.na(s) || s <= 1e-12) {
          lg$warn("Covariate '%s' has zero or undefined SD during zscore transform; setting centered values to 0.", cname)
          data[[cname]] <- rep(0, length(x))
        } else {
          data[[cname]] <- as.vector((x - mean(x, na.rm = TRUE)) / s)
        }
      } else if (ctype == "min") {
        data[[cname]] <- x - min(x, na.rm = TRUE)
      } else if (ctype == "max") {
        data[[cname]] <- x - max(x, na.rm = TRUE)
      }
    }
  }

  # 4. Reference levels (applied before factor_coding so relevel() does not wipe custom contrasts())
  if (!is.null(mobj$reference_level)) {
    for (cc in seq_along(mobj$reference_level)) {
      cname <- names(mobj$reference_level)[cc]
      if (!cname %in% names(data)) next
      if (!is.factor(data[[cname]])) {
        levs <- mobj$factor_levels[[cname]]
        data[[cname]] <- if (!is.null(levs)) factor(as.character(data[[cname]]), levels = as.character(levs)) else factor(data[[cname]])
      }
      ref <- mobj$reference_level[cc]
      if (!ref %in% levels(data[[cname]])) {
        tab <- table(data[[cname]])
        default_lev <- if (length(tab) > 0L && any(tab > 0L)) names(which.max(tab)) else levels(data[[cname]])[1L]
        lg$warn("The requested reference level: %s is not available. Defaulting to %s", ref, default_lev)
        ref <- default_lev
      }
      if (!is.na(ref) && nzchar(ref) && ref %in% levels(data[[cname]])) {
        data[[cname]] <- relevel(data[[cname]], ref = ref)
      }
    }
  }

  # 5. Optional factor contrast coding (treatment, sum, helmert, poly)
  fcoding <- mobj$factor_coding %||% mobj$contrast_spec$factor_coding
  if (!is.null(fcoding) && length(fcoding) > 0L) {
    fac_cols <- names(data)[vapply(data, is.factor, logical(1))]
    for (fc in fac_cols) {
      if (identical(fc, "id")) next
      code_this <- if (!is.null(names(fcoding)) && fc %in% names(fcoding)) {
        fcoding[[fc]]
      } else if (is.null(names(fcoding)) && length(fcoding) == 1L) {
        fcoding[[1L]]
      } else {
        NULL
      }
      if (!is.null(code_this)) {
        data[[fc]] <- apply_factor_coding_to_col(data[[fc]], code_this)
      }
    }
  }

  data
}

#' internal helper function to setup a linear model for a given l1, l2, or l3 model
#' 
#' @param mobj a model object to be populated or modified
#' @param model_formula a character string specifying the formula of the model to be fit
#' @param data a data.frame containing all columns used in model fitting
#' @param id_cols a character vector of column names in \code{data} that identify the observations
#'   and can be used for merging the model against related datasets
#' @param report_alias whether to print a short message when aliased coefficients are detected
#' @return a model object containing the fitted model
#' @keywords internal
#' @importFrom stats alias model.matrix relevel update.formula
#' @importFrom tidyselect any_of
#' @importFrom utils capture.output
mobj_fit_lm <- function(mobj=NULL, model_formula=NULL, data, id_cols=NULL, lg=NULL, report_alias = TRUE) {
  if (is.null(lg)) { lg <- lgr::get_logger() }
  # verify that we have an object of known structure
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec"), null.ok=TRUE)
  if (is.null(mobj)) {
    mobj <- list()
    class(mobj) <- c("list", "hi_model_spec")
  } else {
    # if mobj is passed in, but other elements are NULL, we are refitting a model to a new/updated dataset
    if (is.null(model_formula)) {
      model_formula <- mobj$model_formula
    }

    if (is.null(id_cols)) {
      id_cols <- names(mobj$metadata)
    }
  }

  checkmate::assert_string(model_formula)
  checkmate::assert_data_frame(data)
  checkmate::assert_subset(id_cols, names(data))

  # verify that there is a dummy DV for lm fitting
  if (!"dummy" %in% names(data)) {
    data$dummy <- seq_len(nrow(data))
  }

  model_formula <- update.formula(as.formula(model_formula), "dummy ~ .") # add LHS

  model_vars <- all.vars(model_formula)
  extra_vars <- intersect(mobj$random_effects_vars %||% character(0), names(data))
  keep_vars <- unique(c(model_vars, extra_vars))

  checkmate::assert_subset(model_vars, names(data)) # verify that all predictors are present

  # look for missingness on any predictor variable (including random-effects grouping vars)
  data <- data %>%
    dplyr::select(dplyr::all_of(keep_vars), dplyr::all_of(id_cols)) %>%
    dplyr::mutate(any_miss = rowSums(is.na(dplyr::pick(dplyr::any_of(keep_vars)))) > 0)

  miss_data <- data %>%
    dplyr::filter(.data$any_miss == TRUE) %>%
    dplyr::select(-"any_miss")

  # retain non-missing data
  data <- data %>%
    dplyr::filter(.data$any_miss == FALSE) %>%
    dplyr::select(-"any_miss")

  if (nrow(miss_data) > 0L) {
    lg$warn("Model data contain missing values for one or more covariates.")
    lg$warn("These observations will be dropped from model outputs!")
    lg$warn("%s", capture.output(print(miss_data %>% dplyr::select(-"dummy"))))
  }

  # Add metadata for merging model matrix against identifying columns.
  # Keep `id` as character (unless factor id was explicitly passed without session coercion issues)
  # and keep integer-valued `session` as integer in metadata so downstream joins against
  # integer session columns never fail when `num2fac` converts `session` to a factor in model_data.
  meta_df <- data %>% dplyr::select(dplyr::all_of(id_cols))
  if ("id" %in% names(meta_df) && is.factor(meta_df$id)) {
    meta_df$id <- as.character(meta_df$id)
  }
  if ("session" %in% names(meta_df) && is.factor(meta_df$session)) {
    sess_chr <- as.character(meta_df$session)
    sess_int <- suppressWarnings(as.integer(sess_chr))
    if (!anyNA(sess_int) && all(as.character(sess_int) == sess_chr)) {
      meta_df$session <- sess_int
    } else {
      meta_df$session <- sess_chr
    }
  }
  mobj$metadata <- meta_df

  # Convert character predictors (excluding id unless id is in formula) to factor
  for (vv in setdiff(keep_vars, "dummy")) {
    if (is.character(data[[vv]])) {
      data[[vv]] <- factor(data[[vv]])
    }
  }

  data <- apply_model_covariate_settings(data, mobj, lg)

  # Record canonical factor levels on initial fit so refits preserve factor schema
  if (is.null(mobj$factor_levels)) {
    fac_preds <- setdiff(keep_vars, c("dummy", "id"))
    fac_preds <- fac_preds[vapply(fac_preds, function(v) is.factor(data[[v]]), logical(1))]
    if (length(fac_preds) > 0L) {
      mobj$factor_levels <- stats::setNames(lapply(fac_preds, function(v) levels(data[[v]])), fac_preds)
    }
  } else {
    # Ensure canonical levels reflect post-relevel ordering
    for (vv in intersect(names(mobj$factor_levels), names(data))) {
      if (is.factor(data[[vv]])) {
        mobj$factor_levels[[vv]] <- levels(data[[vv]])
      }
    }
  }

  # keep track of data used for fitting model for refitting in case of missing data
  mobj$model_data <- data %>% dplyr::select(dplyr::all_of(keep_vars))

  # fit model and populate model information
  fit_res <- fit_lm_with_canonical_levels(
    model_formula = model_formula,
    model_data = data,
    canonical_levels = mobj$factor_levels,
    factor_coding = mobj$factor_coding %||% mobj$contrast_spec$factor_coding
  )
  mobj$lmfit <- fit_res$lmfit
  mobj$model_formula <- deparse_formula_str(update.formula(as.formula(model_formula), "NULL ~ ."))
  mm_full <- model.matrix(mobj$lmfit)
  if (isTRUE(fit_res$padded)) {
    mm_sub <- mm_full[fit_res$obs_idx, , drop = FALSE]
    attr(mm_sub, "assign") <- attr(mm_full, "assign")
    attr(mm_sub, "contrasts") <- attr(mm_full, "contrasts")
    mobj$model_matrix <- mm_sub
  } else {
    mobj$model_matrix <- mm_full
  }
  mobj$regressors <- colnames(mobj$model_matrix) # actual regressors after expanding categorical variables

  # handle coefficient aliasing
  al <- tryCatch(alias(mobj$lmfit), error = function(e) list(Complete = NULL))
  coef_na <- names(stats::coef(mobj$lmfit))[is.na(stats::coef(mobj$lmfit))]
  bad_terms <- unique(c(if (!is.null(al$Complete)) rownames(al$Complete) else character(0), coef_na))
  if (length(bad_terms) > 0L) {
    if (isTRUE(report_alias)) {
      cat("Problems with aliased (redundant) terms in model.\n\n")
      cat(paste(bad_terms, collapse = ", "), "\n\n")
    }
    mobj$alias_complete <- al$Complete

    # Retain ALL unaliased columns in model_matrix_noalias (including non-aliased interaction columns like A:B)
    keep_col_names <- setdiff(colnames(mobj$model_matrix), bad_terms)
    if (length(keep_col_names) == 0L) {
      warning(
        "No unaliased (good) terms in this model, suggesting that all covariates are constant or dependent. ",
        "Reverting to intercept-only model."
      )
      mobj$model_matrix_noalias <- matrix(1, nrow = nrow(data), ncol = 1L, dimnames = list(rownames(mobj$model_matrix), "(Intercept)"))
      modeldf <- data.frame(dummy = data$dummy)
      mobj$lmfit_noalias <- stats::lm(dummy ~ 1, data = modeldf)
    } else {
      mobj$model_matrix_noalias <- mobj$model_matrix[, keep_col_names, drop = FALSE]
      modeldf <- as.data.frame(mobj$model_matrix_noalias, check.names = FALSE)
      modeldf$dummy <- data$dummy
      quoted_cols <- paste0("`", gsub("`", "", keep_col_names, fixed = TRUE), "`")
      newf <- stats::as.formula(paste("dummy ~ 0 +", paste(quoted_cols, collapse = " + ")))
      mobj$lmfit_noalias <- stats::lm(newf, data = modeldf)
    }
    mobj$aliased_terms <- bad_terms
  } else {
    mobj$aliased_terms <- NULL
    mobj$alias_complete <- NULL
    mobj$model_matrix_noalias <- NULL
    mobj$lmfit_noalias <- NULL
  }

  return(mobj)

}

respecify_l3_model <- function(mobj, new_data) {
  checkmate::assert_class(mobj, "hi_model_spec") # verify that we have an object of known structure
  checkmate::assert_data_frame(new_data)
  checkmate::assert_subset(c("id", "session"), names(new_data))

  # Combine metadata and model_data without duplicate columns; coerce join keys cleanly
  meta_join <- mobj$metadata
  if ("id" %in% names(meta_join)) meta_join$id <- as.character(meta_join$id)
  if ("session" %in% names(meta_join) && is.factor(meta_join$session)) {
    meta_join$session <- as.character(meta_join$session)
  }

  model_data_subset <- mobj$model_data
  overlap_cols <- intersect(names(meta_join), names(model_data_subset))
  if (length(overlap_cols) > 0L) {
    model_data_subset <- model_data_subset %>% dplyr::select(-dplyr::all_of(overlap_cols))
  }
  full_data <- cbind(meta_join, model_data_subset)

  new_keys <- new_data %>% dplyr::select("id", "session")
  new_keys$id <- as.character(new_keys$id)
  if (!identical(class(new_keys$session), class(full_data$session))) {
    if (is.integer(new_keys$session) || is.numeric(new_keys$session)) {
      sess_int <- suppressWarnings(as.integer(as.character(full_data$session)))
      if (!anyNA(sess_int)) {
        full_data$session <- sess_int
        new_keys$session <- as.integer(new_keys$session)
      } else {
        full_data$session <- as.character(full_data$session)
        new_keys$session <- as.character(new_keys$session)
      }
    } else {
      full_data$session <- as.character(full_data$session)
      new_keys$session <- as.character(new_keys$session)
    }
  }

  joined_data <- new_keys %>%
    dplyr::left_join(full_data, by = c("id", "session"))

  # Auto-inject per-subject indicator EVs for pooled_sessions_subject_ev mode.
  l3_input_mode <- mobj$l3_input_mode
  auto_injected_id <- FALSE
  if (identical(l3_input_mode, "pooled_sessions_subject_ev")) {
    model_vars <- mobj$model_variables
    formula_str <- mobj$model_formula
    id_in_formula <- "id" %in% model_vars ||
      (!is.null(formula_str) && grepl("\\bid\\b", formula_str, perl = TRUE))

    if (!isTRUE(id_in_formula)) {
      auto_injected_id <- TRUE
      lgr::get_logger("glm_pipeline/l3_setup")$info(
        "l3_input_mode='pooled_sessions_subject_ev': auto-injecting 'id' as a factor EV for subject-mean modeling."
      )

      joined_data$id <- as.factor(joined_data$id)

      if (!is.null(formula_str) && nzchar(formula_str)) {
        rhs <- sub("^~\\s*", "", formula_str)
        if (rhs == "1" || rhs == "") {
          augmented_formula <- "~ 0 + id"
        } else {
          augmented_formula <- paste("~ 0 + id +", rhs)
        }
      } else {
        augmented_formula <- "~ 0 + id"
      }

      lg_l3 <- lgr::get_logger("glm_pipeline/l3_setup")
      lg_l3$info("Augmented L3 formula: %s", augmented_formula)
      mobj$model_formula <- augmented_formula
      mobj$model_variables <- unique(c("id", model_vars))
      mobj$suppress_subject_ev_diagonals <- TRUE
    }
  }

  lg_l3 <- lgr::get_logger("glm_pipeline/l3_setup")
  mobj <- mobj_fit_lm(mobj = mobj, data = joined_data, lg = lg_l3, report_alias = !isTRUE(auto_injected_id))

  mobj <- get_contrasts_from_spec(mobj, mobj$lmfit, allow_interactive_dedupe = FALSE)

  # When subject indicator EVs alias pure between-subject main effects (or dropped subjects alias L3 columns),
  # prune non-estimable contrasts and align the contrast matrix onto the full-rank non-aliased design matrix.
  if (!is.null(mobj$aliased_terms) && length(mobj$aliased_terms) > 0L) {
    pruned <- prune_rank_deficient_design_and_contrasts(
      X_obs_full = mobj$model_matrix,
      coef_vec = stats::coef(mobj$lmfit),
      cmat = mobj$contrasts,
      lg = lg_l3,
      context = sprintf("L3 model '%s'", mobj$name %||% "unnamed")
    )
    mobj$model_matrix_noalias <- pruned$dmat
    mobj$contrasts_noalias <- pruned$cmat
    if (identical(l3_input_mode, "pooled_sessions_subject_ev")) {
      cmat_valid <- pruned$cmat
      dropped_all <- pruned$dropped_contrasts
      id_ev_cols <- grep("^id", colnames(pruned$dmat), value = TRUE)
      if (length(id_ev_cols) > 0L && !is.null(cmat_valid) && nrow(cmat_valid) > 0L) {
        custom_names <- rownames(mobj$contrast_list$custom) %||% character(0)
        touches_subj_ev <- apply(cmat_valid[, id_ev_cols, drop = FALSE], 1, function(r) any(abs(r) > 1e-8))
        rnames_cur <- rownames(cmat_valid) %||% paste0("contrast_", seq_len(nrow(cmat_valid)))
        invalid_bs <- touches_subj_ev & !(rnames_cur %in% custom_names)
        if (any(invalid_bs)) {
          dropped_all <- unique(c(dropped_all, rnames_cur[invalid_bs]))
          cmat_valid <- cmat_valid[!invalid_bs, , drop = FALSE]
        }
      }
      if (length(dropped_all) > 0L) {
        warning(
          sprintf(
            "In L3 model '%s' (pooled_sessions_subject_ev), pure between-subject effects are aliased with subject intercepts and cannot be tested against within-subject error; dropped non-estimable/between-subject contrast(s): %s.",
            mobj$name %||% "unnamed",
            paste(dropped_all, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      mobj$contrasts_noalias <- cmat_valid
      mobj$contrasts <- cmat_valid
      mobj$regressors <- colnames(pruned$dmat)
    }
  }

  return(mobj)
}

#' Helper function to obtain the number of volumes in a 4D nifti file
#' 
#' @param nifti a 4D nifti file
#' @return the number of volumes in \code{nifti}
#' @keywords internal
lookup_run_volumes <- function(nifti) {
  if (!file.exists(nifti)) return(NA_integer_)

  # use fast Rcpp function to get dims without reading data
  nii_dim <- get_nifti_dim(nifti)
  stopifnot(length(nii_dim) > 3)

  return(nii_dim[4])
}

#' Helper function to obtain the number of voxels in a 4D file for populating FEAT FSF
#'
#' @param nifti a 4D nifti file
#' @return the number of voxels in \code{nifti} as calculated by x * y * z * t
#' @keywords internal
lookup_nvoxels <- function(nifti) {
  if (!file.exists(nifti)) return(NA_integer_)

  # use fast Rcpp function to get dims without reading data
  nii_dim <- get_nifti_dim(nifti)
  nvoxels <- prod(nii_dim)

  return(nvoxels)
}

#' Helper function to obtain the dimensions of a NIfTI file
#'
#' @param nifti a nifti file
#' @return a named vector of dimensions for the input (dim_x, dim_y, dim_z, dim_t)
#' @keywords internal
lookup_dim <- function(nifti) {
  ret_vec <- c(dim_x = NA_integer_, dim_y = NA_integer_, dim_z = NA_integer_, dim_t = NA_integer_)
  if (!file.exists(nifti)) return(ret_vec)
  
  nii_dim <- get_nifti_dim(nifti)
  ret_vec[seq_along(nii_dim)] <- nii_dim
  
  return(ret_vec)
}

#' small helper function to parse duration syntax of days-hours:minutes:seconds
#'   into lubridate duration object
#'
#' @param str string containing a duration that may include a days specification
#' @importFrom lubridate hms
#' @importFrom checkmate assert_string
#' @keywords internal
dhms <- function(str) {
  checkmate::assert_string(str)
  str <- validate_dhms(str) # force checks
  if (grepl("^\\d+-", str, perl = TRUE)) {
    split_hyphen <- strsplit(str, "-", fixed = TRUE)[[1]]
    days <- as.numeric(split_hyphen[1])
    period <- lubridate::hms(split_hyphen[2:length(split_hyphen)])
    period@day <- days
  } else {
    period <- lubridate::hms(str)
  }
  return(period)
}

#' helper function to convert a period object back into a d-h:m:s string
#'
#' @keywords internal
#' @importFrom checkmate assert_class
#' @importFrom lubridate seconds_to_period period_to_seconds
as_dhms <- function(per) {
  checkmate::assert_class(per, "Period")

  # if any element of the period exceeds intelligent tolerances (e.g., minutes > 60),
  # re-create the period so that it enforces reasonable tolerances.
  per <- seconds_to_period(period_to_seconds(per))
    
  dhms_str <- paste0(
    sprintf("%02d", hour(per)), ":",
    sprintf("%02d", minute(per)), ":",
    sprintf("%02d", second(per))
  )

  if (per@day > 0) {
    dhms_str <- paste0(per@day, "-", dhms_str)
  }
  return(dhms_str)
}

#' helper function to validate format of walltime inputs for HPC submission
#'
#' @param str string containing a duration that may include a days specification
#' @importFrom checkmate assert_string
#' @details this always converts to an hms format, and if days are present, it
#'   converts to dhms. Supported date formats match slurm sbatch:
#'   https://slurm.schedmd.com/sbatch.html
#' @keywords internal
validate_dhms <- function(str) {
  checkmate::assert_string(str)
  if (grepl("^\\d+:\\d+?$", str, perl=T)) { # m:s input
    return(paste0("00:", str)) # add 0 hours prefix
  } else if (grepl("^(\\d+-)?\\d+:\\d+:\\d+?$", str, perl=T)) { # h:m:s or d-h:m:s input
    return(str) 
  } else if (grepl("^(\\d+-)?\\d+:\\d+?$", str, perl=T)) { # d-h:m input
    return(paste0(str, ":00")) # add 0 seconds  
  } else if (grepl("^\\d+-\\d+$", str, perl=T)) { # days-hours input
    return(paste0(str, ":00:00")) # add 0 minutes, 0 seconds
  } else if (grepl("^\\d+$", str, perl=T)) {
    return(paste0("00:", str, ":00")) # minutes only -> hours:minutes:seconds
  } else {
    stop("Invalid duration string: ", str)
  }
}

#' convert a number of hours to a days, hours, minutes, seconds format
#' 
#' @importFrom lubridate day hour minute second seconds_to_period dhours
#' @keywords internal
hours_to_dhms <- function(hours, frac=FALSE) {
  checkmate::assert_number(hours, lower = 0)
  dur <- lubridate::dhours(hours)
  period <- seconds_to_period(dur)

  if (isTRUE(frac)) {
    str <- sprintf("%02d:%02d:%.03f", hour(period), minute(period), second(period))
  } else {
    str <- sprintf("%02d:%02d:%02d", hour(period), minute(period), round(second(period)))
  }

  if (day(period) > 0) {
    str <- paste0(sprintf("%d-", day(period)), str)
  }

  return(str)
}

#' helper function to refresh l3 model status and save gpa object from batch pipeline back to its cache
#' 
#' @param gpa a glm_pipeline_arguments object
#' @param backend_cache_paths optional named character vector of backend cache files to merge
#' @return a refreshed version of the gpa object
#' @importFrom lgr get_logger
#' @export
cleanup_glm_pipeline <- function(gpa, backend_cache_paths = NULL) {
  lg <- lgr::get_logger("glm_pipeline/cleanup_glm_pipeline")
  lg$set_threshold(gpa$lgr_threshold)

  if (!is.null(backend_cache_paths)) {
    checkmate::assert_character(backend_cache_paths)
    if (is.null(names(backend_cache_paths)) || any(names(backend_cache_paths) == "")) {
      warning("backend_cache_paths should be a named character vector keyed by backend.")
    }
    gpa <- merge_backend_caches(gpa, backend_cache_paths, lg = lg)
  }

  gpa <- refresh_glm_status(gpa, level = 1L, lg = lg)
  gpa <- refresh_glm_status(gpa, level = 2L, lg = lg)
  gpa <- refresh_glm_status(gpa, level = 3L, lg = lg)
  res <- tryCatch(saveRDS(gpa, file = gpa$output_locations$object_cache), error = function(e) {
    lg$error("Could not save gpa object to file: %s", gpa$output_locations$object_cache)
    return(NULL)
  })

  # replace copied L1 niftis with symlinks to the corresponding file
  if (isTRUE(gpa$glm_settings$fsl$replace_l1_nifti_symlink)) {
    feat_files <- file.path(
      gpa$l1_model_setup$fsl$feat_dir,
      "filtered_func_data.nii.gz"
    )
    feat_files <- feat_files[file.exists(feat_files)]
    input_files <- gpa$l1_model_setup$fsl$run_nifti[file.exists(feat_files)]
    for (ff in seq_along(feat_files)) {
      unlink(feat_files[ff])
      file.symlink(input_files[ff], feat_files[ff])
    }
  }

  return(gpa)
}

merge_backend_caches <- function(shared_gpa, backend_cache_paths, lg = NULL) {
  checkmate::assert_class(shared_gpa, "glm_pipeline_arguments")
  checkmate::assert_character(backend_cache_paths)
  if (is.null(lg)) lg <- lgr::get_logger()

  for (backend_name in names(backend_cache_paths)) {
    cache_path <- backend_cache_paths[[backend_name]]
    if (!file.exists(cache_path)) {
      lg$warn("Backend cache missing for %s: %s", backend_name, cache_path)
      next
    }

    env <- new.env(parent = emptyenv())
    load(cache_path, envir = env)
    if (!exists("gpa", envir = env)) {
      lg$warn("No gpa object found in backend cache: %s", cache_path)
      next
    }

    backend_gpa <- env$gpa
    if (!inherits(backend_gpa, "glm_pipeline_arguments")) {
      lg$warn("Backend cache %s does not contain a glm_pipeline_arguments object.", cache_path)
      next
    }

    for (level in 1:3) {
      setup_name <- paste0("l", level, "_model_setup")
      setup_class <- paste0("l", level, "_setup")

      backend_setup <- backend_gpa[[setup_name]]
      if (is.null(backend_setup) || !inherits(backend_setup, setup_class)) next

      if (is.null(shared_gpa[[setup_name]]) || !inherits(shared_gpa[[setup_name]], setup_class)) {
        shared_gpa[[setup_name]] <- backend_setup
      } else if (backend_name %in% names(backend_setup)) {
        shared_gpa[[setup_name]][[backend_name]] <- backend_setup[[backend_name]]
      }
    }
  }

  return(shared_gpa)
}

#' helper function to copy any missing fields in target
#' @param target a named list to be populated by defaults if fields are missing
#' @param defaults a named list containing default values
populate_defaults <- function(target = NULL, defaults) {
  if (is.null(target)) target <- list()
  if (is.data.frame(target) && nrow(target) == 1L && is.data.frame(defaults) && nrow(defaults) == 1L) {
    return_df <- TRUE
    target <- as.list(target)
    defaults <- as.list(defaults)
  } else {
    return_df <- FALSE
  }

  miss_fields <- setdiff(names(defaults), names(target))
  if (length(miss_fields) > 0L) {
    for (mm in miss_fields) {
      target[[mm]] <- defaults[[mm]]
    }
  }

  if (isTRUE(return_df)) target <- as.data.frame(target)
  return(target)
}

# little helper function to create named list from objects
named_list <- function(...) {
  vnames <- as.character(match.call())[-1]
  return(stats::setNames(list(...), vnames))
}

# little helper function to create named list from objects
named_vector <- function(...) {
  vnames <- as.character(match.call())[-1]
  vec <- list(...)
  nulls <- sapply(vec, is.null) # NULL will drop from vector
  vec[nulls] <- NULL
  vnames <- vnames[!nulls]

  # for now, all inputs must be length 1 atomic values or 0-length NULL/character(0) etc.
  sapply(vec, checkmate::assert_atomic, len = 1L)
  classes <- sapply(vec, class)
  if (length(unique(classes)) > 1L) {
    stop("All inputs must be of the same data type.")
  }

  return(stats::setNames(do.call(c, vec), vnames))
}

#' Helper function to create named data.frame from a set of objects
#' 
#' @details This function helps with the problem of having several vectors
#'   in the workspace that you want to combine with data.fram
# example:
# 
# named_data_frame <- function(...) {
#   vnames <- as.character(match.call())[-1]
#   return(setNames(data.frame(...), vnames))
# }

#' Enforce that GLM setup objects are complete
#' @param gpa A \code{glm_pipeline_arguments} object.
#' @param level GLM level to validate.
#' @param lg Logger object.
#' @param glm_software Optional backend filter.
#' @noRd
enforce_glms_complete <- function(gpa, level=1L, lg=NULL, glm_software=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower=1L, upper=3L)
  checkmate::assert_class(lg, "Logger")
  checkmate::assert_character(glm_software, null.ok = TRUE)

  obj_name <- glue("l{level}_model_setup")
  expect_class <- glue("l{level}_setup")
  obj <- gpa[[obj_name]]
  backend_names <- if (is.null(glm_software)) normalize_backend_strings(gpa$glm_software) else normalize_backend_strings(glm_software)

  if (is.null(obj) || !inherits(obj, expect_class)) {
    msg <- sprintf("No l%d_model_setup found in the glm pipeline object.", level)
    lg$error(msg)
    stop(msg)
  } else {
    if ("fsl" %in% backend_names && !is.null(obj$fsl)) {
      nmiss <- sum(obj$fsl$feat_complete == FALSE)
      n_feat_runs <- nrow(obj$fsl)
      if (nmiss == n_feat_runs) {
        msg <- sprintf("All feat runs in %s$fsl are incomplete.", obj_name)
        lg$error(msg)
        stop(msg)
      } else if (nmiss > 0) {
        lg$warn(
          "There are %d missing FEAT outputs in %s$fsl. Using complete %d outputs.",
          nmiss, obj_name, n_feat_runs - nmiss
        )
      }
    }

    if ("spm" %in% backend_names) {
      if (!is.null(obj$spm) && "spm_complete" %in% names(obj$spm)) {
        nmiss <- sum(obj$spm$spm_complete == FALSE)
        n_spm_runs <- nrow(obj$spm)
        if (nmiss == n_spm_runs) {
          msg <- sprintf("All SPM runs in %s$spm are incomplete.", obj_name)
          lg$error(msg)
          stop(msg)
        } else if (nmiss > 0) {
          lg$warn(
            "There are %d missing SPM outputs in %s$spm. Using complete %d outputs.",
            nmiss, obj_name, n_spm_runs - nmiss
          )
        }
      }
    }

    if ("afni" %in% backend_names) {
      if (!is.null(obj$afni) && "afni_complete" %in% names(obj$afni)) {
        nmiss <- sum(obj$afni$afni_complete == FALSE)
        n_afni_runs <- nrow(obj$afni)
        if (nmiss == n_afni_runs) {
          msg <- sprintf("All AFNI runs in %s$afni are incomplete.", obj_name)
          lg$error(msg)
          stop(msg)
        } else if (nmiss > 0) {
          lg$warn(
            "There are %d missing AFNI outputs in %s$afni. Using complete %d outputs.",
            nmiss, obj_name, n_afni_runs - nmiss
          )
        }
      }
    }
  }

  return(invisible(NULL))
}


#' helper function to ask user to choose models at a given level for further processing
#' @param gpa a \code{glm_pipeline_arguments} with models specified at a given level
#' @param model_names a user-specified string of models to process/use
#' @param level the level of GLM analysis to be specified (1, 2, or 3)
#' @return a character vector of user-specified model names at this level
#' @keywords internal
choose_glm_models <- function(gpa, model_names, level, lg=NULL) {
  checkmate::assert_integerish(level, lower = 1, upper = 3)
  all_m_names <- names(gpa[[paste0("l", level, "_models")]]$models)
  checkmate::assert_subset(model_names, c("prompt", "all", "none", all_m_names))
  if (is.null(lg)) lg <- lgr::get_logger("")
  checkmate::assert_class(lg, "Logger")

  if (is.null(model_names)) {
    chosen_models <- NULL # happens when user de-selects all models at one level
  } else if (model_names[1L] == "all") {
    chosen_models <- all_m_names
  } else if (model_names[1L] == "none") {
    chosen_models <- NULL
  } else if (model_names[1L] == "prompt") {
    cat("\n")
    chosen_models <- select.list(all_m_names,
      multiple = TRUE,
      title = paste("Choose all level", level, "models to include: ")
    )
    if (identical(chosen_models, character(0))) {
      lg$info(paste("No level", level, "models were selected."))
      chosen_models <- NULL
    }
  } else {
    chosen_models <- model_names # user-specified set
  }
  return(chosen_models)
}


#' helper function to validate numbers in an input string or vector
#' @param inp A string, character vector, or numeric vector
#' @param lower The lowest valid value for each number
#' @param upper The highest valid value for each number
#' @param as_string If TRUE, return numbers as a space-separated string. If FALSE,
#'   return the validated numeric vector itself.
#' @keywords internal
check_nums <- function(inp, lower = 0, upper = 1e10, as_string = TRUE) {
  inp_name <- deparse(substitute(inp)) # get name of object passed in
  if (checkmate::test_string(inp)) {
    inp <- suppressWarnings(as.numeric(strsplit(inp, "[\\s,;]+", perl=TRUE)[[1]]))
  } else if (checkmate::test_character(inp)) {
    inp <- suppressWarnings(as.numeric(inp))
  }

  if (checkmate::test_numeric(inp, lower = lower, upper = upper, any.missing = FALSE)) {
    return(paste(inp, collapse = " "))
  } else {
    stop("Problem with ", inp_name, " specification: ", paste(inp, collapse = " "))
  }
}

#' this helper function replaces matching rows of a current data.frame with rows in a new data.frame
#' 
#' @param current The current data.frame that may contains rows matching the new data.frame
#' @param new A data.frame containing new data that should be added to the current data.frame
#' @param id_cols A character vector of columns names that must exist in both the current
#'   and new data.frames. These are used to determine which rows match in the two datasets.
#' @details The goal here is to keep records in the current data.frame that don't overlap with
#'   the new data.frame and replace overlapping records with those in the new data.frame. This is helpful
#'   when you want to update a master data.frame with new records that may overlap current records or may
#'   be truly new. If there is no overlap in the datasets, this function basically just binds them
#'   together using dplyr::bind_rows.
#' @importFrom dplyr bind_rows across
#' @importFrom tidyselect all_of
#' @importFrom checkmate assert_character assert_data_frame assert_subset test_data_table
#' @keywords internal
update_df <- function(current = NULL, new = NULL, id_cols = NULL, sort = TRUE) {
  checkmate::assert_character(id_cols, any.missing = FALSE)
  checkmate::assert_data_frame(new)

  has_current <- !(is.null(current) || (is.data.frame(current) && nrow(current) == 0L))
  if (has_current) {
    checkmate::assert_data_frame(current)
    checkmate::assert_subset(id_cols, names(current)) # enforce all id columns in current
  } else {
    # no current data.frame (nothing to append to)
    return(new)
  }

  checkmate::assert_subset(id_cols, names(new)) # enforce all id columns in new

  # data.table uses different subsetting semantics; normalize to data.frame
  if (checkmate::test_data_table(current)) current <- as.data.frame(current)
  if (checkmate::test_data_table(new)) new <- as.data.frame(new)

  id_list <- as.list(id_cols) # to make do.call happy

  # I wonder if there's a simple way to do an anti_join approach to the problem, but this string approach is efficient

  # create unique strings for the combination of identifying columns in the current and new data
  current_id <- do.call(paste, c(current[id_cols], sep = "-"))
  new_id <- do.call(paste, c(new[id_cols], sep="-"))

  to_replace <- current_id %in% new_id
  if (any(to_replace)) {
    ret <- current[!to_replace, , drop = F] %>% dplyr::bind_rows(new)
  } else {
    ret <- current %>% dplyr::bind_rows(new) # no overlap, so just rbind
  }

  if (isTRUE(sort)) {
    ret <- ret %>% dplyr::arrange(across(all_of(id_cols)))
  }

  return(ret)
}

# helper function for printing current selections in case of NULL
c_string <- function(vec, null_val="none") {
  if (is.null(vec)) {
    null_val
  } else {
    paste(vec, collapse = ", ")
  }
}

# fast permutation calculator to avoid gtools dependence
# @details This is way faster than permutations in gtools
# and is adapted from: https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
# @param n is the number of elements to permute
# @keys is an optional character vector of keys whose permutations should be calculated
# @return a permutations x keys/elements of n matrix with all possible permutations
compute_permutations <- function(n, keys = NULL) {
  if (!is.null(keys)) {
    if (missing(n) || is.null(n)) {
      # support key input alone
      n <- length(keys)
    }
    checkmate::assert_character(keys, len = n)
  }
  
  if (n == 1) {
    A <- matrix(1)
  } else {
    sp <- Recall(n - 1)
    p <- nrow(sp)
    A <- matrix(nrow = n * p, ncol = n)

    for (i in 1:n) {
      A[(i - 1) * p + 1:p, ] <- cbind(i, sp + (sp >= i))
    }
  }

  # convert to lookup of keys based on indices in A
  if (!is.null(keys)) {
    A <- apply(A, 2, function(x) { keys[x]} )
  }

  return(A)
}

#' internal function to print bidirectional set differences in a list containing vectors to be compared
#' @param l a named list containing vectors to be compared
#' @details All pairwise combinations of vectors in the list will be compared. The setdiff() operation is
#'   run twice for each pair, reflecting the two directions of comparison (s1 %notin% s2 and s2 %notin% s1).
#' @return NULL (invisibly)
#' @keywords internal
#' @importFrom utils combn
setdiff_list_combn <- function(l) {
  combs <- combn(length(l), 2)
  for (i in seq_len(ncol(combs))) {
    v1 <- l[[combs[1, i]]]
    v2 <- l[[combs[2, i]]]
    n1 <- names(l)[combs[1, i]]
    n2 <- names(l)[combs[2, i]]
    s1 <- setdiff(v1, v2)
    s2 <- setdiff(v2, v1)

    if (length(s1) > 0L) {
      cat(glue("Values in {n1} that are not in {n2}: {paste(s1, collapse=', ')}\n", .trim = F))
    }

    if (length(s2) > 0L) {
      cat(glue("Values in {n2} that are not in {n1}: {paste(s2, collapse=', ')}\n", .trim = F))
    }
  }
  return(invisible(NULL))
}

preprocess_all_functional <- function(run_df, output_directory = NULL, cleanup_failed = TRUE) {
  if (is.null(output_directory)) stop("Must specify the root directory for outputs")
  run_df <- run_df %>%
    dplyr::mutate(
      expect_func_dir = glue::glue_data(., "{output_directory}/sub-{id}/func/{task_name}_run-{sprintf('%02d', run_number)}"),
      expect_func_file = glue::glue_data(., "{expect_func_dir}/{basename(run_nifti)}"),
      expect_func_json = sub("\\.nii\\.gz$", ".json", run_nifti),
      expect_complete = glue::glue_data(., "{expect_func_dir}/.preprocessfunctional_complete"),
      expect_mprage_bet = glue::glue_data(., "{output_directory}/sub-{id}/anat/sub-{id}_T1w_bet.nii.gz"),
      expect_mprage_warpcoef = glue::glue_data(., "{output_directory}/sub-{id}/anat/sub-{id}_T1w_warpcoef_withgdc.nii.gz"),
      output_dir_exists = dir.exists(expect_func_dir),
      func_file_exists = file.exists(expect_func_file),
      is_complete = file.exists(expect_complete),
      mprage_bet_exists = file.exists(expect_mprage_bet),
      mprage_warpcoef_exists = file.exists(expect_mprage_warpcoef)
    )

    if (any(!run_df$mprage_bet_exists)) {
      warning("Cannot process some functional data because mprage_bet.nii.gz is not present (run preprocessMprage first?)")
      cat(run_df$expect_mprage_bet[!run_df$mprage_bet_exists], sep = "\n")
      run_df <- run_df %>% dplyr::filter(.data$mprage_bet_exists == TRUE)
    }

    if (any(!run_df$mprage_warpcoef_exists)) {
      warning("Cannot process some functional data because mprage_warpcoef_withgdc.nii.gz is not present (run preprocessMprage first?)")
      cat(run_df$expect_mprage_warpcoef[!run_df$mprage_warpcoef_exists], sep = "\n")
      run_df <- run_df %>% dplyr::filter(.data$mprage_warpcoef_exists == TRUE)
    }

    # cleanup failed runs
    if (isTRUE(cleanup_failed)) {
      to_dump <- with(run_df, output_dir_exists == TRUE & is_complete == FALSE)
      if (any(to_dump)) {
        message(glue("Deleting failed directories: {paste(run_df$expect_func_dir[to_dump], collapse=', ')}"))
        unlink(run_df$expect_func_dir[to_dump], recursive = TRUE)
        run_df$output_dir_exists[to_dump] <- FALSE
        run_df$func_file_exists[to_dump] <- FALSE
      }
    }

    # need to setup output directories before symlinks and calls
    miss_dirs <- with(run_df, output_dir_exists == FALSE)
    if (any(miss_dirs)) {
      sapply(run_df$expect_func_dir[miss_dirs], dir.create, recursive = TRUE)
      run_df$output_dir_exists[miss_dirs] <- TRUE
    }

    # make symlinks
    to_link <- with(run_df, func_file_exists == FALSE)
    if (any(to_link)) {
      file.symlink(run_df$run_nifti[to_link], run_df$expect_func_file[to_link]) # accepts vector inputs
      run_df$func_file_exists[to_link] <- TRUE
    }

    # read and copy through slice times from json
    for (ii in seq_len(nrow(run_df))) {
      jfile <- run_df$expect_func_json[ii]
      if (file.exists(jfile)) {
        ff <- jsonlite::read_json(jfile)
        stimes <- as.numeric(ff$SliceTiming)
        writeLines(paste(as.character(stimes), collapse = ","), con = file.path(run_df$expect_func_dir[ii], ".stimes"))
      }
    }

    to_run <- run_df %>% filter(is_complete == FALSE)
    # don't pass through sbref if it does not exist
    sbref_string <- sapply(to_run$sbref_nifti, function(x) if (file.exists(x)) paste("-func_refimg", x) else "")

    calls <- paste(
      "cd", to_run$expect_func_dir, "&&",
      "preprocessFunctional",
      "-4d", basename(to_run$expect_func_file), sbref_string, "-se_phasepos", to_run$se_pos, "-se_phaseneg", to_run$se_neg,
      "-mprage_bet", to_run$expect_mprage_bet, "-warpcoef", to_run$expect_mprage_warpcoef,
      "-epi_pedir y- -epi_echospacing .00053 -epi_te 30 -tr .635",
      "-hp_filter 120s -rescaling_method 100_voxelmean -template_brain MNI_2.3mm",
      "-func_struc_dof bbr -warp_interpolation spline -constrain_to_template y",
      "-4d_slice_motion -custom_slice_times .stimes",
      "-ica_aroma -motion_censor fd=0.9",
      "-nuisance_file nuisance_regressors.txt -nuisance_compute csf,dcsf,wm,dwm -smoothing_kernel 6 -cleanup",
      ">preprocessFunctional_stdout 2>preprocessFunctional_stderr"
    )

    pre <- c(
      "module use /proj/mnhallqlab/sw/modules",
      "module load afni/23.0.07",
      "module load fsl/6.0.6",
      "module load r/4.2.1",
      "module load c3d/1.1.0",
      "module load freesurfer/6.0.0",
      "module load ants/2.3.1",
      "module load imagemagick/7.1.1-11",
      "source /proj/mnhallqlab/lab_resources/lab_python3/bin/activate"
    )

    scripts_out <- tempdir()
    message("Writing job scripts to: ", scripts_out)

    cluster_submit_shell_jobs(calls,
      commands_per_cpu = 1L, cpus_per_job = 8L, memgb_per_command = 24, time_per_job = "50:00:00",
      pre = pre, debug = FALSE, job_out_dir = scripts_out
    )
}

preprocess_all_mprage <- function(run_df, output_directory = NULL, cleanup_failed = TRUE) {
  if (is.null(output_directory)) stop("Must specify the root directory for outputs")
  run_df <- run_df %>%
    dplyr::distinct(.data$t1w, .keep_all = TRUE) %>% # drop repeated t1 images
    dplyr::mutate(
      expect_mprage_dir = glue::glue_data(., "{output_directory}/sub-{id}/anat"),
      expect_mprage_file = glue::glue_data(., "{output_directory}/sub-{id}/anat/sub-{id}_T1w.nii.gz"),
      expect_complete = glue::glue_data(., "{expect_mprage_dir}/.preprocessmprage_complete"),
      output_dir_exists = dir.exists(expect_mprage_dir),
      mprage_file_exists = file.exists(expect_mprage_file),
      is_complete = file.exists(expect_complete)
    )

  # cleanup failed runs
  if (isTRUE(cleanup_failed)) {
    to_dump <- with(run_df, output_dir_exists == TRUE & is_complete == FALSE)
    if (any(to_dump)) {
      message(glue("Deleting failed directories: {paste(run_df$expect_mprage_dir[to_dump], collapse=', ')}"))
      unlink(run_df$expect_mprage_dir[to_dump], recursive = TRUE)
      run_df$output_dir_exists[to_dump] <- FALSE
      run_df$mprage_file_exists[to_dump] <- FALSE
    }
  }

  # need to setup output directories before symlinks and calls
  miss_dirs <- with(run_df, output_dir_exists == FALSE)
  if (any(miss_dirs)) {
    sapply(run_df$expect_mprage_dir[miss_dirs], dir.create, recursive = TRUE)
    run_df$output_dir_exists[miss_dirs] <- TRUE
  }

  # make symlinks
  to_link <- with(run_df, mprage_file_exists == FALSE)
  if (any(to_link)) {
    file.symlink(run_df$t1w[to_link], run_df$expect_mprage_file[to_link]) # accepts vector inputs
    run_df$mprage_file_exists[to_link] <- TRUE
  }

  to_run <- run_df %>% filter(is_complete == FALSE)

  calls <- paste(
    "cd", to_run$expect_mprage_dir, "&&",
    "preprocessMprage -template_brain MNI_2mm -grad_unwarp prisma.coeff.grad -weak_bias -cleanup",
    "-n", to_run$expect_mprage_file,
    ">preprocessMprage_stdout 2>preprocessMprage_stderr"
  )

  pre <- c(
    "module use /proj/mnhallqlab/sw/modules",
    "module load afni/23.0.07",
    "module load fsl/6.0.6",
    "module load r/4.2.1",
    "module load c3d/1.1.0",
    "module load freesurfer/6.0.0",
    "module load ants/2.3.1",
    "module load imagemagick/7.1.1-11",
    "source /proj/mnhallqlab/lab_resources/lab_python3/bin/activate"
  )

  scripts_out <- tempdir()
  message("Writing job scripts to: ", scripts_out)

  cluster_submit_shell_jobs(calls, time_per_job = "2:00:00",
    commands_per_cpu = 1L, cpus_per_job = 8L, memgb_per_command = 12,
    pre = pre, debug = FALSE, job_out_dir = scripts_out
  )
}
