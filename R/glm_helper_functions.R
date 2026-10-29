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

#' helper function to generate a contrast matrix from an lm() object
#'   and a set of user-specified contrasts using emmeans
#'
#' @param mobj a model object created by build_l<X>_models
#' @param lmfit an optional lm() object used for emmeans calculations. If provided, this object
#'   will be used instead of mobj$lmfit (the parent lm on the overall dataset).
#'
#' @return a modified copy of the model object \code{mobj} with $contrast_list and $contrasts
#'   fully populated 
#' @keywords internal
#' @importFrom emmeans emmeans emtrends
#' @importFrom graphics pairs
#' @importFrom stats coef na.omit
#' @importFrom utils menu
get_contrasts_from_spec <- function(mobj, lmfit=NULL) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "l1_wi_spec", "hi_model_spec")) # verify that we have an object of known structure
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
  c_colnames <- paste0(prefix, c_colnames)

  # If no contrasts were specified in spec or contrast_list for L2/L3 models,
  # default to diagonal contrasts and warn the user.
  if (inherits(mobj, "hi_model_spec")) {
    has_list_contrasts <- any(vapply(contrast_list, function(x) {
      if (is.null(x)) return(FALSE)
      if (is.matrix(x)) return(nrow(x) > 0L)
      length(x) > 0L
    }, logical(1)))

    has_spec_contrasts <- any(c(
      isTRUE(spec$diagonal),
      length(spec$cond_means) > 0L,
      length(spec$pairwise_diffs) > 0L,
      isTRUE(spec$cell_means),
      isTRUE(spec$overall_response),
      length(spec$simple_slopes) > 0L
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
    rownames(diag_mat) <- paste0("EV_", prefix, c_colnames) # simple contrast naming for each individual regressor
    colnames(diag_mat) <- c_colnames # always have columns named by regressor

    c_diagonal <- diag_mat
  }

  ### add condition means, if requested
  c_cond_means <- contrast_list$cond_means
  if (is.null(lmfit) && (length(spec$cond_means) > 0L || length(spec$pairwise_diffs) > 0L ||
    isTRUE(spec$cell_means) || isTRUE(spec$overall_response) || length(spec$simple_slopes) > 0L)) {
    warning(
      "emmeans-based contrasts requested, but no model is available. Dropping emmeans contrasts.",
      call. = FALSE
    )
    spec$cond_means <- character(0)
    spec$pairwise_diffs <- character(0)
    spec$cell_means <- FALSE
    spec$overall_response <- FALSE
    spec$simple_slopes <- list()
  }

  if (length(spec$cond_means) > 0L && is.null(c_cond_means)) {
    for (vv in spec$cond_means) {
      ee <- emmeans(lmfit, as.formula(paste("~", vv)), weights = spec$weights)
      edata <- summary(ee)
      econ <- ee@linfct
      econ[abs(econ) < 1e-8] <- 0 # round tiny contrasts due to floating point
      enames <- unname(apply(edata[strsplit(vv, ":")[[1]]], 1, paste, collapse=".")) # combine multi-factor contrasts, if relevant

      # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
      which_na <- is.na(edata$emmean)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      # add contrast names to matrix
      rownames(econ) <- paste0(prefix, paste(vv, enames, sep = "."))

      # add contrasts to matrix
      c_cond_means <- rbind(c_cond_means, econ)
    }
    colnames(c_cond_means) <- c_colnames
  }

  ### add cell means for all factors, if requested
  c_cell_means <- contrast_list$cell_means
  if (isTRUE(spec$cell_means) && is.null(c_cell_means)) {
    # get model-predicted means for each factor
    ee <- emmeans(lmfit, as.formula(paste("~", paste(spec$cat_vars, collapse = "*"))), weights = spec$weights)
    # pp <- pairs(ee)
    edata <- summary(ee)

    econ <- ee@linfct
    econ[abs(econ) < 1e-8] <- 0 # round tiny contrasts due to floating point
    enames <- apply(edata[, spec$cat_vars, drop = FALSE], 1, function(x) { paste(names(x), x, sep=".", collapse = "_") })

    # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
    which_na <- is.na(edata$emmean)
    if (any(which_na)) {
      econ <- econ[!which_na, , drop = FALSE]
      enames <- enames[!which_na]
    }

    rownames(econ) <- paste0(prefix, enames)
    c_cell_means <- rbind(c_cell_means, econ)
    colnames(c_cell_means) <- c_colnames
  }

  ### add overall response mean
  c_overall <- contrast_list$overall
  if (isTRUE(spec$overall_response) && is.null(c_overall)) {
    ee <- emmeans(lmfit, ~1, weights = spec$weights)
    econ <- ee@linfct
    econ[abs(econ) < 1e-8] <- 0 # round tiny contrasts due to floating point
    rownames(econ) <- paste0(prefix, "overall")
    colnames(econ) <- c_colnames
    c_overall <- econ
  }

  ### add pairwise differences, if requested
  c_pairwise_diffs <- contrast_list$pairwise_diffs
  if (length(spec$pairwise_diffs) > 0L && is.null(c_pairwise_diffs)) {
    for (vv in spec$pairwise_diffs) {
      ee <- emmeans(lmfit, as.formula(paste("~", vv)), weights = spec$weights)
      pp <- pairs(ee)
      edata <- summary(pp)
      econ <- pp@linfct
      econ[abs(econ) < 1e-8] <- 0 # round tiny contrasts due to floating point
      enames <- paste(vv, make.names(sub("\\s+-\\s+", "_M_", edata$contrast, perl = TRUE)), sep=".") # names of contrasts

      # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
      which_na <- is.na(edata$estimate)
      if (any(which_na)) {
        econ <- econ[!which_na, , drop = FALSE]
        enames <- enames[!which_na]
      }

      # add contrast names to matrix
      rownames(econ) <- paste0(prefix, enames)

      # add pairwise differences for this factor to the pairwise matrix
      c_pairwise_diffs <- rbind(c_pairwise_diffs, econ)
    }
    colnames(c_pairwise_diffs) <- c_colnames
  }

  ### add per-cell model-predicted simple slopes
  c_simple_slopes <- contrast_list$simple_slopes
  if (length(spec$simple_slopes) > 0L && is.null(c_simple_slopes)) {
    for (vv in seq_along(spec$simple_slopes)) {
      trend_var <- names(spec$simple_slopes)[vv]
      for (comb in spec$simple_slopes[[vv]]) {
        ee <- emtrends(lmfit,
          specs = as.formula(paste("~", paste(comb, collapse = "*"))),
          var = trend_var, weights = spec$weights
        )
        edata <- summary(ee)
        econ <- ee@linfct
        econ[abs(econ) < 1e-8] <- 0 # round tiny contrasts due to floating point
        enames <- paste(trend_var, apply(edata[, comb, drop = FALSE], 1, paste, collapse = "."), sep = ".")

        # if any emmeans are not estimable, then this is likely due to aliasing. For now, drop
        which_na <- is.na(edata[[paste(trend_var, "trend", sep=".")]])
        if (any(which_na)) {
          econ <- econ[!which_na, , drop = FALSE]
          enames <- enames[!which_na]
        }

        rownames(econ) <- paste0(prefix, enames)
        c_simple_slopes <- rbind(c_simple_slopes, econ)
      }
    }
    colnames(c_simple_slopes) <- c_colnames
  }

  c_custom <- contrast_list$custom

  # walk through all within-subject factor contrasts and return as a single matrix
  c_within <- get_wi_contrast_matrix(mobj, c_colnames)

  #combine each element
  cmat_full <- rbind(
    c_diagonal,
    c_cond_means,
    c_cell_means,
    c_overall,
    c_pairwise_diffs,
    c_simple_slopes,
    c_within,
    c_custom
  )

  if (is.null(cmat_full)) {
    # no contrasts specified! Not good, but for now, just create an empty matrix and let upstream handle bothering the user
    cmat_full <- matrix(numeric(0), nrow=0, ncol=length(c_colnames), dimnames=list(NULL, c_colnames))
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
      warning("Some contrasts to be dropped were not found in the matrix: ", paste(spec$delete[is.na(which_del)], collapse = ", "))
      which_del <- na.omit(which_del)
      cmat <- if (length(which_del) > 0L) cmat_full[-1 * na.omit(which_del), , drop = FALSE] else cmat_full
    }
  } else {
    cmat <- cmat_full
  }

  # check for duplicated contrast rows and resolve naming collisions.
  # In non-interactive contexts (e.g., scripts/tests), retain the first row
  # deterministically and drop the rest.
  dupe_list <- get_dupe_rows(cmat)
  drop_rows <- c()
  if (length(dupe_list) > 0L) {
    is_interactive <- rlang::is_interactive()
    if (is.null(rownames(cmat))) {
      rownames(cmat) <- paste0("contrast_", seq_len(nrow(cmat)))
    }
    if (isTRUE(is_interactive)) {
      message("Duplicate contrasts found in your matrix. This can occur for many benign reasons, but we need you to say what you want
      to call these contrasts in the output so that it is clear to you.\n")
    } else {
      warning(
        "Duplicate contrasts found in non-interactive mode; keeping the first name from each duplicate set.",
        call. = FALSE
      )
    }
    
    for (dd in dupe_list) {
      message("These contrasts are the same: ", paste(rownames(cmat)[dd], collapse=", "))
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
        message(
          "Non-interactive session detected while resolving duplicate contrasts; keeping '",
          rownames(cmat)[dd][whichkeep],
          "' and dropping: ",
          paste(rownames(cmat)[dd][-whichkeep], collapse = ", ")
        )
      }
      drop_rows <- c(drop_rows, dd[-whichkeep])
    }

    # add these back into model specification so that the same decisions are made if the spec is used to repopulate a new matrix
    mobj$contrast_spec$delete <- unique(c(mobj$contrast_spec$delete, rownames(cmat)[drop_rows]))
  } 

  # duplicated does not give user any control over which one to retain (in terms of contrast name) -- now superseded by duplicate resolution above
  # dupes <- duplicated(cmat, MARGIN = 1)
  dupes <- rep(FALSE, nrow(cmat))
  dupes[drop_rows] <- TRUE
  
  cmat <- cmat[!dupes, , drop=FALSE] # drop duplicated contrasts

  # parentheses in contrast names generate problems in running through FEAT (chokes on submission)
  # so far, this comes up only with intercept terms -- may need to do a more general substitution later
  if (!is.null(cmat)) { rownames(cmat) <- gsub("(Intercept)", "Intercept", rownames(cmat), fixed = T) }

  # populate any updates to the contrast_list object based on new calculations
  contrast_list$diagonal <- c_diagonal
  contrast_list$cond_means <- c_cond_means
  contrast_list$cell_means <- c_cell_means
  contrast_list$overall <- c_overall
  contrast_list$pairwise_diffs <- c_pairwise_diffs
  contrast_list$simple_slopes <- c_simple_slopes
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

  cnames <- do.call(c, lapply(clist, rownames))
  cmat <- matrix(0, nrow = length(cnames), ncol = length(c_colnames), dimnames = list(cnames, c_colnames))

  # fill in relevant rows and columns in combined matrix for each wi model contrast matrix
  for (cc in clist) {
    cmat[rownames(cc), colnames(cc)] <- cc
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

  # create a nested data.table object for fitting each grouping unit
  data <- data.table(data, key = split_on)
  data$dummy <- seq_len(nrow(data))
  dsplit <- data[, .(dt = list(.SD)), by = split_on]

  if (!"session" %in% names(dsplit)) {
    dsplit$session <- as.integer(aggregated_session)
  }

  # use model formula from parent object
  model_formula <- terms(mobj$lmfit)

  cope_list <- list()
  contrast_list <- list()
  model_matrix_list <- list()
  n_l2_copes <- rep(NA_integer_, nrow(dsplit))
  for (vv in seq_len(nrow(dsplit))) {
    model_data <- as.data.frame(dsplit[[vv, "dt"]])
    # Match initial specification: remove incomplete observations before computing
    # transformations, then apply settings within the retained subject/session.
    complete <- stats::complete.cases(model_data[, all.vars(model_formula), drop = FALSE])
    model_data <- model_data[complete, , drop = FALSE]
    model_data <- apply_model_covariate_settings(model_data, mobj, lgr::get_logger())
    lmfit <- lm(model_formula, data = model_data)

    mm <- get_contrasts_from_spec(mobj, lmfit)
    n_contrasts <- if (is.null(mm$contrasts)) 0L else nrow(mm$contrasts)
    contrast_names <- if (n_contrasts > 0L) rownames(mm$contrasts) else character(0)
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
    contrast_list[[vv]] <- mm$contrasts
    model_matrix_list[[vv]] <- model.matrix(lmfit)
    n_l2_copes[vv] <- if (is.null(mm$contrasts)) 0L else nrow(mm$contrasts)
  }

  dsplit[, cope_list := cope_list]
  dsplit[, contrasts := contrast_list]
  dsplit[, model_matrix := model_matrix_list]
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

# Apply the same predictor settings during initial specification and L2 refits.
# Callers remove incomplete model rows before computing transformation parameters.
apply_model_covariate_settings <- function(data, mobj, lg) {
  if (!is.null(mobj$covariate_transform)) {
    for (cc in seq_along(mobj$covariate_transform)) {
      cname <- names(mobj$covariate_transform)[cc]
      ctype <- mobj$covariate_transform[cc]
      if (ctype == "mean") {
        data[[cname]] <- data[[cname]] - mean(data[[cname]], na.rm = TRUE)
      } else if (ctype == "zscore") {
        data[[cname]] <- as.vector(scale(data[[cname]]))
      } else if (ctype == "min") {
        data[[cname]] <- data[[cname]] - min(data[[cname]], na.rm = TRUE)
      } else if (ctype == "max") {
        data[[cname]] <- data[[cname]] - max(data[[cname]], na.rm = TRUE)
      }
    }
  }

  if (!is.null(mobj$reference_level)) {
    for (cc in seq_along(mobj$reference_level)) {
      cname <- names(mobj$reference_level)[cc]
      if (!is.factor(data[[cname]])) data[[cname]] <- factor(data[[cname]])
      ref <- mobj$reference_level[cc]
      if (!ref %in% levels(data[[cname]])) {
        default_lev <- names(which.max(table(data[[cname]])))
        lg$warn("The requested reference level: %s is not available. Defaulting to %s", ref, default_lev)
        ref <- default_lev
      }
      data[[cname]] <- relevel(data[[cname]], ref = ref)
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

  # use model formula from parent object
  # model_formula <- terms(mobj$lmfit)
  model_vars <- all.vars(model_formula)

  checkmate::assert_subset(model_vars, names(data)) # verify that all predictors are present

  # look for missingness on any predictor variable
  # https://stackoverflow.com/questions/64287986/create-variable-that-captures-if-there-are-missing-fields-in-4-string-variables
  data <- data %>%
    dplyr::select(dplyr::all_of(model_vars), dplyr::all_of(id_cols)) %>% # just keep model-relevant variables
    dplyr::mutate(any_miss = rowSums(is.na(dplyr::pick(dplyr::any_of(model_vars)))) > 0)
    #dplyr::mutate(across(all_of(model_vars), ~ +(is.na(.))))

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

  # Add metadata, for merging model matrix against identifying columns
  # This occurs before mean centering in case one of the id_cols is also in the model (e.g., run_number)
  mobj$metadata <- data %>% dplyr::select(dplyr::all_of(id_cols))

  data <- apply_model_covariate_settings(data, mobj, lg)

  # keep track of data used for fitting model for refitting in case of missing data
  mobj$model_data <- data %>% dplyr::select(dplyr::all_of(model_vars))

  # fit model and populate model information
  mobj$lmfit <- lm(model_formula, data)
  mobj$model_formula <- as.character(update.formula(as.formula(model_formula), "NULL ~ .")) # remove LHS and convert to string
  mobj$model_matrix <- model.matrix(mobj$lmfit)
  mobj$regressors <- colnames(mobj$model_matrix) # actual regressors after expanding categorical variables

  # handle coefficient aliasing
  al <- alias(mobj$lmfit)
  if (!is.null(al$Complete)) {
    bad_terms <- rownames(al$Complete)
    if (isTRUE(report_alias)) {
      cat("Problems with aliased (redundant) terms in model.\n\n")
      cat(paste(bad_terms, collapse = ", "), "\n\n")
    }
    mobj$alias_complete <- al$Complete

    # find unaliased (good) terms in the model
    good_terms <- colnames(mobj$model_matrix)[!(colnames(mobj$model_matrix) %in% bad_terms)]
    good_terms <- good_terms[!good_terms == "(Intercept)"]

    if (length(good_terms) == 0L) {
      warning(
        "No unaliased (good) terms in this model, suggesting that all covariates are constant or dependent. ",
        "Reverting to intercept-only model."
      )
      good_terms <- "1"
    }

    # build new model formula with only good terms
    newf <- as.formula(paste("dummy ~", paste(good_terms, collapse = " + ")))

    # also generate a model data.frame that expands dummy codes for terms, retaining only good variables
    mobj$model_matrix_noalias <- mobj$model_matrix[, grep(":", good_terms, fixed = TRUE, value = TRUE, invert = TRUE)]
    modeldf <- as.data.frame(mobj$model_matrix_noalias)
    modeldf$dummy <- data$dummy # copy across dummy DV for fitting
    mobj$lmfit_noalias <- lm(newf, modeldf)
    mobj$aliased_terms <- bad_terms

    # N.B. emmeans needs to calculate contrasts on the original design to see the factor structure
    # So, we also need to drop out columns from the emmeans linfct
  } else {
    mobj$aliased_terms <- NULL
    mobj$alias_complete <- NULL
  }

  return(mobj)

}


#OLD version
# respecify_l3_model <- function(mobj, data) {
#   checkmate::assert_multi_class(mobj, c("l1_model_spec", "hi_model_spec")) # verify that we have an object of known structure
#   checkmate::assert_data_frame(data)
#   checkmate::assert_subset(c("id", "session"), names(data))

#   # use model formula from parent object
#   model_formula <- terms(mobj$lmfit)

#   mobj$lmfit <- lm(model_formula, data = data)
#   mobj$regressors <- colnames(mobj$model_matrix)
#   mobj$model_matrix <- model.matrix(mobj$lmfit)
#   mobj$metadata <- data %>% dplyr::select("id", "session")
#   mobj$model_data <- data

#   mobj <- get_contrasts_from_spec(mobj, mobj$lmfit)

#   return(mobj)
# }


# new version, but depends on mobj being specified under new approach
respecify_l3_model <- function(mobj, new_data) {
  checkmate::assert_class(mobj, "hi_model_spec") # verify that we have an object of known structure
  checkmate::assert_data_frame(new_data)
  checkmate::assert_subset(c("id", "session"), names(new_data))

  # use id and session in new_data to subset the full data from initial model specification
  # ensure we combine metadata and model data without creating duplicate columns (e.g., session might be in both)
  model_data_subset <- mobj$model_data
  overlap_cols <- intersect(names(mobj$metadata), names(model_data_subset))
  if (length(overlap_cols) > 0L) {
    model_data_subset <- model_data_subset %>% dplyr::select(-dplyr::all_of(overlap_cols))
  }
  full_data <- cbind(mobj$metadata, model_data_subset)

  new_data <- new_data %>% dplyr::select("id", "session") #just metadata of interest
  new_data <- new_data %>%
    dplyr::left_join(full_data, by=c("id", "session"))

  # Auto-inject per-subject indicator EVs for pooled_sessions_subject_ev mode.
  # FSL's repeated-measures approach requires one EV per subject to model out
  # the subject-specific intercept (mean). If the user's L3 formula does not
  # include 'id', we automatically add it as a factor to ensure correct
  # repeated-measures modeling.
  l3_input_mode <- mobj$l3_input_mode
  if (identical(l3_input_mode, "pooled_sessions_subject_ev")) {
    model_vars <- mobj$model_variables
    formula_str <- mobj$model_formula
    id_in_formula <- "id" %in% model_vars ||
      (!is.null(formula_str) && grepl("\\bid\\b", formula_str, perl = TRUE))

    if (!isTRUE(id_in_formula)) {
      lgr::get_logger("glm_pipeline/l3_setup")$info(
        "l3_input_mode='pooled_sessions_subject_ev': auto-injecting 'id' as a factor EV for subject-mean modeling."
      )

      # ensure id is a factor in the data
      new_data$id <- as.factor(new_data$id)

      # prepend id to the existing model formula
      if (!is.null(formula_str) && nzchar(formula_str)) {
        # strip any leading "~ " or "~" prefix, then rebuild
        rhs <- sub("^~\\s*", "", formula_str)
        if (rhs == "1" || rhs == "") {
          # intercept-only model: replace with id (id already provides intercept per subject)
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
    }
  }

  mobj <- mobj_fit_lm(mobj=mobj, data=new_data)

  mobj <- get_contrasts_from_spec(mobj, mobj$lmfit)

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
