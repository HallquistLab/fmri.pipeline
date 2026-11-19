# SPM design, estimation, and timing helpers

normalize_spm_estimation_method <- function(estimation_method) {
  if (is.null(estimation_method) || !is.character(estimation_method) || length(estimation_method) != 1L) {
    stop("estimation_method must be a single non-empty string.")
  }
  if (!nzchar(estimation_method)) {
    stop("estimation_method must be a single non-empty string.")
  }

  method_norm <- gsub("[[:space:]]+", "", tolower(estimation_method))

  if (method_norm %in% c("classical", "reml")) {
    return("Classical")
  }

  if (method_norm %in% c("bayesian2", "bayes2", "bayesian_2", "bayesian-2", "secondlevel", "l2")) {
    return("Bayesian2")
  }

  if (method_norm %in% c("bayesian", "bayes", "bayesian1", "vb")) {
    stop("Bayesian 1st-level estimation is not yet supported; provide Classical or Bayesian2.")
  }

  stop("estimation_method must be one of: 'Classical' or 'Bayesian2'.")
}

normalize_spm_write_residuals <- function(write_residuals, method_label) {
  if (is.logical(write_residuals)) {
    if (length(write_residuals) != 1L) {
      stop("write_residuals must be a single logical or 0/1.")
    }
  } else if (is.numeric(write_residuals) && length(write_residuals) == 1L && write_residuals %in% c(0, 1)) {
    write_residuals <- as.logical(write_residuals)
  } else {
    stop("write_residuals must be a single logical or 0/1.")
  }

  if (isTRUE(write_residuals) && !identical(method_label, "Classical")) {
    warning("write_residuals is only supported for Classical estimation; disabling.")
    write_residuals <- FALSE
  }

  return(write_residuals)
}

read_slice_timing_from_json <- function(json_path, lg = NULL) {
  if (is.null(json_path) || !nzchar(json_path) || !file.exists(json_path)) {
    return(NULL)
  }

  meta <- tryCatch(
    jsonlite::read_json(json_path, simplifyVector = TRUE),
    error = function(e) {
      if (!is.null(lg)) lg$warn("Failed to read JSON metadata %s: %s", json_path, as.character(e))
      return(NULL)
    }
  )
  if (is.null(meta)) return(NULL)

  st <- meta$SliceTiming
  if (is.null(st)) return(NULL)
  st <- as.numeric(st)
  st <- st[!is.na(st)]
  if (length(st) < 2L) return(NULL)

  ref_fields <- c(
    "SliceTimingReference",
    "SliceTimingRef",
    "SliceReference",
    "ReferenceSlice",
    "RefSlice",
    "ReferenceSliceTime"
  )
  ref_val <- NULL
  for (ff in ref_fields) {
    if (!is.null(meta[[ff]])) {
      ref_val <- meta[[ff]]
      break
    }
  }

  ref_idx <- NULL
  if (!is.null(ref_val)) {
    ref_num <- suppressWarnings(as.numeric(ref_val))
    if (is.finite(ref_num)) {
      if (abs(ref_num - round(ref_num)) < 1e-6 && ref_num >= 1 && ref_num <= length(st)) {
        ref_idx <- as.integer(round(ref_num))
      } else if (ref_num >= min(st) - 1e-6 && ref_num <= max(st) + 1e-6) {
        ref_idx <- which.min(abs(st - ref_num))
      }
    }
  }

  return(list(times = st, ref_idx = ref_idx, source = json_path))
}

read_slice_timing_from_stimes <- function(stimes_path, lg = NULL) {
  if (is.null(stimes_path) || !nzchar(stimes_path) || !file.exists(stimes_path)) {
    return(NULL)
  }

  raw <- tryCatch(
    readLines(stimes_path, warn = FALSE),
    error = function(e) {
      if (!is.null(lg)) lg$warn("Failed to read slice timing file %s: %s", stimes_path, as.character(e))
      return(NULL)
    }
  )
  if (is.null(raw) || length(raw) == 0L) return(NULL)
  st <- as.numeric(unlist(strsplit(paste(raw, collapse = " "), "[,[:space:]]+")))
  st <- st[!is.na(st)]
  if (length(st) < 2L) return(NULL)
  return(list(times = st, ref_idx = NULL, source = stimes_path))
}

infer_spm_microtime <- function(run_nifti, fmri_t = NULL, fmri_t0 = NULL, lg = NULL) {
  run_nifti <- run_nifti[!is.na(run_nifti) & nzchar(run_nifti)]

  st_info <- NULL
  if (length(run_nifti) > 0L) {
    for (rn in run_nifti) {
      json_path <- if (grepl("\\.nii\\.gz$", rn)) {
        sub("\\.nii\\.gz$", ".json", rn)
      } else if (grepl("\\.nii$", rn)) {
        sub("\\.nii$", ".json", rn)
      } else {
        paste0(rn, ".json")
      }

      st_info <- read_slice_timing_from_json(json_path, lg = lg)
      if (!is.null(st_info)) break

      stimes_path <- file.path(dirname(rn), ".stimes")
      st_info <- read_slice_timing_from_stimes(stimes_path, lg = lg)
      if (!is.null(st_info)) break
    }
  }

  if (!is.null(st_info)) {
    n_slices <- length(st_info$times)
    if (is.null(fmri_t)) fmri_t <- n_slices
    if (is.null(fmri_t0) && !is.null(st_info$ref_idx)) fmri_t0 <- st_info$ref_idx
    if (is.null(fmri_t0)) fmri_t0 <- 1L
    if (!is.null(lg)) {
      lg$debug(
        "Inferred SPM microtime from %s: fmri_t=%d fmri_t0=%d",
        st_info$source, fmri_t, fmri_t0
      )
    }
  }

  if (is.null(fmri_t)) {
    if (length(run_nifti) > 0L && file.exists(run_nifti[1L])) {
      fmri_t <- tryCatch(
        {
          hdr <- suppressWarnings(RNifti::niftiHeader(run_nifti[1L]))
          dims <- hdr$dim
          if (length(dims) < 3L) return(NULL)
          as.integer(dims[3])
        },
        error = function(e) {
          if (!is.null(lg)) lg$warn("Could not read NIfTI header for fmri_t: %s", as.character(e))
          return(NULL)
        }
      )
      if (!is.numeric(fmri_t) || length(fmri_t) != 1L || is.na(fmri_t) || fmri_t < 1L) {
        fmri_t <- NULL
      }
    }
  }

  if (is.null(fmri_t)) fmri_t <- 20L
  if (is.null(fmri_t0)) fmri_t0 <- 1L

  if (!is.numeric(fmri_t) || length(fmri_t) != 1L || is.na(fmri_t)) {
    stop("fmri_t must be a single numeric value.")
  }
  if (!is.numeric(fmri_t0) || length(fmri_t0) != 1L || is.na(fmri_t0)) {
    stop("fmri_t0 must be a single numeric value.")
  }
  fmri_t <- as.integer(round(fmri_t))
  fmri_t0 <- as.integer(round(fmri_t0))
  if (fmri_t < 1L) stop("fmri_t must be >= 1.")
  if (fmri_t0 < 1L) stop("fmri_t0 must be >= 1.")
  if (fmri_t0 > fmri_t) {
    stop("fmri_t0 must be <= fmri_t.")
  }

  return(list(
    fmri_t = fmri_t,
    fmri_t0 = fmri_t0,
    st_source = if (!is.null(st_info)) st_info$source else NULL
  ))
}

#' @keywords internal
#' @noRd
spm_hrf <- function(TR, P) {
    ## R port of spm8 functions needed for deconvolution
    ## spm_hrf - Returns a hemodynamic response function as a difference of gammas (canonical)
    ##
    ## USAGE:
    ##   hrf = spm_hrf(TR,[p])
    ##
    ##INPUT:
    ## TR   - scan repetition time (seconds)
    ## P    - 1x7 vector of parameters for the response function (two gamma functions)
    ##                                                 (default value, in seconds)
    ## P[1] - delay of response (relative to onset)    (6)
    ## P[2] - delay of undershoot (relative to onset)  (16)
    ## P[3] - dispersion of response                   (1)
    ## P[4] - dispersion of undershoot                 (1)
    ## P[5] - ratio of response to undershoot          (6)
    ## P[6] - onset (seconds)                          (0)
    ## P[7] - length of kernel (seconds)               (32)
    ##
    ## OUTPUT:
    ##  $hrf    - hemodynamic response function
    ##  $p      - parameters of the response function
    ##_______________________________________________________________________
    ## Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

    ## Karl Friston
    ## $Id: spm_hrf.m 387 2005-12-17 18:31:23Z klaas $

    ## default parameters
    ##-----------------------------------------------------------------------
    fMRI_T <- 16 #microtime resolution is 1/16th of TR

    p <- c(6, 16, 1, 1, 6, 0, 32)

    if (!missing(P)) {
        p[1:length(P)] = P
    }

    ## modelled hemodynamic response function - {mixture of Gammas}
    ##-----------------------------------------------------------------------
    dt    <- TR/fMRI_T
    u     <- 0:(p[7]/dt) - p[6]/dt #sampling grid of HRF in microtime units (e.g., 0:256 for default params and 2s TR)

    ## MH: eliminated use of spm Gpdf in favor of built-in gamma PDF in R. Checked that this yields identical
    ## results, although the third parameter of spm_Gpdf is really rate, not shape.
    hrf   = stats::dgamma(u, shape=p[1]/p[3], rate=dt/p[3]) - stats::dgamma(u, shape=p[2]/p[4], rate=dt/p[4])/p[5]

    hrf   = hrf[0:(p[7]/TR)*fMRI_T + 1] #subsample hrf back onto TR grid
    hrf   = hrf/sum(hrf) #Unit normalize hrf

    #return hrf  explicitly as column vector because other functions cbind additional basis functions
    return(list(hrf=matrix(hrf, ncol=1), p=p))
}

spm_get_bf <- function(xBF) {
    ## fills in basis function structure
    ## FORMAT [xBF] = spm_get_bf(xBF)
    ##
    ## xBF.dt      - time bin length {seconds}
    ## xBF.name    - description of basis functions specified
    ## xBF.length  - window length (seconds)
    ## xBF.order   - order
    ## xBF.bf      - Matrix of basis functions
    ##
    ## xBF.name  'hrf'
    ##       'hrf (with time derivative)'
    ##       'hrf (with time and dispersion derivatives)'
    ##       'Fourier set'
    ##       'Fourier set (Hanning)'
    ##       'Gamma functions'
    ##       'Finite Impulse Response'}
    ##
    ## (any other specification will default to hrf)
    ##__________________________________________________________________________
    ##
    ## spm_get_bf prompts for basis functions to model event or epoch-related
    ## responses.  The basis functions returned are unitary and orthonormal
    ## when defined as a function of peri-stimulus time in time-bins.
    ## It is at this point that the distinction between event and epoch-related
    ## responses enters.
    ##_______________________________________________________________________
    ## Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

    ## Karl Friston
    ## $Id: spm_get_bf.m 3934 2010-06-17 14:58:25Z guillaume $

    ## length of time bin
    ##--------------------------------------------------------------------------
    if (missing(xBF)) {
        xBF <- list()
        while (is.null(xBF$dt) || is.na(xBF$dt)) {
            xBF$dt <- readline('time bin for basis functions in seconds (.0625): ')
            if (nchar(xBF$dt) == 0L) {
                xBF$dt <- .0625  #default
            } else {
                xBF$dt <- as.numeric(xBF$dt)
            }            
        }
    }
    dt   = xBF$dt

    ## assemble basis functions
    ##==========================================================================

    ## model event-related responses
    ##--------------------------------------------------------------------------
    if (!exists('name', where=xBF)) {
        mnames <- c("hrf", "hrf (with time derivative)", "hrf (with time and dispersion derivatives)",
                    "Fourier set", "Fourier set (Hanning)", "Gamma functions", "Finite Impulse Response")
                    
        while (is.null(xBF$name) || is.na(xBF$name) || nchar(xBF$name) == 0L) {
            cat("Hemodynamic Basis functions...\n")
            cat(paste0("  (1) ", mnames[1], "\n  (2) ", mnames[2], "\n  (3) ", mnames[3], 
                       "\n  (4) ", mnames[4], "\n  (5) ", mnames[5], "\n  (6) ", mnames[6], "\n  (7) ", mnames[7], "\n"))
            msel <- as.numeric(readline("Select basis set (1-7): "))
            if (is.na(msel) || msel < 1 || msel > 7) {
                xBF$name <- NULL
            } else {
                xBF$name <- mnames[msel]
            }
        }
    }

    ## get order and length parameters
    ##--------------------------------------------------------------------------
    if (xBF$name %in% c("Fourier set", "Fourier set (Hanning)", "Gamma functions", "Finite Impulse Response")) {
        while (is.null(xBF$length) || is.na(xBF$length) || xBF$length < 1) {
            xBF$length <- readline("window length in seconds (32): ")
            if (nchar(xBF$length) == 0L) { xBF$length <- 32 } else { xBF$length <- as.numeric(xBF$length) }  #default to 32s
        }

        while (is.null(xBF$order) || is.na(xBF$order) || xBF$order < 1) {            
            xBF$order <- readline("order (4): ")
            if (nchar(xBF$order) == 0L) { xBF$order <- 4 } else { xBF$order <- as.numeric(xBF$order) }  #default to fourth-order
        }
    }
    h   = xBF$order
    l   = xBF$length
    
    ## create basis functions
    ##--------------------------------------------------------------------------
    if (xBF$name %in% c("Fourier set", "Fourier set (Hanning)")) {
        pst   = seq(from=0, to=l, by=xBF$dt)
        pst   = pst/max(pst)

        ## hanning window
        ##----------------------------------------------------------------------
        if (xBF$name == "Fourier set (Hanning)") {
            g  = (1 - cos(2*pi*pst))/2
        } else {
            g  = rep(1, length(pst))
        }

        ## zeroth and higher Fourier terms
        ##----------------------------------------------------------------------
        bf    = g
        for (i in 1:h) {
            bf = cbind(bf, g * sin(i*2*pi*pst))
            bf = cbind(bf, g * cos(i*2*pi*pst))
        }
    } else if (xBF$name == "Gamma functions") {
        pst   = seq(from=0, to=l, by=xBF$dt)
        bf    = spm_gamma_bf(pst,h)
    } else if (xBF$name == "Finite Impulse Response") {
        bin   = l/h
        bf    = kronecker(diag(h), rep(1, round(bin/dt)))
    } else if (xBF$name == "NONE") {
        bf = 1
    } else {
        ## canonical hemodynamic response function
        ##----------------------------------------------------------------------
        hrf            = spm_hrf(dt)
        bf             = hrf$hrf
        p              = hrf$p

        ## add time derivative
        ##----------------------------------------------------------------------
        if (grepl("time", xBF$name, ignore.case = TRUE)) {
            dp     = 1
            p[6]   = p[6] + dp
            D      = (bf[,1] - spm_hrf(dt,p)$hrf)/dp
            bf     = cbind(bf, D)
            p[6]   = p[6] - dp

            ## add dispersion derivative
            ##------------------------------------------------------------------
            if (grepl("dispersion", xBF$name, ignore.case = TRUE)) {
                dp    = 0.01
                p[3]  = p[3] + dp
                D     = (bf[,1] - spm_hrf(dt,p)$hrf)/dp
                bf    = cbind(bf, D)
            }
        }

        ## length and order
        ##----------------------------------------------------------------------
        xBF$length = dim(bf)[1]*dt
        xBF$order  = dim(bf)[2]         
    }

    xBF$bf <- spm_orth(bf)
    return(xBF)

}


## compute Gamma functions
##--------------------------------------------------------------------------
spm_gamma_bf <- function(u,h) {
    ## returns basis functions used for Volterra expansion
    ## FORMAT bf = spm_gamma_bf(u,h)
    ## u   - times {seconds}
    ## h   - order
    ## bf  - basis functions (mixture of Gammas)
    ##__________________________________________________________________________
    
    u     = as.vector(u)
    bf    = c()
    for (i in 2:(1 + h)) {
        m   = 2^i
        s   = sqrt(m)
        bf  = cbind(bf, stats::dgamma(u, shape=(m/s)^2, rate=m/s^2)) #replace spm_Gpdf in favor of R built-in dgamma
    }
   
    return(bf)
}

spm_orth <- function(X, OPT) {
    ## Recursive Gram-Schmidt orthogonalisation of basis functions
    ## FORMAT X = spm_orth(X,OPT)
    ##
    ## X   - matrix
    ## OPT - 'norm' for Euclidean normalisation
    ##     - 'pad'  for zero padding of null space [default]
    ##
    ## Serial orthogonalisation starting with the first column
    ##
    ## Reference:
    ## Golub, Gene H. & Van Loan, Charles F. (1996), Matrix Computations (3rd
    ## ed.), Johns Hopkins, ISBN 978-0-8018-5414-9.
    ##__________________________________________________________________________
    ## Copyright (C) 2002-2012 Wellcome Trust Centre for Neuroimaging

    ## Karl Friston
    ## $Id: spm_orth.m 4708 2012-04-02 11:28:14Z guillaume $

    if (missing(OPT)) { OPT = "pad" }

    ##-Recursive Gram-Schmidt orthogonalisation
    ##--------------------------------------------------------------------------
    ## sw = warning('off','all'); ##suppress warnings?
    d     = dim(X)
    n     = d[1]; m <- d[2]
    X     = X[, apply(X, 2, function(col) { any(col != 0.0) }), drop=FALSE] ##Only retain non-zero columns
    rankX = qr(X)$rank

    orth <- tryCatch({
        x <- X[,1, drop=FALSE]
        j <- 1

        if (ncol(X) > 1L) {
            for (i in 2:ncol(X)) {
                D = X[,i, drop=FALSE]
                D = D - x%*%(MASS::ginv(x) %*% D)
                if (norm(D,"1") > exp(-32)) {
                    x <- cbind(x, D)
                    j[length(j) + 1] <- i
                }
                if (length(j) == rankX) { break }
            }
        }
        list(x = x, j = j)
    }, error=function(e) {
        print(e)
        list(x = matrix(0, nrow=n, ncol=0), j = c())
    })
    x <- orth$x
    j <- orth$j

    ## and normalisation, if requested
    ##--------------------------------------------------------------------------
    if (OPT == "pad") {
        X      = matrix(0, nrow=n, ncol=m)
        X[,j]  = x
    } else {
        ## Euclidean norm each column in X        
        X <- apply(x, 2, function(col) { if (any(x != 0.0)) { col/sqrt(sum(col^2)) } else { col } })
    }

    return(X)
}
#' Uses an object from build_design_matrix to generate a corresponding SPM GLM design matrix.
#'
#' @param bdm The object (of class 'bdm') generated by build_design_matrix. In particular, the function
#'            uses the \code{$design} and \code{$run_niftis} elements to extract relevant regressors and
#'            information about the NIfTI files to analyze.
#' @param ts_files A character vector of nuisance text files -- one per run -- to be added as covariates in
#'            SPM's multi_reg field.
#' @param output_dir The path where SPM scripts and outputs should go for this design matrix. The path can be
#'            relative or absolute, but will be converted to absolute internally using normalizePath().
#' @param hpf The high-pass filter cutoff (in seconds) used to remove low frequencies prior to analysis.
#'            Default: 100
#' @param hrf_derivs Whether to include derivatives of the regressors in the model to capture HRF variation.
#'            Options are: "none", "time" (temporal derivative), "time_dispersion" (temporal and dispersion derivatives).
#'            Default: "none"
#' @param cvi Serial correlations model. Options are: "AR(1)", "FAST", or "none". Default: "AR(1)"
#' @param estimation_method Estimation method. Options are: "Classical" or "Bayesian2". Default: "Classical"
#' @param write_residuals Whether to write residual images. Only supported for Classical estimation. Default: FALSE
#' @param fmri_t Microtime resolution (number of time bins per scan). Default: 20
#' @param fmri_t0 Microtime onset (reference bin). Default: 1
#' @param nifti_tmpdir Where to place uncompressed NIfTIs for analysis (SPM doesn't handle .nii.gz)
#' @param cleanup_tmp Whether to remove uncompressed .nii files after completing this function
#' @param condition_contrasts Whether to setup contrasts for each condition. In multi-run data, these contrasts will
#'            represent the condition average across all runs (e.g., [ .5, .5 ] for 2 runs). Default: TRUE
#' @param unit_contrasts Whether to estimate a unit-height contrast for every regressor in the design.
#'            This essentially includes a diagonal matrix in the contrast specification. Default: FALSE
#' @param effects_of_interest_F Whether to include a single F test contrast with all regressors of interest. Useful
#'            in adjusting VOIs for nuisance regressors. Default: TRUE
#' @param spm_execute_setup Whether to run the GLM setup script in MATLAB, creating SPM.mat. Default: FALSE
#' @param spm_execute_glm Whether to run the GLM after creating SPM.mat. This could take a while! Default: FALSE
#' @param spm_execute_contrasts Whether to compute contrasts after GLM is complete. Depends on \code{spm_execute_glm}. Default: FALSE
#' @param concatenate_runs Whether to convert multi-run data into a single concatenated session. This will execute spm_fmri_concatenate
#'            after the GLM setup is complete. This script computes run-specific whitening and high-pass filters while keeping
#'            a single concatenated time series. This is useful as a preamble to DCM, which often uses concatenated time series. Default: FALSE
#' @param spm_path The path to an spm12 directory. This will be included in MATLAB scripts to ensure that spm is found.
#' @param matlab_cmd Command to launch MATLAB/Octave. Default: "matlab"
#' @param matlab_args Arguments passed to MATLAB/Octave. Default: "-batch"
#' @param compute_env Optional character vector of shell commands to prepare the environment
#'   (e.g., module load statements). These are prepended to direct MATLAB execution.
#'
#' @importFrom R.matlab readMat
#'
#' @details
#'
#' This function is intended to setup a first-level (subject) fMRI GLM analysis in SPM12. It accepts an object from build_design_matrix
#' and uses the contents of this object to setup timing and HRF settings (including parametric modulation) for all regressors.
#'
#' At this point, the function is not very flexible in terms of allowing custom tweaks to each field. It is mostly intended to setup
#' the necessary SPM structures to support DCM, which depends on an SPM.mat file. I discovered the spm12r package recently,
#' which has a lot of flexibility for setting up a first-level model. In the future, we may wish to adapt the code below to use
#' this instead. \url{https://cran.rstudio.com/web/packages/spm12r/vignettes/fmri_task_processing.html}
#'
#' @return A list containing the MATLAB syntax/scripts for GLM setup, execution, and contrast estimation
#' @author Michael Hallquist
#' 
#' @export
generate_spm_mat <- function(bdm, ts_files=NULL, output_dir="spm_out",
                             hpf=100, hrf_derivs="none", cvi="AR(1)",
                             estimation_method="Classical", write_residuals=FALSE,
                             fmri_t=20, fmri_t0=1,
                             nifti_tmpdir=tempdir(),
                             cleanup_tmp=FALSE, condition_contrasts=TRUE,
                             unit_contrasts=TRUE, effects_of_interest_F=TRUE,
                             spm_execute_setup=FALSE, spm_execute_glm=FALSE,
                             spm_execute_contrasts=FALSE, concatenate_runs=FALSE,
                             spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12",
                             matlab_cmd="matlab", matlab_args="-batch",
                             compute_env=NULL) {

  #for concatenate runs, need to setup a single SPM mat (one session), run GLM setup, then call spm_fmri_concatenate('SPM.mat', nscans)
  #where the latter is a vector run volumes ($run_volumes in BDM)
  
  spm_syntax <- list()

  if (is.null(spm_path)) { message("No spm_path provided. The scripts will assume that spm is already in the MATLAB path.") }
  if (!is.null(spm_path) && !dir.exists(spm_path)) { stop("Provided spm_path does not exist: ", spm_path) }
  
  if (!inherits(bdm, "bdm")) { stop("generate_spm_mat requires an object of class 'bdm' from build_design_matrix.") }

  if (spm_execute_glm && !spm_execute_setup) {
    message("Because spm_execute_glm is TRUE, I will set spm_execute_setup to TRUE, too.")
    spm_execute_setup <- TRUE #running the GLM depends on first setting up the SPM.mat object
  }

  if (spm_execute_contrasts && !spm_execute_setup) {
    message("Because spm_execute_contrasts is TRUE, I will set spm_execute_setup to TRUE, too.")
    spm_execute_setup <- TRUE #running the contrasts depends on first setting up the SPM.mat object
  }
    
  if (spm_execute_contrasts && !spm_execute_glm) {
    message("Because spm_execute_contrasts is TRUE, I will set spm_execute_glm to TRUE, too.")
    spm_execute_glm <- TRUE #running the contrasts depends on first estimating the model
  }

  method_label <- normalize_spm_estimation_method(estimation_method)
  write_residuals <- normalize_spm_write_residuals(write_residuals, method_label)

  if (is.null(fmri_t) || is.null(fmri_t0)) {
    stop("fmri_t and fmri_t0 must be provided (use spm_l1_model to auto-infer).")
  }
  if (!is.numeric(fmri_t) || length(fmri_t) != 1L || is.na(fmri_t)) {
    stop("fmri_t must be a single numeric value.")
  }
  if (!is.numeric(fmri_t0) || length(fmri_t0) != 1L || is.na(fmri_t0)) {
    stop("fmri_t0 must be a single numeric value.")
  }
  fmri_t <- as.integer(round(fmri_t))
  fmri_t0 <- as.integer(round(fmri_t0))
  if (fmri_t < 1L) stop("fmri_t must be >= 1.")
  if (fmri_t0 < 1L) stop("fmri_t0 must be >= 1.")
  if (fmri_t0 > fmri_t) stop("fmri_t0 must be <= fmri_t.")

  if (is.null(cvi) || !nzchar(cvi)) {
    stop("cvi must be one of: 'AR(1)', 'FAST', or 'none'.")
  }
  cvi_norm <- tolower(cvi)
  if (cvi_norm %in% c("ar(1)", "ar1")) {
    cvi <- "AR(1)"
  } else if (cvi_norm == "fast") {
    cvi <- "FAST"
  } else if (cvi_norm == "none") {
    cvi <- "none"
  } else {
    stop("cvi must be one of: 'AR(1)', 'FAST', or 'none'.")
  }
  
  #extract key ingredients from bdm object  
  if (concatenate_runs) {
    design <- bdm$design_concat #already has events concatenated into one run
  } else {
    design <- bdm$design
  }
  
  run_niftis <- bdm$run_niftis
  tr <- bdm$tr
  
  if (inherits(design, "list")) {
    #in case of 1-run input, convert 1D list to an array with 1 row and nregressors columns
    design <- array(design, dim=c(1, length(design)), dimnames=list(c("run1"), names(design)))
  }

  if (is.null(tr)) { stop("You must specify a TR in seconds") }
  if (!concatenate_runs && length(run_niftis) != nrow(design)) { stop("Length of run_niftis argument is not equal to rows in design") }
  
  stopifnot(all(file.exists(run_niftis)))
  run_niftis <- normalizePath(run_niftis, mustWork = TRUE)

  if (!is.null(ts_files)) {
    stopifnot(all(file.exists(ts_files)))
    if (isTRUE(concatenate_runs) && length(ts_files) != 1L) {
      stop("For concatenated runs, ts_files must have length 1 (a single concatenated file).")
    } else if (!isTRUE(concatenate_runs) && length(ts_files) != length(run_niftis)) {
      stop("Length of ts_files argument is not equal to length of run_niftis.")
    }
  }

  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  output_dir <- normalizePath(output_dir) #convert to absolute path to avoid any ambiguity
  
  nruns <- dim(design)[1]
  nregressors <- dim(design)[2]
  
  aligned <- spm_align_inputs(bdm, run_niftis, ts_files, concatenate_runs, output_dir)
  run_lengths <- aligned$lengths
  scan_offsets <- aligned$offsets
  ts_files <- aligned$ts_files

  run_numbers <- if (!is.null(bdm$runs_to_output)) as.integer(bdm$runs_to_output) else seq_along(run_niftis)
  if (length(run_numbers) != length(run_niftis)) {
    stop("Length of bdm$runs_to_output must match run_niftis for SPM setup.")
  }

  spm_l2_projection <- NULL
  spm_l2_projection_cols <- character(0)
  spm_l2_projection_interactions <- character(0)
  projection_interaction_terms <- character(0)
  if (!is.null(bdm$spm_l2_projection)) {
    if (!is.data.frame(bdm$spm_l2_projection)) {
      stop("bdm$spm_l2_projection must be a data.frame with run_number and regressor columns.")
    }
    if (!"run_number" %in% names(bdm$spm_l2_projection)) {
      stop("bdm$spm_l2_projection must contain a run_number column.")
    }

    spm_l2_projection <- bdm$spm_l2_projection
    if (anyDuplicated(spm_l2_projection$run_number) > 0L) {
      stop("bdm$spm_l2_projection contains duplicated run_number rows.")
    }

    projection_idx <- match(run_numbers, spm_l2_projection$run_number)
    if (anyNA(projection_idx)) {
      missing_runs <- run_numbers[is.na(projection_idx)]
      stop(
        "bdm$spm_l2_projection is missing run_number row(s): ",
        paste(missing_runs, collapse = ", ")
      )
    }
    spm_l2_projection <- spm_l2_projection[projection_idx, , drop = FALSE]
    spm_l2_projection_cols <- setdiff(names(spm_l2_projection), "run_number")

    if (length(spm_l2_projection_cols) > 0L) {
      spm_l2_projection[spm_l2_projection_cols] <- lapply(spm_l2_projection[spm_l2_projection_cols], as.numeric)
      bad_cols <- spm_l2_projection_cols[vapply(
        spm_l2_projection[spm_l2_projection_cols],
        function(x) any(!is.finite(x)),
        logical(1)
      )]
      if (length(bad_cols) > 0L) {
        stop(
          "bdm$spm_l2_projection contains non-finite values in: ",
          paste(bad_cols, collapse = ", ")
        )
      }
    } else {
      spm_l2_projection <- NULL
    }
  }

  if (!is.null(spm_l2_projection) && !is.null(bdm$spm_l2_projection_interactions)) {
    if (is.logical(bdm$spm_l2_projection_interactions)) {
      if (length(bdm$spm_l2_projection_interactions) != 1L || is.na(bdm$spm_l2_projection_interactions)) {
        stop("bdm$spm_l2_projection_interactions must be scalar logical or character.")
      }
      if (isTRUE(bdm$spm_l2_projection_interactions)) {
        spm_l2_projection_interactions <- spm_l2_projection_cols
      }
    } else if (is.character(bdm$spm_l2_projection_interactions)) {
      vals <- unique(bdm$spm_l2_projection_interactions[!is.na(bdm$spm_l2_projection_interactions) & nzchar(bdm$spm_l2_projection_interactions)])
      if (length(vals) > 0L) {
        bad <- setdiff(vals, spm_l2_projection_cols)
        if (length(bad) > 0L) {
          stop(
            "Unknown requested interaction term(s) in bdm$spm_l2_projection_interactions: ",
            paste(bad, collapse = ", "),
            ". Available: ",
            paste(spm_l2_projection_cols, collapse = ", ")
          )
        }
        spm_l2_projection_interactions <- vals
      }
    } else {
      stop("bdm$spm_l2_projection_interactions must be scalar logical or character.")
    }
  }

  skip_projection_main_effects <- FALSE
  if (!is.null(spm_l2_projection) && length(spm_l2_projection_cols) > 0L && !isTRUE(concatenate_runs)) {
    # In non-concatenated SPM, per-run projected covariates are constant within
    # each Sn() block and therefore aliased with session constants.
    skip_projection_main_effects <- TRUE
    warning(
      "Skipping projected L2 main-effect regressors in non-concatenated SPM design because they are aliased with session constants.",
      call. = FALSE
    )
  }

  dropped_constant_interactions <- character(0)

  # spm only works with .nii files, not .nii.gz
  gzipped <- grepl(".nii.gz$", run_niftis)
  gunzip_cmds <- character(0)
  if (any(gzipped)) {
    dir.create(nifti_tmpdir, showWarnings = FALSE, recursive = TRUE)
    nifti_tmpdir <- normalizePath(nifti_tmpdir)
    # Use deterministic names in the nifti_tmpdir so the cluster job can find/create them reliably
    # Isolate temporary images by analysis and run, even for repeated basenames.
    scope_file <- tempfile()
    writeLines(output_dir, scope_file)
    scope <- unname(tools::md5sum(scope_file))
    unlink(scope_file)
    nifti_tmpdir <- file.path(nifti_tmpdir, paste0("spm_", scope))
    dir.create(nifti_tmpdir, showWarnings = FALSE, recursive = TRUE)
    writeLines(nifti_tmpdir, file.path(output_dir, ".nifti_tmpdir"))
    target_niftis <- vapply(seq_along(run_niftis), function(i) {
      x <- run_niftis[i]
      if (grepl("\\.nii\\.gz$", x)) {
        file.path(nifti_tmpdir, paste0(sprintf("run_%04d_", i), sub("\\.gz$", "", basename(x))))
      } else x
    }, character(1))

    # Record the commands to be run either now or on the cluster
    gunzip_cmds <- c(paste("mkdir -p", shQuote(nifti_tmpdir)), paste0("gunzip -c ", shQuote(run_niftis[gzipped]), " > ", shQuote(target_niftis[gzipped])))

    if (spm_execute_setup || spm_execute_glm) {
      for (cmd in gunzip_cmds) {
        message("Executing: ", cmd)
        if (system(cmd) != 0L) stop("SPM image decompression failed.")
      }
    }
    cat(c("#!/bin/bash", "set -euo pipefail", gunzip_cmds), file = file.path(output_dir, "gunzip_commands.sh"), sep = "\n")
    run_niftis <- target_niftis
  } else {
    unlink(file.path(output_dir, c("gunzip_commands.sh", ".nifti_tmpdir")))
  }

  spm_preamble <- c(
    ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
    "spm('defaults', 'fmri');",
    "spm_jobman('initcfg');",
    ""    
  )

  baseobj <- paste0("matlabbatch{1}.spm.stats.fmri_spec")
  escape_spm_string <- function(x) gsub("'", "''", as.character(x), fixed = TRUE)
  format_spm_values <- function(x) {
    paste(format(as.numeric(x), trim = TRUE, scientific = FALSE, digits = 15), collapse = " ")
  }
  sanitize_term_name <- function(x) {
    nm <- make.names(x)
    nm <- gsub("\\.+", ".", nm)
    nm <- gsub("^\\.|\\.$", "", nm)
    nm
  }

  if (hrf_derivs == "none") {
    hrf_string <- paste0(baseobj, ".bases.hrf = struct('derivs', [0 0]);")
  } else if (hrf_derivs == "time") {
    hrf_string <- paste0(baseobj, ".bases.hrf = struct('derivs', [1 0]);")
  } else if (hrf_derivs == "time_dispersion") {
    hrf_string <- paste0(baseobj, ".bases.hrf = struct('derivs', [1 1]);")
  } else { stop("don't understand hrf_derivs argument: ", hrf_derivs) }

  m_string <- c(
    paste0("matlabbatch = []; %initialize empty structure"),
    paste0(baseobj, ".dir = {'", output_dir, "'};"),
    "",
    "% timing",
    paste0(baseobj, ".timing.units = 'secs';"),
    paste0(baseobj, ".timing.RT = ", tr, ";"),
    paste0(baseobj, ".timing.fmri_t = ", fmri_t, "; % microtime resolution"),
    paste0(baseobj, ".timing.fmri_t0 = ", fmri_t0, "; % alignment within TR"),
    #"% factorial",
    paste0(baseobj, ".fact = struct('name', {}, 'levels', {});"),
    "% basis functions",
    hrf_string,
    "% volterra",
    paste0(baseobj, ".volt = 1;"),
    "% global",
    paste0(baseobj, ".global = 'None';"),
    "% mthresh",
    paste0(baseobj, ".mthresh = 0.8000;"),
    "% mask",
    paste0(baseobj, ".mask = {''};"),
    "% cvi",
    paste0(baseobj, ".cvi = '", cvi, "';")
  )

  for (rr in 1:nruns) {
    if (concatenate_runs) {
      run_vol_concat <- cumsum(run_lengths)
      stopifnot(nruns==1) #should be using a concat design that sees only one run
      scans_string <- c()
      for (nn in 1:length(run_niftis)) {
        scans_string <- c(
          scans_string,
          paste0("% Run ", nn, " scans"),
          paste0("Nt = ", run_lengths[nn], ";"),
          paste0("Nt_start = ", ifelse(nn > 1, run_vol_concat[nn-1], 0), ";"),
          "for i = 1:Nt",
          paste0("  ", baseobj, ".sess(", rr, ").scans{i+Nt_start,1} = [ '", run_niftis[nn], "' ',' num2str(i + ", scan_offsets[nn], ") ];"),
          "end",
          ""
        )
      }
    } else {
      scans_string <- c(
        paste0("Nt = ", run_lengths[rr], ";"),
        "for i = 1:Nt",
        paste0("  ", baseobj, ".sess(", rr, ").scans{i,1} = [ '", run_niftis[rr], "' ',' num2str(i + ", scan_offsets[rr], ") ];"),
        "end"
      )
    }

    m_string <- c(m_string,
      paste0(baseobj, ".sess(", rr, ").scans = {};"),
      scans_string,
      paste0(baseobj, ".sess(", rr, ").multi = {''};"),
      paste0(baseobj, ".sess(", rr, ").regress = struct('name', {}, 'val', {});"),
      paste0(baseobj, ".sess(", rr, ").hpf = ", hpf, ";")
    )

    #add additional time series for this run
    if (!is.null(ts_files)) {
      m_string <- c(m_string, paste0(baseobj, ".sess(", rr, ").multi_reg = { '", ts_files[rr], "' };") ) } else {
        m_string <- c(m_string, paste0(baseobj, ".sess(", rr, ").multi_reg = { '' };"))
      }

    if (!isTRUE(skip_projection_main_effects) && !is.null(spm_l2_projection) && length(spm_l2_projection_cols) > 0L) {
      if (isTRUE(concatenate_runs)) {
        projection_values <- lapply(spm_l2_projection_cols, function(cc) {
          unlist(lapply(seq_along(run_lengths), function(run_idx) {
            rep(spm_l2_projection[[cc]][run_idx], run_lengths[run_idx])
          }), use.names = FALSE)
        })
      } else {
        projection_values <- lapply(spm_l2_projection_cols, function(cc) {
          rep(spm_l2_projection[[cc]][rr], run_lengths[rr])
        })
      }

      for (kk in seq_along(spm_l2_projection_cols)) {
        m_string <- c(
          m_string,
          paste0(
            baseobj, ".sess(", rr, ").regress(", kk, ").name = '",
            escape_spm_string(spm_l2_projection_cols[kk]), "';"
          ),
          paste0(
            baseobj, ".sess(", rr, ").regress(", kk, ").val = [ ",
            format_spm_values(projection_values[[kk]]), " ]';"
          )
        )
      }
    }

    #collect regressors that are aligned to the same event
    rmat <- design[rr,]

    event_alignment <- sapply(rmat, function(x) { attr(x, "event") })
    uevents <- unique(event_alignment)

    ulist <- list()

    for (u in uevents) {
      ureg <- rmat[which(event_alignment == u)]

      non_empty <- vapply(ureg, function(x) nrow(x) > 0L, logical(1))
      if (!any(non_empty)) {
        warning("No events found for alignment: ", u, " in run ", rr, ". Skipping.")
        next
      }
      if (!all(non_empty)) {
        warning("Dropping empty regressors for alignment: ", u, " in run ", rr, ".")
        ureg <- ureg[non_empty]
      }

      reg_lengths <- vapply(ureg, nrow, integer(1))
      which_1 <- vapply(ureg, function(x) { all(abs(x[,"value"] - 1.0) < 1e-3, na.rm=TRUE) }, logical(1))

      #SPM does not like NAs or NaNs in parametric modulators.
      #It always generates both the event-ness regressor (based on the onsets vector) and the parametric modulator (modulated by pmod)
      #Therefore, we need to handle NAs in the modulators. If I were using SPM for a GLM, I might just handle the convolution myself and hand off a convolved signal.
      #But because the immediate goal here is to have a design matrix for DCM, I will remove the NA trials first -- for both event-ness and modulator --
      #And put a separate 'eventness' regressor for trials when an event occurred, but there is no modulation value
      
      if (length(ureg) > 1L && (all(which_1) || length(unique(reg_lengths)) > 1L)) {
        # Multiple unit-height regressors or mismatched lengths: treat each as its own condition.
        for (reg_name in names(ureg)) {
          reg <- ureg[[reg_name]]
          reg <- reg[!is.na(reg[,"value"]), , drop = FALSE]
          if (nrow(reg) == 0L) next
          ulist[[reg_name]][["event"]] <- reg
        }
        next
      }

      uvalues <- do.call(cbind, lapply(ureg, function(r) { r[,"value"] }))
      which_na <- apply(uvalues, 1, function(r) { any(is.na(r)) }) # missing values for any event (row)?
      parametric <- !which_1

      if (length(ureg) > 1L) {
        #determine eventness and parametric modulator(s)
        utimes <- do.call(cbind, lapply(ureg, function(r) { r[,"onset"] }))
        stopifnot(all(na.omit(apply(utimes, 1, function(x) { sd(x) }) < 1e-3))) #make sure that all times are the same (within rounding error)

        if (sum(which_1) > 1L) {
          stop("More than one unit-height regressor for event: ", u)
        } else if (sum(which_1) == 0L) {
          stop("No unit-height regressor for event: ", u)
        }

        ulist[[u]][["event"]] <- ureg[[ names(ureg)[which_1] ]][!which_na,,drop=FALSE]
        ulist[[u]][ names(ureg)[!which_1] ] <- lapply(ureg[ names(ureg)[!which_1] ], function(reg) { reg[!which_na,,drop=FALSE] })
      } else {
        ulist[[u]][["event"]] <- ureg[[1]][!which_na,,drop=FALSE] #just copy this across -- only one regressor on this alignment
        if (parametric[1]) { warning("Only one regressor with alignment: ", u, ", but it appears to be parametric. SPM will add the eventness regressor, too.") }
      }

      #add eventness regressor for only missing trials, if needed
      if (any(which_na)) {
        ulist[[paste0(u, "_nopmod")]][["event"]] <- ureg[[ names(ureg)[which_1] ]][which_na,,drop=FALSE]
      }
      
    }

    if (length(spm_l2_projection_interactions) > 0L) {
      run_start_times <- if (isTRUE(concatenate_runs)) {
        c(0, head(cumsum(run_lengths) * tr, -1))
      } else {
        NULL
      }

      for (cc in seq_along(ulist)) {
        cond_name <- names(ulist)[cc]
        event_df <- ulist[[cc]][["event"]]
        if (is.null(event_df) || nrow(event_df) == 0L) next

        if (isTRUE(concatenate_runs)) {
          run_idx <- findInterval(event_df[, "onset"], run_start_times)
          if (any(run_idx < 1L | run_idx > length(run_lengths))) {
            warning("Some concatenated event onsets were outside run boundaries; clamping projected interaction run indices.")
            run_idx <- pmin(pmax(run_idx, 1L), length(run_lengths))
          }
        } else {
          run_idx <- rep(rr, nrow(event_df))
        }

        for (covar in spm_l2_projection_interactions) {
          int_name <- sanitize_term_name(paste0(cond_name, "_x_l2_", covar))
          int_vals <- spm_l2_projection[[covar]][run_idx]
          if (any(!is.finite(int_vals))) {
            stop("Non-finite projected interaction values for term: ", int_name)
          }
          if (length(int_vals) > 0L && all(abs(int_vals - int_vals[1L]) < 1e-8)) {
            dropped_constant_interactions <- unique(c(dropped_constant_interactions, int_name))
            next
          }
          int_df <- event_df
          int_df[, "value"] <- as.numeric(int_vals)
          ulist[[cc]][[int_name]] <- int_df
          projection_interaction_terms <- unique(c(projection_interaction_terms, int_name))
        }
      }
    }

    for (cc in 1:length(ulist)) {
      #event-only lists have length 1
      if (length(ulist[[cc]]) > 1L) {
        pmod_string <- c()
        for (pm in 2:length(ulist[[cc]])) {
          pmod_string <- c(pmod_string, paste0(baseobj, ".sess(", rr, ").cond(", cc, ").pmod(", pm-1, ") = struct('name', {'", names(ulist[[cc]])[pm], "-pmod'}, 'param', { [ ",
            paste(ulist[[cc]][[pm]][,"value"], collapse=", "), " ] }, 'poly', { 1 });"))
        }
      } else {
        pmod_string <- paste0(baseobj, ".sess(", rr, ").cond(", cc, ").pmod = struct('name', {}, 'param', {}, 'poly', {});")
      }
      
      m_string <- c(m_string,
        paste0(baseobj, ".sess(", rr, ").cond(", cc, ").name = '", names(ulist)[cc], "';"),
        paste0(baseobj, ".sess(", rr, ").cond(", cc, ").onset = [ ", paste(ulist[[cc]][["event"]][,"onset"], collapse=", "), " ];"),
        paste0(baseobj, ".sess(", rr, ").cond(", cc, ").duration = [ ", paste(ulist[[cc]][["event"]][,"duration"], collapse=", "), " ];"),
        paste0(baseobj, ".sess(", rr, ").cond(", cc, ").tmod = 0;"),
        pmod_string,
        paste0(baseobj, ".sess(", rr, ").cond(", cc, ").orth = 1;"),
        "",
        ""
      )

    }
  }

  if (length(dropped_constant_interactions) > 0L) {
    message(
      paste0(
        "Skipping projected interaction term(s) with no within-condition variation: ",
        paste(dropped_constant_interactions, collapse = ", "),
        "."
      )
    )
  }

  #Generate syntax for executing SPM GLM setup
  exec_string <- c(
    spm_preamble,
    paste0("if exist('", file.path(output_dir, "SPM.mat"), "', 'file'), delete('", file.path(output_dir, "SPM.mat"), "'); end;"), #remove SPM.mat if it exists
    paste0("run('", file.path(output_dir, "glm_design_batch.m"), "');"), #source the settings for this GLM
    "% RUN DESIGN MATRIX JOB",
    paste0("spm_jobman('run',matlabbatch);")
  )

  if (concatenate_runs) {
    exec_string <- c(exec_string, paste0("spm_fmri_concatenate( [ '", output_dir, "' filesep 'SPM.mat'], [ ", paste(run_lengths, collapse=", "), " ]);"))
  }
  
  #write GLM setup to file
  cat(m_string, file=file.path(output_dir, "glm_design_batch.m"), sep="\n")
  cat(exec_string, file=file.path(output_dir, "setup_glm_design.m"), sep="\n")
  spm_syntax[["glm_design"]] <- m_string
  spm_syntax[["exec_glm_setup"]] <- exec_string

  build_matlab_call <- function(script_path) {
    cmd_str <- paste0(
      "try; run('", script_path, "'); ",
      "catch ME; disp(getReport(ME,'extended')); exit(1); end; exit(0);"
    )
    paste(matlab_cmd, matlab_args, shQuote(cmd_str))
  }
  build_shell_call <- function(script_path) {
    cmd <- build_matlab_call(script_path)
    if (!is.null(compute_env) && length(compute_env) > 0L) {
      return(paste(c(compute_env, cmd), collapse = " && "))
    }
    cmd
  }

  if (spm_execute_setup) {
    if (system(build_shell_call(file.path(output_dir, "setup_glm_design.m"))) != 0L) stop("SPM design setup failed.")
  } else {
    message("To run GLM setup script, execute this command: ", build_shell_call(file.path(output_dir, "setup_glm_design.m")))
  }

  #SPM model estimation setup
  baseobj <- paste0("matlabbatch{1}.spm.stats.fmri_est")
  m_string <- c(
    spm_preamble,
    paste0("matlabbatch = []; %initialize empty structure"),
    "% ESTIMATE MODEL",
    "",
    paste0(baseobj, ".spmmat = { [ '", output_dir, "' filesep 'SPM.mat']};"),
    "% write_residuals",
    paste0(baseobj, ".write_residuals = ", ifelse(isTRUE(write_residuals), 1, 0), ";"),
    "% method",
    paste0(baseobj, ".method.", method_label, " = 1;"),
    "% RUN MODEL ESTIMATION JOB",
    "spm_jobman('run',matlabbatch);"
  )

  #write estimation to file
  cat(m_string, file=file.path(output_dir, "run_glm.m"), sep="\n")
  spm_syntax[["run_glm"]] <- m_string  

  if (spm_execute_glm) {
    if (system(build_shell_call(file.path(output_dir, "run_glm.m"))) != 0L) stop("SPM estimation failed.")
  } else {
    message("To estimate GLM after design matrix setup, execute this command: ", build_shell_call(file.path(output_dir, "run_glm.m")))
  }

  spm_contrast_cmds <- generate_spm_contrasts(output_dir, condition_contrasts, unit_contrasts,
    effects_of_interest_F, spm_path, execute=spm_execute_contrasts, 
    matlab_cmd = matlab_cmd, matlab_args = matlab_args)
  
  if (!spm_execute_contrasts && !is.null(spm_contrast_cmds)) {
    cat(c("#!/bin/bash", "set -euo pipefail", spm_contrast_cmds$extract_cmd, spm_contrast_cmds$setup_cmd),
      file = file.path(output_dir, "setup_spm_contrasts.sh"), sep = "\n")
  }
  
  spm_syntax[["gunzip_cmds"]] <- gunzip_cmds
  spm_syntax[["contrast_cmds"]] <- spm_contrast_cmds
  spm_syntax[["projection_interaction_terms"]] <- projection_interaction_terms
  
  if (cleanup_tmp && (spm_execute_glm || spm_execute_setup)) {
    if (any(gzipped)) {
      sapply(run_niftis[gzipped], function(x) { unlink(x) })
    }
  }
  
  return(spm_syntax)
}
