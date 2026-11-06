#' Function to perform voxelwise deconvolution on an fMRI dataset using the fMRI model arguments object
#'
#' @param niftis A vector of processed fMRI timeseries images (4D files) to be deconvolved
#' @param add_metadata A data.frame with one row per value of \code{niftis}, saved
#'   as atlas-level metadata. The \code{id}, \code{run}, and optional \code{session}
#'   columns are also saved in per-file sidecars for event alignment.
#' @param out_dir Base output directory for deconvolved time series files. Default is \code{getwd()}.
#' @param out_file_expression Expression evaluated to resolve the filename for the deconvolved csv files. Default is
#'   to prefix \code{"decon_"}, replace slashes in \code{niftis} with periods, and add
#'   the atlas image name. Note that the suffix \code{_deconvolved.csv.gz} or \code{_original.csv.gz} will be added
#'   to the expression, so don't pass this piece.
#' @param log_file Name (and path) of log file for any deconvolution errors or messages
#' @param TR the repetition time of the sequence in seconds. Required
#' @param time_offset The number of seconds that will be subtracted or added to the time field. Default: 0. Useful if some number of volumes
#'   have been dropped from the NIfTI data prior to deconvolution.
#' @param atlas_files optional atlas file vector specifying voxels used in deconvolution. If omitted, perform whole-brain deconvolution
#' @param mask an optional character string specifying a mask that should be used to constrain bounds of deconvolution.
#' @param nprocs The number of processors to use simultaneously for deconvolution
#' @param save_original_ts Whether to save the voxelwise BOLD data prior to deconvolution (for comparison/diagnosis). Default: TRUE
#' @param algorithm Which deconvolution algorithm to use for deconvolving voxelwise time series. Default: "bush2011". Alternative is
#'   "bush2015", which implements a resampling approach as well, or "reglin"/"regularized_linear",
#'   which estimates a continuous latent activity time series using regularized linear deconvolution
#'   with optional HRF tuning.
#' @param decon_settings A list of settings passed to the deconvolution algorithm.
#' @param afni_dir Retained for compatibility; ignored. NIfTI extraction is internal
#'   and does not require AFNI.
#'
#' @details The Bush 2011 algorithm uses the package's internal RcppArmadillo
#'   implementation and does not require an external executable.
#'   Atlas, mask, and input images must share a spatial grid. Constant signals
#'   and signals containing nonfinite observations are excluded before fitting;
#'   their output rows are retained with missing estimates. The \code{decon_status}
#'   column identifies \code{"ok"}, \code{"constant"}, or \code{"nonfinite"} signals.
#'   Each deconvolution file has a \code{.metadata.csv} sidecar containing its
#'   source NIfTI path and any \code{id}, \code{run}, and \code{session} columns
#'   in \code{add_metadata}. These identify files during event alignment, including
#'   custom filenames. If absent, alignment falls back to the source path or
#'   requires an explicit \code{file_manifest}.
#'
#'   If you want to use subject metadata to name the output file, use \code{this_subj} in your \code{out_file_expression}, which will give you access to a one-row
#'   data.frame containing the metadata for the current subject in the loop.
#'
#' @examples
#'   \dontrun{
#'
#'     #name outputs according to subject metadata
#'     xx <- voxelwise_deconvolution(
#'       niftis="/proj/mnhallqlab/user/michael/test_nifti.nii.gz",
#'       out_dir="/proj/mnhallqlab/user/michael/decon_outputs",
#'       out_file_expression = expression(
#'         paste0(
#'           this_subj$subid,
#'           "_run",
#'           this_subj$run_num,
#'           "_",
#'           atlas_img_name
#'         )
#'       )
#'     )
#'   }
#'
#' @return Nothing (invisible NULL).
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom checkmate assert_file_exists
#' @importFrom foreach %do% foreach
#' @importFrom dplyr mutate across select left_join
#' @importFrom readr write_csv
#' @importFrom stats setNames
#' @export
voxelwise_deconvolution <- function(
  niftis, add_metadata=NULL, out_dir=getwd(), out_file_expression=NULL,
  log_file=file.path(out_dir, "deconvolve_errors"),
  TR = NULL, time_offset=0, atlas_files=NULL, mask=NULL, nprocs=20, save_original_ts=TRUE, algorithm="bush2011",
  decon_settings=list(nev_lr = .01, #neural events learning rate (default in algorithm)
    epsilon = .005, #convergence criterion (default)
    beta = 60, #best from Bush 2015 update
    kernel = spm_hrf(TR)$hrf, #canonical SPM difference of gammas
    n_resample = 25), afni_dir=NULL) { #for Bush 2015 only

  sapply(niftis, checkmate::assert_file_exists)
  checkmate::assert_data_frame(add_metadata, nrows=length(niftis), null.ok=TRUE)
  sapply(atlas_files, checkmate::assert_file_exists)
  checkmate::assert_numeric(TR, lower=0.01)
  algorithm <- match.arg(
    algorithm,
    c("bush2011", "bush2015", "reglin", "regularized_linear")
  )

  if (!is.null(decon_settings$kernel)) {
    hrf_pad <- length(decon_settings$kernel)
  } else {
    stop("At present no alternative way to set hrf_pad if $kernel not passed in decon_settings")
  }

  if (is.null(atlas_files)) {
    message("Performing whole-brain voxelwise deconvolution")
    if (is.null(mask)) {
      stop("Cannot conduct voxelwise deconvolution without a relevant mask")
    } else {
      checkmate::assert_file_exists(mask)
    }

    atlas_files <- mask
  }
  
  zero_thresh <- 1e-4 #for binarizing/indexing
  atlas_imgs <- lapply(atlas_files, RNifti::readNifti)
  mask_file <- mask
  if (!is.null(mask)) {
    checkmate::assert_file_exists(mask)
    mask <- RNifti::readNifti(mask)
    
    #ensure that mask is binary
    mask <- is.finite(mask) & mask >= zero_thresh
  }
  
  #setup cluster
  if (nprocs > 1) {
    cl <- makeCluster(nprocs)
    registerDoParallel(cl)
    on.exit(try(stopCluster(cl)))
  } else {
    registerDoSEQ()
  }

  #loop over atlas files
  for (ai in seq_along(atlas_files)) {
    cat("Working on atlas: ", atlas_files[ai], "\n")
    aimg <- atlas_imgs[[ai]]

    atlas_geometry <- roi_image_geometry(atlas_files[ai])
    if (atlas_geometry$volumes != 1L) stop("Deconvolution atlases must be single-volume images")
    if (!is.null(mask_file)) {
      validate_roi_geometry(atlas_geometry, roi_image_geometry(mask_file), atlas_files[ai], mask_file)
    }
    
    if (!is.null(mask)) { aimg <- aimg*mask } #multiply against binary mask to apply it
    
    
    # One shared index vector defines both metadata and extraction order.
    # NIfTI and R arrays both store i first, then j, then k within each volume.
    voxel_indices <- which(is.finite(aimg) & aimg > zero_thresh)
    if (length(voxel_indices) == 0L) stop("No voxels remain in atlas after masking: ", atlas_files[ai])
    a_indices <- arrayInd(voxel_indices, dim(aimg)[seq_len(min(3L, length(dim(aimg))))])
    if (ncol(a_indices) == 2L) a_indices <- cbind(a_indices, 1L)
    xyz <- cbind(a_indices - 1L, 1) %*% t(atlas_geometry$affine)
    a_coordinates <- cbind(a_indices, xyz[, 1:3, drop = FALSE])
    
    a_coordinates <- as.data.frame(a_coordinates) %>%
      setNames(c("i", "j", "k", "x", "y", "z")) %>%
      dplyr::mutate(vnum = seq_along(voxel_indices), atlas_value = as.numeric(aimg[voxel_indices]), atlas_name = basename(atlas_files[ai])) %>%
      mutate(across(c("x", "y", "z"), ~round(.x, 2))) %>%
      dplyr::select("vnum", "atlas_value", everything())
    
    #setup output subdirectories for deconvolved files, named according to atlas
    atlas_img_name <- basename(sub(".nii(.gz)*", "", atlas_files[ai], perl=TRUE))
    dir.create(file.path(out_dir, atlas_img_name, "deconvolved"), showWarnings=FALSE, recursive=TRUE)
    if (isTRUE(save_original_ts)) { dir.create(file.path(out_dir, atlas_img_name, "original"), showWarnings=FALSE, recursive=TRUE) }
    
    
    #set a default filename based on the nifti -- hopefully this will not collide with other subjects/runs
    if (is.null(out_file_expression)) {
      out_file_expression <- expression(paste0("decon_", gsub("[/\\]", ".", niftis[si]), "_", atlas_img_name))
    }

    if (!is.null(add_metadata)) {
      # add nifti as a metadata output
      add_metadata$.nifti <- niftis
    }

    #loop over niftis in parallel
    ff <- foreach(si = seq_along(niftis),
      .packages=c("dplyr", "readr", "data.table", "fmri.pipeline", "foreach", "iterators"),
      .export=c("mat2df", "extract_nifti_voxel_ts", "roi_image_geometry",
        "validate_roi_geometry", "decon_signal_status")) %dopar% {

      #get the si-th row of the metadata to match nifti, allow one to use this_subj in out_file_expression
      if (!is.null(add_metadata)) {
        this_subj <- add_metadata %>% dplyr::slice(si)
      }

      out_name <- file.path(out_dir, atlas_img_name, "deconvolved", paste0(eval(out_file_expression), "_deconvolved.csv.gz"))

      if (file.exists(out_name)) {
        message("Deconvolved file already exists: ", out_name)
        return(NULL)
      }
      
      cat("  Deconvolving subject: ", niftis[si], "\n")
      # Only spatial geometry must match; the input may contain many volumes.
      input_geometry <- roi_image_geometry(niftis[si])
      input_geometry$volumes <- 1L
      validate_roi_geometry(atlas_geometry, input_geometry, atlas_files[ai], niftis[si])

      # to_deconvolve is a voxels x time matrix
      to_deconvolve <- extract_nifti_voxel_ts(niftis[si], voxel_indices)

      # Screen before scaling: a constant signal would otherwise become NaN and
      # contaminate global tuning. Keep invalid voxels in the output as NA.
      decon_status <- decon_signal_status(t(to_deconvolve))
      valid_voxels <- decon_status == "ok"
      scaled_input <- matrix(NA_real_, nrow(to_deconvolve), ncol(to_deconvolve))
      if (any(valid_voxels)) {
        scaled_input[valid_voxels, ] <- t(scale(t(to_deconvolve[valid_voxels, , drop = FALSE])))
      }
      to_deconvolve <- scaled_input
      if (any(!valid_voxels)) {
        warning(sum(!valid_voxels), " constant or nonfinite voxel signals excluded from deconvolution: ", niftis[si])
      }

      # just demean, which will rescale to percent signal change around 0 (this matches Bush 2015)
      #to_deconvolve <- t(apply(to_deconvolve, 1, function(x) { scale(x, scale=FALSE) }))

      # pct signal change around 0
      #to_deconvolve <- t(apply(to_deconvolve, 1, function(x) { x/mean(x)*100 - 100 }))

      # Zero pad the tail end to keep the deconvolved signal stable at the
      # boundary and trim the same padding after deconvolution.
      alg_input <- to_deconvolve[valid_voxels, , drop = FALSE] %>%
        cbind(matrix(0, nrow = sum(valid_voxels), ncol = hrf_pad)) %>%
        as.matrix()

      #test1 <- deconvolve_nlreg(to_deconvolve[117,], kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon)
      #test2 <- deconvolve_nlreg(to_deconvolve[118,], kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon)
      
      if (!any(valid_voxels)) {
        deconv_mat <- matrix(NA_real_, nrow(to_deconvolve), ncol(to_deconvolve))
      } else {
        if (algorithm == "bush2015") {
          #use R implementation of Bush 2015 algorithm
          deconv_mat <- foreach(vox_ts=iter(alg_input, by="row"), .combine="rbind", .packages=c("fmri.pipeline")) %do% {
            reg <- tryCatch(deconvolve_nlreg_resample(as.vector(vox_ts), kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon, n_resample=decon_settings$n_resample),
              error=function(e) { cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file=log_file, append=TRUE); return(rep(NA, length(vox_ts))) })

            if (is.list(reg)) { reg <- reg$NEVmean } #just keep the mean resampled events vector
            return(reg)
          }
        } else if (algorithm %in% c("reglin", "regularized_linear")) {
          # Regularized linear deconvolution estimates continuous latent activity without
          # sigmoid bounding. R function is time x voxels, matching the Rcpp path.
          alg_input <- t(alg_input)
          deconv_mat <- tryCatch(
            deconvolve_reglin(
              BOLDobs = alg_input,
              kernel = decon_settings$kernel,
              lambda = decon_settings[["lambda"]],
              lambda_grid = if (is.null(decon_settings$lambda_grid)) 10^seq(-4, 2, length.out = 20) else decon_settings$lambda_grid,
              penalty = if (is.null(decon_settings$penalty)) "diff2" else decon_settings$penalty,
              kernel_grid = decon_settings$kernel_grid,
              hrf_lags = if (is.null(decon_settings$hrf_lags)) 0L else decon_settings$hrf_lags,
              tune_by = if (is.null(decon_settings$tune_by)) "global" else decon_settings$tune_by,
              normalize = if (is.null(decon_settings$normalize)) FALSE else decon_settings$normalize,
              demean = if (is.null(decon_settings$demean)) TRUE else decon_settings$demean,
              prewhiten_gcv = if (is.null(decon_settings$prewhiten_gcv)) FALSE else decon_settings$prewhiten_gcv,
              prewhiten_rho = decon_settings$prewhiten_rho,
              backend = if (is.null(decon_settings$backend)) "cpp" else decon_settings$backend,
              wiener_shrinkage = if (is.null(decon_settings$wiener_shrinkage)) "scalar" else decon_settings$wiener_shrinkage,
              wiener_psd_smooth = if (is.null(decon_settings$wiener_psd_smooth)) 7L else decon_settings$wiener_psd_smooth,
              wiener_psd_floor = if (is.null(decon_settings$wiener_psd_floor)) 1e-8 else decon_settings$wiener_psd_floor,
              gcv_rule = if (is.null(decon_settings$gcv_rule)) "min" else decon_settings$gcv_rule,
              trim_kernel = TRUE
            ),
            error = function(e) {
              cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file = log_file, append = TRUE)
              return(matrix(NA, nrow = nrow(alg_input), ncol = ncol(alg_input)))
            }
          )

          deconv_mat <- t(deconv_mat)
        } else if (algorithm == "bush2011") {
          # This should use the new internal RcppArmadillo function
          # Rcpp function is time x voxels... transpose inputs and outputs to match voxels x time expectations
          alg_input <- t(alg_input)
          deconv_mat <- tryCatch(deconvolve_nlreg(BOLDobs = alg_input, kernel=decon_settings$kernel, nev_lr=decon_settings$nev_lr, epsilon=decon_settings$epsilon, beta=decon_settings$beta),
            error=function(e) {
              cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file=log_file, append=TRUE)
              return(matrix(NA, nrow=nrow(alg_input), ncol=ncol(alg_input)))
            })

          deconv_mat <- t(deconv_mat)
        }
      
        # Trim the HRF end-padding for all variants.
        deconv_mat <- deconv_mat[, seq_len(ncol(to_deconvolve)), drop = FALSE]
        restored <- matrix(NA_real_, nrow(to_deconvolve), ncol(to_deconvolve))
        restored[valid_voxels, ] <- deconv_mat
        deconv_mat <- restored
      }

      # melt this for combination
      deconv_melt <- mat2df(deconv_mat, value_name="decon", varnames=c("vnum", "volume"))
      deconv_melt$time <- (deconv_melt$volume - 1)*TR + time_offset #convert back to seconds; first volume is time 0

      deconv_df <- deconv_melt %>% 
        dplyr::mutate(vnum=as.numeric(.data$vnum)) %>%
        left_join(a_coordinates, by="vnum") %>%
        #mutate(nifti=niftis[si]) %>%
        dplyr::select(-"i", -"j", -"k") #omitting i, j, k for now

      deconv_df$decon_status <- decon_status[deconv_df$vnum]

      #add subject metadata, if relevant
      #if (!is.null(add_metadata)) { deconv_df <- deconv_df %>% cbind(this_subj) }

      readr::write_csv(deconv_df, file = out_name)

      # A per-file sidecar makes custom filenames independent of identity
      # parsing. Preserve identifier strings such as "001" and "participantA".
      identity_metadata <- data.frame(.nifti = niftis[si])
      if (!is.null(add_metadata)) {
        for (column in intersect(c("id", "run", "session"), names(this_subj))) {
          identity_metadata[[column]] <- this_subj[[column]]
        }
      }
      readr::write_csv(identity_metadata, file = paste0(out_name, ".metadata.csv"))

      #handle output of original time series (before deconvolution) -- mostly useful for debugging
      if (isTRUE(save_original_ts)) {
        to_deconvolve_melt <- mat2df(to_deconvolve, value_name="BOLD_z", varnames=c("vnum", "time"))
        to_deconvolve_melt$time <- (to_deconvolve_melt$time - 1)*TR + time_offset #convert back to seconds; first volume is time 0

        orig_df <- to_deconvolve_melt %>% 
          dplyr::mutate(vnum=as.numeric(.data$vnum)) %>%
          left_join(a_coordinates, by="vnum") %>%
          #mutate(nifti=niftis[si]) %>%
          dplyr::select(-"i", -"j", -"k") #omitting i, j, k for now

        #if (!is.null(add_metadata)) { orig_df <- orig_df %>% cbind(this_subj) }

        #update output file for original
        out_name <- file.path(out_dir, atlas_img_name, "original", paste0(eval(out_file_expression), "_original.csv.gz"))
        readr::write_csv(orig_df, file = out_name)
      }
    }
    
    #write metadata as single data.frame that can be merged selectively, cutting down on storage demands in individual files
    if (!is.null(add_metadata)) {
      out_name <- file.path(out_dir, atlas_img_name, paste0(atlas_img_name, "_metadata.csv"))
      if (file.exists(out_name)) {
        readr::write_csv(add_metadata, file = out_name, append = TRUE, col_names = FALSE)
      } else {
        readr::write_csv(add_metadata, file = out_name, append = FALSE, col_names = TRUE)
      }
    }

  }
  
  return(invisible(NULL))
  
}
