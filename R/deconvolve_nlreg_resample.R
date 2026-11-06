#' This function deconvolves the BOLD signal using Bush 2011 method, augmented
#' by the resampling approach of Bush 2015.
#
#' @param bold_obs observed BOLD signal
#' @param kernel assumed kernel of the BOLD signal
#' @param nev_lr learning rate for the assignment of neural events. Default: .01
#' @param epsilon relative error change (termination condition). Default: .005
#' @param beta slope of the sigmoid transfer function. Default: 40
#' @param n_resample number of resampling steps for deconvolution (at least two). Default: 30
#' @param trim_kernel whether to remove the first K time points from the deconvolved vector, corresponding to
#'            kernel leftovers from convolution. Default: TRUE
#'
#' @param normalize whether to z-score the input before both the base fit and
#'   residual resampling. Default: TRUE. FALSE is for inputs already standardized
#'   by the caller; Bush resampling is not calibrated for arbitrary input units.
#'
#' @return
#'   list containing the following fields:
#'     - NEVest  - the base neural event estimate
#'     - NEVmean - the mean neural event estimate
#'     - NEVstd - the std dev. of the neural event estimate
#'     - NEVcupp - the mean (upper limit 95% ci)
#'     - NEVclow - the mean (lower limit 95% ci)
#'     - BLDmean - the mean of the BOLD estimate
#'     - BLDstd - the std dev. of the BOLD estimate
#'     - BLDcupp - the mean BOLD (upper limit 95% ci)
#'     - BLDclow - the mean BOLD (lower limit 95% ci)
#'
#' @details
#' Author: Keith Bush
#' Institution: University of Arkansas at Little Rock
#' Date:        Aug. 12, 2013
#'
#' @author Keith Bush
#' @importFrom stats t.test
#' @importFrom checkmate assert_logical
#' @export
deconvolve_nlreg_resample <- function(bold_obs, kernel, nev_lr=.01, epsilon=.005, beta=40, n_resample=30, trim_kernel=TRUE, normalize=TRUE) {

  checkmate::assert_logical(trim_kernel, len=1L)
  checkmate::assert_flag(normalize)
  checkmate::assert_integerish(n_resample, lower=2L, len=1L)

  #length of HRF
  Kobs <- length(kernel)

  #Scale observe via z-scoring
  bold_obs_scale <- if (normalize) as.vector(scale(bold_obs)) else as.vector(bold_obs)

  #Deconvolve observed BOLD
  nev_dcv <- deconvolve_nlreg(matrix(bold_obs_scale, ncol=1), kernel, nev_lr, epsilon, beta, normalize = FALSE, trim_kernel = FALSE)

  #Reconvolve estimated true BOLD
  bold_dcv_full <- convolve_dcv_hrf(t(nev_dcv), kernel)

  #remove kernel on front end
  bold_dcv_low <- bold_dcv_full[Kobs:length(bold_dcv_full)]

  #Convert to a z-score representation
  bold_dcv_scale <- as.vector(scale(bold_dcv_low))

  #==============================
  # Apply the resampling method
  #==============================
  NEVvariants <- matrix(0, nrow=n_resample, ncol=length(nev_dcv))
  bold_variants <- matrix(0, nrow=n_resample, ncol=length(bold_obs_scale))

  NEVvariants[1, ] <- nev_dcv
  bold_variants[1, ] <- bold_dcv_scale

  # This residual is fixed across surrogate draws.
  residuals <- bold_obs_scale - bold_dcv_scale
  for (z in 2:n_resample) {

    #Randomize the residual order
    RND_residuals <- sample(residuals)

    #Reapply the residual to the filtered bold_obs
    RNDobs <- bold_dcv_scale + RND_residuals

    #Deconvolve the variant
    NEVresidual <- deconvolve_nlreg(BOLDobs = matrix(RNDobs, ncol=1), kernel, nev_lr, epsilon, beta, normalize = FALSE, trim_kernel = FALSE)

    #Store encoding of this variant
    NEVvariants[z, ] <- NEVresidual

    #Reconvolve to form the filtered BOLD of this variant
    bold_dcv_full <- convolve_dcv_hrf(t(NEVresidual), kernel) #need transpose?
    bold_dcv_low <- bold_dcv_full[Kobs:length(bold_dcv_full)]
    BOLDdcv_res <- as.vector(scale(bold_dcv_low))
    bold_variants[z, ] <- BOLDdcv_res

  }

  #========================================================
  # Compute Performance for regular versus precision-guided
  #========================================================

  result <- list()

  ##Compute distribution over NEVvariants
  result[["NEVest"]] <- nev_dcv
  result[["NEVmean"]] <- apply(NEVvariants, 2, mean)
  result[["NEVstd"]] <- apply(NEVvariants, 2, sd)

  result[["NEVcupp"]] <- 0 * result[["NEVest"]]
  result[["NEVclow"]] <- 0 * result[["NEVest"]]

  for (i in seq_along(result[["NEVcupp"]])) {
    t_res <- t.test(NEVvariants[, i])
    result[["NEVclow"]][i] <- t_res$conf.int[1]
    result[["NEVcupp"]][i] <- t_res$conf.int[2]
  }

  #remove the initial timepoints from the deconvolved time series?
  if (isTRUE(trim_kernel)) {
    result[c("NEVest", "NEVmean", "NEVstd", "NEVcupp", "NEVclow")] <-
      lapply(result[c("NEVest", "NEVmean", "NEVstd", "NEVcupp", "NEVclow")], function(x) {
        x[Kobs:length(x)] #trim off the kernel elements
      })
  }

  ##Compute distribution over bold_variants
  result[["BOLDmean"]] <- apply(bold_variants, 2, mean)
  result[["BOLDstd"]] <- apply(bold_variants, 2, sd)
  result[["BOLDcupp"]] <- 0 * result[["BOLDmean"]] #preallocate
  result[["BOLDclow"]] <- 0 * result[["BOLDmean"]]

  for (i in seq_along(result[["BOLDcupp"]])) {
    t_res <- t.test(bold_variants[, i])
    result[["BOLDclow"]][i] <- t_res$conf.int[1]
    result[["BOLDcupp"]][i] <- t_res$conf.int[2]
  }

  return(result)
}

#' This function takes in a matrix of generated neural events and
#' parameters describing the HRF function and convolves the neural
#' events and HRF function neural events
#'
#' @param NEVgen ROIs x time matrix of true neural events generated by the model
#' @param kernel The HRF kernel used for convolution
#'
#' @author Keith Bush
#' @keywords internal
convolve_dcv_hrf <- function(NEVgen, kernel) {
  stopifnot(is.matrix(NEVgen))

  # Compute number of ROIS
  N <- nrow(NEVgen)
  Bn <- ncol(NEVgen)

  # Calc simulation steps related to simulation time
  Bk <- length(kernel)

  # Allocate Memory to store model brfs
  bold_gen <- matrix(0, nrow = N, ncol = Bn + Bk - 1) # zeros(N,Bn+Bk-1);

  # Convert neural events to indices
  for (curr_node in 1:N) {

    # Superimpose all kernels into one time-series
    for (i in seq_along(NEVgen[curr_node, ])) {
      bold_gen[curr_node, i:(i + Bk - 1)] <- NEVgen[curr_node, i] * t(kernel) + bold_gen[curr_node, i:(i + Bk - 1)]
    }
  }

  # Trim the excess Bk-1 time-points from result
  bold_gen <- bold_gen[, 1:Bn]

  return(bold_gen)
}