#' Regularized linear deconvolution of BOLD time series
#'
#' Estimate a continuous latent activity time series by solving a regularized
#' inverse convolution problem for one or more BOLD signals.
#'
#' @param BOLDobs matrix of observed BOLD time series (n_timepoints x n_signals)
#'   or a numeric vector.
#' @param kernel assumed HRF kernel. A vector or single-column matrix.
#' @param lambda regularization strength. If \code{NULL}, chosen from
#'   \code{lambda_grid} by generalized cross-validation.
#' @param lambda_grid candidate regularization strengths used when
#'   \code{lambda} is \code{NULL}.
#' @param penalty regularization penalty: \code{"diff2"} for smooth second
#'   differences, \code{"diff1"} for smooth first differences, \code{"ridge"}
#'   for amplitude shrinkage, \code{"matern1"} for a damped first-order
#'   Matérn/Ornstein-Uhlenbeck penalty (\code{0.10 * I + 0.90 * D1'D1}), or
#'   \code{"sobolev"} for a critically damped second-order Sobolev penalty
#'   (\code{0.05 * I + 0.45 * D1'D1 + 0.50 * D2'D2}) that suppresses negative
#'   sidelobe ringing and stabilizes low frequencies.
#' @param kernel_grid optional list or matrix of candidate HRF kernels. If
#'   supplied, \code{kernel} and \code{hrf_lags} are ignored for tuning.
#' @param hrf_lags integer shifts, in TR bins, applied to \code{kernel}. Negative
#'   values make the candidate HRF peak earlier; positive values make it later.
#' @param tune_by whether to tune \code{lambda}/HRF globally across all signals
#'   or independently for each signal.
#' @param normalize whether to z-score \code{BOLDobs} before deconvolution.
#'   Default: \code{TRUE}.
#' @param demean whether to subtract column means before deconvolution. Ignored
#'   when \code{normalize = TRUE}, because z-scoring already demeans.
#' @param trim_kernel whether to remove the initial \code{length(kernel) - 1}
#'   latent samples, corresponding to pre-observation HRF support.
#' @param ridge_floor small diagonal ridge added for numerical stability.
#' @param return_diagnostics whether to return a list containing diagnostics
#'   rather than only the deconvolved matrix.
#' @param prewhiten_gcv whether to estimate AR(1) autocorrelation and use a
#'   prewhitened objective for lambda/HRF selection and deconvolution.
#' @param prewhiten_rho optional AR(1) coefficient used when
#'   \code{prewhiten_gcv = TRUE}. If \code{NULL}, a pooled median AR(1)
#'   coefficient is estimated from first-pass reconstruction residuals.
#' @param backend implementation backend. \code{"cpp"} uses the RcppArmadillo
#'   regularized inverse implementation; \code{"r"} uses the reference R
#'   implementation; \code{"wiener"} uses a frequency-domain Wiener-style
#'   deconvolution for comparison.
#' @param wiener_shrinkage shrinkage strategy for \code{backend = "wiener"}.
#'   \code{"scalar"} uses \code{Conj(H) / (|H|^2 + lambda)}. \code{"psd"}
#'   estimates smoothed signal and residual noise spectra from a scalar Wiener
#'   first pass, then uses frequency-dependent shrinkage.
#' @param wiener_psd_smooth odd integer smoothing width for the estimated
#'   signal/noise spectra in \code{wiener_shrinkage = "psd"}.
#' @param wiener_psd_floor lower bound for estimated spectra.
#' @param gcv_rule lambda/HRF selection rule. \code{"min"} selects the minimum
#'   mean GCV candidate. \code{"1se"} selects the largest lambda whose mean GCV
#'   is within one standard error of the minimum. The standard error is computed
#'   across signal columns for \code{tune_by = "global"}.
#'
#' @param scaling Explicit input scaling: \code{"zscore"}, \code{"noise_sd"},
#'   \code{"demean"}, \code{"none"}, or \code{"psc"} (percent signal change).
#'   \code{"noise_sd"} centers each signal and divides by a high-frequency
#'   second-difference estimate of observation noise SD (adjusted for AR(1)
#'   autocorrelation when available), standardizing SNR without conflating
#'   neural task variance with noise variance. When supplied,
#'   overrides \code{normalize} and \code{demean}. Otherwise those flags retain
#'   their existing behavior; the default is z-scoring.
#'
#' @return A matrix of deconvolved time series, or a list with \code{activity}
#'   and tuning diagnostics when \code{return_diagnostics = TRUE}.
#' @details Constant signals and signals containing nonfinite observations are
#'   excluded from fitting and tuning. Their activity columns are returned as
#'   \code{NA}; diagnostics include a per-column \code{signal_status} vector
#'   (\code{"ok"}, \code{"constant"}, or \code{"nonfinite"}). Partial missingness
#'   excludes the entire signal; observations are not interpolated or deleted.
#'   Z-scoring remains the default and is recommended for standardized within-run
#'   trial responses. For fixed HRF/regularization, the linear inverse is homogeneous
#'   in its input: rescaling the observations rescales activity without changing the
#'   conditioning of the inverse matrix. Without z-scoring, global GCV gives signals
#'   with larger variance more weight. Demeaning retains input units; \code{"none"}
#'   retains offsets too and is mainly useful for already-centered data or simulations.
#'   Percent signal change uses \code{100 * (y - mean(y)) / mean(y)} and requires a
#'   positive baseline well away from zero. It is unsuitable for zero-mean residuals.
#'   No scaling choice removes HRF uncertainty or regularization-induced attenuation.
#'
#'   Scaling is applied once. The matrix attributes \code{scaled:center},
#'   \code{scaled:scale}, and \code{scaling} record the actual BOLD transformation
#'   \code{(y - center) / scale}; diagnostics expose the same values. With
#'   \code{scaling = "none"}, the recorded center is zero. The activity scale depends
#'   on this transformation and the unit-sum HRF; it is not a calibrated physiological
#'   neural-amplitude scale. The removed BOLD mean is not an identified neural baseline.
#'
#'   Scalar Wiener GCV uses the exact trace of the finite observed smoother, including
#'   zero extension, latent truncation, and the implemented filtered HRF when whitening.
#'   This trace is conditional on preprocessing and the AR(1) coefficient; it does not
#'   include derivatives of estimated preprocessing/noise parameters. Adaptive PSD
#'   shrinkage requires a fixed lambda and a single HRF. Its diagnostics label the
#'   filter trace \code{conditional_fixed_psd} and return missing GCV values, because
#'   the trace holding spectra fixed is not the degrees of freedom of the full
#'   nonlinear, data-adaptive estimator. Scalar diagnostics use
#'   \code{df_method = "finite_observation_trace"}.
#' @importFrom stats fft median mvfft sd
#' @export
deconvolve_reglin <- function(
  BOLDobs,
  kernel,
  lambda = NULL,
  lambda_grid = 10^seq(-4, 2, length.out = 20),
  penalty = c("diff2", "diff1", "ridge", "matern1", "sobolev"),
  kernel_grid = NULL,
  hrf_lags = 0L,
  tune_by = c("global", "signal"),
  normalize = TRUE,
  demean = TRUE,
  trim_kernel = TRUE,
  ridge_floor = 1e-8,
  return_diagnostics = FALSE,
  prewhiten_gcv = FALSE,
  prewhiten_rho = NULL,
  backend = c("cpp", "r", "wiener"),
  wiener_shrinkage = c("scalar", "psd"),
  wiener_psd_smooth = 7L,
  wiener_psd_floor = 1e-8,
  gcv_rule = c("min", "1se"),
  scaling = NULL
) {
  penalty <- match.arg(penalty)
  tune_by <- match.arg(tune_by)
  backend <- match.arg(backend)
  wiener_shrinkage <- match.arg(wiener_shrinkage)
  gcv_rule <- match.arg(gcv_rule)

  if (is.null(dim(BOLDobs))) {
    BOLDobs <- matrix(BOLDobs, ncol = 1L)
  } else {
    BOLDobs <- as.matrix(BOLDobs)
  }

  checkmate::assert_matrix(BOLDobs, mode = "numeric", min.rows = 3L, min.cols = 1L)
  checkmate::assert_numeric(kernel, any.missing = FALSE, min.len = 2L)
  checkmate::assert_number(ridge_floor, lower = 0)
  checkmate::assert_logical(normalize, len = 1L)
  checkmate::assert_logical(demean, len = 1L)
  checkmate::assert_logical(trim_kernel, len = 1L)
  checkmate::assert_logical(return_diagnostics, len = 1L)
  checkmate::assert_logical(prewhiten_gcv, len = 1L)
  checkmate::assert_number(prewhiten_rho, lower = -0.99, upper = 0.99, null.ok = TRUE)
  checkmate::assert_integerish(wiener_psd_smooth, lower = 1L, len = 1L)
  checkmate::assert_number(wiener_psd_floor, lower = .Machine$double.eps)
  wiener_psd_smooth <- as.integer(wiener_psd_smooth)
  scaling <- resolve_decon_scaling(scaling, normalize, demean)

  if (is.null(lambda)) {
    checkmate::assert_numeric(lambda_grid, lower = 0, any.missing = FALSE, min.len = 1L)
    lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  } else {
    checkmate::assert_number(lambda, lower = 0)
    lambda_grid <- lambda
  }

  kernels <- build_decon_kernel_grid(kernel, kernel_grid = kernel_grid, hrf_lags = hrf_lags)
  kernel_len <- length(kernels[[1L]])
  validate_wiener_tuning(backend, wiener_shrinkage, length(lambda_grid), length(kernels))
  if (nrow(BOLDobs) < 3L) {
    stop("Need at least 3 timepoints for regularized linear deconvolution")
  }

  signal_status <- decon_signal_status(BOLDobs)
  valid <- signal_status == "ok"
  if (!all(valid)) {
    warning(sum(!valid), " constant or nonfinite signals excluded from deconvolution", call. = FALSE)
    n_output <- nrow(BOLDobs) + if (trim_kernel) 0L else kernel_len - 1L
    activity <- matrix(NA_real_, n_output, ncol(BOLDobs))
    centers <- scales <- rep(NA_real_, ncol(BOLDobs))

    if (any(valid)) {
      # Evaluate the public arguments before recursing so all backends and the
      # first-pass noise estimate see only valid signals with identical settings.
      fit_args <- mget(names(formals(deconvolve_reglin)), envir = environment())
      fit_args$BOLDobs <- BOLDobs[, valid, drop = FALSE]
      fit_args$return_diagnostics <- TRUE
      result <- do.call(deconvolve_reglin, fit_args)
      activity[, valid] <- result$activity
      centers[valid] <- as.numeric(result$scaled_center)
      scales[valid] <- as.numeric(result$scaled_scale)
      if (tune_by == "signal") {
        lambdas <- rep(NA_real_, ncol(BOLDobs))
        kernel_indices <- rep(NA_integer_, ncol(BOLDobs))
        lambdas[valid] <- result$lambda
        kernel_indices[valid] <- result$kernel_index
        result$lambda <- lambdas
        result$kernel_index <- kernel_indices
      }
    } else {
      n_selected <- if (tune_by == "global") 1L else ncol(BOLDobs)
      result <- list(
        lambda = rep(NA_real_, n_selected),
        kernel_index = rep(NA_integer_, n_selected),
        kernels = kernels,
        tuning = data.frame(kernel = integer(), lambda = numeric(), df = numeric(),
          rss = numeric(), gcv = numeric(), gcv_se = numeric()),
        penalty = if (backend == "wiener") "wiener" else penalty,
        tune_by = tune_by, gcv_rule = gcv_rule, prewhiten_gcv = prewhiten_gcv,
        ar1_rho = if (prewhiten_gcv) {
          if (is.null(prewhiten_rho)) NA_real_ else prewhiten_rho
        } else 0
      )
      if (backend == "wiener") {
        result$wiener_shrinkage <- wiener_shrinkage
        result$df_method <- if (wiener_shrinkage == "psd") "conditional_fixed_psd" else "finite_observation_trace"
      }
    }

    attr(activity, "scaled:center") <- centers
    attr(activity, "scaled:scale") <- scales
    result$scaled_center <- centers
    result$scaled_scale <- scales
    result$signal_status <- signal_status
    result$scaling <- scaling
    attr(activity, "scaling") <- scaling
    result$activity <- activity
    return(if (return_diagnostics) result else activity)
  }

  # Keep native scaling for the usual C++ modes, avoiding an R column loop
  # and extra matrix copies. R scaling is needed for PSC, noise_sd, or noise estimation.
  scaled <- if (backend != "cpp" || scaling %in% c("psc", "noise_sd") ||
      (prewhiten_gcv && is.null(prewhiten_rho))) {
    scale_decon_input(BOLDobs, scaling = scaling, noise_rho = prewhiten_rho)
  } else NULL

  # Two-pass AR(1) prewhitened GCV:
  # Estimating AR(1) autocorrelation directly from raw BOLD would conflate slow
  # hemodynamic task signal with colored noise. Instead, we run an unwhitened
  # first pass at max(lambda_grid) (strong smoothing), subtract the reconvolved
  # task fit, and estimate a pooled median AR(1) rho from the reconstruction
  # residuals before running prewhitened GCV over the full lambda_grid.
  if (isTRUE(prewhiten_gcv) && is.null(prewhiten_rho)) {
    first_pass <- deconvolve_reglin(
      BOLDobs = BOLDobs,
      kernel = kernel,
      lambda = max(lambda_grid),
      lambda_grid = max(lambda_grid),
      penalty = penalty,
      kernel_grid = kernel_grid,
      hrf_lags = hrf_lags,
      tune_by = tune_by,
      normalize = normalize,
      demean = demean,
      trim_kernel = FALSE,
      ridge_floor = ridge_floor,
      return_diagnostics = TRUE,
      prewhiten_gcv = FALSE,
      backend = backend,
      wiener_shrinkage = wiener_shrinkage,
      wiener_psd_smooth = wiener_psd_smooth,
      wiener_psd_floor = wiener_psd_floor,
      gcv_rule = gcv_rule,
      scaling = scaling
    )
    first_resid <- decon_first_pass_residuals(
      y = scaled$y,
      activity = first_pass$activity,
      kernels = kernels,
      kernel_index = first_pass$kernel_index,
      tune_by = tune_by
    )
    prewhiten_rho <- resolve_decon_ar1(first_resid, prewhiten_gcv = TRUE, prewhiten_rho = NULL)
  }

  if (backend == "wiener") {
    result <- deconvolve_reglin_wiener(
      BOLDobs = scaled$y,
      kernels = kernels,
      lambda_grid = lambda_grid,
      tune_by = tune_by,
      normalize = FALSE,
      demean = FALSE,
      trim_kernel = trim_kernel,
      return_diagnostics = return_diagnostics,
      prewhiten_gcv = prewhiten_gcv,
      prewhiten_rho = prewhiten_rho,
      wiener_shrinkage = wiener_shrinkage,
      wiener_psd_smooth = wiener_psd_smooth,
      wiener_psd_floor = wiener_psd_floor,
      gcv_rule = gcv_rule
    )
    return(set_decon_scaling(result, scaled, scaling))
  }

  if (backend == "cpp") {
    kernel_mat <- do.call(cbind, kernels)
    cpp_res <- deconvolve_reglin_cpp(
      BOLDobs = if (is.null(scaled)) BOLDobs else scaled$y,
      kernels = kernel_mat,
      lambda_grid = lambda_grid,
      penalty = penalty,
      tune_by = tune_by,
      normalize = is.null(scaled) && scaling == "zscore",
      demean = is.null(scaled) && scaling == "demean",
      trim_kernel = trim_kernel,
      ridge_floor = ridge_floor,
      return_diagnostics = return_diagnostics,
      prewhiten_gcv = prewhiten_gcv,
      prewhiten_rho = if (is.null(prewhiten_rho)) NA_real_ else prewhiten_rho,
      gcv_rule = gcv_rule
    )

    if (isTRUE(return_diagnostics)) {
      cpp_res$kernels <- kernels
      cpp_res$signal_status <- signal_status
    }
    if (is.null(scaled)) {
      activity <- if (is.list(cpp_res)) cpp_res$activity else cpp_res
      scaled <- list(center = as.numeric(attr(activity, "scaled:center")),
        scale = as.numeric(attr(activity, "scaled:scale")))
    }
    return(set_decon_scaling(cpp_res, scaled, scaling))
  }

  y <- scaled$y
  ar1_rho <- resolve_decon_ar1(y, prewhiten_gcv = prewhiten_gcv, prewhiten_rho = prewhiten_rho)
  y_fit <- if (isTRUE(prewhiten_gcv)) apply_ar1_whitening(y, ar1_rho) else y
  n_time <- nrow(y)
  n_signals <- ncol(y)
  n_latent <- n_time + kernel_len - 1L
  penalty_crossprod <- decon_penalty_crossprod(n_latent, penalty)

  best_score <- if (tune_by == "global") Inf else rep(Inf, n_signals)
  best_activity <- matrix(NA_real_, nrow = n_latent, ncol = n_signals)
  best_lambda <- if (tune_by == "global") NA_real_ else rep(NA_real_, n_signals)
  best_kernel <- if (tune_by == "global") NA_integer_ else rep(NA_integer_, n_signals)
  tuning <- vector("list", length(kernels) * length(lambda_grid))
  tuning_i <- 0L

  for (ki in seq_along(kernels)) {
    hrf <- kernels[[ki]]
    conv_mat <- decon_convolution_matrix(n_time, hrf)
    fit_mat <- if (isTRUE(prewhiten_gcv)) apply_ar1_whitening(conv_mat, ar1_rho) else conv_mat
    hth <- crossprod(fit_mat)
    hty <- crossprod(fit_mat, y_fit)

    for (lam in lambda_grid) {
      lhs <- hth + lam * penalty_crossprod
      diag(lhs) <- diag(lhs) + ridge_floor
      activity <- solve(lhs, hty)
      fitted <- fit_mat %*% activity
      residual <- y_fit - fitted
      rss_by_signal <- colSums(residual^2)
      df <- sum(diag(solve(lhs, hth)))
      gcv_denom <- max((n_time - df)^2, .Machine$double.eps)
      gcv_by_signal <- rss_by_signal / gcv_denom
      gcv_se <- decon_gcv_se(gcv_by_signal)

      tuning_i <- tuning_i + 1L
      tuning[[tuning_i]] <- data.frame(
        kernel = ki,
        lambda = lam,
        df = df,
        rss = sum(rss_by_signal),
        gcv = mean(gcv_by_signal),
        gcv_se = gcv_se,
        stringsAsFactors = FALSE
      )

      if (tune_by == "global") {
        score <- mean(gcv_by_signal)
        if (score < best_score) {
          best_score <- score
          best_activity <- activity
          best_lambda <- lam
          best_kernel <- ki
        }
      } else {
        improved <- gcv_by_signal < best_score
        if (any(improved)) {
          best_score[improved] <- gcv_by_signal[improved]
          best_activity[, improved] <- activity[, improved, drop = FALSE]
          best_lambda[improved] <- lam
          best_kernel[improved] <- ki
        }
      }
    }
  }

  tuning <- do.call(rbind, tuning)
  if (gcv_rule == "1se" && tune_by == "global") {
    selected <- select_1se_candidate(tuning)
    hrf <- kernels[[selected$kernel]]
    conv_mat <- decon_convolution_matrix(n_time, hrf)
    fit_mat <- if (isTRUE(prewhiten_gcv)) apply_ar1_whitening(conv_mat, ar1_rho) else conv_mat
    hth <- crossprod(fit_mat)
    hty <- crossprod(fit_mat, y_fit)
    lhs <- hth + selected$lambda * penalty_crossprod
    diag(lhs) <- diag(lhs) + ridge_floor
    best_activity <- solve(lhs, hty)
    best_lambda <- selected$lambda
    best_kernel <- selected$kernel
  }

  if (isTRUE(trim_kernel)) {
    best_activity <- best_activity[kernel_len:nrow(best_activity), , drop = FALSE]
  }

  attr(best_activity, "scaled:center") <- scaled$center
  attr(best_activity, "scaled:scale") <- scaled$scale
  attr(best_activity, "scaling") <- scaling

  if (!isTRUE(return_diagnostics)) {
    return(best_activity)
  }

  list(
    activity = best_activity,
    lambda = best_lambda,
    kernel_index = best_kernel,
    kernels = kernels,
    tuning = tuning,
    penalty = penalty,
    tune_by = tune_by,
    gcv_rule = gcv_rule,
    prewhiten_gcv = prewhiten_gcv,
    ar1_rho = ar1_rho,
    scaled_center = scaled$center,
    scaled_scale = scaled$scale,
    signal_status = signal_status,
    scaling = scaling
  )
}

# Shared screening for standalone and voxelwise deconvolution. A nonfinite
# observation excludes its entire signal, preserving the acquisition clock.
decon_signal_status <- function(y) {
  vapply(seq_len(ncol(y)), function(i) {
    x <- y[, i]
    if (any(!is.finite(x))) return("nonfinite")
    spread <- stats::sd(x)
    if (!is.finite(spread)) return("nonfinite")
    if (spread <= .Machine$double.eps) "constant" else "ok"
  }, character(1L))
}

build_decon_kernel_grid <- function(kernel, kernel_grid = NULL, hrf_lags = 0L) {
  if (!is.null(kernel_grid)) {
    if (is.list(kernel_grid)) {
      kernels <- lapply(kernel_grid, as.numeric)
    } else {
      kernel_grid <- as.matrix(kernel_grid)
      kernels <- lapply(seq_len(ncol(kernel_grid)), function(i) as.numeric(kernel_grid[, i]))
    }
  } else {
    checkmate::assert_integerish(hrf_lags, any.missing = FALSE)
    base_kernel <- as.numeric(kernel)
    kernels <- lapply(as.integer(hrf_lags), shift_decon_kernel, kernel = base_kernel)
  }

  lens <- vapply(kernels, length, integer(1L))
  if (length(unique(lens)) != 1L) {
    stop("All candidate HRF kernels must have the same length")
  }

  lapply(kernels, normalize_decon_kernel)
}

normalize_decon_kernel <- function(kernel) {
  kernel <- as.numeric(kernel)
  if (any(!is.finite(kernel))) {
    stop("HRF kernels must contain only finite values")
  }
  ksum <- sum(kernel)
  if (isTRUE(all.equal(ksum, 0))) {
    stop("HRF kernel sum is zero; cannot normalize")
  }
  kernel / ksum
}

shift_decon_kernel <- function(lag, kernel) {
  lag <- as.integer(lag)
  n <- length(kernel)

  if (lag == 0L) {
    return(kernel)
  }

  if (abs(lag) >= n) {
    stop("HRF lag magnitude must be smaller than the kernel length")
  }

  if (lag > 0L) {
    c(rep(0, lag), kernel[seq_len(n - lag)])
  } else {
    lag <- abs(lag)
    c(kernel[(lag + 1L):n], rep(0, lag))
  }
}

validate_wiener_tuning <- function(backend, shrinkage, n_lambda, n_kernels) {
  if (backend == "wiener" && shrinkage == "psd" && (n_lambda > 1L || n_kernels > 1L)) {
    stop("Automatic GCV tuning is not supported for data-adaptive PSD Wiener filters; ",
      "supply a fixed lambda and one kernel, or use wiener_shrinkage = 'scalar'.")
  }
}

resolve_decon_scaling <- function(scaling = NULL, normalize = TRUE, demean = TRUE) {
  checkmate::assert_flag(normalize)
  checkmate::assert_flag(demean)
  if (is.null(scaling)) {
    return(if (isTRUE(normalize)) "zscore" else if (isTRUE(demean)) "demean" else "none")
  }
  match.arg(scaling, c("zscore", "demean", "none", "psc", "noise_sd"))
}

# Estimate marginal observation noise SD from second differences, where smooth
# HRF-convolved task signal has negligible energy (< 0.5% of ||h||_2^2).
# Under AR(1) noise with lag-1 correlation rho, Var(Delta^2 e_t) = 2*(1-rho)*(3-rho)*sigma_e^2.
estimate_decon_noise_sd <- function(y, noise_rho = NULL) {
  y <- as.matrix(y)
  n <- nrow(y)
  p <- ncol(y)
  if (n < 4L) {
    divisor <- apply(y, 2L, stats::sd)
    divisor[!is.finite(divisor) | divisor <= .Machine$double.eps] <- 1
    return(as.numeric(divisor))
  }
  d1 <- y[2:n, , drop = FALSE] - y[1:(n - 1L), , drop = FALSE]
  d2 <- y[3:n, , drop = FALSE] - 2 * y[2:(n - 1L), , drop = FALSE] + y[1:(n - 2L), , drop = FALSE]
  rob_scale <- function(mat) {
    vapply(seq_len(ncol(mat)), function(j) {
      col <- mat[, j]
      col <- col[is.finite(col)]
      if (length(col) < 2L) return(NA_real_)
      s_sd <- stats::sd(col)
      s_mad <- stats::mad(col, center = 0)
      if (is.finite(s_mad) && s_mad > .Machine$double.eps) 0.5 * s_sd + 0.5 * s_mad else s_sd
    }, numeric(1L))
  }
  s2 <- rob_scale(d2)
  if (!is.null(noise_rho)) {
    rho_vec <- rep_len(as.numeric(noise_rho), p)
  } else {
    rho_vec <- rep(NA_real_, p)
  }
  need_rho <- !is.finite(rho_vec)
  if (any(need_rho)) {
    # Estimate AR(1) rho purely from the lag-1 autocorrelation r1 of Delta^2 y,
    # where smooth HRF task signal is annihilated: r1 = -(4 - 3*rho + rho^2) / (2*(3 - rho)).
    rho_est <- vapply(which(need_rho), function(j) {
      col <- d2[, j]
      if (length(col) < 4L || any(!is.finite(col)) || stats::sd(col) <= .Machine$double.eps) return(0)
      r1 <- stats::cor(col[-length(col)], col[-1L])
      if (!is.finite(r1)) return(0)
      r1 <- min(max(r1, -2 / 3), -0.46)
      b <- 3 + 2 * r1
      c_val <- 4 + 6 * r1
      disc <- max(0, b * b - 4 * c_val)
      0.5 * (b - sqrt(disc))
    }, numeric(1L))
    rho_vec[need_rho] <- pmin(pmax(rho_est, 0), 0.90)
  }
  rho_vec <- pmin(pmax(rho_vec, -0.90), 0.95)
  denom <- sqrt(pmax(2 * (1 - rho_vec) * (3 - rho_vec), .Machine$double.eps))
  divisor <- s2 / denom
  divisor[!is.finite(divisor) | divisor <= .Machine$double.eps] <- 1
  as.numeric(divisor)
}

# Record the actual transform y_scaled = (y - center) / scale. Scaling is
# independent of the inverse operator; it does not remove regularization bias.
scale_decon_input <- function(y, normalize = TRUE, demean = TRUE, scaling = NULL, noise_rho = NULL) {
  scaling <- resolve_decon_scaling(scaling, normalize, demean)
  center <- if (scaling == "none") rep(0, ncol(y)) else colMeans(y)
  divisor <- rep(1, ncol(y))
  if (scaling == "zscore") {
    divisor <- apply(y, 2L, stats::sd)
    divisor[!is.finite(divisor) | divisor <= .Machine$double.eps] <- 1
  } else if (scaling == "noise_sd") {
    divisor <- estimate_decon_noise_sd(y, noise_rho = noise_rho)
  } else if (scaling == "psc") {
    spread <- apply(y, 2L, stats::sd)
    if (any(!is.finite(center) | center <= sqrt(.Machine$double.eps) * spread)) {
      stop("Percent-signal-change scaling requires a positive baseline mean well away from zero; ",
        "use zscore or demean for residualized/zero-mean signals.")
    }
    divisor <- center / 100
  }
  if (scaling != "none") y <- sweep(y, 2L, center, "-")
  if (scaling %in% c("zscore", "psc", "noise_sd")) y <- sweep(y, 2L, divisor, "/")
  list(y = y, center = center, scale = divisor)
}

set_decon_scaling <- function(result, scaled, scaling) {
  activity <- if (is.list(result)) result$activity else result
  attr(activity, "scaled:center") <- scaled$center
  attr(activity, "scaled:scale") <- scaled$scale
  attr(activity, "scaling") <- scaling
  if (!is.list(result)) return(activity)
  result$activity <- activity
  result$scaled_center <- scaled$center
  result$scaled_scale <- scaled$scale
  result$scaling <- scaling
  result
}

resolve_decon_ar1 <- function(y, prewhiten_gcv = FALSE, prewhiten_rho = NULL) {
  if (!isTRUE(prewhiten_gcv)) {
    return(0)
  }

  if (!is.null(prewhiten_rho)) {
    return(prewhiten_rho)
  }

  rhos <- apply(y, 2L, function(x) {
    if (length(x) < 3L || stats::sd(x) <= .Machine$double.eps) {
      return(NA_real_)
    }
    stats::cor(x[-length(x)], x[-1L], use = "complete.obs")
  })
  rhos <- rhos[is.finite(rhos)]
  if (length(rhos) == 0L) {
    return(0)
  }
  max(min(stats::median(rhos), 0.95), 0)
}

apply_ar1_whitening <- function(x, rho) {
  x <- as.matrix(x)
  if (!is.finite(rho) || abs(rho) <= .Machine$double.eps) {
    return(x)
  }

  out <- x
  out[1L, ] <- sqrt(1 - rho^2) * x[1L, ]
  out[-1L, ] <- x[-1L, , drop = FALSE] - rho * x[-nrow(x), , drop = FALSE]
  out
}

decon_gcv_se <- function(gcv_by_signal) {
  gcv_by_signal <- as.numeric(gcv_by_signal)
  gcv_by_signal <- gcv_by_signal[is.finite(gcv_by_signal)]
  if (length(gcv_by_signal) <= 1L) {
    return(0)
  }
  stats::sd(gcv_by_signal) / sqrt(length(gcv_by_signal))
}

select_1se_candidate <- function(tuning) {
  min_idx <- which.min(tuning$gcv)
  threshold <- tuning$gcv[min_idx] + tuning$gcv_se[min_idx]
  eligible <- tuning[tuning$gcv <= threshold, , drop = FALSE]
  eligible <- eligible[order(-eligible$lambda, eligible$gcv, eligible$kernel), , drop = FALSE]
  eligible[1L, , drop = FALSE]
}

deconvolve_reglin_wiener <- function(
  BOLDobs,
  kernels,
  lambda_grid,
  tune_by = c("global", "signal"),
  normalize = TRUE,
  demean = TRUE,
  trim_kernel = TRUE,
  return_diagnostics = FALSE,
  prewhiten_gcv = FALSE,
  prewhiten_rho = NULL,
  wiener_shrinkage = c("scalar", "psd"),
  wiener_psd_smooth = 7L,
  wiener_psd_floor = 1e-8,
  gcv_rule = c("min", "1se")
) {
  tune_by <- match.arg(tune_by)
  wiener_shrinkage <- match.arg(wiener_shrinkage)
  gcv_rule <- match.arg(gcv_rule)
  scaled <- scale_decon_input(BOLDobs, normalize = normalize, demean = demean)
  y <- scaled$y
  n_time <- nrow(y)
  n_signals <- ncol(y)
  kernel_len <- length(kernels[[1L]])
  n_latent <- n_time + kernel_len - 1L
  ar1_rho <- resolve_decon_ar1(y, prewhiten_gcv = prewhiten_gcv, prewhiten_rho = prewhiten_rho)
  y_fit <- if (isTRUE(prewhiten_gcv)) apply_ar1_whitening(y, ar1_rho) else y

  best_score <- if (tune_by == "global") Inf else rep(Inf, n_signals)
  best_activity <- matrix(NA_real_, nrow = n_latent, ncol = n_signals)
  best_lambda <- if (tune_by == "global") NA_real_ else rep(NA_real_, n_signals)
  best_kernel <- if (tune_by == "global") NA_integer_ else rep(NA_integer_, n_signals)
  tuning <- vector("list", length(kernels) * length(lambda_grid))
  tuning_i <- 0L
  needs_scores <- return_diagnostics || length(kernels) * length(lambda_grid) > 1L

  for (ki in seq_along(kernels)) {
    fit_kernel <- if (isTRUE(prewhiten_gcv)) ar1_filter_kernel(kernels[[ki]], ar1_rho) else kernels[[ki]]
    plan <- wiener_fft_plan(n_time, n_latent, fit_kernel)
    y_fft <- mvfft(rbind(y_fit, matrix(0, plan$n_fft - n_time, n_signals)))

    for (lam in lambda_grid) {
      if (wiener_shrinkage == "psd") {
        psd <- estimate_wiener_psd(
          y = y_fit,
          kernel = fit_kernel,
          lambda = lam,
          n_latent = n_latent,
          smooth_width = wiener_psd_smooth,
          psd_floor = wiener_psd_floor,
          fft_plan = plan, y_fft = y_fft
        )
        activity <- wiener_deconvolve_matrix(
          y_fit,
          fit_kernel,
          lambda = lam,
          n_latent = n_latent,
          signal_psd = psd$signal_psd,
          noise_psd = psd$noise_psd, fft_plan = plan, y_fft = y_fft
        )
        filter <- wiener_inverse_filter(plan$h_fft, lam, psd$signal_psd, psd$noise_psd)
      } else {
        activity <- wiener_deconvolve_matrix(y_fit, fit_kernel, lambda = lam,
          n_latent = n_latent, fft_plan = plan, y_fft = y_fft)
        filter <- wiener_inverse_filter(plan$h_fft, lam)
      }
      if (!needs_scores) {
        best_activity <- activity
        best_lambda <- lam
        best_kernel <- ki
        if (tune_by == "signal") {
          best_lambda <- rep(lam, n_signals)
          best_kernel <- rep(ki, n_signals)
        }
        next
      }
      df <- wiener_degrees_freedom(filter, fit_kernel, n_time)
      fitted <- wiener_reconvolve_matrix(activity, fit_kernel, n_time = n_time, fft_plan = plan)
      residual <- y_fit - fitted
      rss_by_signal <- colSums(residual^2)
      gcv_denom <- max((n_time - df)^2, .Machine$double.eps)
      gcv_by_signal <- if (wiener_shrinkage == "psd") rep(NA_real_, n_signals) else rss_by_signal / gcv_denom
      gcv_se <- if (wiener_shrinkage == "psd") NA_real_ else decon_gcv_se(gcv_by_signal)

      tuning_i <- tuning_i + 1L
      tuning[[tuning_i]] <- data.frame(
        kernel = ki,
        lambda = lam,
        df = df,
        rss = sum(rss_by_signal),
        gcv = mean(gcv_by_signal),
        gcv_se = gcv_se,
        stringsAsFactors = FALSE
      )

      if (tune_by == "global") {
        score <- mean(gcv_by_signal)
        if (wiener_shrinkage == "psd" || score < best_score) {
          best_score <- score
          best_activity <- activity
          best_lambda <- lam
          best_kernel <- ki
        }
      } else {
        improved <- if (wiener_shrinkage == "psd") rep(TRUE, n_signals) else gcv_by_signal < best_score
        if (any(improved)) {
          best_score[improved] <- gcv_by_signal[improved]
          best_activity[, improved] <- activity[, improved, drop = FALSE]
          best_lambda[improved] <- lam
          best_kernel[improved] <- ki
        }
      }
    }
  }

  tuning <- do.call(rbind, tuning)
  if (needs_scores && wiener_shrinkage == "scalar" && gcv_rule == "1se" && tune_by == "global") {
    selected <- select_1se_candidate(tuning)
    fit_kernel <- if (isTRUE(prewhiten_gcv)) ar1_filter_kernel(kernels[[selected$kernel]], ar1_rho) else kernels[[selected$kernel]]
    plan <- wiener_fft_plan(n_time, n_latent, fit_kernel)
    best_activity <- wiener_deconvolve_matrix(y_fit, fit_kernel, lambda = selected$lambda,
      n_latent = n_latent, fft_plan = plan, y_fft = y_fft)
    best_lambda <- selected$lambda
    best_kernel <- selected$kernel
  }

  best_activity <- align_wiener_activity(
    activity = best_activity,
    n_time = n_time,
    kernel_len = kernel_len,
    trim_kernel = trim_kernel
  )

  attr(best_activity, "scaled:center") <- scaled$center
  attr(best_activity, "scaled:scale") <- scaled$scale

  if (!isTRUE(return_diagnostics)) {
    return(best_activity)
  }

  list(
    activity = best_activity,
    lambda = best_lambda,
    kernel_index = best_kernel,
    kernels = kernels,
    tuning = tuning,
    penalty = "wiener",
    wiener_shrinkage = wiener_shrinkage,
    df_method = if (wiener_shrinkage == "psd") "conditional_fixed_psd" else "finite_observation_trace",
    tune_by = tune_by,
    gcv_rule = gcv_rule,
    prewhiten_gcv = prewhiten_gcv,
    ar1_rho = ar1_rho,
    scaled_center = scaled$center,
    scaled_scale = scaled$scale,
    signal_status = rep("ok", n_signals)
  )
}

ar1_filter_kernel <- function(kernel, rho) {
  if (!is.finite(rho) || abs(rho) <= .Machine$double.eps) {
    return(kernel)
  }
  c(kernel, 0) - rho * c(0, kernel)
}

wiener_fft_length <- function(n_time, n_latent, kernel_len) {
  2^ceiling(log2(n_latent + kernel_len - 1L))
}

wiener_fft_plan <- function(n_time, n_latent, kernel) {
  n_fft <- wiener_fft_length(n_time, n_latent, length(kernel))
  list(n_fft = n_fft, h_fft = fft(c(kernel, rep(0, n_fft - length(kernel)))))
}

wiener_inverse_filter <- function(h_fft, lambda, signal_psd = NULL, noise_psd = NULL) {
  h_power <- Mod(h_fft)^2
  if (is.null(signal_psd) || is.null(noise_psd)) {
    denom <- h_power + lambda
    numerator <- Conj(h_fft)
  } else {
    denom <- h_power * signal_psd + lambda * noise_psd
    numerator <- Conj(h_fft) * signal_psd
  }
  denom <- pmax(denom, .Machine$double.eps)
  numerator / denom
}

wiener_deconvolve_matrix <- function(y, kernel, lambda, n_latent, signal_psd = NULL,
    noise_psd = NULL, fft_plan = NULL, y_fft = NULL) {
  y <- as.matrix(y)
  if (is.null(fft_plan)) fft_plan <- wiener_fft_plan(nrow(y), n_latent, kernel)
  n_fft <- fft_plan$n_fft
  filt <- wiener_inverse_filter(fft_plan$h_fft, lambda, signal_psd, noise_psd)
  if (is.null(y_fft)) {
    y_fft <- mvfft(rbind(y, matrix(0, nrow = n_fft - nrow(y), ncol = ncol(y))))
  }
  x <- Re(mvfft(y_fft * filt, inverse = TRUE) / n_fft)
  x[seq_len(n_latent), , drop = FALSE]
}

align_wiener_activity <- function(activity, n_time, kernel_len, trim_kernel = TRUE) {
  trimmed <- activity[seq_len(n_time), , drop = FALSE]
  if (isTRUE(trim_kernel)) {
    return(trimmed)
  }
  rbind(matrix(0, nrow = kernel_len - 1L, ncol = ncol(activity)), trimmed)
}

wiener_reconvolve_matrix <- function(activity, kernel, n_time, fft_plan = NULL) {
  activity <- as.matrix(activity)
  n_latent <- nrow(activity)
  if (is.null(fft_plan)) fft_plan <- wiener_fft_plan(n_time, n_latent, kernel)
  n_fft <- fft_plan$n_fft
  h_fft <- fft_plan$h_fft

  activity_pad <- rbind(activity, matrix(0, nrow = n_fft - n_latent, ncol = ncol(activity)))
  fitted_fft <- mvfft(activity_pad) * h_fft
  fitted <- Re(mvfft(fitted_fft, inverse = TRUE) / n_fft)
  fitted[seq_len(n_time), , drop = FALSE]
}

# For S = R C_h P C_f E, the observed diagonal at time t is
# sum_{lag=0}^{min(t,K-1)} h[lag] f[-lag]. Each lag occurs n_time-lag
# times. This is the exact trace of the implemented zero-extended, truncated
# smoother, rather than the trace over all padded Fourier modes.
# With PSD shrinkage it is conditional on the fitted spectra, not a Jacobian df.
wiener_degrees_freedom <- function(filter, kernel, n_time) {
  impulse <- Re(fft(filter, inverse = TRUE) / length(filter))
  lags <- seq_len(min(length(kernel), n_time)) - 1L
  sum((n_time - lags) * kernel[lags + 1L] * impulse[(-lags) %% length(filter) + 1L])
}

estimate_wiener_psd <- function(y, kernel, lambda, n_latent, smooth_width = 7L, psd_floor = 1e-8, fft_plan = NULL, y_fft = NULL) {
  scalar_activity <- wiener_deconvolve_matrix(y, kernel, lambda = lambda, n_latent = n_latent,
    fft_plan = fft_plan, y_fft = y_fft)
  scalar_fitted <- wiener_reconvolve_matrix(scalar_activity, kernel, n_time = nrow(y), fft_plan = fft_plan)
  residual <- y - scalar_fitted
  n_fft <- wiener_fft_length(n_time = nrow(y), n_latent = n_latent, kernel_len = length(kernel))
  signal_psd <- average_matrix_periodogram(scalar_activity, n_fft = n_fft)
  noise_psd <- average_matrix_periodogram(residual, n_fft = n_fft)
  signal_psd <- smooth_periodogram(signal_psd, width = smooth_width)
  noise_psd <- smooth_periodogram(noise_psd, width = smooth_width)
  list(
    signal_psd = pmax(signal_psd, psd_floor),
    noise_psd = pmax(noise_psd, psd_floor)
  )
}

average_matrix_periodogram <- function(x, n_fft) {
  x <- as.matrix(x)
  if (nrow(x) < n_fft) {
    x <- rbind(x, matrix(0, nrow = n_fft - nrow(x), ncol = ncol(x)))
  } else if (nrow(x) > n_fft) {
    x <- x[seq_len(n_fft), , drop = FALSE]
  }
  rowMeans(Mod(mvfft(x))^2) / nrow(x)
}

smooth_periodogram <- function(x, width = 7L) {
  width <- as.integer(width)
  if (width <= 1L) {
    return(as.numeric(x))
  }
  if (width %% 2L == 0L) {
    width <- width + 1L
  }
  pad <- (width - 1L) / 2L
  padded <- c(tail(x, pad), x, head(x, pad))
  kernel <- rep(1 / width, width)
  as.numeric(stats::filter(padded, filter = kernel, sides = 2L))[(pad + 1L):(pad + length(x))]
}

# Construct the m x m quadratic penalty matrix P for latent activity x:
# - "ridge":   I_m (pure L2 amplitude shrinkage)
# - "diff1":   D1' D1 (first-difference roughness penalty)
# - "diff2":   D2' D2 (second-difference curvature penalty)
# - "matern1": 0.10 * I_m + 0.90 * D1' D1 (damped Ornstein-Uhlenbeck / Matérn-1/2
#              prior that eliminates the zero-frequency nullspace of pure D1' D1)
# - "sobolev": 0.05 * I_m + 0.45 * D1' D1 + 0.50 * D2' D2 (critically damped
#              second-order Sobolev prior that suppresses oscillatory sidelobe
#              ringing around sharp impulses while stabilizing low frequencies)
decon_penalty_crossprod <- function(n_latent, penalty) {
  if (penalty == "ridge") {
    return(diag(n_latent))
  }
  if (penalty == "matern1") {
    d1 <- diff(diag(n_latent), differences = 1L)
    return(0.10 * diag(n_latent) + 0.90 * crossprod(d1))
  }
  if (penalty == "sobolev") {
    d1 <- diff(diag(n_latent), differences = 1L)
    d2 <- diff(diag(n_latent), differences = 2L)
    return(0.05 * diag(n_latent) + 0.45 * crossprod(d1) + 0.50 * crossprod(d2))
  }

  differences <- switch(
    penalty,
    diff1 = 1L,
    diff2 = 2L
  )
  dmat <- diff(diag(n_latent), differences = differences)
  crossprod(dmat)
}

decon_convolution_matrix <- function(n_time, kernel) {
  kernel <- as.numeric(kernel)
  kernel_len <- length(kernel)
  n_latent <- n_time + kernel_len - 1L
  conv_mat <- matrix(0, nrow = n_time, ncol = n_latent)

  for (ki in seq_len(kernel_len)) {
    lag <- ki - 1L
    cols <- seq_len(n_time) + kernel_len - 1L - lag
    conv_mat[cbind(seq_len(n_time), cols)] <- kernel[ki]
  }

  conv_mat
}

decon_first_pass_residuals <- function(y, activity, kernels, kernel_index, tune_by) {
  n_time <- nrow(y)
  residuals <- matrix(NA_real_, nrow = n_time, ncol = ncol(y))

  if (tune_by == "global") {
    conv_mat <- decon_convolution_matrix(n_time, kernels[[kernel_index]])
    return(y - conv_mat %*% activity)
  }

  for (si in seq_len(ncol(y))) {
    conv_mat <- decon_convolution_matrix(n_time, kernels[[kernel_index[si]]])
    residuals[, si] <- y[, si] - as.vector(conv_mat %*% activity[, si, drop = FALSE])
  }

  residuals
}
