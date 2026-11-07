# Event-blind continuous waveform validation. Forward generation is independent
# of package HRF/event/convolution routines. TR=1 s; no upsampled-resolution claim.
wv_sources <- function() file.path("inst", "validation", c("waveform_solver.cpp",
  "waveform_recovery.R", "run_waveform_recovery.R"))
wv_load_solver <- function(cache_dir = file.path(tempdir(), "waveform-sourceCpp")) {
  if (!exists("wave_solve_cpp", mode = "function", inherits = TRUE))
    Rcpp::sourceCpp(file.path("inst", "validation", "waveform_solver.cpp"),
      env = environment(wv_load_solver), cacheDir = cache_dir)
  invisible(TRUE)
}
wv_scenarios <- function() {
  base <- list(n = 256L, voxels = 8L, runs = 3L, family = "brief", noise_sd = .04,
    noise_rho = .5, mask = FALSE, peaks = c(5, 7), heterogeneous = FALSE)
  changes <- list(brief = list(), blocks = list(family = "block"),
    smooth_heldout = list(family = "smooth"), biphasic_heldout = list(family = "biphasic"),
    neural_overlap = list(family = "overlap"),
    hrf_mismatch = list(family = "smooth", peaks = c(4.3, 8.3)),
    censored_heldout = list(family = "smooth", mask = TRUE),
    spatial_heterogeneity = list(family = "smooth", heterogeneous = TRUE))
  lapply(changes, function(x) utils::modifyList(base, x))
}
wv_grid <- function() {
  out <- rbind(data.frame(method = "diff2", lambda = c(.001, .01, .1, 1), spatial = 0),
    data.frame(method = "tv", lambda = c(.003, .01, .03, .1), spatial = 0),
    data.frame(method = "mixture", lambda = c(.003, .01, .03, .1), spatial = 0),
    data.frame(method = "tv_spatial", lambda = c(.003, .01, .03, .1), spatial = .01))
  out$id <- paste(out$method, format(out$lambda, scientific = FALSE, trim = TRUE), sep = "/")
  out
}
wv_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
wv_cor <- function(a, b) {
  a <- as.matrix(a); b <- as.matrix(b)
  vapply(seq_len(ncol(a)), function(j) if (stats::sd(a[,j]) > 1e-10 && stats::sd(b[,j]) > 1e-10)
    stats::cor(a[,j], b[,j]) else NA_real_, numeric(1))
}
wv_hrf <- function(t, peak) {
  value <- (stats::dgamma(pmax(t, 0), peak) - stats::dgamma(pmax(t, 0), 16) / 6) /
    (stats::pgamma(32, peak) - stats::pgamma(32, 16) / 6)
  value[t < 0 | t > 32] <- 0
  value
}
wv_ar <- function(n, rho, sd) as.numeric(stats::filter(stats::rnorm(n + 200) * sd * sqrt(1-rho^2),
  rho, method = "recursive"))[-seq_len(200)]

# Smooth compact pulses, blocks and biphasic shapes at a 0.1-second forward
# grid. Per-trial/per-voxel width and latency vary independently of RPE.
wv_shape <- function(t, duration, family) {
  u <- t / duration; inside <- u >= 0 & u < 1
  y <- numeric(length(t))
  if (family %in% c("brief", "block", "overlap")) y[inside] <- 1
  if (family == "smooth") y[inside] <- sin(pi * u[inside])^2
  if (family == "biphasic") y[inside] <- sin(2*pi*u[inside]) * (1 - .3*u[inside])
  y
}
wv_forward <- function(neural, peak, step = .1) {
  h <- wv_hrf(seq(0, 32, by = step), peak)
  # stats::convolve is an independent FFT forward implementation, with an
  # explicit reversal for convolution and sufficient zero padding.
  full <- stats::convolve(neural, rev(h), type = "open") * step
  full[seq(1L, length(neural), by = as.integer(round(1 / step)))]
}
wv_generate <- function(config, seed) {
  stopifnot(config$n >= 160L, config$voxels %% 2L == 0L, config$runs >= 3L)
  set.seed(seed); n <- config$n; v <- config$voxels; step <- .1
  time <- seq(0, n - 1, by = step); sampled <- seq(1L, length(time), by = 10L)
  groups <- rep(1:2, each = v / 2L); ids <- paste0("v", seq_len(v))
  gain <- stats::runif(v, .8, 1.2)
  if (config$heterogeneous) gain[seq(2L, v, by = 2L)] <- -gain[seq(2L, v, by = 2L)]
  kernels <- vapply(config$peaks[groups], function(p) { h <- wv_hrf(0:32, p); h / sum(h) }, numeric(33))
  out <- list(config = config, gain = gain, kernels = kernels,
    groups = stats::setNames(paste0("parcel", groups), ids))
  fields <- c("bold", "without_previous", "truth", "previous", "target", "background",
    "signal", "noise", "events", "pairs", "masks")
  for (field in fields) out[[field]] <- stats::setNames(vector("list", config$runs), paste0("run", seq_len(config$runs)))
  for (r in seq_len(config$runs)) {
    bases <- seq(36, n - 50, by = 40) + stats::runif(1, -.4, .4)
    gap <- rep(c(3, 6, 10), length.out = length(bases))
    gap <- gap[(seq_along(gap) + r - 2L) %% length(gap) + 1L]
    m <- length(bases); null <- seq_len(m) %% 2L == 1L
    rpe <- stats::rnorm(m)
    previous <- target <- matrix(0, length(time), v)
    # This background and an unmarked event are retained in neural truth.
    background <- vapply(seq_len(v), function(j) .04*sin(time/13 + stats::runif(1, 0, 6)) +
      .12 * wv_shape(time - (bases[1L] + 22), 4.7, "smooth"), numeric(length(time)))
    for (i in seq_len(m)) for (j in seq_len(v)) {
      family <- config$family
      duration <- switch(family, brief = stats::runif(1, .8, 2), block = stats::runif(1, 3, 6),
        smooth = stats::runif(1, 2.5, 6.5), biphasic = stats::runif(1, 4, 8), overlap = stats::runif(1, 11, 14))
      latency <- stats::runif(1, 0, .8)
      previous[,j] <- previous[,j] + gain[j] * stats::runif(1, 1.5, 2.5) *
        wv_shape(time - bases[i] - latency, duration, family)
      target_duration <- if (family == "overlap") stats::runif(1, 3, 6) else duration * stats::runif(1, .7, 1.3)
      amplitude <- if (null[i]) 0 else gain[j] * (1 + .5*rpe[i] + stats::rnorm(1, 0, .3))
      target[,j] <- target[,j] + amplitude * wv_shape(time - bases[i] - gap[i] -
        stats::runif(1, 0, .8), target_duration, if (family == "overlap") "block" else family)
    }
    signal <- previous_signal <- matrix(0, n, v)
    for (j in seq_len(v)) {
      signal[,j] <- wv_forward(previous[,j] + target[,j] + background[,j], config$peaks[groups[j]])
      previous_signal[,j] <- wv_forward(previous[,j], config$peaks[groups[j]])
    }
    common <- wv_ar(n, config$noise_rho, config$noise_sd)
    noise <- vapply(seq_len(v), function(j) .35 * common + sqrt(1-.35^2) *
      wv_ar(n, config$noise_rho, config$noise_sd), numeric(n))
    y <- 100 + signal + noise
    mask <- rep(TRUE, n)
    if (config$mask) { missing <- c(82:88,139:145,185:190); mask[missing[missing<=n]] <- FALSE }
    colnames(y) <- colnames(kernels) <- ids
    out$bold[[r]] <- y
    out$without_previous[[r]] <- y - previous_signal
    out$truth[[r]] <- (previous + target + background)[sampled,,drop=FALSE]
    out$previous[[r]] <- previous[sampled,,drop=FALSE]
    out$target[[r]] <- target[sampled,,drop=FALSE]
    out$background[[r]] <- background[sampled,,drop=FALSE]
    out$signal[[r]] <- signal; out$noise[[r]] <- noise
    out$events[[r]] <- data.frame(onset = c(bases, bases + gap), condition = "event")
    out$pairs[[r]] <- data.frame(previous = bases, onset = bases + gap, gap = gap, null = null, rpe = rpe)
    out$masks[[r]] <- mask
  }
  out
}

# Independently constructed finite convolution, including 32 pre-run samples.
wv_design <- function(n, h, rows, rho, nuisance = cbind(1, seq(-1, 1, length.out = n))) {
  stopifnot(n>=3L,length(h)>=2L,all(is.finite(h)),length(rows)>=3L,
    all(rows>=1L & rows<=n),!is.unsorted(rows,strictly=TRUE),is.finite(rho),abs(rho)<1,
    nrow(nuisance)==n,ncol(nuisance)>=1L,all(is.finite(nuisance)))
  p <- n + length(h) - 1L
  lag <- outer(0:(n-1L), -(length(h)-1L):(n-1L), "-")
  forward <- matrix(0, n, p); valid <- lag >= 0 & lag < length(h)
  forward[valid] <- h[lag[valid] + 1L]
  w <- diag(length(rows)); w[1L,1L] <- sqrt(1-rho^2)
  for (i in seq.int(2L, length(rows))) {
    a <- rho^(rows[i]-rows[i-1L]); b <- sqrt((1-rho^2)/(1-a^2))
    w[i,i] <- b; w[i,i-1L] <- -a*b
  }
  q <- qr.Q(qr(w %*% nuisance[rows,,drop=FALSE]))
  transform <- w - q %*% crossprod(q, w)
  list(a = transform %*% forward[rows,,drop=FALSE], transform = transform, forward = forward)
}
wv_laplacian <- function(v) diag(v) - matrix(1/v, v, v)

# All preprocessing, including counterfactual fits, uses the observed-data scale
# once. No event timestamps, RPEs, or neural truth enter reconstruction.
wv_reconstruct <- function(y, h, rho, mask, setting, groups, cache = new.env(parent = emptyenv()),
    scale_reference = y, tolerance = 1e-5, max_iterations = 12000L) {
  n <- nrow(y); v <- ncol(y); rows <- which(mask)
  scales <- apply(scale_reference[rows,,drop=FALSE], 2L, stats::sd)
  if (any(!is.finite(scales) | scales <= 0)) stop("Nonpositive BOLD scale")
  scaled <- sweep(y, 2L, colMeans(y[rows,,drop=FALSE]), "-")
  scaled <- sweep(scaled, 2L, scales, "/")
  # Spatial regularization is in standardized units, within supplied parcels
  # and exact HRF groups; no sign-based regrouping is allowed.
  keys <- paste(apply(h, 2L, function(z) paste(sprintf("%.17g", z), collapse = ",")),
    sprintf("%.17g", rho), if (setting$spatial > 0) groups else "all")
  activity <- matrix(NA_real_, n, v); diagnostic <- list()
  started <- proc.time()[["elapsed"]]
  for (key in unique(keys)) {
    j <- which(keys == key); first <- j[1L]
    method <- if (setting$method == "tv_spatial") "tv" else setting$method
    plan_key <- paste(key, paste(rows, collapse = ","), n, length(j), method, setting$spatial,
      if (method == "diff2") setting$lambda else "common", sep = "|")
    if (!exists(plan_key, envir = cache, inherits = FALSE)) {
      design <- wv_design(n, h[,first], rows, rho[first])
      pointer <- wave_prepare_cpp(design$a, method, wv_laplacian(length(j)),
        setting$spatial, setting$lambda)
      assign(plan_key, list(pointer = pointer, transform = design$transform), cache)
    }
    plan <- get(plan_key, cache, inherits = FALSE)
    answer <- wave_solve_cpp(plan$pointer, plan$transform %*% scaled[rows,j,drop=FALSE],
      setting$lambda, max_iterations, tolerance / 10, tolerance)
    activity[,j] <- sweep(answer$activity[length(h[,first])-1L+seq_len(n),,drop=FALSE], 2L, scales[j], "*")
    diagnostic[[length(diagnostic)+1L]] <- data.frame(voxel = j, iterations = answer$iterations,
      converged = as.logical(answer$converged), primal = as.numeric(answer$primal),
      dual = as.numeric(answer$dual), primal_limit = as.numeric(answer$primal_limit), dual_limit = as.numeric(answer$dual_limit))
  }
  list(activity = activity, diagnostic = do.call(rbind, diagnostic), elapsed = proc.time()[["elapsed"]]-started)
}

wv_interpolate <- function(x, times) {
  if (any(times < 0 | times > nrow(x)-1L)) stop("Epoch outside acquisition grid")
  lo <- floor(times); frac <- times-lo
  (1-frac)*x[lo+1L,,drop=FALSE] + frac*x[pmin(lo+2L,nrow(x)),,drop=FALSE]
}
wv_align_offset <- function(estimate, truth, rows) sweep(estimate, 2L,
  colMeans(estimate[rows,,drop=FALSE] - truth[rows,,drop=FALSE]), "-")

wv_metrics <- function(answer, counterfactual, data, r) {
  n <- data$config$n; v <- data$config$voxels; truth <- data$truth[[r]]
  rows <- 33:(n-33L); rows <- rows[data$masks[[r]][rows]]
  estimate <- wv_align_offset(answer$activity, truth, rows)
  difference <- wv_align_offset(answer$activity-counterfactual$activity, data$previous[[r]], rows)
  pair <- data$pairs[[r]]; lags <- -2:14
  epoch <- leakage <- list()
  epochs_estimate <- epochs_truth <- array(NA_real_,c(length(lags),nrow(pair),v))
  for (i in seq_len(nrow(pair))) {
    at <- pair$onset[i] + lags
    e <- wv_interpolate(estimate, at); t <- wv_interpolate(truth, at)
    valid <- data$masks[[r]][floor(at)+1L] & data$masks[[r]][pmin(floor(at)+2L,n)]
    epochs_estimate[valid,i,] <- e[valid,,drop=FALSE]
    epochs_truth[valid,i,] <- t[valid,,drop=FALSE]
    e <- e[valid,,drop=FALSE]; t <- t[valid,,drop=FALSE]
    denom <- colSums(t^2)
    epoch[[i]] <- c(waveform_r = wv_mean(wv_cor(e,t)),
      waveform_nmse = mean(colSums((e-t)^2)/pmax(denom,.01)),
      peak_absolute_error = mean(abs(apply(abs(e),2L,max)-apply(abs(t),2L,max))),
      area_absolute_error = mean(abs(colSums(e)-colSums(t))),
      peak_latency_error = mean(abs(vapply(seq_len(v),function(j)
        lags[valid][which.max(abs(e[,j]))]-lags[valid][which.max(abs(t[,j]))],numeric(1)))),
      effective_width_error = mean(abs(colSums(abs(e))/pmax(apply(abs(e),2L,max),1e-10) -
        colSums(abs(t))/pmax(apply(abs(t),2L,max),1e-10))))
    at <- pair$onset[i] + 0:3
    delta <- wv_interpolate(difference, at); previous <- wv_interpolate(data$previous[[r]], at)
    valid <- data$masks[[r]][floor(at)+1L] & data$masks[[r]][pmin(floor(at)+2L,n)]
    scale <- apply(abs(data$previous[[r]]),2L,max)
    leakage[[i]] <- data.frame(gap = pair$gap[i], null = pair$null[i],
      excess_leakage = mean(sqrt(colMeans((delta[valid,,drop=FALSE]-previous[valid,,drop=FALSE])^2))/scale),
      previous_activity = mean(sqrt(colMeans(previous[valid,,drop=FALSE]^2))/scale),
      recovered_previous_activity = mean(sqrt(colMeans(delta[valid,,drop=FALSE]^2))/scale),
      target_error = mean(sqrt(colMeans((wv_interpolate(estimate,at)[valid,,drop=FALSE] -
        wv_interpolate(truth,at)[valid,,drop=FALSE])^2))/scale))
  }
  leaks <- do.call(rbind, leakage)
  # Remove the across-trial mean separately at each lag. This checks whether
  # recovery retains individual-trial departures from an average waveform.
  deviations <- function(x) {
    for(j in seq_len(v)) x[,,j] <- x[,,j] - rowMeans(x[,,j],na.rm=TRUE)
    x
  }
  de <- deviations(epochs_estimate); dt <- deviations(epochs_truth)
  deviation_r <- deviation_nmse <- numeric(v)
  for(j in seq_len(v)) {
    valid <- is.finite(de[,,j]) & is.finite(dt[,,j])
    a <- de[,,j][valid]; b <- dt[,,j][valid]
    deviation_r[j] <- if(stats::sd(a)>1e-10 && stats::sd(b)>1e-10)stats::cor(a,b)else NA_real_
    deviation_nmse[j] <- sum((a-b)^2)/max(sum(b^2),1e-10)
  }
  metrics <- c(neural_r = wv_mean(wv_cor(estimate[rows,,drop=FALSE],truth[rows,,drop=FALSE])),
    neural_nmse = mean(colSums((estimate[rows,,drop=FALSE]-truth[rows,,drop=FALSE])^2)/colSums(truth[rows,,drop=FALSE]^2)),
    unaligned_neural_nmse = mean(colSums((answer$activity[rows,,drop=FALSE]-truth[rows,,drop=FALSE])^2)/colSums(truth[rows,,drop=FALSE]^2)),
    trial_deviation_r = wv_mean(deviation_r), trial_deviation_nmse = mean(deviation_nmse),
    rowMeans(do.call(cbind,epoch)), excess_leakage = mean(leaks$excess_leakage),
    previous_activity = mean(leaks$previous_activity),
    recovered_previous_activity = mean(leaks$recovered_previous_activity),
    null_target_excess = mean(leaks$excess_leakage[leaks$null]),
    null_target_error = mean(leaks$target_error[leaks$null]),
    convergence = mean(c(answer$diagnostic$converged,counterfactual$diagnostic$converged)),
    iterations = max(c(answer$diagnostic$iterations,counterfactual$diagnostic$iterations)),
    elapsed = answer$elapsed+counterfactual$elapsed)
  list(metrics = metrics, leakage = leaks)
}

wv_hrf_fit <- function(data) {
  fmri.pipeline::estimate_crossrun_hrf(data$bold, data$events, tr = 1,
    kernel_grid = fmri.pipeline::hrf_library(1), valid_rows = data$masks,
    spatial_groups = data$groups, amplitude_model = "trial", spatial_strength = .02)
}
wv_kernels <- function(data, source, fit = NULL) {
  lapply(seq_along(data$bold), function(r) {
    if (source == "oracle") return(list(h = data$kernels, rho = rep(data$config$noise_rho,data$config$voxels)))
    if (source == "canonical") return(list(h = matrix(rep(fmri.pipeline::hrf_library(1)[,1L],data$config$voxels),
      ncol=data$config$voxels), rho = rep(data$config$noise_rho,data$config$voxels)))
    list(h = fit$kernels[,fit$selections[[r]]$kernel_index,drop=FALSE], rho = fit$noise[[r]]$rho)
  })
}
