# Repeated, independently generated simulations for cross-run HRF validation.
# Run with this PR installed, one BLAS/OpenMP thread:
# Rscript --vanilla inst/validation/hrf_recovery.R /tmp/hrf-validation 40
# Sourcing defines helpers without running the experiment.

hrfv_scenarios <- function() {
  base <- list(n = 220L, voxels = 24L, runs = 3L, noise_sd = .04, rho = .5,
    hrf = "shifted", spacing = c(9, 13), feedback_mean = 1, rpe_beta = .6,
    heterogeneous = FALSE, drift = FALSE, neural_duration = 0, masked = FALSE,
    null = FALSE)
  changes <- list(
    canonical = list(hrf = "canonical"),
    canonical_low_snr = list(hrf = "canonical", noise_sd = .15),
    shifted = list(),
    off_library = list(hrf = "off_library"),
    signed_rpe = list(feedback_mean = 0),
    dense_overlap = list(spacing = c(4, 6)),
    heterogeneous_parcels = list(heterogeneous = TRUE),
    two_runs = list(runs = 2L),
    run_drift = list(runs = 2L, hrf = "canonical", drift = TRUE),
    sustained_neural = list(hrf = "canonical", neural_duration = 2.5),
    masked_training = list(masked = TRUE),
    null = list(hrf = "canonical", null = TRUE, rpe_beta = 0)
  )
  lapply(changes, function(change) utils::modifyList(base, change))
}

# Continuous, truncated difference-of-gammas generating model. This does not
# call the package's HRF/event-design/convolution helpers. Duration responses use
# analytic gamma CDF integrals, not the inverse's discrete convolution matrix.
hrfv_cdf <- function(t, peak, dispersion, ratio = 6, support = 32) {
  x <- pmax(0, pmin(support, t))
  mass <- stats::pgamma(support, peak / dispersion, scale = dispersion) -
    stats::pgamma(support, 16, scale = 1) / ratio
  (stats::pgamma(x, peak / dispersion, scale = dispersion) -
    stats::pgamma(x, 16, scale = 1) / ratio) / mass
}

hrfv_response <- function(lag, peak, dispersion, ratio = 6, duration = 0) {
  if (duration > 0) return((hrfv_cdf(lag, peak, dispersion, ratio) -
    hrfv_cdf(lag - duration, peak, dispersion, ratio)) / duration)
  mass <- stats::pgamma(32, peak / dispersion, scale = dispersion) -
    stats::pgamma(32, 16, scale = 1) / ratio
  out <- (stats::dgamma(pmax(lag, 0), peak / dispersion, scale = dispersion) -
    stats::dgamma(pmax(lag, 0), 16, scale = 1) / ratio) / mass
  out[lag < 0 | lag > 32] <- 0
  out
}

hrfv_project_truth <- function(n, onsets, amplitudes, durations) {
  # Resolution-limited latent truth: linearly deposit event mass on the TR=1
  # grid. For sustained neural activity, integrate deposits at 0.025-s steps.
  # Continuous generating BOLD does NOT use this projection.
  truth <- numeric(n)
  for (e in seq_along(onsets)) {
    pieces <- if (durations[e] > 0) ceiling(durations[e] / .025) else 1L
    times <- onsets[e] + if (durations[e] > 0) ((seq_len(pieces) - .5) / pieces) * durations[e] else 0
    for (at in times) {
      lo <- floor(at); fraction <- at - lo
      if (lo >= 0 && lo < n) truth[lo + 1L] <- truth[lo + 1L] + amplitudes[e] * (1 - fraction) / pieces
      if (lo + 1 >= 0 && lo + 1 < n) truth[lo + 2L] <- truth[lo + 2L] + amplitudes[e] * fraction / pieces
    }
  }
  truth
}

hrfv_ar1 <- function(n, rho, sd) {
  as.numeric(stats::filter(stats::rnorm(n + 100L, sd = sd * sqrt(1 - rho^2)),
    filter = rho, method = "recursive"))[-seq_len(100L)]
}

hrfv_generate <- function(config, seed) {
  set.seed(seed)
  n <- config$n; v <- config$voxels
  ids <- paste0("v", seq_len(v)); runs <- paste0("run", seq_len(config$runs))
  groups <- rep(1:2, each = v / 2)
  territory <- groups
  if (config$heterogeneous) {
    flipped <- c(seq_len(v / 8), v / 2 + seq_len(v / 8))
    territory[flipped] <- 3L - territory[flipped]
  }
  gain <- stats::runif(v, .75, 1.25) * rep(c(1, 1, 1, -1), length.out = v)
  bold <- events <- truth <- trials <- kernels <- masks <- vector("list", length(runs))
  names(bold) <- names(events) <- names(truth) <- names(trials) <- names(kernels) <- names(masks) <- runs
  for (r in seq_along(runs)) {
    feedback <- 20 + cumsum(stats::runif(100, config$spacing[1], config$spacing[2]))
    feedback <- feedback[feedback < n - 35]
    rpe <- as.numeric(scale(hrfv_ar1(length(feedback), .3, 1)))
    cue <- feedback - stats::runif(length(feedback), 3, 5)
    response <- feedback - stats::runif(length(feedback), .7, 1.5)
    events[[r]] <- rbind(data.frame(onset = cue, condition = "cue"),
      data.frame(onset = response, condition = "response"),
      data.frame(onset = feedback, condition = "feedback"))
    durations <- c(rep(0, 2 * length(feedback)), rep(config$neural_duration, length(feedback)))
    amplitudes <- c(.3 * (1 + stats::rnorm(length(feedback), sd = .2)),
      -.2 * (1 + stats::rnorm(length(feedback), sd = .2)),
      config$feedback_mean + config$rpe_beta * rpe)
    if (config$feedback_mean == 0) amplitudes[seq_len(2 * length(feedback))] <- 0
    if (config$null) amplitudes[] <- 0
    peaks <- switch(config$hrf, canonical = rep(6, v), shifted = c(5, 7)[territory],
      off_library = c(4.3, 8.3)[territory])
    dispersions <- if (config$hrf == "off_library") c(1.35, .9)[territory] else rep(1, v)
    ratios <- if (config$hrf == "off_library") rep(8, v) else rep(6, v)
    if (config$drift) peaks[] <- c(4.5, 7.5)[r]
    kernels[[r]] <- vapply(seq_len(v), function(j) {
      h <- hrfv_response(0:32, peaks[j], dispersions[j], ratios[j])
      h / sum(h)
    }, numeric(33))
    shared_noise <- hrfv_ar1(n, config$rho, config$noise_sd)
    y <- matrix(0, n, v); x <- matrix(0, n, v)
    projected <- hrfv_project_truth(n, events[[r]]$onset, amplitudes, durations)
    # Reuse continuous generating responses only for identical true shapes.
    shape_ids <- paste(peaks, dispersions, ratios, sep = ":")
    signal <- list()
    for (shape in unique(shape_ids)) {
      j <- match(shape, shape_ids)
      design <- vapply(seq_len(nrow(events[[r]])), function(e) {
        hrfv_response(0:(n - 1L) - events[[r]]$onset[e], peaks[j],
          dispersions[j], ratios[j], durations[e])
      }, numeric(n))
      signal[[shape]] <- as.numeric(design %*% amplitudes)
    }
    for (j in seq_len(v)) {
      noise <- .35 * shared_noise + sqrt(1 - .35^2) * hrfv_ar1(n, config$rho, config$noise_sd)
      y[, j] <- 100 + gain[j] * signal[[shape_ids[j]]] + noise
      x[, j] <- gain[j] * projected
    }
    colnames(y) <- colnames(x) <- colnames(kernels[[r]]) <- ids
    bold[[r]] <- y; truth[[r]] <- x
    trials[[r]] <- data.frame(onset = feedback, rpe = rpe)
    masks[[r]] <- if (config$masked) stats::runif(n) > .2 else rep(TRUE, n)
  }
  list(bold = bold, events = events, truth = truth, trials = trials,
    kernels = kernels, groups = stats::setNames(paste0("parcel", groups), ids), masks = masks)
}

hrfv_sample <- function(y, at) {
  lo <- floor(at); f <- at - lo
  stopifnot(all(lo >= 0), all(lo + 2L <= nrow(y)))
  y[lo + 1L, , drop = FALSE] * (1 - f) + y[lo + 2L, , drop = FALSE] * f
}

hrfv_center <- function(y) sweep(y, 2L, colMeans(y), "-")
hrfv_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
hrfv_cor <- function(x, y) {
  x <- hrfv_center(x); y <- hrfv_center(y)
  denominator <- sqrt(colSums(x^2) * colSums(y^2))
  out <- colSums(x * y) / denominator
  out[denominator < 1e-14] <- NA_real_
  out
}

hrfv_metrics <- function(activity, truth, trials, selected_hrf, true_hrf) {
  interior <- 33:(nrow(activity) - 32L)
  x <- hrfv_center(truth[interior, , drop = FALSE])
  y <- hrfv_center(activity[interior, , drop = FALSE])
  true_sd <- sqrt(colMeans(x^2))
  nrmse <- sqrt(colMeans((x - y)^2)) / true_sd
  nrmse[true_sd < 1e-14] <- NA_real_
  rpe <- trials$rpe - mean(trials$rpe)
  beta <- function(a) as.numeric(crossprod(rpe, a) / sum(rpe^2))
  tx <- hrfv_sample(truth, trials$onset); ty <- hrfv_sample(activity, trials$onset)
  bx <- beta(tx); by <- beta(ty)
  direction <- sign(bx); gain <- by / bx; gain[abs(bx) < 1e-10] <- NA_real_
  rpe_matrix <- matrix(rpe, nrow(tx), ncol(tx))
  area <- function(a) Reduce(`+`, lapply(-1:3, function(lag) hrfv_sample(a, trials$onset + lag)))
  ax <- area(truth); ay <- area(activity); abx <- beta(ax)
  area_gain <- beta(ay) / abx; area_gain[abs(abx) < 1e-10] <- NA_real_
  lag_beta <- function(a) vapply(-2:5, function(lag) beta(hrfv_sample(a, trials$onset + lag)), numeric(ncol(a)))
  peak_x <- max.col(abs(lag_beta(truth)), ties.method = "first")
  peak_y <- max.col(abs(lag_beta(activity)), ties.method = "first")
  peak_error <- abs(peak_x - peak_y); peak_error[abs(bx) < 1e-10] <- NA_real_
  lagged <- cbind(1, rpe[-1], rpe[-length(rpe)])
  neighbor_beta <- function(a) stats::lm.fit(lagged, a[-1, , drop = FALSE])$coefficients[3, ]
  leakage <- abs(neighbor_beta(ty) - neighbor_beta(tx)) / abs(bx)
  leakage[abs(bx) < 1e-10] <- NA_real_
  rpe_cor <- hrfv_cor(rpe_matrix, ty)
  c(hrf_nmse = mean(colSums((selected_hrf - true_hrf)^2) / colSums(true_hrf^2)),
    hrf_peak_error = mean(abs(max.col(t(selected_hrf), ties.method = "first") -
      max.col(t(true_hrf), ties.method = "first"))),
    neural_r = hrfv_mean(hrfv_cor(x, y)), neural_nrmse = hrfv_mean(nrmse),
    rpe_r = if (any(abs(bx) > 1e-10)) hrfv_mean(direction * rpe_cor) else NA_real_,
    rpe_gain = hrfv_mean(gain), rpe_gain_error = hrfv_mean(abs(gain - 1)),
    rpe_area_r = if (any(abs(abx) > 1e-10))
      hrfv_mean(sign(abx) * hrfv_cor(rpe_matrix, ay)) else NA_real_,
    rpe_area_gain = hrfv_mean(area_gain), rpe_area_gain_error = hrfv_mean(abs(area_gain - 1)),
    epoch_peak_error = hrfv_mean(peak_error),
    neighbor_leakage = hrfv_mean(leakage), noise_rms = mean(sqrt(colMeans(y^2))),
    null_abs_rpe_r = if (all(abs(bx) < 1e-10)) hrfv_mean(abs(rpe_cor)) else NA_real_)
}

hrfv_operator <- function(kernel, n, lambda, cache) {
  key <- paste(n, lambda, paste(sprintf("%.17g", kernel), collapse = ","), sep = ":")
  if (!exists(key, envir = cache, inherits = FALSE)) {
    operator <- fmri.pipeline::deconvolve_reglin(diag(n), kernel, lambda = lambda,
      scaling = "none", backend = "cpp", penalty = "diff2")
    assign(key, operator, envir = cache)
  }
  get(key, envir = cache, inherits = FALSE)
}

hrfv_apply <- function(bold, kernels, lambda, cache) {
  out <- matrix(NA_real_, nrow(bold), ncol(bold))
  keys <- apply(kernels, 2, function(h) paste(sprintf("%.17g", h), collapse = ","))
  y <- hrfv_center(bold)
  for (key in unique(keys)) {
    indices <- which(keys == key)
    out[, indices] <- hrfv_operator(kernels[, indices[1]], nrow(bold), lambda, cache) %*%
      y[, indices, drop = FALSE]
  }
  out
}

hrfv_summary <- function(metrics) {
  identifiers <- c("scenario", "replicate", "method", "lambda")
  measures <- setdiff(names(metrics), identifiers)
  groups <- split(metrics, interaction(metrics$scenario, metrics$method, metrics$lambda, drop = TRUE))
  rows <- lapply(groups, function(d) {
    key <- d[1, c("scenario", "method", "lambda")]
    do.call(rbind, lapply(measures, function(measure) {
      x <- d[[measure]]; x <- x[is.finite(x)]; n <- length(x)
      se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
      data.frame(key, metric = measure, n = n, mean = if (n) mean(x) else NA_real_,
        mc_se = se, lower = if (n > 1) mean(x) - stats::qt(.975, n - 1) * se else NA_real_,
        upper = if (n > 1) mean(x) + stats::qt(.975, n - 1) * se else NA_real_)
    }))
  })
  summary <- do.call(rbind, rows); rownames(summary) <- NULL
  baseline <- metrics[metrics$method == "canonical", c("scenario", "replicate", "lambda", measures)]
  paired <- merge(metrics[metrics$method != "canonical", ], baseline,
    by = c("scenario", "replicate", "lambda"), suffixes = c("", "_canonical"))
  delta <- paired[identifiers]
  for (measure in measures) delta[[measure]] <- paired[[measure]] - paired[[paste0(measure, "_canonical")]]
  # Reuse the aggregation without recursion for the paired table.
  delta_groups <- split(delta, interaction(delta$scenario, delta$method, delta$lambda, drop = TRUE))
  paired_rows <- lapply(delta_groups, function(d) {
    key <- d[1, c("scenario", "method", "lambda")]
    do.call(rbind, lapply(measures, function(measure) {
      x <- d[[measure]]; x <- x[is.finite(x)]; n <- length(x)
      se <- if (n > 1) stats::sd(x) / sqrt(n) else NA_real_
      larger_better <- measure %in% c("neural_r", "rpe_r", "rpe_area_r")
      data.frame(key, metric = measure, n = n, delta = if (n) mean(x) else NA_real_,
        mc_se = se, lower = if (n > 1) mean(x) - stats::qt(.975, n - 1) * se else NA_real_,
        upper = if (n > 1) mean(x) + stats::qt(.975, n - 1) * se else NA_real_,
        improved_fraction = if (n && measure %in% c("neural_r", "rpe_r", "rpe_area_r",
          "neural_nrmse", "hrf_nmse", "epoch_peak_error", "neighbor_leakage",
          "rpe_gain_error", "rpe_area_gain_error"))
          mean(if (larger_better) x > 1e-10 else x < -1e-10) else NA_real_)
    }))
  })
  list(summary = summary, paired = do.call(rbind, paired_rows))
}

run_hrf_validation <- function(out_dir, repetitions = 40L, scenario_names = names(hrfv_scenarios())) {
  stopifnot(repetitions >= 1L, repetitions == as.integer(repetitions))
  if (!requireNamespace("fmri.pipeline", quietly = TRUE)) stop("Install this PR's fmri.pipeline first")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- hrfv_scenarios()[scenario_names]
  stopifnot(!any(vapply(scenarios, is.null, logical(1))))
  library <- fmri.pipeline::hrf_library(1)
  lambdas <- c(.001, .01, .1) # .01 is primary; all are reported, none tuned.
  cache <- new.env(parent = emptyenv())
  methods <- list(condition_voxel = list(amplitude_model = "condition"),
    trial_voxel = list(amplitude_model = "trial"),
    trial_parcel = list(amplitude_model = "trial", pooling = "parcel"),
    trial_spatial = list(amplitude_model = "trial", spatial_strength = .02))
  all_metrics <- all_runs <- list(); checked_public <- FALSE
  start <- proc.time()[["elapsed"]]
  for (s in seq_along(scenarios)) {
    scenario <- names(scenarios)[s]; config <- scenarios[[s]]
    message("Scenario ", scenario, ": ", repetitions, " independent datasets")
    for (replicate in seq_len(repetitions)) {
      seed <- 8320000L + match(scenario, names(hrfv_scenarios())) * 10000L + replicate
      data <- hrfv_generate(config, seed)
      fitted <- lapply(names(methods), function(method) {
        settings <- c(list(bold = data$bold, events = data$events, tr = 1, kernel_grid = library,
          valid_rows = data$masks), methods[[method]])
        if (method %in% c("trial_parcel", "trial_spatial")) settings$spatial_groups <- data$groups
        do.call(fmri.pipeline::estimate_crossrun_hrf, settings)
      })
      names(fitted) <- names(methods)
      selections <- list(canonical = lapply(data$bold, function(y) library[, rep(1L, ncol(y)), drop = FALSE]),
        oracle_hrf = data$kernels)
      for (method in names(fitted)) selections[[method]] <- lapply(fitted[[method]]$selections,
        function(d) library[, d$kernel_index, drop = FALSE])
      for (lambda in lambdas) for (method in names(selections)) {
        rows <- vector("list", length(data$bold))
        for (r in seq_along(data$bold)) {
          activity <- hrfv_apply(data$bold[[r]], selections[[method]][[r]], lambda, cache)
          values <- hrfv_metrics(activity, data$truth[[r]], data$trials[[r]],
            selections[[method]][[r]], data$kernels[[r]])
          values <- c(values, fallback_fraction = if (method %in% names(fitted))
            mean(grepl("^fallback", fitted[[method]]$selections[[r]]$status)) else NA_real_)
          rows[[r]] <- data.frame(scenario = scenario, replicate = replicate, method = method,
            lambda = lambda, run = names(data$bold)[r], as.list(values))
        }
        run_metrics <- do.call(rbind, rows)
        all_runs[[length(all_runs) + 1L]] <- run_metrics
        measures <- setdiff(names(run_metrics), c("scenario", "replicate", "method", "lambda", "run"))
        all_metrics[[length(all_metrics) + 1L]] <- data.frame(scenario = scenario,
          replicate = replicate, method = method, lambda = lambda,
          as.list(vapply(run_metrics[measures], hrfv_mean, numeric(1))))
      }
      if (!checked_public) {
        fit <- fitted$trial_spatial
        direct <- fmri.pipeline::deconvolve_crossrun_hrf(data$bold, fit, .01, scaling = "demean")
        differences <- vapply(seq_along(data$bold), function(r) {
          expected <- hrfv_apply(data$bold[[r]], selections$trial_spatial[[r]], .01, cache)
          max(abs(expected - direct$activity[[r]]))
        }, numeric(1))
        stopifnot(max(differences) < 1e-9)
        writeLines(paste("Cached evaluation versus public application max abs difference:",
          max(differences)), file.path(out_dir, "public_application_check.txt"))
        checked_public <- TRUE
      }
      if (replicate %% 10L == 0L || replicate == repetitions) {
        write.csv(do.call(rbind, all_metrics), file.path(out_dir, "replicate_metrics.csv"), row.names = FALSE)
        message("  completed ", replicate, "; elapsed ", round(proc.time()[["elapsed"]] - start, 1), " s")
      }
    }
  }
  metrics <- do.call(rbind, all_metrics)
  write.csv(do.call(rbind, all_runs), file.path(out_dir, "run_metrics.csv"), row.names = FALSE)
  tables <- hrfv_summary(metrics)
  write.csv(tables$summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(tables$paired, file.path(out_dir, "paired_vs_canonical.csv"), row.names = FALSE)
  saveRDS(list(scenarios = scenarios, repetitions = repetitions, lambdas = lambdas,
    primary_lambda = .01, methods = methods, library = library,
    seed_formula = "8320000 + scenario_index * 10000 + replicate", tr = 1,
    nuisance = "Measured-run centering; training also fits intercept and linear trend",
    inference_unit = "independently generated multirun dataset, not voxel or run",
    package_version = as.character(utils::packageVersion("fmri.pipeline")),
    package_path = find.package("fmri.pipeline"),
    script_md5 = if (file.exists("inst/validation/hrf_recovery.R"))
      unname(tools::md5sum("inst/validation/hrf_recovery.R")) else NA_character_),
    file.path(out_dir, "settings.rds"))
  writeLines(c(capture.output(sessionInfo()), paste("Elapsed seconds:", proc.time()[["elapsed"]] - start),
    paste("OPENBLAS_NUM_THREADS", Sys.getenv("OPENBLAS_NUM_THREADS")),
    paste("OMP_NUM_THREADS", Sys.getenv("OMP_NUM_THREADS"))), file.path(out_dir, "environment.txt"))
  invisible(tables)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_hrf_validation(if (length(args)) args[1] else "/tmp/hrf-validation",
    if (length(args) > 1L) as.integer(args[2]) else 40L)
}
