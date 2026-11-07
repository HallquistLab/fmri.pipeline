# Run from the repository root, with the current package installed.
# FMRI_VALIDATION_LIBRARY=/path/to/library OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
#   Rscript --vanilla inst/validation/run_parametric_recovery.R /tmp/parametric evaluation
source(file.path("inst", "validation", "parametric_recovery.R"), local = TRUE)

pr_methods <- function() c("canonical_estimated_noise", "learned_estimated_noise",
  "learned_frozen_known_noise", "oracle_hrf_estimated_noise", "oracle_hrf_known_noise",
  "learned_frozen_noiseless", "oracle_noiseless_continuous", "oracle_noiseless_discrete",
  "oracle_noisy_discrete", "latent_truth")

pr_fit <- function(data) {
  fmri.pipeline::estimate_crossrun_hrf(data$bold, data$events, tr = 1,
    kernel_grid = fmri.pipeline::hrf_library(1), spatial_groups = data$groups,
    spatial_strength = .02, canonical_strength = .01, amplitude_model = "trial")
}

# All references share schedules, latent amplitudes, and noise innovations.
# Noiseless references freeze the HRFs/covariance selected from noisy data.
pr_recover <- function(data, fit, method, lambda, cache) {
  lapply(seq_along(data$bold), function(r) {
    if (method == "latent_truth") return(list(activity = data$truth[[r]], eta = data$eta_truth[[r]]))
    selected <- fit$kernels[, fit$selections[[r]]$kernel_index, drop = FALSE]
    rho <- fit$noise[[r]]$rho
    h <- if (grepl("^oracle", method)) data$kernels else selected
    if (method == "canonical_estimated_noise") h <- fit$kernels[, rep(1L, ncol(selected)), drop = FALSE]
    if (method %in% c("learned_frozen_known_noise", "oracle_hrf_known_noise",
        "oracle_noiseless_continuous", "oracle_noiseless_discrete", "oracle_noisy_discrete"))
      rho[] <- data$config$noise_rho
    discrete <- method %in% c("oracle_noiseless_discrete", "oracle_noisy_discrete")
    noiseless <- grepl("noiseless", method)
    signal <- if (discrete) data$discrete_signal[[r]] else data$signal[[r]]
    eta <- if (discrete) data$discrete_eta_signal[[r]] else data$eta_signal[[r]]
    y <- signal + if (noiseless) 0 else data$noise[[r]]
    pr_apply(y, eta, h, rho, lambda, cache)
  })
}

pr_contrasts <- function() list(
  shape_estimation = c("learned_estimated_noise", "oracle_hrf_estimated_noise"),
  covariance_application = c("learned_estimated_noise", "learned_frozen_known_noise"),
  covariance_with_true_shape = c("oracle_hrf_estimated_noise", "oracle_hrf_known_noise"),
  adaptive_vs_canonical = c("learned_estimated_noise", "canonical_estimated_noise"),
  measurement_noise = c("oracle_hrf_known_noise", "oracle_noiseless_continuous"),
  measurement_noise_fitted = c("learned_estimated_noise", "learned_frozen_noiseless"),
  sampling_mismatch = c("oracle_noiseless_continuous", "oracle_noiseless_discrete"),
  sampling_mismatch_noisy = c("oracle_hrf_known_noise", "oracle_noisy_discrete"))

pr_paired <- function(metrics) {
  ids <- c("scenario", "replicate", "method", "lambda")
  measures <- setdiff(names(metrics), ids)
  do.call(rbind, lapply(names(pr_contrasts()), function(contrast) {
    pair <- pr_contrasts()[[contrast]]
    a <- metrics[metrics$method == pair[1L], ]; b <- metrics[metrics$method == pair[2L], ]
    d <- merge(a, b, by = c("scenario", "replicate", "lambda"), suffixes = c("_a", "_b"))
    if (!nrow(d)) return(NULL)
    delta <- d[c("scenario", "replicate", "lambda")]; delta$method <- contrast
    for (measure in measures) delta[[measure]] <- d[[paste0(measure, "_a")]] - d[[paste0(measure, "_b")]]
    pr_summarize(delta, probability_bounds = FALSE)
  }))
}

run_parametric_recovery <- function(out_dir, phase = c("evaluation", "development"),
    repetitions = NULL, null_repetitions = NULL, scenario_names = names(pr_scenarios())) {
  phase <- match.arg(phase)
  if (is.null(repetitions)) repetitions <- if (phase == "evaluation") 40L else 3L
  if (is.null(null_repetitions)) null_repetitions <- if (phase == "evaluation") 200L else repetitions
  stopifnot(repetitions >= 1L, null_repetitions >= 1L,
    repetitions == as.integer(repetitions), null_repetitions == as.integer(null_repetitions))
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
    stop("Use an empty output directory; completed/frozen experiments are never overwritten")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- pr_scenarios()[scenario_names]
  stopifnot(length(scenarios) > 0L, !any(vapply(scenarios, is.null, logical(1))))
  lambdas <- c(0, 10^seq(-6, 0))
  # Fixed, disjoint seed blocks; scenario IDs are stable under subsetting.
  seed_base <- if (phase == "evaluation") 18000000L else 16000000L
  counts <- vapply(scenarios, function(x) if (x$beta == 0) null_repetitions else repetitions, integer(1))
  stopifnot(max(counts) < 10000L)
  files <- file.path("inst", "validation", c("parametric_recovery.R", "run_parametric_recovery.R"))
  settings <- list(phase = phase, scenarios = scenarios, counts = counts, lambdas = lambdas,
    primary_lambda = .01, methods = pr_methods(), contrasts = pr_contrasts(),
    primary_metrics = c("amplitude_rmse", "eta_observed_nrmse", "adjusted_beta_rmse",
      "adjusted_coverage", "adjusted_false_positive"),
    seed_base = seed_base, seed_formula = "base + original_scenario_index * 10000 + replicate",
    script_md5 = tools::md5sum(files), inference = "two-sided run-slope t interval, df=runs-1; alpha=.05",
    gain_convention = "unit signed HRF mass, demeaned BOLD; fractional onset geometry only",
    development_use = "implementation/debugging only; no data-dependent method or lambda selection",
    package_description = utils::packageDescription("fmri.pipeline")[c("Package", "Version", "Built")])
  # Freeze design, seed split, methods, and primary endpoint before generating data.
  saveRDS(settings, file.path(out_dir, "settings.rds"))
  metrics <- diagnostics <- list(); start <- proc.time()[["elapsed"]]
  # Fixed n/library/rho strata; reuse native factorizations across datasets.
  cache <- new.env(parent = emptyenv())
  for (scenario in names(scenarios)) {
    message(scenario, ": ", counts[[scenario]], " independent datasets")
    for (replicate in seq_len(counts[[scenario]])) {
      seed <- seed_base + match(scenario, names(pr_scenarios())) * 10000L + replicate
      data <- pr_generate(scenarios[[scenario]], seed)
      fit <- pr_fit(data)
      latent_metrics <- pr_dataset_metrics(lapply(seq_along(data$truth), function(r)
        list(activity = data$truth[[r]], eta = data$eta_truth[[r]])), data)
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(scenario = scenario, replicate = replicate,
        seed = seed, trials_per_run = mean(vapply(data$trials, nrow, integer(1))),
        rho_bias = mean(vapply(fit$noise, function(x) mean(x$rho - data$config$noise_rho), numeric(1))),
        rho_mae = mean(vapply(fit$noise, function(x) mean(abs(x$rho - data$config$noise_rho)), numeric(1))),
        noncanonical_fraction = mean(vapply(fit$selections, function(x) mean(x$kernel_index != 1L), numeric(1))),
        hrf_nmse = mean(vapply(fit$selections, function(x) mean(colSums((fit$kernels[, x$kernel_index,
          drop = FALSE] - data$kernels)^2) / colSums(data$kernels^2)), numeric(1))))
      for (lambda in lambdas) for (method in pr_methods()) {
        measures <- if (method == "latent_truth") latent_metrics else
          pr_dataset_metrics(pr_recover(data, fit, method, lambda, cache), data)
        metrics[[length(metrics) + 1L]] <- data.frame(scenario = scenario, replicate = replicate,
          method = method, lambda = lambda, as.list(measures))
      }
      if (replicate %% 10L == 0L || replicate == counts[[scenario]]) {
        write.csv(do.call(rbind, metrics), file.path(out_dir, "replicate_metrics.csv"), row.names = FALSE)
        write.csv(do.call(rbind, diagnostics), file.path(out_dir, "fit_diagnostics.csv"), row.names = FALSE)
        message("  ", replicate, "; elapsed ", round(proc.time()[["elapsed"]] - start, 1), " s; ",
          length(ls(cache)), " cached operators")
      }
    }
  }
  metrics <- do.call(rbind, metrics)
  write.csv(pr_summarize(metrics), file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(pr_paired(metrics), file.path(out_dir, "paired.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()), paste("Elapsed seconds:", proc.time()[["elapsed"]] - start),
    paste("OPENBLAS_NUM_THREADS:", Sys.getenv("OPENBLAS_NUM_THREADS")),
    paste("OMP_NUM_THREADS:", Sys.getenv("OMP_NUM_THREADS")),
    paste("Source MD5:", names(settings$script_md5), settings$script_md5)), file.path(out_dir, "environment.txt"))
  invisible(metrics)
}

if (sys.nframe() == 0L) {
  .libPaths(c(Sys.getenv("FMRI_VALIDATION_LIBRARY", unset = .libPaths()[1L]), .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  run_parametric_recovery(if (length(args)) args[1L] else "/tmp/parametric-recovery",
    if (length(args) > 1L) args[2L] else "evaluation",
    if (length(args) > 2L) as.integer(args[3L]) else NULL,
    if (length(args) > 3L) as.integer(args[4L]) else NULL,
    if (length(args) > 4L) strsplit(args[5L], ",", fixed = TRUE)[[1L]] else names(pr_scenarios()))
}
