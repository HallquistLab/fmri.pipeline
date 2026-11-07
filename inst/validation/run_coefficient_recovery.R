# Run from the repository root with the PR #46 native package installed.
source(file.path("inst", "validation", "coefficient_recovery.R"), local = TRUE)

run_coefficient_recovery <- function(out_dir, phase = c("evaluation", "development"),
    repetitions = NULL, null_repetitions = NULL, scenario_names = names(cr_scenarios())) {
  phase <- match.arg(phase)
  if (is.null(repetitions)) repetitions <- if (phase == "evaluation") 40L else 3L
  if (is.null(null_repetitions)) null_repetitions <- if (phase == "evaluation") 200L else repetitions
  stopifnot(repetitions >= 1L, null_repetitions >= 1L, repetitions == as.integer(repetitions),
    null_repetitions == as.integer(null_repetitions), max(repetitions, null_repetitions) < 10000L)
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) stop("Use an empty output directory")
  scenarios <- cr_scenarios()[scenario_names]
  stopifnot(length(scenarios) > 0L, !any(vapply(scenarios, is.null, logical(1))))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  counts <- vapply(scenarios, function(x) if (x$beta == 0) null_repetitions else repetitions, integer(1))
  schedule <- cr_schedule()
  files <- file.path("inst", "validation", c("parametric_recovery.R", "run_parametric_recovery.R",
    "coefficient_recovery.R", "run_coefficient_recovery.R"))
  seed_base <- if (phase == "evaluation") 22000000L else 20000000L
  settings <- list(phase = phase, scenarios = scenarios, counts = counts, schedule = schedule,
    primary_lambda = .01, seed_base = seed_base,
    seed_formula = "base + original_scenario_index * 10000 + replicate",
    script_md5 = tools::md5sum(files),
    primary_metrics = c("beta_bias", "beta_rmse", "coefficient_coverage", "coefficient_false_positive", "amplitude_r"),
    interval = "t over six run coefficients, df=5; empirical coverage audit",
    readouts = list(onset = 0, area = -1:3, template = -3:5),
    regression = "ordinary adjusted feedback design versus full transformed event-source design",
    contrast_rank_tolerance = 1e-8, minimum_relative_target_information = 1e-10,
    rng_kind = RNGkind(), package = utils::packageDescription("fmri.pipeline")[c("Package", "Version", "Built")])
  saveRDS(settings, file.path(out_dir, "settings.rds"))
  metrics <- diagnostics <- list(); cache <- new.env(parent = emptyenv())
  start <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    message(scenario, ": ", counts[[scenario]], " independent datasets")
    for (replicate in seq_len(counts[[scenario]])) {
      seed <- seed_base + match(scenario, names(cr_scenarios())) * 10000L + replicate
      data <- cr_generate(scenarios[[scenario]], seed); fit <- pr_fit(data)
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(scenario = scenario, replicate = replicate, seed = seed,
        hrf_nmse = mean(vapply(fit$selections, function(x)
          mean(colSums((fit$kernels[, x$kernel_index, drop = FALSE] - data$kernels)^2) /
            colSums(data$kernels^2)), numeric(1))))
      for (source in unique(schedule$source)) for (lambda in unique(schedule$lambda[schedule$source == source])) {
        runs <- lapply(seq_along(data$bold), function(r) cr_run(data, fit, r, source, lambda, cache))
        selected <- schedule[schedule$source == source & schedule$lambda == lambda, ]
        for (i in seq_len(nrow(selected))) {
          value <- if (selected$readout[i] == "forward") cr_forward(data, fit, source) else
            cr_metrics(runs, data, selected$readout[i], selected$estimator[i])
          metrics[[length(metrics) + 1L]] <- data.frame(scenario = scenario, replicate = replicate,
            method = selected$method[i], lambda = lambda, as.list(value))
        }
      }
      if (replicate %% 10L == 0L || replicate == counts[[scenario]]) {
        utils::write.csv(do.call(rbind, metrics), file.path(out_dir, "replicate_metrics.csv"), row.names = FALSE)
        utils::write.csv(do.call(rbind, diagnostics), file.path(out_dir, "fit_diagnostics.csv"), row.names = FALSE)
        message("  ", replicate, "; elapsed ", round(proc.time()[["elapsed"]] - start, 1), " s")
      }
    }
  }
  writeLines(c(capture.output(sessionInfo()), paste("Elapsed seconds:", proc.time()[["elapsed"]] - start),
    paste("OPENBLAS_NUM_THREADS", Sys.getenv("OPENBLAS_NUM_THREADS")),
    paste("OMP_NUM_THREADS", Sys.getenv("OMP_NUM_THREADS"))), file.path(out_dir, "environment.txt"))
  invisible(do.call(rbind, metrics))
}

if (sys.nframe() == 0L) {
  .libPaths(c(Sys.getenv("FMRI_VALIDATION_LIBRARY", unset = .libPaths()[1L]), .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  run_coefficient_recovery(if (length(args)) args[1L] else "/tmp/coefficient-recovery",
    if (length(args) > 1L) args[2L] else "evaluation",
    if (length(args) > 2L) as.integer(args[3L]) else NULL,
    if (length(args) > 3L) as.integer(args[4L]) else NULL,
    if (length(args) > 4L) strsplit(args[5L], ",", fixed = TRUE)[[1L]] else names(cr_scenarios()))
}
