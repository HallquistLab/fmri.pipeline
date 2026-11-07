# Compare frozen legacy recovery, covariance estimation, and calibrated selection.
# Install this checkout first. One BLAS/OpenMP thread is recommended.
# Rscript --vanilla inst/validation/hrf_noise_recovery.R /tmp/hrf-noise 40
# Optional third argument is a comma-separated scenario subset; fourth is seed base.
source(file.path("inst", "validation", "hrf_recovery.R"))

run_hrf_noise_validation <- function(out_dir, repetitions = 40L,
    scenario_names = names(hrfv_scenarios()), seed_base = 9420000L) {
  stopifnot(repetitions >= 1L, repetitions == as.integer(repetitions))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- hrfv_scenarios()[scenario_names]
  stopifnot(!any(vapply(scenarios, is.null, logical(1))))
  kernels <- fmri.pipeline::hrf_library(1)
  lambdas <- c(.001, .01, .1)
  settings <- list(repetitions = repetitions, scenarios = scenarios,
    seed_base = seed_base, seed_formula = "seed_base + original_scenario_index * 10000 + replicate",
    lambdas = lambdas, primary_lambda = .01, amplitude_model = "trial",
    spatial_strength = .02, canonical_strength = .01, evidence_alpha = .05,
    nuisance_mode = "ignore (prepared inverse inputs); joint in HRF training",
    script_md5 = unname(tools::md5sum(c("inst/validation/hrf_noise_recovery.R",
      "inst/validation/hrf_recovery.R"))))
  # Freeze settings BEFORE generating any datasets or inspecting results.
  saveRDS(settings, file.path(out_dir, "settings.rds"))
  all_metrics <- all_runs <- list(); started <- proc.time()[["elapsed"]]
  for (scenario in names(scenarios)) {
    message(scenario, ": ", repetitions, " independent datasets")
    for (replicate in seq_len(repetitions)) {
      data <- hrfv_generate(scenarios[[scenario]], seed_base +
        match(scenario, names(hrfv_scenarios())) * 10000L + replicate)
      args <- list(bold = data$bold, events = data$events, tr = 1, kernel_grid = kernels,
        valid_rows = data$masks, spatial_groups = data$groups, spatial_strength = .02,
        amplitude_model = "trial")
      modern <- do.call(fmri.pipeline::estimate_crossrun_hrf, args)
      legacy <- do.call(fmri.pipeline::estimate_crossrun_hrf,
        c(args, list(noise_model = "ols", evidence = "legacy")))
      heuristic <- do.call(fmri.pipeline::estimate_crossrun_hrf,
        c(args, list(evidence = "legacy")))
      canonical <- legacy
      canonical_gls <- modern
      for (r in seq_along(data$bold)) {
        canonical$selections[[r]]$kernel_index[] <- 1L
        canonical_gls$selections[[r]]$kernel_index[] <- 1L
      }
      fits <- list(canonical = canonical, canonical_gls = canonical_gls,
        legacy = legacy, gls_heuristic = heuristic, gls_calibrated = modern)
      for (lambda in lambdas) for (method in names(fits)) {
        fit <- fits[[method]]
        result <- fmri.pipeline::deconvolve_crossrun_hrf(data$bold, fit, lambda,
          scaling = "demean", backend = "cpp")
        rows <- lapply(seq_along(data$bold), function(r) {
          selection <- fit$selections[[r]]
          h <- fit$kernels[, selection$kernel_index, drop = FALSE]
          metrics <- hrfv_metrics(result$activity[[r]], data$truth[[r]], data$trials[[r]], h, data$kernels[[r]])
          learned <- method %in% c("legacy", "gls_heuristic", "gls_calibrated")
          extra <- c(noncanonical_fraction = mean(selection$kernel_index != 1L),
            fallback_fraction = if (learned) mean(grepl("^fallback", selection$status)) else NA_real_,
            unidentified_fraction = if (method == "gls_calibrated")
              mean(selection$voxel_evidence %in% c("unidentified", "insufficient_training")) else NA_real_,
            incompatible_fraction = if (method == "gls_calibrated")
              mean(selection$n_compatible_shapes == 0L) else NA_real_,
            task_support_fraction = if (method == "gls_calibrated")
              mean(selection$selected_task_p <= fit$settings$evidence_alpha) else NA_real_,
            noise_rho = mean(fit$noise[[r]]$rho),
            noise_rho_error = mean(abs(fit$noise[[r]]$rho - scenarios[[scenario]]$rho)))
          data.frame(scenario = scenario, replicate = replicate, method = method,
            lambda = lambda, run = names(data$bold)[r], as.list(c(metrics, extra)))
        })
        runs <- do.call(rbind, rows)
        all_runs[[length(all_runs) + 1L]] <- runs
        measures <- setdiff(names(runs), c("scenario", "replicate", "method", "lambda", "run"))
        all_metrics[[length(all_metrics) + 1L]] <- data.frame(scenario = scenario,
          replicate = replicate, method = method, lambda = lambda,
          as.list(vapply(runs[measures], hrfv_mean, numeric(1))))
      }
      if (replicate %% 10L == 0L || replicate == repetitions) {
        write.csv(do.call(rbind, all_metrics), file.path(out_dir, "replicate_metrics.csv"), row.names = FALSE)
        message("  ", replicate, "; elapsed ", round(proc.time()[["elapsed"]] - started, 1), " s")
      }
    }
  }
  metrics <- do.call(rbind, all_metrics)
  tables <- hrfv_summary(metrics)
  # A matched-noise baseline isolates HRF selection from whitening and nuisance.
  matched <- metrics[metrics$method %in% c("canonical_gls", "gls_heuristic", "gls_calibrated"), ]
  matched$method[matched$method == "canonical_gls"] <- "canonical"
  paired <- hrfv_summary(matched)$paired
  write.csv(do.call(rbind, all_runs), file.path(out_dir, "run_metrics.csv"), row.names = FALSE)
  write.csv(tables$summary, file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(tables$paired, file.path(out_dir, "paired_vs_legacy_canonical.csv"), row.names = FALSE)
  write.csv(paired, file.path(out_dir, "paired_vs_matched_canonical.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()),
    paste("Elapsed seconds:", proc.time()[["elapsed"]] - started),
    paste("OPENBLAS_NUM_THREADS", Sys.getenv("OPENBLAS_NUM_THREADS")),
    paste("OMP_NUM_THREADS", Sys.getenv("OMP_NUM_THREADS"))), file.path(out_dir, "environment.txt"))
  invisible(list(summary = tables$summary, paired = paired))
}

if (sys.nframe() == 0L) {
  .libPaths(c(Sys.getenv("FMRI_VALIDATION_LIBRARY", unset = .libPaths()[1L]), .libPaths()))
  args <- commandArgs(trailingOnly = TRUE)
  run_hrf_noise_validation(if (length(args)) args[1] else "/tmp/hrf-noise",
    if (length(args) > 1L) as.integer(args[2]) else 40L,
    if (length(args) > 2L) strsplit(args[3], ",", fixed = TRUE)[[1L]] else names(hrfv_scenarios()),
    if (length(args) > 3L) as.integer(args[4]) else 9420000L)
}
