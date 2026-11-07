# Preserve the original illustrative comparison. New defaults and repeated
# validation are documented in NOISE_AWARE_HRF_RECOVERY.md.
hrfe_estimate_legacy <- function(...) {
  fmri.pipeline::estimate_crossrun_hrf(..., noise_model = "ols", evidence = "legacy")
}

# Reproducible synthetic review; no participant data or observed RPE optimization.
# Run with the installed PR version:
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
#   inst/examples/crossrun_hrf_review.R /tmp/crossrun-hrf-review
library(fmri.pipeline)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else tempfile("crossrun-hrf-review-")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(7206)
tr <- 1; n <- 260L; n_voxels <- 40L
kernels <- hrf_library(tr, peak_delays = c(5, 6, 7), dispersions = 1)
truth <- rep(c(2L, 3L), each = n_voxels / 2)
voxel_ids <- paste0("v", seq_len(n_voxels))
parcels <- stats::setNames(rep(c("early", "late"), each = n_voxels / 2), voxel_ids)
events <- trials <- bold <- neural_inputs <- list()
for (r in 1:4) {
  run <- paste0("run", r)
  cue <- seq(5, 225, by = 10) + runif(23, -.3, .3)
  rpe <- as.numeric(scale(rnorm(length(cue))))
  trials[[run]] <- data.frame(run = run, trial = seq_along(cue), feedback = cue + 4, rpe = rpe)
  ev <- rbind(data.frame(onset = cue, condition = "cue", amplitude = .3),
    data.frame(onset = cue + 2, condition = "response", amplitude = -.2),
    data.frame(onset = cue + 4, condition = "feedback", amplitude = 1 + .4 * rpe))
  events[[run]] <- ev[c("onset", "condition")]
  # Independent generating convolution, with known graded feedback amplitudes.
  neural <- numeric(n)
  for (e in seq_len(nrow(ev))) {
    position <- ev$onset[e] / tr
    lo <- floor(position); f <- position - lo
    neural[lo + 1L] <- neural[lo + 1L] + ev$amplitude[e] * (1 - f)
    if (f > 0) neural[lo + 2L] <- neural[lo + 2L] + ev$amplitude[e] * f
  }
  neural_inputs[[run]] <- neural
  y <- vapply(truth, function(k) {
    as.numeric(convolve(neural, rev(kernels[, k]), type = "open"))[seq_len(n)] +
      100 + seq(-.2, .2, length.out = n) + rnorm(n, sd = .015)
  }, numeric(n))
  colnames(y) <- voxel_ids
  bold[[run]] <- y
}
parcel_fit <- hrfe_estimate_legacy(bold, events, tr, kernels,
  spatial_groups = parcels, pooling = "parcel", return_scores = TRUE)
voxel_fit <- hrfe_estimate_legacy(bold, events, tr, kernels,
  spatial_groups = parcels, spatial_strength = .02, return_scores = TRUE,
  amplitude_model = "trial")
# Prespecified lambda and strength: neither is optimized against the RPE result.
recovered <- deconvolve_crossrun_hrf(bold, voxel_fit, lambda = .01, scaling = "demean")
canonical_fit <- voxel_fit
for (r in names(bold)) canonical_fit$selections[[r]]$kernel_index[] <- 1L
canonical <- deconvolve_crossrun_hrf(bold, canonical_fit, lambda = .01, scaling = "demean")
selection <- do.call(rbind, voxel_fit$selections)
selection$truth <- rep(truth, length(bold))
write.csv(selection, file.path(out_dir, "voxel_selections.csv"), row.names = FALSE)
write.csv(do.call(rbind, parcel_fit$parcels), file.path(out_dir, "parcel_anchors.csv"), row.names = FALSE)
write.csv(voxel_fit$designs, file.path(out_dir, "design_manifest.csv"), row.names = FALSE)
write.csv(voxel_fit$signal_status, file.path(out_dir, "training_signal_status.csv"), row.names = FALSE)
write.csv(data.frame(time = (seq_len(nrow(kernels)) - 1) * tr, kernels),
  file.path(out_dir, "hrf_library.csv"), row.names = FALSE)
summary <- list(); epochs <- list()
for (r in names(bold)) for (g in unique(parcels)) for (method in c("crossrun", "canonical")) {
  y <- if (method == "crossrun") recovered$activity[[r]] else canonical$activity[[r]]
  region <- rowMeans(y[, parcels == g, drop = FALSE])
  at_onset <- approx((seq_len(n) - 1) * tr, region, xout = trials[[r]]$feedback)$y
  association <- lm(at_onset ~ trials[[r]]$rpe)
  summary[[length(summary) + 1L]] <- data.frame(run = r, parcel = g, method = method,
    rpe_slope = unname(coef(association)[2]), rpe_correlation = cor(at_onset, trials[[r]]$rpe),
    correct_hrf_fraction = mean(selection$kernel_index[selection$target_run == r & selection$parcel == g] ==
      selection$truth[selection$target_run == r & selection$parcel == g]))
  for (i in seq_len(nrow(trials[[r]]))) {
    lags <- seq(-3, 6, .5)
    epochs[[length(epochs) + 1L]] <- data.frame(run = r, parcel = g, method = method,
      trial = i, lag = lags, rpe = trials[[r]]$rpe[i],
      activity = approx((seq_len(n) - 1) * tr, region,
        xout = trials[[r]]$feedback[i] + lags)$y)
  }
}
summary <- do.call(rbind, summary)
write.csv(summary, file.path(out_dir, "graded_recovery.csv"), row.names = FALSE)
write.csv(do.call(rbind, epochs), file.path(out_dir, "aligned_simulated_epochs.csv"), row.names = FALSE)
saveRDS(list(hrf_fit = voxel_fit, parcel_fit = parcel_fit, inverse = recovered$settings,
  generating_kernel_index = truth, seed = 7206), file.path(out_dir, "frozen_models.rds"))
png(file.path(out_dir, "crossrun_hrf_review.png"), width = 1300, height = 650, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
matplot((seq_len(nrow(kernels)) - 1) * tr, kernels, type = "l", lty = 1,
  xlab = "Lag (seconds)", ylab = "Unit-sum HRF", main = "Prespecified HRF library")
legend("topright", legend = c("canonical", "early", "late"), col = 1:3, lty = 1, bty = "n")
boxplot(rpe_slope ~ interaction(method, parcel), summary, las = 1, ylim = c(0, .42),
  names = c("Canonical\nearly", "Cross-run\nearly", "Canonical\nlate", "Cross-run\nlate"),
  xlab = "", ylab = "Recovered onset slope versus RPE", main = "Known injected slope = 0.4")
abline(h = .4, lty = 2, col = "gray40")
dev.off()
print(summary, row.names = FALSE)

# Stress the assumptions without selecting new settings: generating HRFs lie
# outside the candidate library, noise is AR(1), and one case mixes HRFs inside
# the supplied parcels. This is descriptive sensitivity analysis, not a claim
# of universal improvement or a full parameter-selection evaluation.
generating <- hrf_library(tr, peak_delays = c(5.5, 7.5), dispersions = 1,
  undershoot_ratios = 8)
stress_models <- stress_rows <- list()
for (scenario in c("off_library_ar1", "heterogeneous_parcels_ar1")) {
  generating_index <- rep(1:2, each = n_voxels / 2)
  if (scenario == "heterogeneous_parcels_ar1") generating_index[c(5:9, 25:29)] <-
    3L - generating_index[c(5:9, 25:29)]
  stress_bold <- lapply(names(bold), function(r) {
    y <- vapply(seq_len(n_voxels), function(v) {
      noise <- as.numeric(stats::filter(rnorm(n + 100, sd = .015 * sqrt(1 - .6^2)),
        filter = .6, method = "recursive"))[-seq_len(100)]
      as.numeric(convolve(neural_inputs[[r]], rev(generating[, generating_index[v]]),
        type = "open"))[seq_len(n)] + 100 + seq(-.2, .2, length.out = n) + noise
    }, numeric(n))
    colnames(y) <- voxel_ids
    y
  })
  names(stress_bold) <- names(bold)
  for (method in c("parcel", "voxel_spatial")) {
    model <- hrfe_estimate_legacy(stress_bold, events, tr, kernels,
      spatial_groups = parcels, pooling = if (method == "parcel") "parcel" else "voxel",
      spatial_strength = if (method == "parcel") 0 else .02, amplitude_model = "trial")
    estimate <- deconvolve_crossrun_hrf(stress_bold, model, .01, scaling = "demean")
    stress_models[[paste(scenario, method, sep = ":")]] <- model
    for (r in names(bold)) for (g in unique(parcels)) {
      voxels <- which(parcels == g)
      region <- rowMeans(estimate$activity[[r]][, voxels, drop = FALSE])
      response <- approx(0:(n - 1), region, xout = trials[[r]]$feedback)$y
      indices <- model$selections[[r]]$kernel_index[voxels]
      shape_error <- mean((kernels[, indices] - generating[, generating_index[voxels]])^2)
      stress_rows[[length(stress_rows) + 1L]] <- data.frame(scenario = scenario, method = method,
        run = r, parcel = g, shape_mse = shape_error,
        rpe_slope = unname(coef(lm(response ~ trials[[r]]$rpe))[2]),
        rpe_correlation = cor(response, trials[[r]]$rpe),
        prespecified_fallback_fraction = mean(model$selections[[r]]$status[voxels] == "fallback_prespecified"))
    }
  }
}
write.csv(do.call(rbind, stress_rows), file.path(out_dir, "sensitivity_recovery.csv"), row.names = FALSE)
saveRDS(list(models = stress_models, generating_kernels = generating),
  file.path(out_dir, "sensitivity_models.rds"))
cat("Artifacts:", normalizePath(out_dir), "\n")
