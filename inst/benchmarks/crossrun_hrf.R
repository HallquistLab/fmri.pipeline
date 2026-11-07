# Run against an installed version of this checkout, excluding installation:
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
#   inst/benchmarks/crossrun_hrf.R /tmp/hrf-benchmark 10000 3
library(fmri.pipeline)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else tempfile("hrf-benchmark-")
n_voxels <- if (length(args) > 1) as.integer(args[2]) else 10000L
repetitions <- if (length(args) > 2) as.integer(args[3]) else 3L
stopifnot(is.finite(n_voxels), n_voxels >= 10, is.finite(repetitions), repetitions >= 1)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(716)
n <- 300L
kernel <- hrf_library(1)
events <- lapply(1:3, function(r) data.frame(onset = seq(5, 280, by = 5) + runif(56),
  condition = rep(c("cue", "response", "feedback", "cue"), 14)))
names(events) <- paste0("run", 1:3)
designs <- lapply(events, function(ev) fmri.pipeline:::hrf_event_design_cpp(n, 1, 0,
  ev$onset, rep(0, nrow(ev)), match(ev$condition, unique(ev$condition)), 3L, kernel))
groups <- stats::setNames(rep(paste0("p", 1:20), length.out = n_voxels), paste0("v", seq_len(n_voxels)))
truth <- rep(seq_len(ncol(kernel)), length.out = n_voxels)
bold <- lapply(designs, function(design) {
  y <- matrix(rnorm(n * n_voxels, sd = .03), n, n_voxels)
  for (k in seq_len(ncol(kernel))) {
    columns <- which(truth == k)
    y[, columns] <- y[, columns] + as.numeric(design[[k]] %*% c(1, -.5, 2))
  }
  colnames(y) <- names(groups)
  y
})
z <- cbind(1, seq(-1, 1, length.out = n))
reference_voxels <- seq_len(min(300L, n_voxels))
reference_y <- bold[[1]][, reference_voxels, drop = FALSE]
# Deliberately conventional per-voxel lm.fit baseline, not an optimized R
# projection implementation. Both solve the same unregularized training models.
reference <- function() {
  out <- matrix(NA_real_, ncol(kernel), length(reference_voxels))
  for (v in seq_along(reference_voxels)) {
    base <- lm.fit(z, reference_y[, v])
    for (k in seq_len(ncol(kernel))) {
      full <- lm.fit(cbind(z, designs[[1]][[k]]), reference_y[, v])
      out[k, v] <- 1 - (sum(full$residuals^2) / (n - full$rank)) /
        (sum(base$residuals^2) / (n - base$rank))
    }
  }
  out
}
native <- function() fmri.pipeline:::hrf_score_run_cpp(reference_y, designs[[1]], z, seq_len(n))$scores
expected <- reference(); observed <- native()
stopifnot(max(abs(expected - observed)) < 1e-9)
median_time <- function(fun) median(replicate(repetitions, system.time(fun())[["elapsed"]]))
rows <- list(
  data.frame(operation = "per_voxel_lm_fit_one_run", voxels = length(reference_voxels), seconds = median_time(reference)),
  data.frame(operation = "native_batched_one_run", voxels = length(reference_voxels), seconds = median_time(native)))
fit <- estimate_crossrun_hrf(bold, events, 1, kernel, spatial_groups = groups,
  spatial_strength = .02, batch_size = 2048)
rows[[3]] <- data.frame(operation = "all_three_crossrun_selections", voxels = n_voxels,
  seconds = median_time(function() estimate_crossrun_hrf(bold, events, 1, kernel,
    spatial_groups = groups, spatial_strength = .02, batch_size = 2048)))
rows[[4]] <- data.frame(operation = "all_three_frozen_deconvolutions", voxels = n_voxels,
  seconds = median_time(function() deconvolve_crossrun_hrf(bold, fit, lambda = .01)))
rows[[5]] <- data.frame(operation = "all_three_frozen_wiener", voxels = n_voxels,
  seconds = median_time(function() deconvolve_crossrun_hrf(bold, fit, lambda = .01, backend = "wiener")))
rows[[6]] <- data.frame(operation = "all_three_trial_amplitude_selections", voxels = n_voxels,
  seconds = median_time(function() estimate_crossrun_hrf(bold, events, 1, kernel,
    spatial_groups = groups, spatial_strength = .02, amplitude_model = "trial")))
timings <- do.call(rbind, rows)
timings$timepoints <- n; timings$runs <- c(1L, 1L, rep(3L, length(rows) - 2L))
timings$candidates <- ncol(kernel); timings$repetitions <- repetitions
write.csv(timings, file.path(out_dir, "timings.csv"), row.names = FALSE)
writeLines(c(paste("max_abs_score_difference", max(abs(expected - observed))),
  capture.output(sessionInfo()), paste("OPENBLAS_NUM_THREADS", Sys.getenv("OPENBLAS_NUM_THREADS")),
  paste("OMP_NUM_THREADS", Sys.getenv("OMP_NUM_THREADS"))), file.path(out_dir, "environment.txt"))
print(timings, row.names = FALSE)
cat("Maximum absolute score difference:", max(abs(expected - observed)), "\n")
