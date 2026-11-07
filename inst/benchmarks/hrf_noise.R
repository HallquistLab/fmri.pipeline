# Install this checkout first. Excludes package loading and data generation.
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
#   inst/benchmarks/hrf_noise.R /tmp/hrf-noise-benchmark 50000 3
.libPaths(c(Sys.getenv("FMRI_VALIDATION_LIBRARY", unset = .libPaths()[1L]), .libPaths()))
library(fmri.pipeline)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1L] else "/tmp/hrf-noise-benchmark"
voxels <- if (length(args) > 1L) as.integer(args[2L]) else 50000L
repetitions <- if (length(args) > 2L) as.integer(args[3L]) else 3L
modern_only <- length(args) > 3L && args[4L] == "modern"
stopifnot(voxels >= 20L, repetitions >= 1L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(81920)
n <- 300L; h <- hrf_library(1)
ids <- paste0("v", seq_len(voxels))
groups <- setNames(rep(paste0("p", 1:20), length.out = voxels), ids)
events <- setNames(lapply(1:3, function(r) data.frame(onset = seq(5, 280, by = 5) + runif(56),
  condition = rep(c("cue", "feedback"), 28))), paste0("run", 1:3))
bold <- lapply(events, function(ev) {
  designs <- fmri.pipeline:::hrf_event_design_cpp(n, 1, 0, ev$onset, rep(0, nrow(ev)),
    1:nrow(ev), nrow(ev), h)
  y <- matrix(rnorm(n * voxels, sd = .03 * sqrt(1 - .5^2)), n, voxels)
  y[1, ] <- y[1, ] / sqrt(1 - .5^2)
  for (t in 2:n) y[t, ] <- y[t, ] + .5 * y[t - 1L, ]
  truth <- rep(seq_len(ncol(h)), length.out = voxels)
  for (k in seq_len(ncol(h))) {
    index <- which(truth == k)
    y[, index] <- y[, index] + as.numeric(designs[[k]] %*% rep(c(1, -.5, 2, .5), 14))
  }
  colnames(y) <- ids; y + 100
})
args_fit <- list(bold = bold, events = events, tr = 1, kernel_grid = h,
  spatial_groups = groups, spatial_strength = .02, amplitude_model = "trial")
median_time <- function(fun) median(replicate(repetitions, system.time(fun())[["elapsed"]]))
variants <- list(legacy_ols = list(noise_model = "ols", evidence = "legacy"),
  calibrated_ols = list(noise_model = "ols"),
  ar1_heuristic = list(evidence = "legacy"), ar1_calibrated = list())
if (modern_only) variants <- variants["ar1_calibrated"]
rows <- lapply(names(variants), function(name) {
  message("Timing ", name)
  seconds <- median_time(function() do.call(estimate_crossrun_hrf, c(args_fit, variants[[name]])))
  message("Completed ", name, ": ", seconds, " s")
  data.frame(operation = paste0("fit_", name), seconds = seconds)
})
fit <- do.call(estimate_crossrun_hrf, args_fit)
message("Timing grouped inverse")
rows[[length(rows) + 1L]] <- data.frame(operation = "apply_ar1_calibrated",
  seconds = median_time(function() deconvolve_crossrun_hrf(bold, fit, .01, scaling = "demean")))
result <- deconvolve_crossrun_hrf(bold, fit, .01, scaling = "demean")
message("Timing diagnostics")
rows[[length(rows) + 1L]] <- data.frame(operation = "diagnose_onset_and_leakage",
  seconds = median_time(function() diagnose_crossrun_recovery(result)))
timings <- do.call(rbind, rows)
timings$voxels <- voxels; timings$runs <- 3L; timings$timepoints <- n
timings$candidates <- ncol(h); timings$repetitions <- repetitions
write.csv(timings, file.path(out_dir, "timings.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()),
  paste("OPENBLAS_NUM_THREADS", Sys.getenv("OPENBLAS_NUM_THREADS")),
  paste("OMP_NUM_THREADS", Sys.getenv("OMP_NUM_THREADS"))), file.path(out_dir, "environment.txt"))
print(timings, row.names = FALSE)
