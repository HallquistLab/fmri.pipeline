# Compare the current source package with a Git revision in this checkout.
# Run from the repository root, with one BLAS/OpenMP thread, for example:
# OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla \
#   inst/benchmarks/deconvolution.R bfbffcb 5
# Compilation/loading are excluded. Times include preprocessing and construction.
args <- commandArgs(trailingOnly = TRUE)
revision <- if (length(args)) args[1L] else "bfbffcb"
repetitions <- if (length(args) > 1L) as.integer(args[2L]) else 5L
stopifnot(length(repetitions) == 1L, is.finite(repetitions), repetitions >= 1L,
  grepl("^[A-Za-z0-9_][A-Za-z0-9_./-]*$", revision))
pkgload::load_all(".", quiet = TRUE)
original <- new.env(parent = asNamespace("fmri.pipeline"))
source_dir <- tempfile("deconvolution_baseline_")
dir.create(source_dir)
for (path in c("R/deconvolve_reglin.R", "src/deconvolve_reglin.cpp", "src/fmri.pipeline.h")) {
  status <- system2("git", c("show", shQuote(paste0(revision, ":", path))),
    stdout = file.path(source_dir, basename(path)))
  if (status != 0L) stop("Cannot read baseline source: ", path)
}
sys.source(file.path(source_dir, "deconvolve_reglin.R"), original)
Rcpp::sourceCpp(file.path(source_dir, "deconvolve_reglin.cpp"), env = original)
set.seed(6106)
kernel <- fmri.pipeline:::spm_hrf(1)$hrf
results <- list()
for (n_time in c(300L, 1200L)) {
  y <- matrix(rnorm(n_time * 200L), nrow = n_time)
  for (backend in c("cpp", "wiener")) {
    for (tuned in c(FALSE, TRUE)) {
      settings <- list(BOLDobs = y, kernel = kernel, backend = backend)
      settings <- c(settings, if (tuned) {
        list(lambda_grid = 10^seq(-4, 2, length.out = 10L))
      } else list(lambda = 0.01))
      before <- do.call(original$deconvolve_reglin, settings)
      after <- do.call(deconvolve_reglin, settings)
      timings <- vapply(seq_len(repetitions), function(i) {
        # Alternate ordering to reduce systematic warm-cache/order effects.
        run_before <- function() unname(system.time(do.call(original$deconvolve_reglin, settings))["elapsed"])
        run_after <- function() unname(system.time(do.call(deconvolve_reglin, settings))["elapsed"])
        if (i %% 2L) c(before = run_before(), after = run_after()) else {
          current <- run_after()
          c(before = run_before(), after = current)
        }
      }, numeric(2L))
      row <- data.frame(timepoints = n_time, signals = ncol(y), backend = backend,
        tuned = tuned, before = median(timings[1L, ]), after = median(timings[2L, ]),
        max_abs_difference = max(abs(after - before)))
      row$speedup <- row$before / row$after
      results[[length(results) + 1L]] <- row
    }
  }
}
print(do.call(rbind, results), row.names = FALSE)
# Wiener tuning intentionally changes when corrected GCV selects another lambda;
# fixed-parameter estimates should agree to numerical precision.
unlink(source_dir, recursive = TRUE)
