# Rscript --vanilla inst/validation/report_parametric_recovery.R INPUT_DIR OUTPUT_DIR
# INPUT_DIR may contain one complete run, or disjoint scenario subdirectories.
source(file.path("inst", "validation", "run_parametric_recovery.R"), local = TRUE)

pr_collect <- function(input_dir) {
  paths <- if (file.exists(file.path(input_dir, "settings.rds"))) file.path(input_dir, "settings.rds") else
    list.files(input_dir, pattern = "^settings[.]rds$", recursive = TRUE, full.names = TRUE)
  if (!length(paths)) stop("No frozen study settings found")
  settings <- lapply(paths, readRDS)
  invariant <- c("phase", "lambdas", "primary_lambda", "methods", "contrasts", "primary_metrics",
    "seed_base", "seed_formula", "script_md5", "inference", "gain_convention", "package_description")
  for (s in settings) {
    if (!identical(s[invariant], settings[[1L]][invariant])) stop("Incompatible study protocols")
    if (!identical(s$script_md5, tools::md5sum(names(s$script_md5)))) stop("Study sources changed since freezing")
  }
  all_scenarios <- unlist(lapply(settings, function(s) names(s$scenarios)), use.names = FALSE)
  if (anyDuplicated(all_scenarios)) stop("Scenario shards overlap")
  for (i in seq_along(settings)) {
    s <- settings[[i]]; path <- dirname(paths[i])
    d <- utils::read.csv(file.path(path, "replicate_metrics.csv"))
    expected <- expand.grid(scenario = names(s$scenarios), method = s$methods, lambda = s$lambdas,
      stringsAsFactors = FALSE)
    expected <- do.call(rbind, lapply(seq_len(nrow(expected)), function(j) {
      count <- s$counts[[expected$scenario[j]]]
      x <- expected[rep(j, count), , drop = FALSE]; x$replicate <- seq_len(count)
      x
    }))
    key <- function(x) paste(x$scenario, x$replicate, x$method, sprintf("%.10g", x$lambda), sep = ":")
    if (anyDuplicated(key(d)) || !setequal(key(d), key(expected))) stop("Incomplete or duplicate dataset results")
    diag <- utils::read.csv(file.path(path, "fit_diagnostics.csv"))
    original_index <- match(diag$scenario, names(pr_scenarios()))
    if (anyNA(original_index) || any(diag$seed != s$seed_base + 10000L * original_index + diag$replicate))
      stop("Incorrect dataset seed provenance")
    expected_datasets <- unique(expected[c("scenario", "replicate")])
    dataset_key <- function(x) paste(x$scenario, x$replicate, sep = ":")
    if (anyDuplicated(diag$seed) || anyDuplicated(dataset_key(diag)) ||
        !setequal(dataset_key(diag), dataset_key(expected_datasets))) stop("Incomplete dataset diagnostics")
  }
  list(settings = settings, paths = paths,
    metrics = do.call(rbind, lapply(dirname(paths), function(p) utils::read.csv(file.path(p, "replicate_metrics.csv")))),
    diagnostics = do.call(rbind, lapply(dirname(paths), function(p) utils::read.csv(file.path(p, "fit_diagnostics.csv")))))
}

pr_write_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con))
  utils::write.csv(x, con, row.names = FALSE)
}

pr_write_provenance <- function(paths, out_dir) {
  lines <- c(paste("Report source MD5:", unname(tools::md5sum(
    file.path("inst", "validation", "report_parametric_recovery.R")))),
    unlist(lapply(paths, function(path) c(paste("Shard:", dirname(path)),
      readLines(file.path(dirname(path), "environment.txt")), ""))))
  lines <- sub("[ \\t]+$", "", lines)
  while (length(lines) && !nzchar(tail(lines, 1L))) lines <- head(lines, -1L)
  writeLines(lines, file.path(out_dir, "provenance.txt"))
}

pr_write_checksums <- function(out_dir) {
  files <- list.files(out_dir, full.names = TRUE)
  files <- files[basename(files) != "checksums.md5"]
  checksums <- tools::md5sum(files)
  writeLines(paste(unname(checksums), basename(names(checksums)), sep = "  "), file.path(out_dir, "checksums.md5"))
}

report_parametric_recovery <- function(input_dir, out_dir) {
  data <- pr_collect(input_dir)
  if (data$settings[[1L]]$phase != "evaluation" ||
      !setequal(unique(data$metrics$scenario), names(pr_scenarios())))
    stop("Final report requires every evaluation scenario")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  settings <- data$settings[[1L]]; primary_lambda <- settings$primary_lambda
  summary <- pr_summarize(data$metrics); paired <- pr_paired(data$metrics)
  saveRDS(data$settings, file.path(out_dir, "settings.rds"))
  saveRDS(data$metrics, file.path(out_dir, "replicate_metrics.rds"), compress = "xz")
  utils::write.csv(data$diagnostics, file.path(out_dir, "fit_diagnostics.csv"), row.names = FALSE)
  pr_write_gz(summary, file.path(out_dir, "summary.csv.gz"))
  pr_write_gz(paired, file.path(out_dir, "paired.csv.gz"))
  primary <- subset(summary, lambda == primary_lambda)
  wide <- stats::reshape(primary[c("scenario", "method", "lambda", "metric", "mean")],
    idvar = c("scenario", "method", "lambda"), timevar = "metric", direction = "wide")
  names(wide) <- sub("^mean[.]", "", names(wide))
  utils::write.csv(wide, file.path(out_dir, "primary.csv"), row.names = FALSE)
  utils::write.csv(primary, file.path(out_dir, "primary_intervals.csv"), row.names = FALSE)
  utils::write.csv(subset(paired, lambda == primary_lambda), file.path(out_dir, "paired_primary.csv"), row.names = FALSE)
  pr_write_provenance(data$paths, out_dir)

  methods <- c("canonical_estimated_noise", "learned_estimated_noise", "oracle_hrf_estimated_noise",
    "oracle_hrf_known_noise", "oracle_noiseless_continuous", "oracle_noiseless_discrete")
  labels <- c("Canonical / estimated noise", "Learned / estimated noise", "True HRF / estimated noise",
    "True HRF / true noise", "Continuous / noiseless / true HRF", "Discrete / noiseless / true HRF")
  colors <- c("#b35806", "#2166ac", "#762a83", "#1b7837", "#d73027", "#555555")
  grDevices::png(file.path(out_dir, "regularization.png"), width = 2000, height = 1800, res = 180)
  graphics::layout(matrix(c(1:9, 10, 10, 10), ncol = 3, byrow = TRUE), heights = c(1, 1, 1, .32))
  graphics::par(mar = c(3.6, 4.2, 2.1, .6), mgp = c(2.2, .6, 0))
  for (scenario in c("moderate_trial_variation", "dependent_trial_variation", "confounded_null")) {
    for (metric in c("amplitude_rmse", "adjusted_beta_rmse", "eta_observed_nrmse")) {
      d <- summary[summary$scenario == scenario & summary$metric == metric &
        summary$method %in% methods & summary$lambda > 0, ]
      graphics::plot(range(d$lambda), range(d$mean, finite = TRUE), type = "n", log = "xy",
        xlab = "lambda", ylab = metric, main = gsub("_", " ", scenario), cex.main = .78)
      graphics::abline(v = primary_lambda, lty = 3, col = "gray60")
      for (i in seq_along(methods)) {
        x <- d[d$method == methods[i], ]; x <- x[order(x$lambda), ]
        graphics::lines(x$lambda, x$mean, col = colors[i], lwd = 1.5, lty = i)
      }
    }
  }
  graphics::par(mar = rep(0, 4)); graphics::plot.new()
  graphics::legend("center", legend = labels, col = colors, lty = seq_along(methods),
    lwd = 1.5, ncol = 2, bty = "n", cex = .85)
  grDevices::dev.off()

  null_methods <- c("latent_truth", "canonical_estimated_noise", "learned_estimated_noise",
    "learned_frozen_known_noise", "oracle_hrf_estimated_noise", "oracle_hrf_known_noise", "oracle_noiseless_discrete")
  null_labels <- c("Latent truth", "Canonical", "Learned", "Learned / true noise",
    "True HRF", "True HRF / true rho", "Discrete / noiseless")
  grDevices::png(file.path(out_dir, "null_calibration.png"), width = 2000, height = 1400, res = 170)
  graphics::par(mfrow = c(2, 2), mar = c(8, 4, 2, 1), oma = c(0, 0, 2, 0))
  for (scenario in c("active_null", "confounded_null")) for (metric in c("false_positive", "coverage")) {
    d <- primary[primary$scenario == scenario & primary$method %in% null_methods &
      primary$metric %in% paste0(c("adjusted_", "unadjusted_"), metric), ]
    graphics::plot(c(.5, length(null_methods) + .5), range(c(d$lower, d$upper), finite = TRUE), type = "n",
      xaxt = "n", xlab = "", ylab = metric, main = gsub("_", " ", scenario))
    graphics::axis(1, at = seq_along(null_methods), labels = null_labels, las = 2, cex.axis = .8)
    graphics::abline(h = if (metric == "coverage") .95 else .05, lty = 3, col = "gray50")
    for (adjusted in c(FALSE, TRUE)) {
      x <- d[d$metric == paste0(if (adjusted) "adjusted_" else "unadjusted_", metric), ]
      at <- match(x$method, null_methods) + if (adjusted) .1 else -.1
      color <- if (adjusted) "#2166ac" else "#b35806"
      graphics::segments(at, x$lower, at, x$upper, col = color)
      graphics::points(at, x$mean, pch = if (adjusted) 19 else 17, col = color)
    }
  }
  graphics::mtext("Blue: adjusted trial regression; orange: RPE-only regression. Bars: Monte Carlo intervals.", outer = TRUE, cex = .85)
  grDevices::dev.off()

  grDevices::png(file.path(out_dir, "trial_variation.png"), width = 2000, height = 650, res = 170)
  graphics::par(mfrow = c(1, 3), mar = c(5, 4, 2, 1))
  scenarios <- c("no_trial_variation", "moderate_trial_variation", "large_trial_variation", "dependent_trial_variation")
  chosen <- c("latent_truth", "learned_estimated_noise", "oracle_hrf_known_noise")
  for (metric in c("rpe_r", "amplitude_r", "eta_observed_r")) {
    d <- primary[primary$scenario %in% scenarios & primary$method %in% chosen & primary$metric == metric, ]
    graphics::plot(c(1, 4), c(0, 1), type = "n", xaxt = "n", xlab = "Neural trial variation",
      ylab = metric, main = gsub("_", " ", metric))
    graphics::axis(1, 1:4, c("None", "Moderate", "Large", "Large AR"), cex.axis = .8)
    for (i in seq_along(chosen)) {
      x <- d[d$method == chosen[i], ]; x <- x[match(scenarios, x$scenario), ]
      graphics::lines(1:4, x$mean, col = c("black", "#2166ac", "#1b7837")[i], type = "b", pch = i)
    }
    if (metric == "rpe_r") graphics::legend("bottomleft", c("Latent truth", "Learned", "True HRF / true rho"),
      col = c("black", "#2166ac", "#1b7837"), lty = 1, pch = 1:3, bty = "n", cex = .65)
  }
  grDevices::dev.off()
  pr_write_checksums(out_dir)
  invisible(list(summary = summary, paired = paired, primary = wide))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Supply input directory and output directory")
  report_parametric_recovery(args[1L], args[2L])
}
