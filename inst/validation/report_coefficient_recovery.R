# Run from repository root: Rscript --vanilla ... INPUT_DIR OUTPUT_DIR
source(file.path("inst", "validation", "run_coefficient_recovery.R"), local = TRUE)

cr_collect <- function(input_dir) {
  paths <- if (file.exists(file.path(input_dir, "settings.rds"))) file.path(input_dir, "settings.rds") else
    list.files(input_dir, pattern = "^settings[.]rds$", recursive = TRUE, full.names = TRUE)
  if (!length(paths)) stop("No frozen study settings found")
  settings <- lapply(paths, readRDS)
  invariant <- setdiff(names(settings[[1L]]), c("scenarios", "counts"))
  scenarios <- unlist(lapply(settings, function(x) names(x$scenarios)), use.names = FALSE)
  if (anyDuplicated(scenarios)) stop("Scenario shards overlap")
  key <- function(d) paste(d$scenario, d$replicate, d$method, sprintf("%.10g", d$lambda), sep = ":")
  metrics <- diagnostics <- list()
  for (i in seq_along(paths)) {
    s <- settings[[i]]; dir <- dirname(paths[i])
    if (!identical(s[invariant], settings[[1L]][invariant])) stop("Incompatible study protocols")
    if (!identical(s$script_md5, tools::md5sum(names(s$script_md5)))) stop("Study sources changed since freezing")
    expected <- do.call(rbind, lapply(names(s$scenarios), function(scenario) {
      d <- s$schedule[rep(seq_len(nrow(s$schedule)), each = s$counts[[scenario]]), ]
      d$scenario <- scenario; d$replicate <- rep(seq_len(s$counts[[scenario]]), nrow(s$schedule)); d
    }))
    d <- utils::read.csv(file.path(dir, "replicate_metrics.csv"))
    if (anyDuplicated(key(d)) || !setequal(key(d), key(expected))) stop("Incomplete or duplicate results")
    diagnostic <- utils::read.csv(file.path(dir, "fit_diagnostics.csv"))
    index <- match(diagnostic$scenario, names(cr_scenarios()))
    dataset_key <- function(x) paste(x$scenario, x$replicate, sep = ":")
    if (anyNA(index) || anyDuplicated(diagnostic$seed) ||
        any(diagnostic$seed != s$seed_base + index * 10000L + diagnostic$replicate) ||
        anyDuplicated(dataset_key(diagnostic)) ||
        !setequal(dataset_key(diagnostic), dataset_key(unique(expected[c("scenario", "replicate")]))))
      stop("Incorrect dataset seed provenance or diagnostics")
    metrics[[i]] <- d; diagnostics[[i]] <- diagnostic
  }
  list(paths = paths, settings = settings, metrics = do.call(rbind, metrics), diagnostics = do.call(rbind, diagnostics))
}

cr_pairs <- function() list(
  area_readout = c("learned/area/ordinary", "learned/onset/ordinary"),
  template_readout = c("learned/template/ordinary", "learned/onset/ordinary"),
  onset_regression = c("learned/onset/transformed", "learned/onset/ordinary"),
  area_regression = c("learned/area/transformed", "learned/area/ordinary"),
  template_regression = c("learned/template/transformed", "learned/template/ordinary"),
  combined_area = c("learned/area/transformed", "learned/onset/ordinary"),
  combined_template = c("learned/template/transformed", "learned/onset/ordinary"),
  area_vs_forward = c("learned/area/transformed", "learned/forward/gls"),
  area_oracle_gap = c("learned/area/transformed", "oracle/area/transformed"))

cr_paired <- function(metrics) {
  measures <- setdiff(names(metrics), c("scenario", "replicate", "method", "lambda"))
  do.call(rbind, lapply(names(cr_pairs()), function(name) {
    pair <- cr_pairs()[[name]]
    d <- merge(metrics[metrics$method == pair[1L], ], metrics[metrics$method == pair[2L], ],
      by = c("scenario", "replicate", "lambda"), suffixes = c("_a", "_b"))
    if (!nrow(d)) return(NULL)
    delta <- d[c("scenario", "replicate", "lambda")]; delta$method <- name
    for (metric in measures) delta[[metric]] <- d[[paste0(metric, "_a")]] - d[[paste0(metric, "_b")]]
    pr_summarize(delta, probability_bounds = FALSE)
  }))
}

cr_write_gz <- function(x, path) {
  con <- gzfile(path, "wt"); on.exit(close(con)); utils::write.csv(x, con, row.names = FALSE)
}
cr_write_metadata <- function(paths, out_dir) {
  lines <- c(paste("Report source MD5:", unname(tools::md5sum("inst/validation/report_coefficient_recovery.R"))),
    unlist(lapply(paths, function(path) c(paste("Shard:", dirname(path)),
      readLines(file.path(dirname(path), "environment.txt")), ""))))
  lines <- sub("[[:blank:]]+$", "", lines)
  while (length(lines) && !nzchar(tail(lines, 1L))) lines <- head(lines, -1L)
  writeLines(lines, file.path(out_dir, "provenance.txt"))
  files <- list.files(out_dir, full.names = TRUE)
  files <- files[basename(files) != "checksums.md5"]
  checksums <- tools::md5sum(files)
  writeLines(paste(unname(checksums), basename(files), sep = "  "), file.path(out_dir, "checksums.md5"))
}

report_coefficient_recovery <- function(input_dir, out_dir) {
  data <- cr_collect(input_dir)
  if (data$settings[[1L]]$phase != "evaluation" || !setequal(unique(data$metrics$scenario), names(cr_scenarios())))
    stop("Final report requires every evaluation scenario")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- pr_summarize(data$metrics); paired <- cr_paired(data$metrics)
  primary <- subset(summary, lambda == .01)
  wide <- stats::reshape(primary[c("scenario", "method", "lambda", "metric", "mean")],
    idvar = c("scenario", "method", "lambda"), timevar = "metric", direction = "wide")
  names(wide) <- sub("^mean[.]", "", names(wide))
  utils::write.csv(wide, file.path(out_dir, "primary.csv"), row.names = FALSE)
  utils::write.csv(primary, file.path(out_dir, "primary_intervals.csv"), row.names = FALSE)
  utils::write.csv(subset(paired, lambda == .01), file.path(out_dir, "paired_primary.csv"), row.names = FALSE)
  utils::write.csv(data$diagnostics, file.path(out_dir, "fit_diagnostics.csv"), row.names = FALSE)
  saveRDS(data$settings, file.path(out_dir, "settings.rds"))
  saveRDS(data$metrics, file.path(out_dir, "replicate_metrics.rds"), compress = "xz")
  cr_write_gz(summary, file.path(out_dir, "summary.csv.gz"))
  cr_write_gz(paired, file.path(out_dir, "paired.csv.gz"))
  methods <- c("learned/onset/ordinary", "learned/area/ordinary", "learned/template/ordinary",
    "learned/onset/transformed", "learned/area/transformed", "learned/template/transformed", "learned/forward/gls")
  labels <- c("Onset", "Area", "Template", "Onset + R", "Area + R", "Template + R", "Forward GLS")
  colors <- c("#b35806", "#1b7837", "#762a83", "#e08214", "#2166ac", "#b2182b", "black")
  point_plot <- function(scenario, metric, reference) {
    d <- primary[primary$scenario == scenario & primary$method %in% methods & primary$metric == metric, ]
    d <- d[match(methods, d$method), ]
    limits <- range(c(d$lower, d$upper, reference), finite = TRUE)
    graphics::plot(c(.5, 7.5), limits, type = "n", xaxt = "n", xlab = "", ylab = metric,
      main = gsub("_", " ", scenario), cex.main = .85)
    graphics::axis(1, 1:7, labels, las = 2, cex.axis = .8)
    graphics::abline(h = reference, lty = 3, col = "gray50")
    graphics::segments(1:7, d$lower, 1:7, d$upper, col = colors)
    graphics::points(1:7, d$mean, col = colors, pch = 19)
  }
  grDevices::png(file.path(out_dir, "coefficient_recovery.png"), width = 1800, height = 1300, res = 160)
  graphics::par(mfrow = c(2, 2), mar = c(6.5, 4, 2.1, 1))
  for (scenario in c("moderate_trial_variation", "dependent_trial_variation", "confounded_effect", "off_library_effect"))
    point_plot(scenario, "beta", .6)
  grDevices::dev.off()
  grDevices::png(file.path(out_dir, "null_calibration.png"), width = 2200, height = 850, res = 160)
  graphics::par(mfrow = c(1, 3), mar = c(6.5, 4, 2.1, 1))
  for (scenario in c("active_null", "confounded_null", "off_library_null"))
    point_plot(scenario, "coefficient_false_positive", .05)
  grDevices::dev.off()
  grDevices::png(file.path(out_dir, "regularization.png"), width = 2000, height = 1350, res = 170)
  graphics::layout(matrix(c(1:6, 7, 7, 7), ncol = 3, byrow = TRUE), heights = c(1, 1, .3))
  graphics::par(mar = c(3.6, 4.2, 2.1, .6), mgp = c(2.2, .6, 0))
  for (scenario in c("moderate_trial_variation", "off_library_effect")) for (metric in c("beta_bias", "beta_rmse", "coefficient_coverage")) {
    d <- summary[summary$scenario == scenario & summary$metric == metric & summary$method %in% methods, ]
    graphics::plot(range(d$lambda), range(d$mean, finite = TRUE), log = "x", type = "n", xlab = "lambda",
      ylab = metric, main = gsub("_", " ", scenario), cex.main = .8)
    graphics::abline(v = .01, col = "gray60", lty = 3)
    for (i in 1:7) {
      x <- d[d$method == methods[i], ]; x <- x[order(x$lambda), ]
      if (i == 7) graphics::points(x$lambda, x$mean, pch = 4) else
        graphics::lines(x$lambda, x$mean, col = colors[i], lty = i, lwd = 1.5)
    }
  }
  graphics::par(mar = rep(0, 4)); graphics::plot.new()
  graphics::legend("center", labels, col = colors, lty = c(1:6, NA), pch = c(rep(NA, 6), 4), ncol = 4, bty = "n", cex = .8)
  grDevices::dev.off()
  cr_write_metadata(data$paths, out_dir)
  invisible(list(summary = summary, paired = paired))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Supply input and output directories")
  report_coefficient_recovery(args[1L], args[2L])
}
