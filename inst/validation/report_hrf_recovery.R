# Summarize an already completed validation run. No model/metric selection.
# Rscript --vanilla inst/validation/report_hrf_recovery.R /tmp/hrf-validation
report_hrf_validation <- function(out_dir) {
  paired <- read.csv(file.path(out_dir, "paired_vs_canonical.csv"))
  summary <- read.csv(file.path(out_dir, "summary.csv"))
  settings <- readRDS(file.path(out_dir, "settings.rds"))
  metrics <- read.csv(file.path(out_dir, "replicate_metrics.csv"))
  stopifnot(nrow(metrics) == length(settings$scenarios) * settings$repetitions * 6L * length(settings$lambdas),
    !anyDuplicated(metrics[c("scenario", "replicate", "method", "lambda")]))
  methods <- c("condition_voxel", "trial_voxel", "trial_parcel", "trial_spatial")
  labels <- c("Condition / voxel", "Event / voxel", "Event / parcel", "Event / spatial")
  colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  scenarios <- setdiff(names(settings$scenarios), "null")
  forest <- function(metric, title, larger_better) {
    d <- paired[paired$lambda == .01 & paired$metric == metric & paired$method %in% methods &
      paired$scenario %in% scenarios, ]
    span <- range(c(0, d$lower, d$upper), na.rm = TRUE)
    padding <- max(diff(span) * .08, .005)
    plot(NA, xlim = span + c(-padding, padding), ylim = c(.5, length(scenarios) + .6),
      yaxt = "n", ylab = "", xlab = if (larger_better) "Change vs canonical; higher is better" else
        "Change vs canonical; lower is better", main = title)
    axis(2, at = seq_along(scenarios), labels = gsub("_", " ", rev(scenarios)), las = 1, cex.axis = .9)
    abline(v = 0, col = "gray55", lty = 2)
    abline(h = seq_along(scenarios), col = "gray92")
    for (m in seq_along(methods)) {
      sub <- d[d$method == methods[m], ]
      at <- match(sub$scenario, rev(scenarios)) + (m - 2.5) * .13
      segments(sub$lower, at, sub$upper, at, col = colors[m], lwd = 1.4)
      points(sub$delta, at, col = colors[m], pch = 16, cex = .7)
    }
  }
  draw <- function() {
    par(mfrow = c(1, 3), mar = c(5, 10, 4, 1), oma = c(3, 0, 2, 0))
    forest("neural_r", "Continuous neural recovery", TRUE)
    forest("rpe_r", "Trial-onset RPE tracking", TRUE)
    forest("rpe_gain_error", "Error in RPE gain (ideal gain = 1)", FALSE)
    mtext("Fixed lambda = 0.01; paired means and 95% Monte Carlo intervals across independent datasets",
      outer = TRUE, side = 3, line = .2, cex = 1)
    par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4), oma = rep(0, 4))
    plot.new()
    legend("bottom", legend = labels, col = colors, pch = 16, horiz = TRUE, bty = "n", inset = .01)
  }
  png(file.path(out_dir, "hrf_recovery_comparison.png"), width = 2400, height = 1150, res = 140)
  draw(); dev.off()
  pdf(file.path(out_dir, "hrf_recovery_comparison.pdf"), width = 17, height = 8.2)
  draw(); dev.off()
  heatmap <- function(metric, title) {
    d <- paired[paired$method == "trial_spatial" & paired$metric == metric, ]
    limit <- max(abs(d$delta), na.rm = TRUE)
    palette <- grDevices::colorRampPalette(c("#B2182B", "white", "#2166AC"))(101)
    plot(NA, xlim = c(.5, 3.5), ylim = c(.5, length(scenarios) + .5), axes = FALSE,
      xlab = "Fixed lambda", ylab = "", main = title)
    axis(1, at = 1:3, labels = settings$lambdas)
    axis(2, at = seq_along(scenarios), labels = gsub("_", " ", rev(scenarios)), las = 1)
    for (i in seq_len(nrow(d))) {
      y <- match(d$scenario[i], rev(scenarios)); x <- match(d$lambda[i], settings$lambdas)
      if (!is.na(y) && is.finite(d$delta[i])) {
        color <- palette[pmax(1, pmin(101, round(51 + 50 * d$delta[i] / limit)))]
        rect(x - .5, y - .5, x + .5, y + .5, col = color, border = "gray90")
        text(x, y, sprintf("%+.3f", d$delta[i]), cex = .85,
          col = if (abs(d$delta[i]) > .65 * limit) "white" else "black")
      }
    }
  }
  png(file.path(out_dir, "hrf_lambda_sensitivity.png"), width = 1500, height = 1150, res = 140)
  par(mfrow = c(1, 2), mar = c(5, 10, 4, 1))
  heatmap("neural_r", "Spatial HRF\nChange in neural r")
  heatmap("rpe_r", "Spatial HRF\nChange in onset RPE r")
  dev.off()
  # Compact committed tables retain all lambdas and Monte Carlo intervals.
  wide <- function(d, fields) {
    d <- d[c("scenario", "method", "lambda", "metric", fields)]
    reshape(d, idvar = c("scenario", "method", "lambda"), timevar = "metric", direction = "wide")
  }
  write.csv(wide(summary, c("mean", "lower", "upper", "n")),
    file.path(out_dir, "summary_wide.csv"), row.names = FALSE)
  write.csv(wide(paired, c("delta", "lower", "upper", "n", "improved_fraction")),
    file.path(out_dir, "paired_wide.csv"), row.names = FALSE)
  number <- function(x) ifelse(is.finite(x), sprintf("%.3f", x), "NA")
  value <- function(scenario, method, metric) {
    row <- summary[summary$scenario == scenario & summary$method == method &
      summary$metric == metric & summary$lambda == .01, ]
    number(row$mean)
  }
  interval <- function(scenario, method, metric) {
    row <- paired[paired$scenario == scenario & paired$method == method &
      paired$metric == metric & paired$lambda == .01, ]
    paste0(number(row$delta), " [", number(row$lower), ", ", number(row$upper), "]")
  }
  table <- c("| Scenario | Canonical neural r | Spatial change in neural r | Spatial change in onset RPE r | Spatial change in neural NRMSE |",
    "| --- | ---: | ---: | ---: | ---: |")
  for (scenario in scenarios) table <- c(table, paste0("| ", gsub("_", " ", scenario), " | ",
    value(scenario, "canonical", "neural_r"), " | ", interval(scenario, "trial_spatial", "neural_r"),
    " | ", interval(scenario, "trial_spatial", "rpe_r"), " | ",
    interval(scenario, "trial_spatial", "neural_nrmse"), " |"))
  writeLines(table, file.path(out_dir, "primary_table.md"))
  invisible(list(summary = summary, paired = paired))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  report_hrf_validation(if (length(args)) args[1] else "/tmp/hrf-validation")
}
