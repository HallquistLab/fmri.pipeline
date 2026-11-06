# Run from the checkout: Rscript --vanilla inst/examples/event_overlap_review.R /tmp/overlap-review
# Or run against an installed package. Writes review tables and a PNG.
if (file.exists('DESCRIPTION') && requireNamespace('pkgload', quietly = TRUE)) {
  pkgload::load_all('.', quiet = TRUE)
} else library(fmri.pipeline)
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1L] else tempfile('event-overlap-review-')
dir.create(out, recursive = TRUE, showWarnings = FALSE)

samples <- data.frame(time = 0:45, value = sin((0:45) / 4), usable = TRUE)
samples$usable[samples$time %in% 19:21] <- FALSE
events <- data.frame(id = 'demo', run = 1L, trial = c('early', 'middle', 'neighbor', 'edge'),
  onset = c(2, 20, 23.5, 43), response = c(1, 19, 22, 42))
object <- fmri_ts$new(samples, events, tr = 1)
results <- lapply(c('allow', 'truncate', 'exclude'), function(policy) {
  get_medusa_interpolated_ts(object, event = 'onset', time_before = -3, time_after = 3,
    output_resolution = 0.5, valid_col = 'usable', neighbor_events = c('onset', 'response'),
    overlap_policy = policy, collide_before = 'onset', collide_after = 'onset', time_audit = TRUE)
})
names(results) <- c('allow', 'truncate', 'exclude')
for (policy in names(results)) write.csv(results[[policy]], file.path(out, paste0('alignment_', policy, '.csv')), row.names = FALSE)
write.csv(attr(results$allow, 'source_clock'), file.path(out, 'source_clock.csv'), row.names = FALSE)

# The review uses fixed synthetic settings, not data-derived scientific defaults.
kernel <- c(0, 0.05, 0.2, 0.35, 0.25, 0.1, 0.05)
recovery <- deconvolution_impulse_diagnostics(events, kernel, n_time = 46, tr = 1,
  lambda = 0.01, time_before = -3, time_after = 3, output_resolution = 0.5,
  input_demean = TRUE, tail_padding = length(kernel), return_timeseries = TRUE)
uncentered <- deconvolution_impulse_diagnostics(events, kernel, n_time = 46, tr = 1,
  lambda = 0.01, time_before = -3, time_after = 3, output_resolution = 0.5,
  input_demean = FALSE, tail_padding = length(kernel), return_timeseries = TRUE)
write.csv(uncentered$impulses, file.path(out, 'impulse_recovery_uncentered.csv'), row.names = FALSE)
write.csv(uncentered$windows, file.path(out, 'source_target_windows_uncentered.csv'), row.names = FALSE)
write.csv(recovery$impulses, file.path(out, 'impulse_recovery.csv'), row.names = FALSE)
write.csv(recovery$windows, file.path(out, 'source_target_windows.csv'), row.names = FALSE)
write.csv(recovery$time_series, file.path(out, 'impulse_timeseries.csv'), row.names = FALSE)
saveRDS(recovery$settings, file.path(out, 'frozen_settings.rds'))

png(file.path(out, 'overlap_review.png'), width = 1300, height = 1350, res = 140)
par(mfrow = c(3, 1), mar = c(4.5, 4.8, 3, 1))
d <- results$allow[results$allow$trial == 'middle', ]
plot(d$evt_time, d$value_mean, type = 'b', ylim = range(samples$value),
  xlab = 'Seconds from middle event', ylab = 'Aligned value', main = 'Censored samples stay missing')
abline(v = c(-1, 1), col = 'grey75', lty = 3)
d <- recovery$time_series[recovery$time_series$source_trial == 'middle', ]
plot(d$relative_time, d$activity, type = 'l', xlim = c(-8, 12), ylim = range(c(d$activity, d$truth)),
  xlab = 'Seconds from injected event', ylab = 'Recovered activity', main = 'One injected event spreads into neighboring windows')
u <- uncentered$time_series[uncentered$time_series$source_trial == 'middle', ]
lines(u$relative_time, u$activity, col = '#0072B2', lty = 2)
lines(d$relative_time, d$truth, col = '#555555', lty = 3)
abline(v = 3.5, col = '#D55E00', lty = 2)
legend('topright', c('Centered input', 'Uncentered input', 'Injected grid impulse', 'Neighbor onset'),
  col = c('black', '#0072B2', '#555555', '#D55E00'), lty = c(1, 2, 3, 2), bty = 'n')
w <- recovery$windows[recovery$windows$source_trial == 'middle', ]
wu <- uncentered$windows[uncentered$windows$source_trial == 'middle', ]
barplot(rbind(w$window_l1_ratio, wu$window_l1_ratio), beside = TRUE,
  names.arg = w$target_trial, col = c('grey35', '#0072B2'),
  ylab = 'Absolute-area ratio', main = 'Target-window area relative to source window (not attribution)')
abline(h = 1, lty = 3)
legend('topright', c('Centered input', 'Uncentered input'), fill = c('grey35', '#0072B2'), bty = 'n')
invisible(dev.off())
cat('Review artifacts:', normalizePath(out), '\n')
