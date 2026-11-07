cr_test_env <- function() {
  withr::local_dir(pkg_root)
  skip_if_not(file.exists("inst/validation/run_coefficient_recovery.R"), "Validation scripts require repository sources")
  e <- new.env(parent = globalenv())
  sys.source(file.path("inst", "validation", "run_coefficient_recovery.R"), e)
  e
}

test_that("contrast-specific estimation matches independent regression and detects nonidentification", {
  e <- cr_test_env(); set.seed(7731L)
  x <- cbind(1, matrix(rnorm(100), 20, 5)); y <- matrix(rnorm(60), 20, 3)
  contrast <- e$cr_contrast(x, 3L)
  expect_true(contrast$identified)
  expect_equal(as.numeric(crossprod(contrast$weights, y)),
    as.numeric(stats::lm.fit(x, y)$coefficients[3L, ]), tolerance = 1e-10)
  expect_equal(as.numeric(crossprod(contrast$weights, x)), c(0, 0, 1, 0, 0, 0), tolerance = 1e-10)
  repeated <- e$cr_contrast(cbind(x, x[, 1L]), 3L)
  expect_equal(repeated$weights, contrast$weights, tolerance = 1e-10)
  unidentified <- e$cr_contrast(cbind(x, x[, 3L]), 3L)
  expect_false(unidentified$identified); expect_true(all(is.na(unidentified$weights)))
})

test_that("event model represents cue, response and feedback effects without using true amplitudes", {
  e <- cr_test_env(); config <- e$cr_scenarios()$confounded_effect
  config$eta_sd <- 0
  d <- e$cr_generate(config, 7732L)
  coefficients <- c(.7, .5, -.35, .4, 1, .6, 0, 0, .35, .25)
  for (r in seq_along(d$bold)) {
    x <- e$cr_event_model(d$trials[[r]])
    expect_equal(sweep(matrix(as.numeric(x %*% coefficients), nrow(x), length(d$gain)),
      2L, d$gain, "*"), d$amplitudes[[r]], ignore_attr = TRUE, tolerance = 1e-12)
  }
})

test_that("off-library stress changes HRFs while retaining neural truth and noise", {
  e <- cr_test_env(); config <- e$cr_scenarios()$off_library_effect
  base <- config; base$peaks <- NULL
  shifted <- e$cr_generate(config, 7733L); original <- e$cr_generate(base, 7733L)
  expect_equal(shifted$noise, original$noise)
  expect_equal(shifted$truth, original$truth)
  expect_equal(shifted$trials, original$trials)
  expect_equal(shifted$amplitudes, original$amplitudes)
  expect_gt(max(abs(shifted$signal[[1L]] - original$signal[[1L]])), .01)
  expect_equal(colSums(shifted$kernels), rep(1, config$voxels), ignore_attr = TRUE)
})

test_that("readout geometry and fitted templates have the intended unit gains", {
  e <- cr_test_env(); n <- 180L; feedback <- c(30.2, 60.8, 90.5, 120.1)
  times <- c(feedback - 5, feedback - 2.5, feedback); b <- e$pr_deposit(n, times)
  h <- hrf_library(1)[, 1L]
  cache <- new.env(parent = emptyenv())
  d <- e$pr_operator(n, h, .5, .01, cache)
  forward <- e$pr_discrete_forward(n, h) %*% b
  g <- d %*% e$pr_center(forward)
  readouts <- e$cr_readouts(n, feedback, g)
  for (name in c("onset", "area")) expect_equal(
    diag(readouts[[name]] %*% b[, 9:12]), rep(1, 4), tolerance = 1e-12)
  expect_equal(diag(readouts$template %*% g[, 9:12]), rep(1, 4), tolerance = 1e-12)
  expect_error(e$cr_readouts(n, feedback, 0 * g), "Unidentified response template")
  expect_error(e$cr_interpolate(n, 0, -1:3), "outside")
})

test_that("transformed designs recover exact matched mean-model coefficients including signed voxels", {
  e <- cr_test_env(); config <- e$cr_scenarios()$confounded_effect
  config$eta_sd <- 0; config$voxels <- 4L; config$runs <- 3L
  d <- e$cr_generate(config, 7734L); fit <- e$pr_fit(d)
  cache <- new.env(parent = emptyenv())
  runs <- lapply(seq_along(d$bold), function(r) e$cr_run(d, fit, r, "matched_noiseless", .01, cache))
  for (name in c("onset", "area", "template")) {
    for (r in seq_along(runs)) expect_equal(runs[[r]]$readouts[[name]]$transformed,
      .6 * d$gain, tolerance = 1e-7)
    metrics <- e$cr_metrics(runs, d, name, "transformed")
    expect_equal(metrics[["beta"]], .6, tolerance = 1e-7)
    expect_equal(metrics[["identified_fraction"]], 1)
    expect_true(is.na(metrics[["coefficient_coverage"]]))
    ordinary <- e$cr_metrics(runs, d, name, "ordinary")
    expect_equal(metrics[c("amplitude_r", "amplitude_rmse", "neural_r")],
      ordinary[c("amplitude_r", "amplitude_rmse", "neural_r")])
  }
  config$beta <- 0
  null <- e$cr_generate(config, 7734L)
  result <- e$cr_run(null, fit, 1L, "matched_noiseless", .01, cache)
  for (name in c("onset", "area", "template")) expect_equal(result$readouts[[name]]$transformed,
    rep(0, config$voxels), tolerance = 1e-7)
  expect_equal(e$cr_forward(null, fit, "matched_noiseless")[["beta"]], 0, tolerance = 1e-8)
  expect_true(is.na(e$cr_forward(null, fit, "matched_noiseless")[["amplitude_r"]]))
})

test_that("readouts remain RPE independent and match the public neural inverse", {
  e <- cr_test_env(); config <- e$cr_scenarios()$moderate_trial_variation
  config$voxels <- 4L; config$runs <- 3L
  d <- e$cr_generate(config, 7735L); fit <- e$pr_fit(d)
  cache <- new.env(parent = emptyenv())
  result <- e$cr_run(d, fit, 1L, "learned", .01, cache)
  public <- deconvolve_crossrun_hrf(d$bold, fit, .01, scaling = "demean")
  expect_equal(result$activity, public$activity[[1L]], ignore_attr = TRUE, tolerance = 1e-8)
  changed <- d; changed$trials[[1L]]$rpe <- 2 * changed$trials[[1L]]$rpe
  second <- e$cr_run(changed, fit, 1L, "learned", .01, cache)
  expect_equal(second$activity, result$activity)
  for (name in c("onset", "area", "template")) {
    expect_equal(second$readouts[[name]]$amplitude, result$readouts[[name]]$amplitude)
    expect_equal(second$readouts[[name]]$transformed, .5 * result$readouts[[name]]$transformed, tolerance = 1e-8)
  }
  runs <- lapply(seq_along(d$bold), function(r) e$cr_run(d, fit, r, "learned", .01, cache))
  previous <- lapply(runs, function(x) list(activity = x$activity, eta = matrix(0, config$n, config$voxels)))
  reference <- e$pr_dataset_metrics(previous, d)
  current <- e$cr_metrics(runs, d, "onset", "ordinary")
  expect_equal(current[["beta_bias"]], reference[["adjusted_beta_bias"]], tolerance = 1e-10)
  expect_equal(current[["amplitude_r"]], reference[["amplitude_r"]], tolerance = 1e-10)
  expect_error(e$cr_run(d, fit, 1L, "invalid", .01, cache))
})

test_that("coefficient study schedule keeps references at primary lambda and protects outputs", {
  e <- cr_test_env(); schedule <- e$cr_schedule()
  expect_equal(nrow(schedule), 39L)
  expect_true(all(schedule$lambda[schedule$source != "learned"] == .01))
  expect_equal(length(unique(schedule$lambda[schedule$source == "learned"])), 4L)
  expect_false(anyDuplicated(schedule[c("method", "lambda")]) > 0L)
  tmp <- tempfile(); dir.create(tmp); writeLines("frozen", file.path(tmp, "settings.rds"))
  expect_error(e$run_coefficient_recovery(tmp, "development", 1L, 1L), "empty output")
})

test_that("coefficient report preserves pairing and rejects incomplete or changed protocols", {
  withr::local_dir(pkg_root)
  skip_if_not(file.exists("inst/validation/report_coefficient_recovery.R"), "Validation scripts require repository sources")
  e <- new.env(parent = globalenv())
  sys.source("inst/validation/report_coefficient_recovery.R", e)
  d <- data.frame(scenario = "test", replicate = rep(1:4, 2), lambda = .01,
    method = rep(c("learned/area/transformed", "learned/onset/ordinary"), each = 4),
    beta_rmse = c(.05, .1, .15, .2, .1, .2, .3, .4), coefficient_coverage = c(rep(1, 4), rep(.5, 4)))
  paired <- subset(e$cr_paired(d), method == "combined_area")
  expect_equal(paired$n, c(4L, 4L))
  expect_equal(paired$mean, c(-.125, .5))
  expect_equal(paired$mc_se[1L], sd(c(-.05, -.1, -.15, -.2)) / 2)
  tmp <- tempfile(); dir.create(tmp)
  files <- file.path("inst", "validation", c("parametric_recovery.R", "run_parametric_recovery.R",
    "coefficient_recovery.R", "run_coefficient_recovery.R"))
  settings <- list(phase = "evaluation", scenarios = e$cr_scenarios()["active_null"], counts = c(active_null = 1L),
    schedule = e$cr_schedule(), script_md5 = tools::md5sum(files), seed_base = 22000000L)
  saveRDS(settings, file.path(tmp, "settings.rds"))
  metrics <- settings$schedule[c("method", "lambda")]; metrics$scenario <- "active_null"; metrics$replicate <- 1L
  metrics$beta_rmse <- .1
  write.csv(metrics, file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  diag <- data.frame(scenario = "active_null", replicate = 1L, seed = 22050001L)
  write.csv(diag, file.path(tmp, "fit_diagnostics.csv"), row.names = FALSE)
  expect_equal(nrow(e$cr_collect(tmp)$metrics), 39L)
  write.csv(metrics[-1L, ], file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  expect_error(e$cr_collect(tmp), "Incomplete or duplicate")
  write.csv(rbind(metrics, metrics[1L, ]), file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  expect_error(e$cr_collect(tmp), "Incomplete or duplicate")
  write.csv(metrics, file.path(tmp, "replicate_metrics.csv"), row.names = FALSE)
  diag$seed <- 22050002L
  write.csv(diag, file.path(tmp, "fit_diagnostics.csv"), row.names = FALSE)
  expect_error(e$cr_collect(tmp), "Incorrect dataset seed")
  settings$script_md5[] <- "changed"
  saveRDS(settings, file.path(tmp, "settings.rds"))
  expect_error(e$cr_collect(tmp), "sources changed")
})
