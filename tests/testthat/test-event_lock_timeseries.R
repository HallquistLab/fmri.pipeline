make_epoch_fixture <- function(onsets = c(4.5, 10.5, 16.5)) {
  fmri_ts$new(ts_data = data.frame(time = 0:20, value = 0:20), tr = 1,
    event_data = data.frame(id = "A", run = 1L, trial = seq_along(onsets), onset = onsets))
}

test_that("off-grid event alignment retains amplitude and audit times", {
  result <- get_medusa_interpolated_ts(make_epoch_fixture(), event = "onset",
    time_before = -2, time_after = 2, time_audit = TRUE)
  expect_equal(nrow(result), 15L)
  expect_equal(result$value_mean, result$evt_onset + result$evt_time)
  expect_equal(result$time, result$value_mean)
  expect_true(all(result$epoch_status == "ok"))
})

test_that("zero, one-sided, and sub-sample windows return unique requested times", {
  cases <- list(c(0, 0), c(-0.2, 0.2), c(0, 2), c(-2, 0), c(-1.1, 1.2))
  expected <- list(0, 0, 0:2, -2:0, -1:1)
  for (i in seq_along(cases)) {
    result <- get_medusa_interpolated_ts(make_epoch_fixture(10), event = "onset",
      time_before = cases[[i]][1L], time_after = cases[[i]][2L], output_resolution = 1)
    expect_equal(result$evt_time, as.numeric(expected[[i]]))
    expect_false(anyDuplicated(result$evt_time) > 0L)
    expect_true(all(result$evt_time >= cases[[i]][1L] & result$evt_time <= cases[[i]][2L]))
    expect_equal(result$value_mean, 10 + result$evt_time)
  }
  # A single observed sample is usable at exactly its own time without padding.
  result <- get_medusa_interpolated_ts(make_epoch_fixture(10), event = "onset",
    time_before = 0, time_after = 0, pad_before = 0, pad_after = 0)
  expect_equal(result$value_mean, 10)
})

test_that("missing first, intermediate, final, and all onsets retain trial rows", {
  for (onsets in list(c(NA, 10, 16), c(4, NA, 16), c(4, 10, NA),
      c(NA_real_, NaN, Inf))) {
    for (collisions in list(NULL, "onset")) {
      result <- get_medusa_interpolated_ts(make_epoch_fixture(onsets), event = "onset",
        collide_before = collisions, collide_after = collisions,
        time_before = -1, time_after = 1, time_audit = TRUE)
      expect_equal(sort(unique(result$trial)), 1:3)
      expect_equal(nrow(result), 9L)
      missing <- !is.finite(onsets[result$trial])
      expect_true(all(is.na(result$value_mean[missing])))
      expect_true(all(result$epoch_status[missing] == "missing_onset"))
      expect_true(all(is.na(result$time[missing])))
      expect_true(all(is.finite(result$value_mean[!missing])))
    }
  }
})

test_that("fully censored, partial, and missing-signal epochs have explicit support statuses", {
  result <- get_medusa_interpolated_ts(make_epoch_fixture(c(10.2, 10.5, 10.8)),
    event = "onset", collide_before = "onset", collide_after = "onset",
    time_before = -1, time_after = 1, pad_before = 0, pad_after = 0)
  expect_true(all(result$epoch_status[result$trial == 2L] == "censored"))
  expect_true(all(is.na(result$value_mean[result$trial == 2L])))
  partial <- get_medusa_interpolated_ts(make_epoch_fixture(0), event = "onset",
    time_before = -2, time_after = 2)
  expect_true(all(partial$epoch_status == "partial"))
  object <- make_epoch_fixture(10)
  object$ts_data$value <- NA_real_
  missing <- get_medusa_interpolated_ts(object, event = "onset")
  expect_true(all(missing$epoch_status == "insufficient_data"))
})

test_that("nonfinite collision candidates are ignored and finite bounds remain strict", {
  object <- make_epoch_fixture(c(4, 10, 16))
  object$event_data$boundary <- c(2, NA_real_, 18)
  object$event_data$other <- c(Inf, NaN, -Inf)
  raw <- event_lock_ts(object, event = "onset", time_before = -6, time_after = 6,
    collide_before = c("boundary", "other"), collide_after = c("boundary", "other"),
    time_audit = TRUE)$get_ts(orig_names = TRUE)
  expect_true(all(raw$time > raw$collide_before))
  expect_true(all(raw$time < raw$collide_after))
  expect_equal(unique(raw$collide_before), 2)
  expect_equal(unique(raw$collide_after), 18)
})

test_that("unsupported epochs remain explicit missing grids for each group", {
  object <- fmri_ts$new(
    ts_data = transform(expand.grid(time = 0:20, region = c("a", "b")), value = time),
    event_data = data.frame(id = 1, run = 1, trial = 1:2, onset = c(100, NA_real_)),
    vm = list(key = "region"), tr = 1)
  result <- get_medusa_interpolated_ts(object, event = "onset", time_before = -1,
    time_after = 1, group_by = c("trial", "region"))
  expect_equal(nrow(result), 12L)
  expect_true(all(is.na(result$value_mean)))
  expect_equal(unique(result$epoch_status[result$trial == 1]), "outside_run")
  expect_equal(unique(result$epoch_status[result$trial == 2]), "missing_onset")
})

test_that("original trial mappings survive event copying and missing epochs", {
  object <- fmri_ts$new(ts_data = data.frame(seconds = 0:20, signal = 0:20), tr = 1,
    event_data = data.frame(subject = "A", block = 1, attempt = c("first", "last"), onset = c(10, NA)),
    vm = list(id = "subject", run = "block", trial = "attempt", time = "seconds", value = "signal"))
  result <- get_medusa_interpolated_ts(object, event = "onset", time_before = 0,
    time_after = 0, group_by = "attempt")
  expect_equal(result$attempt, c("first", "last"))
  expect_equal(result$signal_mean, c(10, NA_real_))
  expect_equal(names(object$get_events(orig_names = TRUE)), c("subject", "block", "attempt", "onset"))
})

test_that("compress_mts_pca retains zero-mean columns and standardizes by SD", {
  t_seq <- seq(0, 2 * pi, length.out = 25L)
  base_mat <- cbind(
    v1 = sin(t_seq),
    v2 = sin(t_seq) + 0.25 * cos(t_seq),
    v3 = cos(t_seq),
    v4 = rep(0, length(t_seq)) # constant column should be screened out
  )
  # Zero-mean and negative-mean columns with disparate scales must be retained
  # and produce identical correlation-PCA compression to unit-scaled positive-mean columns.
  shifted_scaled <- cbind(
    v1 = 10 * (base_mat[, 1L] - mean(base_mat[, 1L])),
    v2 = 0.05 * (base_mat[, 2L] - mean(base_mat[, 2L])) - 5,
    v3 = 100 * (base_mat[, 3L] - mean(base_mat[, 3L])) - 50,
    v4 = base_mat[, 4L]
  )
  positive_unit <- cbind(
    v1 = as.numeric(scale(base_mat[, 1L])) + 10,
    v2 = as.numeric(scale(base_mat[, 2L])) + 10,
    v3 = as.numeric(scale(base_mat[, 3L])) + 10,
    v4 = base_mat[, 4L]
  )
  res_shifted <- fmri.pipeline:::compress_mts_pca(shifted_scaled, pexp_target = 0.8, scale_columns = TRUE)
  res_pos <- fmri.pipeline:::compress_mts_pca(positive_unit, pexp_target = 0.8, scale_columns = TRUE)
  expect_equal(res_shifted$compress_nt, 3L)
  expect_equal(res_shifted, res_pos, tolerance = 1e-12)
})

test_that("epoch kernel readouts match unit oracle and respect censor and gap boundaries", {
  set.seed(6101)
  ts_dt <- data.table::CJ(id = 1L, run = 1L, vnum = 1:3, time = seq(0, 30, by = 1.0))
  ts_dt[, signal := cos(time / 2.5) + 0.1 * vnum]
  ts_dt[, valid := !(time %in% c(9, 21))]
  ev_dt <- data.table::data.table(id = 1L, run = 1L, trial = 1:4, onset = c(4.4, 8.8, 15.2, 24.6))
  f_obj <- fmri_ts$new(ts_data = ts_dt, event_data = ev_dt, tr = 1.0,
    vm = list(value = "signal", key = "vnum"))

  locked <- event_lock_ts(f_obj, event = "onset", time_before = -2, time_after = 2,
    pad_before = -1.5, pad_after = 1.5, valid_col = "valid")
  melted <- data.table::melt(locked$get_ts(orig_names = TRUE),
    measure.vars = "signal", variable.name = ".vkey", value.name = "value")
  grid <- fmri.pipeline:::epoch_output_grid(-2, 2, 0.5)

  for (rk in c("linear", "boxcar", "triangle", "gaussian")) {
    mat_res <- melted[, fmri.pipeline:::interpolate_epoch_group(.SD, "vnum", "evt_time", grid, 1.0,
      trial_grouped = TRUE, use_matrix = TRUE, readout_kernel = rk, readout_width = 1.5), by = c(".vkey", "trial")]
    unit_res <- melted[, fmri.pipeline:::interpolate_epoch_group(.SD, "vnum", "evt_time", grid, 1.0,
      trial_grouped = TRUE, use_matrix = FALSE, readout_kernel = rk, readout_width = 1.5), by = c(".vkey", "trial")]
    expect_equal(mat_res$value_mean, unit_res$value_mean, tolerance = 1e-12)
    expect_equal(mat_res$value_n_censored, unit_res$value_n_censored)

    medusa <- get_medusa_interpolated_ts(f_obj, event = "onset", time_before = -2, time_after = 2,
      output_resolution = 0.5, valid_col = "valid", readout_kernel = rk, readout_width = 1.5)
    # Trial 2 (onset 8.8) at evt_time = 0 is bracketed by time = 8 and censored time = 9 -> must be NA
    t2_zero <- medusa[medusa$trial == 2 & abs(medusa$evt_time) < 1e-8, ]
    expect_true(is.na(t2_zero$signal_mean))
    expect_equal(t2_zero$signal_n_censored, 3L)
  }
})
