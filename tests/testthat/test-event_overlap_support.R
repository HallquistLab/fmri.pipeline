overlap_fixture <- function(data = data.frame(time = 0:12, value = 0:12), onsets = 6, vm = NULL) {
  fmri_ts$new(ts_data = data, tr = 1, vm = vm,
    event_data = data.frame(id = "A", run = 1, trial = seq_along(onsets), onset = onsets))
}

test_that("explicit missing samples fence interpolation at and across gaps", {
  d <- data.frame(time = 0:12, value = 0:12)
  d$value[d$time %in% 5:7] <- NA_real_
  result <- get_medusa_interpolated_ts(overlap_fixture(d), event = "onset",
    time_before = -3, time_after = 3, output_resolution = 0.5)
  missing <- abs(result$evt_time) < 2
  expect_true(all(is.na(result$value_mean[missing])))
  expect_true(all(result$value_n_missing[missing] == 1))
  expect_true(all(result$value_source_ids[missing] == ""))
  expect_equal(result$value_mean[!missing], 6 + result$evt_time[!missing])
  expect_true(all(result$epoch_status == "partial"))
})

test_that("deleted rows are detected on the acquisition clock, independently of output resolution", {
  d <- data.frame(time = c(0:4, 8:12), value = c(0:4, 8:12))
  object <- overlap_fixture(d)
  for (resolution in c(0.5, 1)) {
    result <- get_medusa_interpolated_ts(object, event = "onset", time_before = -3,
      time_after = 3, output_resolution = resolution)
    missing <- abs(result$evt_time) < 2
    expect_true(all(is.na(result$value_mean[missing])))
    expect_true(all(result$value_n_gap[missing] == 1))
  }
  allowed <- get_medusa_interpolated_ts(object, event = "onset", time_before = -3,
    time_after = 3, max_gap = 4)
  expect_equal(allowed$value_mean, 6 + allowed$evt_time)
  d <- data.frame(time = 0:12, value = 0:12)
  d$value[7] <- NA
  explicit <- get_medusa_interpolated_ts(overlap_fixture(d), event = "onset",
    max_gap = 4, time_before = 0, time_after = 0)
  expect_true(is.na(explicit$value_mean))
})

test_that("logical validity masks preserve censored endpoints and brackets", {
  d <- data.frame(time = 0:12, value = 0:12, usable = TRUE)
  d$usable[d$time == 6] <- FALSE
  d$usable[d$time == 8] <- NA
  result <- get_medusa_interpolated_ts(overlap_fixture(d), event = "onset",
    time_before = -1, time_after = 3, output_resolution = 0.5, valid_col = "usable")
  censored <- result$evt_time %in% c(-0.5, 0, 0.5, 1.5, 2, 2.5)
  expect_true(all(is.na(result$value_mean[censored])))
  expect_true(all(result$value_n_censored[censored] == 1))
  expect_equal(result$value_mean[!censored], 6 + result$evt_time[!censored])
  expect_error(get_medusa_interpolated_ts(overlap_fixture(d), event = "onset", valid_col = "value"),
    "logical column")
})

test_that("each unit must supply its own valid interpolation brackets before aggregation", {
  d <- data.frame(time = rep(0:1, 2), unit = rep(c("a", "b"), each = 2), value = c(0, 10, 100, NA))
  result <- get_medusa_interpolated_ts(overlap_fixture(d, 0.5, list(key = "unit")),
    event = "onset", time_before = 0, time_after = 0)
  expect_equal(result$value_mean, 5)
  expect_equal(result$value_median, 5)
  expect_true(is.na(result$value_sd))
  expect_equal(result$value_n, 1)
  expect_equal(result$value_n_total, 2L)
  expect_equal(result$value_n_missing, 1)
  expect_equal(result$value_source_ids, "1;2")
})

test_that("off-grid median and SD describe the interpolated unit distribution", {
  d <- data.frame(time = rep(0:1, 3), unit = rep(1:3, each = 2), value = c(0, 100, 10, 10, 100, 0))
  result <- get_medusa_interpolated_ts(overlap_fixture(d, 0.5, list(key = "unit")),
    event = "onset", time_before = 0, time_after = 0)
  expect_equal(result$value_mean, mean(c(50, 10, 50)))
  expect_equal(result$value_median, 50)
  expect_equal(result$value_sd, sd(c(50, 10, 50)))
})

test_that("shared source IDs count distinct epochs rather than repeated output lags", {
  result <- get_medusa_interpolated_ts(overlap_fixture(onsets = c(4, 6)), event = "onset",
    time_before = -2, time_after = 2, output_resolution = 0.5, time_audit = TRUE)
  overlap <- result$time >= 4 & result$time <= 6
  expect_true(all(result$value_source_max_epochs[overlap] == 2L))
  expect_true(all(result$align_window_overlap_n == 1L))
  single <- get_medusa_interpolated_ts(overlap_fixture(), event = "onset", output_resolution = 0.25)
  expect_true(all(single$value_source_max_epochs == 1L))
  same <- get_medusa_interpolated_ts(overlap_fixture(onsets = c(4, 4)), event = "onset")
  expect_true(all(same$align_simultaneous_n == 1L))
  expect_true(all(same$value_source_max_epochs == 2L))
  expect_equal(same$value_mean[same$trial == 1], same$value_mean[same$trial == 2])
})

test_that("allow, legacy truncation and whole-window exclusion have explicit distinct behavior", {
  object <- overlap_fixture(onsets = c(4, 6, 12, NA))
  legacy <- event_lock_ts(object, event = "onset", time_before = -3, time_after = 3,
    collide_before = "onset", collide_after = "onset", time_audit = TRUE)$get_ts(orig_names = TRUE)
  expect_true(all(legacy$time[legacy$trial == 1] < 6))
  expect_true(all(legacy$align_overlap_policy == "truncate"))
  allowed <- get_medusa_interpolated_ts(object, event = "onset", time_before = -3, time_after = 3,
    collide_before = "onset", collide_after = "onset", overlap_policy = "allow")
  expect_equal(allowed$value_mean[allowed$trial == 1 & allowed$evt_time == 3], 7)
  excluded <- get_medusa_interpolated_ts(object, event = "onset", time_before = -1, time_after = 1,
    overlap_policy = "exclude")
  expect_true(all(excluded$epoch_status[excluded$trial %in% 1:2] == "overlap_excluded"))
  expect_true(all(is.na(excluded$value_mean[excluded$trial %in% 1:2])))
  expect_true(all(is.finite(excluded$value_mean[excluded$trial == 3 & excluded$evt_time <= 0])))
  expect_true(all(excluded$epoch_status[excluded$trial == 4] == "missing_onset"))
})

test_that("neighbor event annotations include other event types and run boundaries", {
  object <- overlap_fixture(onsets = c(4, 6, NA))
  object$event_data$action <- c(4, 4, NA)
  result <- get_medusa_interpolated_ts(object, event = "onset", time_before = 0, time_after = 0,
    neighbor_events = c("onset", "action"))
  first <- result[result$trial == 1, ]
  expect_equal(first$align_simultaneous_n, 2L)
  expect_equal(first$align_next_event, 6)
  expect_true(is.na(first$align_previous_event))
  expect_equal(first$align_run_start_distance, 4)
  expect_equal(first$align_run_end_distance, 8)
  expect_true(is.na(result$align_window_overlap_n[result$trial == 3]))
})

test_that("multivariate support and averaging across trials retain contributor counts", {
  d <- data.frame(time = 0:12, a = 0:12, b = NA_real_)
  result <- get_medusa_interpolated_ts(overlap_fixture(d, c(4, 6), list(value = c("a", "b"))),
    event = "onset", time_before = 0, time_after = 0, group_by = NULL)
  expect_equal(result$a_mean, 5)
  expect_equal(result$a_n, 2)
  expect_equal(result$b_n, 0)
  expect_equal(result$b_n_missing, 2)
})

test_that("duplicate acquisition times within a unit are errors", {
  d <- data.frame(time = c(0:12, 6), value = c(0:12, 6))
  expect_error(get_medusa_interpolated_ts(overlap_fixture(d), event = "onset"), "Duplicate source times")
})

test_that("matrix interpolation agrees with independent per-unit interpolation on masked grids", {
  set.seed(20261106)
  data <- data.table::as.data.table(expand.grid(time = 0:12, unit = 1:5))
  data$value <- rnorm(nrow(data))
  data$.align_source_id <- data$time + 1L
  data$.align_censored <- sample(c(TRUE, FALSE), nrow(data), replace = TRUE, prob = c(0.2, 0.8))
  data$value[sample(nrow(data), 10)] <- NA_real_
  data$.align_epoch_id <- 1L
  grid <- seq(-0.5, 13, 0.5)
  actual <- fmri.pipeline:::interpolate_epoch_group(data, "unit", "time", grid, 1, TRUE)
  oracle <- fmri.pipeline:::interpolate_epoch_group(data, "unit", "time", grid, 1, TRUE, use_matrix = FALSE)
  expect_equal(actual, oracle, ignore_attr = TRUE)
})
