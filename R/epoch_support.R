# Shared grid and support primitives for event alignment and impulse diagnostics.
epoch_output_grid <- function(time_before, time_after, output_resolution) {
  checkmate::assert_number(time_before, upper = 0, finite = TRUE)
  checkmate::assert_number(time_after, lower = 0, finite = TRUE)
  checkmate::assert_number(output_resolution, lower = .Machine$double.eps, finite = TRUE)
  grid <- c(-rev(seq_len(floor(-time_before / output_resolution + 1e-10))) * output_resolution,
    0, seq_len(floor(time_after / output_resolution + 1e-10)) * output_resolution)
  grid <- pmax(time_before, pmin(time_after, grid))
  if (anyDuplicated(grid)) stop("Duplicate epoch grid times")
  grid
}

resolve_overlap_policy <- function(overlap_policy, collide_before, collide_after) {
  if (is.null(overlap_policy)) {
    return(if (length(c(collide_before, collide_after))) "truncate" else "allow")
  }
  match.arg(overlap_policy, c("allow", "truncate", "exclude"))
}

# Build 9-point quadrature offsets and normalized weights around each target
# output time for local window-integrated readouts ("boxcar", "triangle",
# "gaussian"). Unlike two-point linear interpolation (which only evaluates a
# point value), local kernel quadrature integrates over a window of width
# `readout_width` (default 1.5 * max_gap) within each contiguous uncensored
# segment, improving SNR for slower or duration-modulated epoch signals while
# strictly refusing to integrate across censored observations or TR gaps.
epoch_kernel_quadrature <- function(readout_kernel = "linear", readout_width = NULL, max_gap = 1.0) {
  readout_kernel <- match.arg(readout_kernel, c("linear", "boxcar", "triangle", "gaussian"))
  if (readout_kernel == "linear") {
    return(list(kernel = "linear", width = 0, offsets = 0, weights = 1))
  }
  if (is.null(readout_width)) readout_width <- 1.5 * max_gap
  checkmate::assert_number(readout_width, lower = .Machine$double.eps, finite = TRUE)
  h <- 0.5 * readout_width
  offsets <- seq(-0.9 * h, 0.9 * h, length.out = 9L)
  weights <- switch(
    readout_kernel,
    boxcar = rep(1, length(offsets)),
    triangle = pmax(0, 1 - abs(offsets) / h),
    gaussian = exp(-0.5 * (offsets / (0.5 * h))^2)
  )
  weights <- weights / sum(weights)
  list(kernel = readout_kernel, width = readout_width, offsets = offsets, weights = weights)
}

# Never remove invalid rows before locating interpolation brackets. A missing
# row blocks its adjacent intervals; deleted rows are detected by elapsed time.
interpolate_epoch_unit <- function(time, value, source_id, censored, tout, max_gap,
    readout_kernel = "linear", readout_width = NULL) {
  quad <- epoch_kernel_quadrature(readout_kernel, readout_width, max_gap)
  keep <- is.finite(time)
  order <- order(time[keep])
  time <- time[keep][order]
  value <- value[keep][order]
  source_id <- source_id[keep][order]
  censored <- censored[keep][order]
  out <- data.table::data.table(evt_time = tout, .estimate = NA_real_,
    .support = if (length(time)) "outside_support" else "missing_epoch",
    .left = NA_integer_, .right = NA_integer_)
  if (!length(time)) return(out)
  if (anyDuplicated(time)) stop("Duplicate source times within an epoch unit")
  tol <- 1e-8 * max(1, max_gap)
  left <- findInterval(tout, time)
  near_left <- left > 0L & abs(tout - time[pmax(1L, left)]) <= tol
  near_right <- left < length(time) &
    abs(tout - time[pmin(length(time), left + 1L)]) <= tol
  exact <- near_left | near_right
  index <- ifelse(near_left, left, left + 1L)
  if (any(exact)) {
    rows <- which(exact)
    ii <- index[rows]
    valid <- !censored[ii] & is.finite(value[ii])
    out$.support[rows] <- ifelse(censored[ii], "censored", "missing_value")
    out$.estimate[rows[valid]] <- value[ii[valid]]
    out$.support[rows[valid]] <- "observed"
    out$.left[rows[valid]] <- source_id[ii[valid]]
  }
  rows <- which(!exact & left > 0L & left < length(time))
  if (length(rows)) {
    lo <- left[rows]
    hi <- lo + 1L
    blocked <- censored[lo] | censored[hi]
    missing <- !is.finite(value[lo]) | !is.finite(value[hi])
    gap <- time[hi] - time[lo] > max_gap + tol
    valid <- !(blocked | missing | gap)
    out$.support[rows] <- ifelse(blocked, "censored",
      ifelse(missing, "missing_value", ifelse(gap, "gap", "interpolated")))
    good <- rows[valid]
    weight <- (tout[good] - time[lo[valid]]) / (time[hi[valid]] - time[lo[valid]])
    out$.estimate[good] <- (1 - weight) * value[lo[valid]] + weight * value[hi[valid]]
    out$.left[good] <- source_id[lo[valid]]
    out$.right[good] <- source_id[hi[valid]]
  }
  if (quad$kernel != "linear") {
    good_rows <- which(out$.support %in% c("observed", "interpolated"))
    if (length(good_rows)) {
      usable <- !censored & is.finite(value)
      bridge <- if (length(time) > 1L) {
        usable[-length(time)] & usable[-1L] & (diff(time) <= max_gap + tol)
      } else logical(0L)
      seg_id <- c(1L, 1L + cumsum(!bridge))
      center_lo <- ifelse(exact[good_rows], index[good_rows], left[good_rows])
      for (idx in seq_along(good_rows)) {
        r <- good_rows[idx]
        s_id <- seg_id[center_lo[idx]]
        in_seg <- which(usable & seg_id == s_id)
        if (length(in_seg) == 1L) {
          out$.estimate[r] <- value[in_seg]
        } else {
          t_seg <- time[in_seg]
          v_seg <- value[in_seg]
          q_times <- tout[r] + quad$offsets
          in_bounds <- q_times >= min(t_seg) - tol & q_times <= max(t_seg) + tol
          if (any(in_bounds)) {
            q_clamped <- pmin(max(t_seg), pmax(min(t_seg), q_times[in_bounds]))
            q_vals <- stats::approx(t_seg, v_seg, xout = q_clamped, rule = 2L, ties = "ordered")$y
            w_sub <- quad$weights[in_bounds]
            out$.estimate[r] <- sum(w_sub * q_vals) / sum(w_sub)
          }
        }
      }
    }
  }
  out
}

# These are temporal annotations, not estimates of neural attribution. Windows
# are closed intervals, so touching endpoints count as overlap.
epoch_event_annotations <- function(events, event, trial, neighbor_events, before, after,
    run_start, run_end, policy) {
  onsets <- events[[event]]
  neighbors <- data.table::rbindlist(lapply(neighbor_events, function(column) {
    data.table::data.table(.row = seq_len(nrow(events)), .column = column,
      .onset = events[[column]])
  }))
  neighbors <- neighbors[is.finite(.onset)]
  data.table::rbindlist(lapply(seq_along(onsets), function(i) {
    onset <- onsets[i]
    other <- neighbors[!(.row == i & .column == event)]$.onset
    if (!is.finite(onset)) {
      return(data.table::data.table(align_onset = NA_real_, align_previous_event = NA_real_,
        align_next_event = NA_real_, align_simultaneous_n = NA_integer_,
        align_window_overlap_n = NA_integer_, align_run_start_distance = NA_real_,
        align_run_end_distance = NA_real_, align_overlap_policy = policy))
    }
    past <- other[other < onset]
    future <- other[other > onset]
    overlap <- is.finite(onsets) & seq_along(onsets) != i &
      onsets + before <= onset + after & onsets + after >= onset + before
    data.table::data.table(align_onset = onset,
      align_previous_event = if (length(past)) max(past) else NA_real_,
      align_next_event = if (length(future)) min(future) else NA_real_,
      align_simultaneous_n = sum(other == onset), align_window_overlap_n = sum(overlap),
      align_run_start_distance = onset - run_start,
      align_run_end_distance = run_end - onset, align_overlap_policy = policy)
  }))
}

# Return one row of statistics and exact epoch/acquisition memberships. Source
# IDs refer to acquisition times; memberships count distinct aligned epochs.
epoch_stat_row <- function(value, support, left, right, epoch, time) {
  valid <- is.finite(value)
  n <- sum(valid)
  ids <- sort(unique(c(left[valid], right[valid])))
  ids <- ids[!is.na(ids)]
  members <- c(paste(epoch[valid], left[valid], sep = ":"),
    paste(epoch[valid], right[valid], sep = ":"))
  members <- members[!grepl(":NA$", members)]
  list(evt_time = time, value_mean = if (n) mean(value[valid]) else NA_real_,
    value_median = if (n) stats::median(value[valid]) else NA_real_,
    value_sd = if (n > 1L) stats::sd(value[valid]) else NA_real_,
    value_n = n, value_n_total = length(value),
    value_n_censored = sum(support == "censored"),
    value_n_missing = sum(support %in% c("missing_value", "missing_epoch")),
    value_n_gap = sum(support == "gap"), value_n_outside = sum(support == "outside_support"),
    value_source_ids = paste(ids, collapse = ";"), .usage = list(unique(members)))
}

# Process one output group at a time to avoid materializing voxel x trial x
# output-grid arrays for the entire run. Common rectangular acquisition grids
# use shared brackets and matrix arithmetic; ragged grids use the unit oracle.
interpolate_epoch_group <- function(data, unit_keys, time_column, tout, max_gap,
    trial_grouped, use_matrix = TRUE, readout_kernel = "linear", readout_width = NULL) {
  quad <- epoch_kernel_quadrature(readout_kernel, readout_width, max_gap)
  data <- data.table::copy(data)
  if (length(unit_keys)) data[, .unit := .GRP, by = unit_keys] else data[, .unit := 1L]
  n_unit <- max(data$.unit)
  times <- sort(unique(data[[time_column]][is.finite(data[[time_column]])]))
  rectangular <- trial_grouped && length(times) && nrow(data) == n_unit * length(times) &&
    all(is.finite(data[[time_column]])) &&
    (length(times) < 2L || all(diff(times) <= max_gap + 1e-8 * max(1, max_gap)))
  if (anyDuplicated(data[, c(unit_keys, time_column), with = FALSE])) {
    stop("Duplicate source times within an epoch unit")
  }
  if (use_matrix && rectangular) {
    row <- match(data[[time_column]], times)
    index <- row + length(times) * (data$.unit - 1L)
    values <- censored <- matrix(NA_real_, length(times), n_unit)
    values[index] <- data$value
    censored[index] <- data$.align_censored
    ids <- data$.align_source_id[match(times, data[[time_column]])]
    template <- interpolate_epoch_unit(times, times, ids, rep(FALSE, length(times)), tout, max_gap)
    supported <- !is.na(template$.left)
    lo <- pmax(1L, match(template$.left, ids), na.rm = TRUE)
    hi <- match(template$.right, ids)
    hi[is.na(hi)] <- lo[is.na(hi)]
    weight <- rep(0, length(tout))
    between <- supported & !is.na(template$.right)
    weight[between] <- (tout[between] - times[lo[between]]) / (times[hi[between]] - times[lo[between]])
    blocked <- censored[lo, , drop = FALSE] | censored[hi, , drop = FALSE]
    missing <- !is.finite(values[lo, , drop = FALSE]) | !is.finite(values[hi, , drop = FALSE])
    estimates <- values[lo, , drop = FALSE] * (1 - weight) + values[hi, , drop = FALSE] * weight
    valid_center <- !(blocked | missing) & matrix(supported, nrow = length(tout), ncol = n_unit)
    if (quad$kernel != "linear" && length(times) > 1L) {
      bad_sample <- censored | !is.finite(values)
      num <- matrix(0, nrow = length(tout), ncol = n_unit)
      den <- matrix(0, nrow = length(tout), ncol = n_unit)
      for (m in seq_along(quad$offsets)) {
        tm <- interpolate_epoch_unit(times, times, ids, rep(FALSE, length(times)),
          tout + quad$offsets[m], max_gap)
        sup_m <- !is.na(tm$.left)
        lo_m <- pmax(1L, match(tm$.left, ids), na.rm = TRUE)
        hi_m <- match(tm$.right, ids)
        hi_m[is.na(hi_m)] <- lo_m[is.na(hi_m)]
        w_m <- rep(0, length(tout))
        btw_m <- sup_m & !is.na(tm$.right)
        w_m[btw_m] <- (tout[btw_m] + quad$offsets[m] - times[lo_m[btw_m]]) /
          (times[hi_m[btw_m]] - times[lo_m[btw_m]])
        est_m <- values[lo_m, , drop = FALSE] * (1 - w_m) + values[hi_m, , drop = FALSE] * w_m
        for (i in which(supported & sup_m)) {
          span_idx <- min(lo[i], lo_m[i]):max(hi[i], hi_m[i])
          ok_u <- valid_center[i, ] & (colSums(bad_sample[span_idx, , drop = FALSE]) == 0L)
          if (any(ok_u)) {
            num[i, ok_u] <- num[i, ok_u] + quad$weights[m] * est_m[i, ok_u]
            den[i, ok_u] <- den[i, ok_u] + quad$weights[m]
          }
        }
      }
      has_den <- den > 0
      estimates[has_den] <- num[has_den] / den[has_den]
    }
    estimates[!valid_center] <- NA_real_
    return(data.table::rbindlist(lapply(seq_along(tout), function(i) {
      support <- if (!supported[i]) rep("outside_support", n_unit) else
        ifelse(blocked[i, ], "censored", ifelse(missing[i, ], "missing_value", template$.support[i]))
      epoch_stat_row(estimates[i, ], support, rep(template$.left[i], n_unit),
        rep(template$.right[i], n_unit), rep(data$.align_epoch_id[1L], n_unit), tout[i])
    })))
  }
  units <- data[, {
    result <- interpolate_epoch_unit(get(..time_column), value, .align_source_id,
      .align_censored, tout, max_gap, readout_kernel = quad$kernel, readout_width = quad$width)
    result$.epoch <- .align_epoch_id[1L]
    result
  }, by = unit_keys]
  units[, {
    row <- epoch_stat_row(.estimate, .support, .left, .right, .epoch, evt_time[1L])
    row$evt_time <- NULL
    row
  }, by = "evt_time"]
}
