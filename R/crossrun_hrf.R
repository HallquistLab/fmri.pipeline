#' Build a small, constrained library of causal HRFs
#'
#' @param tr Acquisition interval in seconds.
#' @param peak_delays Positive response-delay parameters in seconds.
#' @param dispersions Positive response-dispersion parameters.
#' @param undershoot_ratios Positive response-to-undershoot ratios (\code{> 1};
#'   \code{Inf} gives a pure single-gamma HRF with zero post-stimulus undershoot).
#' @param duration Kernel duration in seconds (at least two acquisition intervals).
#' @param preset Optional preset grid configuration: \code{"custom"} (default, uses
#'   supplied \code{peak_delays}, \code{dispersions}, \code{undershoot_ratios}),
#'   \code{"standard"} (delays 5, 6, 7s; dispersions 0.8, 1, 1.2),
#'   \code{"extended"} (delays 4--8s; dispersions 0.8, 1, 1.25),
#'   \code{"canonical_27"} (fixed-undershoot 0.5s canonical grid: delays 4--8s by
#'   0.5s; dispersions 0.9, 1.0, 1.1; 27 candidates), \code{"dense"} (delays
#'   4--8.5s by 0.5s; dispersions 0.8, 1, 1.2, 1.4), \code{"manifold"} or
#'   \code{"medusa"} (pruned 3D physiological manifold spanning
#'   fast/narrow/minimal-undershoot subcortical through slow cortical shapes;
#'   ~44 candidates), or \code{"enet"} (maximin cosine epsilon-net covering the
#'   expanded 3D parameter cube with \code{n_enet} candidates).
#' @param normalize Normalization convention across candidate columns: \code{"auto"}
#'   (uses \code{"positive_lobe"} whenever undershoot varies across candidates,
#'   includes \code{Inf}, differs from the canonical fixed ratio 6, or uses
#'   \code{"manifold"}, \code{"medusa"}, or \code{"enet"} so varying undershoot
#'   ratios do not inflate deconvolved neural amplitudes; uses \code{"signed_sum"}
#'   for legacy fixed-ratio-6 grids), \code{"signed_sum"} (\code{sum(h) == 1}),
#'   \code{"positive_lobe"} (\code{sum(pmax(h, 0)) == 1}), or \code{"peak"}
#'   (\code{max(h) == 1}).
#' @param n_enet Number of candidates retained when \code{preset = "enet"} (default 45).
#' @return A normalized kernel matrix, with time rows, candidate columns, and
#'   \code{parameters}, \code{tr}, and \code{normalize} attributes. The canonical
#'   candidate \code{(6, 1, 6)} is first if present.
#' @details Candidates use the SPM difference-of-gammas parameterization with
#'   undershoot delay 16 seconds and undershoot dispersion 1. Setting
#'   \code{undershoot_ratios = Inf} omits the negative gamma lobe (single-gamma HRF).
#' @export
hrf_library <- function(tr, peak_delays = c(5, 6, 7), dispersions = c(0.8, 1, 1.2),
    undershoot_ratios = 6, duration = 32,
    preset = c("custom", "standard", "extended", "canonical_27", "dense", "manifold", "enet", "medusa"),
    normalize = c("auto", "signed_sum", "positive_lobe", "peak"),
    n_enet = 45L) {
  checkmate::assert_number(tr, lower = 0.01, finite = TRUE)
  checkmate::assert_number(duration, lower = 2 * tr, finite = TRUE)
  preset <- match.arg(preset)
  normalize <- match.arg(normalize)

  disp_min <- if (tr <= 1.2) 0.60 else min(0.80, max(0.60, round(0.45 * tr, 2L)))

  if (preset %in% c("medusa", "manifold", "enet")) {
    peak_delays <- seq(2.5, 8.0, by = 0.5)
    dispersions <- sort(unique(c(disp_min, 0.8, 1.0, 1.2, 1.4)))
    undershoot_ratios <- c(3, 6, 12, Inf)
  } else if (preset == "standard") {
    peak_delays <- c(5, 6, 7)
    dispersions <- c(0.8, 1, 1.2)
    undershoot_ratios <- 6
  } else if (preset == "extended") {
    peak_delays <- c(4, 5, 6, 7, 8)
    dispersions <- c(0.8, 1, 1.25)
    undershoot_ratios <- 6
  } else if (preset == "canonical_27") {
    peak_delays <- seq(4.0, 8.0, by = 0.5)
    dispersions <- c(0.9, 1.0, 1.1)
    undershoot_ratios <- 6
  } else if (preset == "dense") {
    peak_delays <- seq(4, 8.5, by = 0.5)
    dispersions <- c(0.8, 1, 1.2, 1.4)
    undershoot_ratios <- 6
  }

  for (x in list(peak_delays, dispersions)) {
    checkmate::assert_numeric(x, lower = .Machine$double.eps, min.len = 1L,
      any.missing = FALSE, finite = TRUE, unique = TRUE)
  }
  checkmate::assert_numeric(undershoot_ratios, lower = .Machine$double.eps,
    min.len = 1L, any.missing = FALSE, unique = TRUE)
  if (any(is.nan(undershoot_ratios))) {
    stop("undershoot_ratios must not contain NaN")
  }
  if (any(undershoot_ratios <= 1)) {
    stop("undershoot_ratios must be strictly greater than 1 (use Inf for zero undershoot)")
  }

  norm_mode <- if (normalize == "auto") {
    if (preset %in% c("medusa", "manifold", "enet") ||
        length(undershoot_ratios) > 1L ||
        any(is.infinite(undershoot_ratios)) ||
        !isTRUE(all.equal(as.numeric(undershoot_ratios), 6.0))) {
      "positive_lobe"
    } else {
      "signed_sum"
    }
  } else {
    normalize
  }

  pars <- expand.grid(peak_delay = peak_delays, dispersion = dispersions,
    undershoot_ratio = undershoot_ratios, KEEP.OUT.ATTRS = FALSE)
  if (preset %in% c("medusa", "manifold")) {
    keep_phys <- with(pars,
      !(peak_delay <= 2.5 & (dispersion > 0.85 | undershoot_ratio < 10)) &
      !(peak_delay <= 3.5 & (dispersion > 1.00 | undershoot_ratio < 6)) &
      !(peak_delay >= 6.5 & dispersion < 0.80) &
      !(peak_delay >= 7.5 & is.infinite(undershoot_ratio))
    )
    pars <- pars[keep_phys, , drop = FALSE]
  }
  canonical <- with(pars, peak_delay == 6 & dispersion == 1 & undershoot_ratio == 6)
  pars <- pars[order(!canonical, seq_len(nrow(pars))), , drop = FALSE]
  rownames(pars) <- NULL

  u_tr <- seq(0, floor(duration / tr)) * 16
  dt_hrf <- tr / 16
  raw_kernels <- vapply(seq_len(nrow(pars)), function(i) {
    h <- as.numeric(spm_hrf(tr, c(pars$peak_delay[i], 16, pars$dispersion[i], 1,
      pars$undershoot_ratio[i], 0, duration))$hrf)
    g1_tr <- stats::dgamma(u_tr, shape = pars$peak_delay[i] / pars$dispersion[i],
      rate = dt_hrf / pars$dispersion[i])
    idx_peak <- which.max(g1_tr)
    if (!all(is.finite(h)) || !length(idx_peak) || h[idx_peak] <= 0) {
      stop("Candidate HRF has non-positive net response on the TR grid; increase undershoot_ratios or duration")
    }
    normalize_decon_kernel(h, mode = norm_mode)
  }, numeric(floor(duration / tr) + 1L))

  if (preset %in% c("medusa", "manifold", "enet")) {
    n_target <- if (preset %in% c("medusa", "manifold")) {
      min(44L, ncol(raw_kernels))
    } else {
      checkmate::assert_integerish(n_enet, lower = 2L, len = 1L, any.missing = FALSE)
      min(as.integer(n_enet), ncol(raw_kernels))
    }
    keep_idx <- hrf_select_maximin_subset(raw_kernels, n_target, seed_index = 1L)
    raw_kernels <- raw_kernels[, keep_idx, drop = FALSE]
    pars <- pars[keep_idx, , drop = FALSE]
    rownames(pars) <- NULL
  }

  colnames(raw_kernels) <- paste0("hrf_", seq_len(ncol(raw_kernels)))
  attr(raw_kernels, "parameters") <- pars
  attr(raw_kernels, "tr") <- tr
  attr(raw_kernels, "normalize") <- norm_mode
  raw_kernels
}

# Deterministic greedy farthest-point subset selection in L2-normalized waveform
# angle (maximin cosine epsilon-net). Never touches .Random.seed.
hrf_select_maximin_subset <- function(kernels, n_keep, seed_index = 1L) {
  k_tot <- ncol(kernels)
  n_keep <- min(max(1L, as.integer(n_keep)), k_tot)
  if (n_keep >= k_tot) return(seq_len(k_tot))
  seed_index <- min(max(1L, as.integer(seed_index)), k_tot)
  norms <- sqrt(colSums(kernels^2))
  norms[norms <= 1e-12] <- 1.0
  u_mat <- sweep(kernels, 2L, norms, "/")
  cos_mat <- crossprod(u_mat)
  selected <- integer(n_keep)
  selected[1L] <- seed_index
  max_cos <- cos_mat[, seed_index]
  max_cos[seed_index] <- Inf
  for (step in seq.int(2L, n_keep)) {
    nxt <- which.min(max_cos)
    selected[step] <- nxt
    max_cos <- pmax(max_cos, cos_mat[, nxt])
    max_cos[nxt] <- Inf
  }
  selected
}

# Deterministic AR-stratum quantization for Step-1 candidate scoring when K is
# large and many parcels/voxels carry distinct AR parameters. Keeps exact strata
# whenever n_kernels <= 15 or unique strata per design <= max_strata, and never
# modifies fit$noise used in Step-2 deconvolution.
hrf_partition_ar_strata <- function(des_id, rho, psi1, psi2, psi3, psi4,
    n_kernels, max_strata = 6L) {
  n_vox <- length(des_id)
  ar_mat <- cbind(rho, psi1, psi2, psi3, psi4)
  raw_keys <- paste(des_id, rho, psi1, psi2, psi3, psi4, sep = ":")
  if (n_kernels <= 15L || length(unique(raw_keys)) <= max_strata) {
    return(list(keys = raw_keys, ar_mat = ar_mat))
  }
  keys_out <- raw_keys
  ar_out <- ar_mat
  for (d_val in unique(des_id)) {
    idx_d <- which(des_id == d_val)
    ok_d <- idx_d[rowSums(!is.finite(ar_mat[idx_d, , drop = FALSE])) == 0L]
    if (!length(ok_d)) next
    sub_keys <- raw_keys[ok_d]
    if (length(unique(sub_keys)) <= max_strata) next
    feat <- cbind(2.0 * rho[ok_d], psi1[ok_d], psi2[ok_d], psi3[ok_d], psi4[ok_d])
    ord <- order(rho[ok_d], psi1[ok_d], seq_along(ok_d))
    seed_pos <- unique(round(seq(1, length(ok_d), length.out = max_strata)))
    ctrs <- feat[ord[seed_pos], , drop = FALSE]
    cl <- rep(1L, length(ok_d))
    for (iter in seq_len(12L)) {
      d2 <- matrix(0.0, nrow = length(ok_d), ncol = nrow(ctrs))
      for (ci in seq_len(nrow(ctrs))) {
        diff_c <- sweep(feat, 2L, ctrs[ci, ], "-")
        d2[, ci] <- rowSums(diff_c * diff_c)
      }
      new_cl <- max.col(-d2, ties.method = "first")
      if (iter > 1L && identical(new_cl, cl)) break
      cl <- new_cl
      for (ci in seq_len(nrow(ctrs))) {
        m_ci <- which(cl == ci)
        if (length(m_ci)) ctrs[ci, ] <- colMeans(feat[m_ci, , drop = FALSE])
      }
    }
    for (ci in sort(unique(cl))) {
      m_vox <- ok_d[cl == ci]
      rep_ar <- colMeans(ar_mat[m_vox, , drop = FALSE])
      rep_ar[1L] <- round(rep_ar[1L], 4L)
      rep_ar[-1L] <- round(rep_ar[-1L], 4L)
      ar_out[m_vox, ] <- rep(rep_ar, each = length(m_vox))
      keys_out[m_vox] <- paste(c(d_val, rep_ar), collapse = ":")
    }
  }
  list(keys = keys_out, ar_mat = ar_out)
}

hrf_validate_bold <- function(bold) {
  checkmate::assert_list(bold, min.len = 1L, names = "unique")
  if (any(!nzchar(names(bold)))) stop("BOLD runs must have nonempty unique names")
  for (y in bold) {
    checkmate::assert_matrix(y, mode = "numeric", min.rows = 3L, min.cols = 1L)
    if (is.null(colnames(y)) || anyNA(colnames(y)) || any(!nzchar(colnames(y))) ||
        anyDuplicated(colnames(y))) stop("BOLD matrices require unique voxel column names")
    if (!identical(colnames(y), colnames(bold[[1L]])))
      stop("All BOLD runs must have identical voxel IDs in the same order")
  }
  invisible(TRUE)
}

hrf_run_list <- function(x, bold, argument, default) {
  if (is.null(x)) return(lapply(bold, default))
  if (!is.list(x) || !identical(names(x), names(bold)))
    stop(argument, " must be a list with the same run names and order as bold")
  x
}

hrf_run_numeric <- function(x, bold, argument) {
  checkmate::assert_numeric(x, any.missing = FALSE, finite = TRUE)
  if (length(x) == 1L) return(stats::setNames(rep(x, length(bold)), names(bold)))
  if (length(x) != length(bold) || !identical(names(x), names(bold)))
    stop(argument, " must be scalar or named in BOLD run order")
  x
}

# Build the per-run event-convolved design matrices across all candidate HRFs.
# In "condition" mode, events of the same condition share one regressor per
# candidate HRF. In "trial" mode, each event gets its own free amplitude column
# so zero-mean parametric signals (e.g., reward prediction errors) can still
# inform HRF shape selection. However, when multiple conditions are closely
# coupled within trials (median cross-condition onset separation <= 8 s) or
# event count exceeds residual degrees of freedom, simultaneous free columns
# across coupled conditions act like a local Taylor/derivative basis that can
# absorb HRF timing differences into trial amplitudes. In that case, we screen
# candidate designs that grant free trial columns to one focal condition at a
# time (pooling the remaining conditions) and select the focal condition with
# the strongest mean maximum partial R^2 across usable voxels.
hrf_build_run_design <- function(y, ev, tr, offset, kernels, z, mask,
    amplitude_model, rho_ref, batch_size, fallback_kernel, spatial_groups = NULL) {
  labels <- as.character(ev$condition)
  durations <- if (is.null(ev$duration)) rep(0, nrow(ev)) else ev$duration
  u_cond <- unique(labels)
  n_cond <- length(u_cond)
  if (amplitude_model == "condition") {
    design_columns <- match(labels, u_cond)
    return(hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
      design_columns, max(design_columns), kernels))
  }
  if (n_cond <= 1L || sum(mask) < 3L) {
    design_columns <- seq_len(nrow(ev))
    return(hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
      design_columns, max(design_columns), kernels))
  }
  # Check if multiple conditions are closely coupled within trials (forming a
  # local derivative basis when given simultaneous free trial parameters).
  cross_dists <- vapply(seq_len(nrow(ev)), function(e) {
    other <- ev$onset[labels != labels[e]]
    if (length(other)) min(abs(ev$onset[e] - other)) else Inf
  }, numeric(1))
  coupled <- stats::median(cross_dists) <= 8.0 || nrow(ev) > sum(mask) / 5.0
  if (!coupled) {
    design_columns <- seq_len(nrow(ev))
    return(hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
      design_columns, max(design_columns), kernels))
  }
  # Screen focal-condition trial-free designs using the reference AR(1) rho.
  # Allow spatial parcels (or strongly selective ungrouped voxels) to retain
  # their own focal condition rather than forcing a single whole-brain winner
  # that would collapse zero-mean trial-varying conditions (e.g., RPEs) in
  # focal regions when another condition dominates cortical voxel counts.
  cond_cols <- match(labels, u_cond)
  screen_all <- ncol(kernels) <= 12L
  screen_idx <- if (screen_all) seq_len(ncol(kernels)) else
    hrf_select_maximin_subset(kernels, 12L, seed_index = fallback_kernel)
  kernels_screen <- if (screen_all) kernels else kernels[, screen_idx, drop = FALSE]
  fb_screen <- if (screen_all) as.integer(fallback_kernel) else 1L
  cols_list <- vector("list", n_cond)
  des_list <- vector("list", n_cond)
  voxel_val_mat <- matrix(NA_real_, nrow = n_cond, ncol = ncol(y))
  comp_ic <- logical(n_cond)
  best_ic <- NA_integer_
  best_val <- -Inf
  for (ic in seq_len(n_cond)) {
    idx_e <- which(cond_cols == ic)
    cols_ic <- cond_cols
    cols_ic[idx_e] <- n_cond + seq_along(idx_e)
    cols_ic <- match(cols_ic, unique(cols_ic))
    cols_list[[ic]] <- cols_ic
    des_ic <- hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
      cols_ic, max(cols_ic), kernels_screen)
    if (screen_all) des_list[[ic]] <- des_ic
    part_ic <- hrf_score_run_cpp(y, des_ic, z, which(mask), as.integer(batch_size),
      rho = rho_ref, evidence = FALSE, fallback_index = fb_screen)
    if (isTRUE(part_ic$comparable)) {
      comp_ic[ic] <- TRUE
      max_scores <- apply(part_ic$scores, 2L, function(col) {
        if (any(is.finite(col))) max(col, na.rm = TRUE) else NA_real_
      })
      voxel_val <- max_scores
      voxel_val_mat[ic, ] <- voxel_val
      val <- if (any(is.finite(voxel_val))) mean(voxel_val, na.rm = TRUE) else -Inf
      if (is.finite(val) && val > best_val) {
        best_val <- val
        best_ic <- ic
      }
    }
  }
  if (!is.na(best_ic)) {
    vox_ic <- rep(best_ic, ncol(y))
    if (sum(comp_ic) > 1L && ncol(y) > 1L) {
      valid_ic <- which(comp_ic)
      if (!is.null(spatial_groups) && length(spatial_groups) == ncol(y)) {
        for (g in unique(spatial_groups)) {
          g_idx <- which(spatial_groups == g)
          if (!length(g_idx)) next
          g_vals <- rowMeans(voxel_val_mat[valid_ic, g_idx, drop = FALSE], na.rm = TRUE)
          if (any(is.finite(g_vals))) {
            ic_g <- valid_ic[which.max(ifelse(is.finite(g_vals), g_vals, -Inf))]
            base_g <- voxel_val_mat[best_ic, g_idx]
            base_g_mean <- if (any(is.finite(base_g))) mean(base_g, na.rm = TRUE) else -Inf
            best_g_mean <- g_vals[match(ic_g, valid_ic)]
            thr_floor <- if (length(g_idx) >= 3L) 0.025 else 0.05
            thr_gain <- if (length(g_idx) >= 3L) 0.025 else 0.05
            if (ic_g != best_ic && is.finite(best_g_mean) && best_g_mean > thr_floor &&
                (!is.finite(base_g_mean) || best_g_mean > base_g_mean + thr_gain)) {
              vox_ic[g_idx] <- ic_g
            }
          }
        }
      } else {
        for (v in seq_len(ncol(y))) {
          v_vals <- voxel_val_mat[valid_ic, v]
          if (any(is.finite(v_vals))) {
            ic_v <- valid_ic[which.max(ifelse(is.finite(v_vals), v_vals, -Inf))]
            base_v <- voxel_val_mat[best_ic, v]
            best_v <- voxel_val_mat[ic_v, v]
            if (ic_v != best_ic && is.finite(best_v) && best_v > 0.05 &&
                (!is.finite(base_v) || best_v > base_v + 0.03)) {
              vox_ic[v] <- ic_v
            }
          }
        }
      }
    }
    if (!screen_all) {
      for (ic in unique(c(best_ic, vox_ic))) {
        cols_ic <- cols_list[[ic]]
        des_list[[ic]] <- hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
          cols_ic, max(cols_ic), kernels)
      }
    }
    best_des <- des_list[[best_ic]]
    if (any(vox_ic != best_ic)) {
      attr(best_des, "voxel_design_id") <- vox_ic
      attr(best_des, "candidate_designs") <- des_list
    }
    return(best_des)
  }
  design_columns <- seq_len(nrow(ev))
  hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
    design_columns, max(design_columns), kernels)
}

# Partition a single run's usable observations into temporal sub-blocks
# separated by a short causal guard gap and score candidate HRFs on each
# sub-block using the full causal event design and full-run AR(1) noise
# estimates. Boundary carryover events whose HRF tail enters a sub-block are
# retained as a pooled carryover column so active BOLD tails are never omitted.
hrf_single_run_subblocks <- function(run_name, y, ev, tr, offset, kernels,
    candidate_ids, z, mask, noise_r, amplitude_model, batch_size,
    fallback_kernel, evidence, noise_model, spatial_groups = NULL) {
  usable_rows <- which(mask)
  n_obs <- length(usable_rows)
  n_kernels <- ncol(kernels)
  n_vox <- ncol(y)
  empty_out <- list(results = list(), scores = list(), names = character(0))
  if (n_obs < 24L) return(empty_out)

  labels <- as.character(ev$condition)
  durations <- if (is.null(ev$duration)) rep(0, nrow(ev)) else ev$duration
  u_cond <- unique(labels)
  cond_cols <- match(labels, u_cond)
  des_cond_full <- hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
    cond_cols, max(cond_cols), kernels)
  des_trial_full <- if (amplitude_model == "trial") {
    hrf_event_design_cpp(nrow(y), tr, offset, ev$onset, durations,
      seq_len(nrow(ev)), nrow(ev), kernels)
  } else NULL

  rho_ref <- if (any(is.finite(noise_r$rho))) stats::median(noise_r$rho, na.rm = TRUE) else 0
  block_candidates <- if (amplitude_model == "condition") {
    if (n_obs >= 150L) c(3L, 2L) else 2L
  } else {
    if (n_obs >= 240L && nrow(ev) >= 18L) c(3L, 2L) else 2L
  }

  append_carryover <- function(des_obj, carry_idx) {
    if (is.null(des_obj) || !length(carry_idx)) return(des_obj)
    add_col <- function(d_list) {
      Map(function(x, full_x) {
        cbind(x, rowSums(full_x[, carry_idx, drop = FALSE]))
      }, d_list, des_trial_full)
    }
    vox_id <- attr(des_obj, "voxel_design_id")
    cand_d <- attr(des_obj, "candidate_designs")
    out_d <- add_col(des_obj)
    if (!is.null(vox_id) && !is.null(cand_d)) {
      attr(out_d, "voxel_design_id") <- vox_id
      attr(out_d, "candidate_designs") <- lapply(cand_d, function(dl) if (is.null(dl)) NULL else add_col(dl))
    }
    out_d
  }

  best_out <- empty_out
  for (n_blocks in block_candidates) {
    block_id <- ceiling(seq_along(usable_rows) * n_blocks / n_obs)
    # Trim a short interior boundary buffer when sufficient frames exist so
    # AR(p) innovations across adjacent sub-blocks are weakly correlated.
    half_gap <- if (n_obs >= 72L * n_blocks / 2L) {
      min(max(1L, as.integer(floor(2.0 / max(tr, 0.25)))), as.integer(floor(n_obs / (10L * n_blocks))))
    } else 0L
    sub_results <- sub_scores <- vector("list", n_blocks)
    sub_names <- paste0(run_name, "_block", seq_len(n_blocks))
    names(sub_results) <- names(sub_scores) <- sub_names
    n_comparable_blocks <- 0L

    for (b in seq_len(n_blocks)) {
      idx_b <- which(block_id == b)
      if (half_gap > 0L && length(idx_b) > 2L * half_gap + 16L) {
        lo_trim <- if (b > 1L) half_gap else 0L
        hi_trim <- if (b < n_blocks) half_gap else 0L
        idx_b <- idx_b[(1L + lo_trim):(length(idx_b) - hi_trim)]
      }
      rows_b <- usable_rows[idx_b]
      mask_b <- rep(FALSE, nrow(y))
      mask_b[rows_b] <- TRUE

      design_b <- NULL
      if (amplitude_model == "condition") {
        col_ok <- vapply(seq_len(ncol(des_cond_full[[1L]])), function(c_idx) {
          norms_b <- vapply(des_cond_full, function(x) sqrt(sum(x[rows_b, c_idx]^2)), numeric(1))
          norms_all <- vapply(des_cond_full, function(x) sqrt(sum(x[, c_idx]^2)), numeric(1))
          min(norms_b) > 1e-3 * max(norms_all, 1e-8)
        }, logical(1))
        if (any(col_ok)) {
          design_b <- if (all(col_ok)) des_cond_full else
            lapply(des_cond_full, function(x) x[, col_ok, drop = FALSE])
        }
      } else {
        ev_ratios <- vapply(seq_len(nrow(ev)), function(e_idx) {
          norms_b <- vapply(des_trial_full, function(x) sqrt(sum(x[rows_b, e_idx]^2)), numeric(1))
          norms_all <- vapply(des_trial_full, function(x) sqrt(sum(x[, e_idx]^2)), numeric(1))
          c(min(norms_b) / max(norms_all, 1e-8), max(norms_b) / max(norms_all, 1e-8))
        }, numeric(2))
        ev_ok <- ev_ratios[1L, ] >= 0.12
        ev_carry <- (!ev_ok) & (ev_ratios[2L, ] > 1e-3)
        if (any(ev_ok)) {
          design_b <- hrf_build_run_design(y, ev[ev_ok, , drop = FALSE], tr, offset,
            kernels, z, mask_b, amplitude_model, rho_ref, batch_size, fallback_kernel,
            spatial_groups = spatial_groups)
          if (any(ev_carry)) {
            design_b <- append_carryover(design_b, which(ev_carry))
          }
        }
      }

      res_b <- list(
        scores = matrix(NA_real_, n_kernels, n_vox),
        signal_status = rep("insufficient_observations", n_vox)
      )
      if (evidence == "calibrated") {
        for (field in c("task_p", "fallback_p", "incompatible_p")) {
          res_b[[field]] <- matrix(NA_real_, n_kernels, n_vox)
        }
      }

      block_comp <- FALSE
      if (!is.null(design_b) && length(rows_b) >= 6L) {
        psi1_r <- if ("psi1" %in% names(noise_r)) noise_r$psi1 else rep(0, nrow(noise_r))
        psi2_r <- if ("psi2" %in% names(noise_r)) noise_r$psi2 else rep(0, nrow(noise_r))
        psi3_r <- if ("psi3" %in% names(noise_r)) noise_r$psi3 else rep(0, nrow(noise_r))
        psi4_r <- if ("psi4" %in% names(noise_r)) noise_r$psi4 else rep(0, nrow(noise_r))
        des_id_b <- attr(design_b, "voxel_design_id")
        cand_des_b <- attr(design_b, "candidate_designs")
        if (is.null(des_id_b) || is.null(cand_des_b)) {
          des_id_b <- rep(1L, n_vox)
          cand_des_b <- list(design_b)
        }
        ar_part_b <- hrf_partition_ar_strata(
          des_id_b, noise_r$rho, psi1_r, psi2_r, psi3_r, psi4_r, n_kernels
        )
        ar_keys_b <- ar_part_b$keys
        for (ak in unique(ar_keys_b)) {
          indices <- which(ar_keys_b == ak)
          idx1 <- indices[1L]
          ar_vec <- as.numeric(ar_part_b$ar_mat[idx1, ])
          des_use_b <- cand_des_b[[des_id_b[idx1]]]
          part <- hrf_score_run_ar(y[, indices, drop = FALSE], des_use_b, z, rows_b,
            as.integer(batch_size), ar_vec = ar_vec, evidence = evidence == "calibrated",
            fallback_index = as.integer(fallback_kernel))
          res_b$scores[, indices] <- part$scores
          res_b$signal_status[indices] <- part$signal_status
          if (evidence == "calibrated") {
            for (field in c("task_p", "fallback_p", "incompatible_p")) {
              res_b[[field]][, indices] <- part[[field]]
            }
          }
          if (isTRUE(part$comparable)) block_comp <- TRUE
        }
        if (noise_model == "ar1" && evidence == "calibrated") {
          unknown_noise <- grepl("^insufficient_noise", noise_r$method) &
            !res_b$signal_status %in% c("nonfinite", "constant", "nuisance_only")
          res_b$signal_status[unknown_noise] <- "noise_unidentified"
          for (field in c("scores", "task_p", "fallback_p", "incompatible_p")) {
            res_b[[field]][, unknown_noise] <- NA_real_
          }
        }
      }
      if (block_comp) n_comparable_blocks <- n_comparable_blocks + 1L
      sub_results[[b]] <- res_b
      sub_scores[[b]] <- res_b$scores
    }

    best_out <- list(results = sub_results, scores = sub_scores, names = sub_names)
    if (n_comparable_blocks >= 2L) break
  }
  best_out
}

# Compute a parcel-masked Gaussian spatial affinity matrix W_ij = exp(-d_ij^2 / (2 * sigma^2))
# from 3D voxel coordinates. When spatial_fwhm is NULL, sigma defaults to 1.5x
# the median nearest-neighbor distance. Cross-parcel weights are strictly zeroed
# so anatomical boundaries are respected.
hrf_compute_spatial_weights <- function(spatial_coords, spatial_groups, spatial_fwhm) {
  if (is.null(spatial_coords)) return(NULL)
  v <- nrow(spatial_coords)
  dmat <- as.matrix(stats::dist(spatial_coords))
  if (is.null(spatial_fwhm)) {
    pos_d <- apply(dmat, 1L, function(row) {
      p <- row[row > 0]
      if (length(p)) min(p) else NA_real_
    })
    sigma <- if (any(is.finite(pos_d))) stats::median(pos_d, na.rm = TRUE) * 1.5 else 1.0
    sigma <- max(sigma, 1e-6)
  } else {
    sigma <- max(spatial_fwhm / sqrt(8 * log(2)), 1e-6)
  }
  w <- exp(-(dmat^2) / (2 * sigma^2))
  if (!is.null(spatial_groups)) {
    same_group <- outer(spatial_groups, spatial_groups, "==")
    w[!same_group] <- 0
  }
  w
}

#' Estimate HRF shapes from other runs, with optional spatial regularization
#'
#' @param bold Named list of time-by-voxel numeric matrices for one participant.
#'   Voxel column IDs must be unique and identical in all runs, in the same order.
#' @param events Named list of data.frames with onset (seconds), condition, and
#'   optional duration (seconds, default zero). Names/order must match bold.
#' @param tr Common acquisition interval in seconds. All runs use the same grid.
#' @param kernel_grid Matrix or list of equal-length candidate causal HRFs, sampled
#'   at tr starting at lag zero. Use \code{hrf_library()} or supply explicit
#'   candidates; honors \code{attr(kernel_grid, "normalize")} (\code{"signed_sum"},
#'   \code{"positive_lobe"}, or \code{"peak"}), automatically switching to
#'   \code{"positive_lobe"} when unattributed candidates span heterogeneous
#'   undershoot ratios (and defaulting to \code{"signed_sum"} otherwise).
#' @param nuisance Optional matching list of time-by-regressor matrices. An
#'   intercept and linear trend are always included separately for each run.
#' @param valid_rows Optional matching list of logical observation masks; FALSE
#'   and NA omit rows from training fits without deleting their acquisition times.
#' @param time_offset First acquisition time, scalar or named numeric run vector.
#' @param spatial_groups Optional voxel parcel labels, named in voxel column order.
#'   Labels must be nonmissing and nonempty. Parcels must be defined independently
#'   of target-run responses; no implicit whole-brain parcel is constructed.
#' @param pooling Select HRFs per voxel or share one HRF per supplied parcel.
#' @param spatial_strength Nonnegative fixed penalty on voxel candidate distance
#'   from the training-derived parcel HRF; requires spatial_groups. Zero disables
#'   the soft penalty. Only used with voxel pooling; not automatically tuned.
#' @param fallback_kernel Index of the prespecified fallback candidate, default 1.
#' @param min_r2 Minimum mean adjusted partial R-squared for a supported selection,
#'   default 0.01. This is a configurable heuristic, not a significance threshold.
#' @param min_training_runs Minimum usable other runs per voxel, default 1.
#' @param min_parcel_voxels Minimum eligible training voxels per parcel, default 5.
#' @param batch_size Maximum voxel columns per native scoring block.
#' @param return_scores Retain candidate-by-voxel training score matrices per fold.
#' @param amplitude_model Training amplitudes shared within condition (fast
#'   default), or free for each event record (trial). Trial amplitudes can support
#'   HRF estimation when signed trial responses average to zero, but cost more and
#'   require a design with sufficient residual degrees of freedom.
#' @param noise_model Default ar1 uses pooled, grid-REML covariance estimates;
#'   ols disables whitening. Pooling uses supplied spatial_groups or the run.
#' @param noise_rho Optional fixed AR(1) coefficient, scalar or named in run
#'   order. Bypasses noise estimation. Requires noise_model = "ar1".
#' @param noise_grid Ordered AR(1) grid within [-0.95, 0.95], default -0.8 to 0.9
#'   in steps of 0.1. Boundary estimates are explicitly labeled.
#' @param evidence Default calibrated uses conditional GLS F tests with library
#'   corrections, prespecified-shape comparison, and training-run compatibility.
#'   Legacy retains the adjusted-R-squared heuristic. Conditional tests are exact
#'   for known Gaussian AR covariance and a correct mean model; estimated noise,
#'   non-Gaussian errors, and model misspecification require empirical validation.
#' @param evidence_alpha Per-voxel evidence threshold, default 0.05. Corrects the
#'   candidate search and compatibility comparisons across training runs, not
#'   multiple testing across the brain. Not a posterior probability.
#' @param canonical_strength Fixed shape-distance penalty toward fallback_kernel
#'   in calibrated mode, default 0.01. Zero removes this additional shrinkage.
#' @param spatial_coords Optional numeric matrix of 3D voxel coordinates (N x 3)
#'   used for local Gaussian-weighted anchor pooling.
#' @param spatial_fwhm Optional numeric spatial Gaussian FWHM (in coordinate
#'   units) for local anchor pooling.
#' @param selection_mode Candidate selection rule: `"hard"` (default) selects a
#'   single discrete winning HRF column index per voxel/parcel; `"soft_average"`
#'   constructs a continuous Empirical-Bayes / Boltzmann evidence-weighted
#'   composite HRF kernel across compatible neighborhood candidates blended with
#'   the prespecified fallback kernel.
#' @param noise_groups Optional voxel labels defining noise covariance pooling
#'   neighborhoods, in voxel column order. Defaults to \code{spatial_groups} when
#'   omitted. Allows noise neighborhoods (e.g., tissue- or sub-parcel-informed
#'   groups) to differ from HRF selection parcels.
#' @return An object of class crossrun_hrf with normalized kernels, per-target-run
#'   selections and training manifests, parcel anchors, run status/design summaries,
#'   optional candidate scores, voxel/run IDs, and frozen settings (including
#'   \code{settings$kernel_normalize}). When \code{selection_mode = "soft_average"},
#'   also includes per-target-run candidate weight matrices (\code{kernel_weights})
#'   and normalized composite HRF kernel matrices (\code{selected_kernels}).
#' @details When multiple runs are supplied, each target run is excluded from
#'   its own HRF selection, including parcel anchors and fallback decisions.
#'   When only a single run is supplied (`length(bold) == 1L`), usable
#'   observations are partitioned into contiguous temporal sub-blocks (3 blocks
#'   for longer runs, 2 halves for shorter runs) that are scored against the
#'   full causal event design using full-run AR(1) covariance estimates and
#'   pooled across sub-blocks with cross-block compatibility testing, while
#'   returning a single-run `crossrun_hrf` object for full-run deconvolution.
#'   Training fits have run-specific
#'   condition or event amplitudes and nuisance coefficients, with one shared
#'   candidate HRF across conditions. They do not fit unrestricted trial
#'   timecourses or use RPE.
#'   Native SVD projections are shared across voxel blocks. Scores are adjusted
#'   partial R-squared after removing nuisance, averaged equally across usable
#'   training runs. Parcel scores average these voxel scores equally, retaining
#'   voxel-specific amplitudes rather than averaging BOLD first. In voxel mode,
#'   selection maximizes score minus spatial_strength times normalized squared
#'   distance to the reliable parcel anchor. Weak/insufficient voxel fits fall
#'   back to a reliable parcel anchor, then the prespecified fallback kernel.
#'   This is a discrete shape penalty, not Bayesian posterior inference, a local
#'   adjacency prior, or spatial smoothing of neural activity. Library, thresholds,
#'   parcels, and penalty must be fixed externally or tuned within training data.
#'   Calibrated mode additionally penalizes distance from the prespecified shape.
#'   The legacy gate accepted null free-event fits in colored-noise simulations.
#'   Calibrated mode tests task support, departures from the prespecified shape
#'   using pairwise encompassing models, and compatibility with every usable
#'   training run. Fisher combination uses independent runs; Bonferroni corrections
#'   account for candidate comparisons. Unsupported or inconsistent shapes fall
#'   back. Parcel anchors are admitted either when at least min_parcel_voxels
#'   voxels are individually supported or when the parcel-level Cauchy-combined
#'   p-value (hrf_cauchy_p(), valid under arbitrary spatial correlation within a
#'   parcel) passes the candidate-corrected evidence threshold.
#'   Run compatibility cannot be assessed with a single training run. There is
#'   no access to target-run vascular shape or guarantee of transfer to that run.
#'   AR covariance is estimated per run using a fixed event-amplitude mean model
#'   and, when residual degrees of freedom allow, its first two temporal
#'   derivatives. Noise likelihood profiles are pooled, not the BOLD signals.
#'   All whitening uses the original acquisition clock, including gaps. No
#'   automatic nested tuning or physiological amplitude calibration is included.
#'   Nonfinite/constant/nuisance-only voxel runs are omitted; observation masks are
#'   common across voxels. Candidate task ranks must match within a run and leave
#'   residual degrees of freedom, otherwise that run cannot inform HRF selection.
#' @export
estimate_crossrun_hrf <- function(bold, events, tr, kernel_grid, nuisance = NULL,
    valid_rows = NULL, time_offset = 0, spatial_groups = NULL,
    pooling = c("voxel", "parcel"), spatial_strength = 0, fallback_kernel = 1L,
    min_r2 = 0.01, min_training_runs = 1L, min_parcel_voxels = 5L,
    batch_size = 2048L, return_scores = FALSE, amplitude_model = c("condition", "trial"),
    noise_model = c("ar1", "ols"), noise_rho = NULL,
    noise_grid = seq(-0.8, 0.9, by = 0.1), evidence = c("calibrated", "legacy"),
    evidence_alpha = 0.05, canonical_strength = 0.01,
    spatial_coords = NULL, spatial_fwhm = NULL,
    selection_mode = c("hard", "soft_average"),
    noise_groups = NULL) {
  hrf_validate_bold(bold)
  pooling <- match.arg(pooling)
  amplitude_model <- match.arg(amplitude_model)
  noise_model <- match.arg(noise_model)
  evidence <- match.arg(evidence)
  selection_mode <- match.arg(selection_mode)
  checkmate::assert_number(evidence_alpha, lower = .Machine$double.eps, upper = 0.5, finite = TRUE)
  checkmate::assert_number(canonical_strength, lower = 0, finite = TRUE)
  checkmate::assert_numeric(noise_grid, lower = -0.95, upper = 0.95, min.len = 2L,
    any.missing = FALSE, finite = TRUE, unique = TRUE, sorted = TRUE)
  if (noise_model == "ols" && !is.null(noise_rho)) stop("noise_rho requires noise_model = 'ar1'")
  if (!is.null(noise_rho)) {
    noise_rho <- hrf_run_numeric(noise_rho, bold, "noise_rho")
    checkmate::assert_numeric(noise_rho, lower = -0.95, upper = 0.95)
  }
  checkmate::assert_number(tr, lower = 0.01, finite = TRUE)
  if (!is.null(attr(kernel_grid, "tr")) && !isTRUE(all.equal(attr(kernel_grid, "tr"), tr)))
    stop("kernel_grid sampling TR does not match tr")
  parameters <- attr(kernel_grid, "parameters")
  explicit_norm <- if (!is.null(attr(kernel_grid, "normalize"))) {
    match.arg(attr(kernel_grid, "normalize"), c("signed_sum", "positive_lobe", "peak"))
  } else {
    NULL
  }
  kernels_list <- build_decon_kernel_grid(c(0, 1), kernel_grid, normalize = explicit_norm)
  kernel_normalize <- attr(kernels_list, "normalize")
  if (kernel_normalize == "signed_sum" &&
      any(vapply(kernels_list, function(h) sum(abs(h)) > 100, logical(1))))
    stop("Candidate signed sums are too close to zero for stable unit-sum normalization")
  kernels <- do.call(cbind, kernels_list)
  attr(kernels, "normalize") <- kernel_normalize
  candidate_ids <- if (is.matrix(kernel_grid)) colnames(kernel_grid) else names(kernel_grid)
  if (is.null(candidate_ids)) candidate_ids <- paste0("hrf_", seq_len(ncol(kernels)))
  if (anyNA(candidate_ids) || any(!nzchar(candidate_ids)) || anyDuplicated(candidate_ids))
    stop("Candidate names must be unique and nonempty")
  colnames(kernels) <- candidate_ids
  is_single_run <- length(bold) == 1L
  checkmate::assert_integerish(fallback_kernel, lower = 1, upper = ncol(kernels), len = 1L)
  checkmate::assert_number(spatial_strength, lower = 0, finite = TRUE)
  checkmate::assert_number(min_r2, lower = 0, upper = 1, finite = TRUE)
  checkmate::assert_integerish(min_training_runs, lower = 1,
    upper = if (is_single_run) 3L else length(bold) - 1L, len = 1L)
  checkmate::assert_integerish(min_parcel_voxels, lower = 1, len = 1L)
  checkmate::assert_integerish(batch_size, lower = 1, len = 1L)
  checkmate::assert_flag(return_scores)
  voxel_ids <- colnames(bold[[1L]])
  if (!is.null(spatial_coords)) {
    if (is.data.frame(spatial_coords)) spatial_coords <- as.matrix(spatial_coords)
    if (is.atomic(spatial_coords) && is.null(dim(spatial_coords)))
      spatial_coords <- matrix(spatial_coords, ncol = 1L, dimnames = list(names(spatial_coords), NULL))
    checkmate::assert_matrix(spatial_coords, mode = "numeric", nrows = length(voxel_ids),
      min.cols = 1L, any.missing = FALSE)
    if (any(!is.finite(spatial_coords))) stop("spatial_coords must be finite")
    if (!is.null(rownames(spatial_coords)) && !identical(rownames(spatial_coords), voxel_ids))
      stop("spatial_coords rownames must match voxel column order")
    if (is.null(spatial_groups))
      spatial_groups <- stats::setNames(rep("volume", length(voxel_ids)), voxel_ids)
  }
  if (!is.null(spatial_fwhm))
    checkmate::assert_number(spatial_fwhm, lower = .Machine$double.eps, finite = TRUE)
  if (!is.null(spatial_groups)) {
    if (!is.atomic(spatial_groups) || !identical(names(spatial_groups), voxel_ids) ||
        anyNA(spatial_groups) || any(!nzchar(as.character(spatial_groups))))
      stop("spatial_groups must be named nonmissing labels in voxel column order")
    spatial_groups <- as.character(spatial_groups)
  } else if (pooling == "parcel" || spatial_strength > 0) {
    stop("Parcel pooling/spatial regularization requires explicit spatial_groups")
  }
  if (pooling == "parcel" && spatial_strength != 0)
    stop("spatial_strength applies only to voxel pooling")
  if (!is.null(noise_groups)) {
    if (!is.atomic(noise_groups) || length(noise_groups) != length(voxel_ids) ||
        (!is.null(names(noise_groups)) && !identical(names(noise_groups), voxel_ids)) ||
        anyNA(noise_groups) || any(!nzchar(as.character(noise_groups))))
      stop("noise_groups must be nonmissing labels in voxel column order")
    noise_groups <- stats::setNames(as.character(noise_groups), voxel_ids)
  }
  spatial_weights <- hrf_compute_spatial_weights(spatial_coords, spatial_groups, spatial_fwhm)
  offsets <- hrf_run_numeric(time_offset, bold, "time_offset")
  events <- hrf_run_list(events, bold, "events", function(y) NULL)
  nuisance <- hrf_run_list(nuisance, bold, "nuisance", function(y) matrix(numeric(), nrow(y), 0))
  valid_rows <- hrf_run_list(valid_rows, bold, "valid_rows", function(y) rep(TRUE, nrow(y)))
  per_run <- run_results <- noise <- design_manifest <- event_manifest <- vector("list", length(bold))
  names(per_run) <- names(run_results) <- names(noise) <- names(design_manifest) <- names(event_manifest) <- names(bold)
  if (!is.null(spatial_groups)) names(spatial_groups) <- voxel_ids
  eff_noise_groups <- if (!is.null(noise_groups)) noise_groups else spatial_groups
  if (!any(vapply(events, function(e) is.data.frame(e) && nrow(e) >= 1L, logical(1L))))
    stop("At least one run must contain >= 1 event")
  for (r in seq_along(bold)) {
    y <- bold[[r]]; ev <- events[[r]]
    checkmate::assert_data_frame(ev, min.rows = 0L)
    checkmate::assert_subset(c("onset", "condition"), names(ev))
    checkmate::assert_numeric(ev$onset, any.missing = FALSE, finite = TRUE)
    labels <- as.character(ev$condition)
    if (nrow(ev) > 0L && (anyNA(labels) || any(!nzchar(labels))))
      stop("Event conditions must be nonmissing and nonempty")
    durations <- if (is.null(ev$duration)) rep(0, nrow(ev)) else ev$duration
    checkmate::assert_numeric(durations, lower = 0, any.missing = FALSE, finite = TRUE)
    checkmate::assert_matrix(nuisance[[r]], mode = "numeric", nrows = nrow(y))
    mask <- valid_rows[[r]]
    checkmate::assert_logical(mask, len = nrow(y))
    mask <- !is.na(mask) & mask
    z <- cbind(intercept = 1, trend = seq(-1, 1, length.out = nrow(y)), nuisance[[r]])
    if (any(!is.finite(z[mask, , drop = FALSE]))) stop("Nuisance values must be finite on usable rows")
    noise[[r]] <- hrf_estimate_noise(y, ev, tr, offsets[r], kernels, fallback_kernel,
      z, mask, eff_noise_groups, noise_grid, batch_size,
      if (noise_model == "ols") 0 else if (!is.null(noise_rho)) noise_rho[r] else NULL,
      spatial_coords = spatial_coords, spatial_fwhm = spatial_fwhm)
    rho_ref <- if (any(is.finite(noise[[r]]$rho))) stats::median(noise[[r]]$rho, na.rm = TRUE) else 0
    if (nrow(ev) == 0L) {
      result <- list(scores = matrix(NA_real_, ncol(kernels), ncol(y)),
        signal_status = rep("no_events", ncol(y)))
      if (evidence == "calibrated") for (field in c("task_p", "fallback_p", "incompatible_p"))
        result[[field]] <- matrix(NA_real_, ncol(kernels), ncol(y))
      z_obs <- z[mask, , drop = FALSE]
      sv_z <- if (nrow(z_obs) > 0L && ncol(z_obs) > 0L) svd(z_obs, nv = 0)$d else numeric(0L)
      rk_z <- if (length(sv_z) && sv_z[1L] > 0) sum(sv_z > 1e-8 * sv_z[1L]) else 0L
      run_results[[r]] <- result
      per_run[[r]] <- result$scores
      design_manifest[[r]] <- data.frame(run = names(bold)[r], candidate = candidate_ids,
        noise_rho = rho_ref, task_rank = rep(0L, ncol(kernels)), nuisance_rank = as.integer(rk_z),
        n_observations = sum(mask), comparable = FALSE)
      event_manifest[[r]] <- data.frame(run = names(bold)[r], voxel = voxel_ids,
        status = result$signal_status)
      next
    }
    design <- hrf_build_run_design(y, ev, tr, offsets[r], kernels, z, mask,
      amplitude_model, rho_ref, batch_size, fallback_kernel, spatial_groups = spatial_groups)
    result <- list(scores = matrix(NA_real_, ncol(kernels), ncol(y)),
      signal_status = rep("insufficient_observations", ncol(y)))
    if (evidence == "calibrated") for (field in c("task_p", "fallback_p", "incompatible_p"))
      result[[field]] <- matrix(NA_real_, ncol(kernels), ncol(y))
    manifest <- list()
    psi1_r <- if ("psi1" %in% names(noise[[r]])) noise[[r]]$psi1 else rep(0, nrow(noise[[r]]))
    psi2_r <- if ("psi2" %in% names(noise[[r]])) noise[[r]]$psi2 else rep(0, nrow(noise[[r]]))
    psi3_r <- if ("psi3" %in% names(noise[[r]])) noise[[r]]$psi3 else rep(0, nrow(noise[[r]]))
    psi4_r <- if ("psi4" %in% names(noise[[r]])) noise[[r]]$psi4 else rep(0, nrow(noise[[r]]))
    des_id_r <- attr(design, "voxel_design_id")
    cand_des_r <- attr(design, "candidate_designs")
    if (is.null(des_id_r) || is.null(cand_des_r)) {
      des_id_r <- rep(1L, ncol(y))
      cand_des_r <- list(design)
    }
    ar_part_r <- hrf_partition_ar_strata(
      des_id_r, noise[[r]]$rho, psi1_r, psi2_r, psi3_r, psi4_r, ncol(kernels)
    )
    ar_keys_r <- ar_part_r$keys
    for (ak in unique(ar_keys_r)) {
      indices <- which(ar_keys_r == ak)
      idx1 <- indices[1L]
      ar_vec <- as.numeric(ar_part_r$ar_mat[idx1, ])
      rho <- ar_vec[1L]
      des_use <- cand_des_r[[des_id_r[idx1]]]
      part <- if (sum(mask) >= 3L) {
        hrf_score_run_ar(y[, indices, drop = FALSE], des_use, z, which(mask),
          as.integer(batch_size), ar_vec = ar_vec, evidence = evidence == "calibrated",
          fallback_index = as.integer(fallback_kernel))
      } else list(scores = matrix(NA_real_, ncol(kernels), length(indices)),
        signal_status = rep("insufficient_observations", length(indices)),
        task_rank = rep(NA_integer_, ncol(kernels)), nuisance_rank = NA_integer_,
        n_observations = sum(mask), comparable = FALSE)
      result$scores[, indices] <- part$scores
      result$signal_status[indices] <- part$signal_status
      if (evidence == "calibrated" && sum(mask) >= 3L)
        for (field in c("task_p", "fallback_p", "incompatible_p")) result[[field]][, indices] <- part[[field]]
      manifest[[length(manifest) + 1L]] <- data.frame(run = names(bold)[r], candidate = candidate_ids,
        noise_rho = rho, task_rank = part$task_rank, nuisance_rank = part$nuisance_rank,
        n_observations = part$n_observations, comparable = part$comparable)
    }
    if (noise_model == "ar1" && evidence == "calibrated" && sum(mask) >= 3L) {
      unknown_noise <- grepl("^insufficient_noise", noise[[r]]$method) &
        !result$signal_status %in% c("nonfinite", "constant", "nuisance_only")
      result$signal_status[unknown_noise] <- "noise_unidentified"
      for (field in c("scores", "task_p", "fallback_p", "incompatible_p"))
        result[[field]][, unknown_noise] <- NA_real_
    }
    run_results[[r]] <- result
    per_run[[r]] <- result$scores
    design_manifest[[r]] <- do.call(rbind, manifest)
    event_manifest[[r]] <- data.frame(run = names(bold)[r], voxel = voxel_ids,
      status = result$signal_status)
  }
  subblock_fit <- NULL
  if (is_single_run) {
    mask_1 <- !is.na(valid_rows[[1L]]) & valid_rows[[1L]]
    z_1 <- cbind(intercept = 1, trend = seq(-1, 1, length.out = nrow(bold[[1L]])), nuisance[[1L]])
    subblock_fit <- hrf_single_run_subblocks(names(bold)[1L], bold[[1L]], events[[1L]],
      tr, offsets[1L], kernels, candidate_ids, z_1, mask_1, noise[[1L]],
      amplitude_model, batch_size, fallback_kernel, evidence, noise_model,
      spatial_groups = spatial_groups)
    if (return_scores && evidence == "calibrated") {
      attr(run_results[[1L]], "subblocks") <- subblock_fit$results
    }
  }
  # Relative squared shape distance is bounded by two; candidate indices have
  # no metric meaning. Penalizing an index difference would be arbitrary.
  distance <- outer(seq_len(ncol(kernels)), seq_len(ncol(kernels)), Vectorize(function(i, j) {
    sum((kernels[, i] - kernels[, j])^2) / sum(kernels[, i]^2 + kernels[, j]^2)
  }))
  k_norms <- sqrt(colSums(kernels^2))
  k_norms[k_norms <= 0] <- 1
  h_norm <- sweep(kernels, 2L, k_norms, "/")
  kernel_cos <- crossprod(h_norm)
  selections <- anchors <- candidate_scores <- vector("list", length(bold))
  kernel_weights <- selected_kernels <- vector("list", length(bold))
  training_runs <- vector("list", length(bold))
  names(selections) <- names(anchors) <- names(candidate_scores) <- names(training_runs) <- names(bold)
  names(kernel_weights) <- names(selected_kernels) <- names(bold)
  for (target in seq_along(bold)) {
    if (is_single_run) {
      fold_results <- subblock_fit$results
      fold_scores <- subblock_fit$scores
      train <- seq_along(fold_scores)
      training_runs[[target]] <- subblock_fit$names
    } else {
      fold_results <- run_results
      fold_scores <- per_run
      train <- setdiff(seq_along(bold), target)
      training_runs[[target]] <- names(bold)[train]
    }
    # Sum only training matrices, not total-minus-target: even extreme target
    # values cannot enter a fold through cancellation or missingness decisions.
    total <- matrix(0, ncol(kernels), length(voxel_ids))
    counts <- integer(length(voxel_ids))
    for (r in train) {
      usable <- is.finite(fold_scores[[r]][1L, ])
      total[, usable] <- total[, usable, drop = FALSE] + fold_scores[[r]][, usable, drop = FALSE]
      counts <- counts + usable
    }
    scores <- sweep(total, 2L, pmax(counts, 1L), "/")
    scores[, counts < min_training_runs] <- NA_real_
    dimnames(scores) <- list(candidate_ids, voxel_ids)
    gate <- if (evidence == "calibrated") hrf_fold_evidence(fold_results, train, counts,
      scores, fallback_kernel, evidence_alpha, min_r2, kernels = kernels,
      subblock_mode = is_single_run) else NULL
    anchor <- rep(NA_integer_, length(voxel_ids))
    parcel_n <- integer(length(voxel_ids))
    parcel_rows <- list()
    w_fold <- matrix(0, ncol(kernels), length(voxel_ids), dimnames = list(candidate_ids, voxel_ids))
    w_fold[fallback_kernel, ] <- 1.0
    w_parcel_map <- list()
    if (!is.null(spatial_groups)) for (g in unique(spatial_groups)) {
      voxels <- which(spatial_groups == g)
      eligible <- voxels[counts[voxels] >= min_training_runs]
      pooled <- if (length(eligible)) rowMeans(scores[, eligible, drop = FALSE]) else rep(NA_real_, ncol(kernels))
      p_ev <- hrf_parcel_pooled_evidence(gate, scores, eligible, fallback_kernel,
        evidence_alpha, min_r2, min_parcel_voxels, h_norm, kernel_cos,
        distance[, fallback_kernel], canonical_strength)
      admissible <- p_ev$admissible
      pooled_objective <- pooled
      if (!is.null(gate)) pooled_objective <- pooled_objective - canonical_strength * distance[, fallback_kernel]
      pooled_objective[!admissible] <- -Inf
      best <- if (any(is.finite(pooled_objective))) which.max(pooled_objective) else NA_integer_
      reliable <- length(eligible) >= min_parcel_voxels && !is.na(best) && pooled[best] >= min_r2
      if (reliable) anchor[voxels] <- best
      w_p <- rep(0, ncol(kernels))
      w_p[if (reliable) best else fallback_kernel] <- 1.0
      if (selection_mode == "soft_average" && length(eligible) >= min_parcel_voxels &&
          any(is.finite(pooled)) && max(pooled, na.rm = TRUE) >= min_r2) {
        omega_p <- if (is.null(gate)) 1.0 else if (reliable && best != fallback_kernel)
          max(0.85, p_ev$omega) else p_ev$omega
        w_p <- hrf_soft_boltzmann_weights(
          scores_vec = pooled,
          kernel_cos = kernel_cos,
          dist_fb = distance[, fallback_kernel],
          fallback_kernel = fallback_kernel,
          comp_ok = p_ev$comp_ok,
          omega = omega_p,
          canonical_strength = canonical_strength
        )
      }
      w_parcel_map[[g]] <- w_p
      if (pooling == "parcel") {
        w_fold[, voxels] <- w_p
      }
      parcel_n[voxels] <- length(eligible)
      parcel_rows[[length(parcel_rows) + 1L]] <- data.frame(target_run = names(bold)[target],
        parcel = g, n_voxels = length(eligible), kernel_index = if (reliable) best else fallback_kernel,
        adjusted_partial_r2 = if (!is.na(best)) pooled[best] else NA_real_,
        status = if (reliable) "supported" else "fallback_prespecified")
    }
    anchors[[target]] <- if (length(parcel_rows)) do.call(rbind, parcel_rows) else data.frame()
    chosen <- rep(as.integer(fallback_kernel), length(voxel_ids))
    reason <- rep("fallback_prespecified", length(voxel_ids))
    best <- rep(NA_integer_, length(voxel_ids)); margin <- rep(NA_real_, length(voxel_ids))
    selected_score <- best_score <- rep(NA_real_, length(voxel_ids))
    supported <- which(!is.na(anchor))
    if (!is.null(gate) && length(supported)) supported <- supported[
      gate$compatible[cbind(anchor[supported], supported)]]
    chosen[supported] <- anchor[supported]
    reason[supported] <- if (pooling == "parcel") "pooled_parcel" else "fallback_parcel"
    eligible <- which(counts >= min_training_runs)
    if (length(eligible)) {
      best[eligible] <- max.col(t(scores[, eligible, drop = FALSE]), ties.method = "first")
      best_score[eligible] <- scores[cbind(best[eligible], eligible)]
      if (ncol(kernels) > 1L) {
        second <- rep(-Inf, length(eligible))
        for (k in seq_len(ncol(kernels)))
          second <- pmax(second, ifelse(best[eligible] == k, -Inf, scores[k, eligible]))
        margin[eligible] <- best_score[eligible] - second
      }
      if (pooling == "voxel") {
        objective <- scores[, eligible, drop = FALSE]
        eff_anchor <- anchor[eligible]
        neigh_scores <- objective
        concordance <- rep(1.0, length(eligible))
        # Empirical-Bayes concordance-weighted spatial pooling:
        # Rather than blindly forcing every voxel toward the parcel anchor (which
        # would bias focal sub-regions or edge voxels that genuinely differ from
        # their neighbors), we measure the correlation (concordance c_v) between
        # voxel v's candidate score profile and its leave-one-voxel-out
        # neighborhood profile (optionally weighted by 3D Gaussian distance via
        # spatial_coords / spatial_fwhm). Shrinkage weight scales with c_v^2 so
        # discordant focal voxels retain their own supported HRF while noisy but
        # concordant voxels borrow strength from the parcel neighborhood.
        use_eb_spatial <- spatial_strength > 0 && (!is.null(gate) || !is.null(spatial_weights)) &&
          !is.null(spatial_groups)
        if (use_eb_spatial) {
          for (idx_e in seq_along(eligible)) {
            v <- eligible[idx_e]
            same_p <- which(spatial_groups[eligible] == spatial_groups[v])
            others_e <- setdiff(same_p, idx_e)
            if (!is.null(spatial_weights) && length(others_e) >= 1L) {
              w_vec <- spatial_weights[v, eligible[others_e]]
              if (sum(w_vec) > 1e-12) {
                w_norm <- w_vec / sum(w_vec)
                nv <- as.numeric(objective[, others_e, drop = FALSE] %*% w_norm)
                neigh_scores[, idx_e] <- nv
                # Local coordinate-weighted anchor within parcel
                if (length(same_p) >= min_parcel_voxels) {
                  loc_obj <- nv
                  if (!is.null(gate)) {
                    loc_sup <- as.numeric(gate$supported[, eligible[others_e], drop = FALSE] %*% w_norm) >= 0.35
                    loc_obj <- loc_obj - canonical_strength * distance[, fallback_kernel]
                    loc_obj[!loc_sup] <- -Inf
                  }
                  if (any(is.finite(loc_obj)) && max(nv, na.rm = TRUE) >= min_r2)
                    eff_anchor[idx_e] <- which.max(loc_obj)
                }
              }
            } else if (length(others_e) >= 1L) {
              nv <- rowMeans(objective[, others_e, drop = FALSE])
              neigh_scores[, idx_e] <- nv
            }
            sv <- objective[, idx_e]
            nv <- neigh_scores[, idx_e]
            if (ncol(kernels) > 2L) {
              concordance[idx_e] <- if (stats::sd(sv) > 1e-8 && stats::sd(nv) > 1e-8)
                max(0.0, stats::cor(sv, nv)) else 0.0
            } else if (ncol(kernels) == 2L) {
              concordance[idx_e] <- if ((sv[2L] - sv[1L]) * (nv[2L] - nv[1L]) >= 0) 1.0 else 0.0
            }
          }
          eb_base <- min(0.85, spatial_strength / (spatial_strength + 0.02))
          w_eb <- if (spatial_strength >= 1.0) rep(eb_base, length(eligible)) else eb_base * (concordance^2)
          objective <- sweep(objective, 2L, 1.0 - w_eb, "*") + sweep(neigh_scores, 2L, w_eb, "*")
        }
        penalized <- which(!is.na(eff_anchor))
        admissible_voxel <- if (is.null(gate)) NULL else gate$supported[, eligible, drop = FALSE]
        if (!is.null(gate)) {
          objective <- objective - canonical_strength * distance[, fallback_kernel]
          if (spatial_strength > 0 && length(penalized)) {
            has_own_depart <- if (ncol(kernels) > 1L)
              colSums(gate$supported[-fallback_kernel, eligible[penalized], drop = FALSE]) > 0L
            else rep(FALSE, length(penalized))
            borrow <- penalized[!has_own_depart | concordance[penalized] > 0.25 | spatial_strength >= 1.0]
            if (length(borrow)) {
              anchor_idx <- cbind(eff_anchor[borrow], borrow)
              admissible_voxel[anchor_idx] <- admissible_voxel[anchor_idx] |
                gate$task_supported[cbind(eff_anchor[borrow], eligible[borrow])]
            }
          }
          objective[!admissible_voxel] <- -Inf
        }
        if (spatial_strength > 0 && length(penalized)) {
          if (!is.null(gate) && spatial_strength < 1.0) {
            for (idx_p in penalized) {
              has_dep <- ncol(kernels) > 1L && any(gate$supported[-fallback_kernel, eligible[idx_p]])
              pen_w <- if (has_dep) concordance[idx_p]^2 else 1.0
              objective[, idx_p] <- objective[, idx_p] -
                spatial_strength * pen_w * distance[, eff_anchor[idx_p]]
            }
          } else {
            objective[, penalized] <- objective[, penalized, drop = FALSE] -
              spatial_strength * distance[, eff_anchor[penalized], drop = FALSE]
          }
        }
        candidates <- max.col(t(objective), ties.method = "first")
        keep <- scores[cbind(candidates, eligible)] >= min_r2
        if (!is.null(gate)) {
          keep <- keep & admissible_voxel[cbind(candidates, seq_along(eligible))]
          if (spatial_strength == 0 && length(penalized)) {
            no_depart <- colSums(gate$supported[-fallback_kernel, eligible[penalized], drop = FALSE]) == 0L
            anchor_task <- gate$task_supported[cbind(eff_anchor[penalized], eligible[penalized])]
            anchor_obj <- scores[cbind(eff_anchor[penalized], eligible[penalized])] -
              canonical_strength * distance[eff_anchor[penalized], fallback_kernel]
            fallback_obj <- scores[fallback_kernel, eligible[penalized]]
            keep[penalized[no_depart & anchor_task & anchor_obj > fallback_obj]] <- FALSE
          }
        }
        chosen[eligible[keep]] <- candidates[keep]
        reason[eligible[keep]] <- ifelse(spatial_strength > 0 & !is.na(eff_anchor[keep]),
          "voxel_spatial", "voxel")
        anchor[eligible] <- eff_anchor
        if (selection_mode == "soft_average") {
          if (!is.null(spatial_groups)) {
            for (g in names(w_parcel_map)) {
              g_sup <- intersect(which(spatial_groups == g), supported)
              if (length(g_sup)) w_fold[, g_sup] <- w_parcel_map[[g]]
            }
          }
          for (idx_e in seq_along(eligible)) {
            v <- eligible[idx_e]
            sv <- scores[, v]
            if (!any(is.finite(sv)) || max(sv, na.rm = TRUE) < min_r2) {
              if (!is.null(spatial_groups) && (v %in% supported) &&
                  !is.null(w_parcel_map[[spatial_groups[v]]])) {
                w_fold[, v] <- w_parcel_map[[spatial_groups[v]]]
                chosen[v] <- which.max(w_fold[, v])
                reason[v] <- "fallback_parcel"
              } else {
                w_fold[, v] <- 0
                w_fold[fallback_kernel, v] <- 1.0
                chosen[v] <- as.integer(fallback_kernel)
                reason[v] <- "fallback_prespecified"
              }
              next
            }
            compat_v <- if (is.null(gate)) rep(TRUE, ncol(kernels)) else gate$compatible[, v]
            if (!any(compat_v)) {
              w_fold[, v] <- 0
              w_fold[fallback_kernel, v] <- 1.0
              chosen[v] <- as.integer(fallback_kernel)
              reason[v] <- "fallback_prespecified"
              next
            }
            omega_v <- if (is.null(gate)) {
              1.0
            } else if (keep[idx_e] && candidates[idx_e] != fallback_kernel) {
              1.0
            } else {
              obj_v <- sv - canonical_strength * distance[, fallback_kernel]
              obj_v[!compat_v] <- -Inf
              k_v <- which.max(obj_v)
              if (length(k_v) == 1L && is.finite(obj_v[k_v]) && k_v != fallback_kernel &&
                  sv[k_v] >= sv[fallback_kernel] && isTRUE(gate$task_supported[k_v, v])) {
                p_sh_v <- gate$shape_p_raw[k_v, v]
                p_tk_v <- gate$task_p_raw[k_v, v]
                if (is.finite(p_sh_v) && is.finite(p_tk_v) && p_tk_v <= 0.05 && p_sh_v <= 0.05) {
                  max(0.0, min(1.0, (-log10(max(p_sh_v, 1e-12)) - 1.3) / 1.2))
                } else 0.0
              } else 0.0
            }
            w_v_own <- hrf_soft_boltzmann_weights(
              scores_vec = sv,
              kernel_cos = kernel_cos,
              dist_fb = distance[, fallback_kernel],
              fallback_kernel = fallback_kernel,
              comp_ok = compat_v,
              omega = omega_v,
              canonical_strength = canonical_strength
            )
            if (!is.null(spatial_groups)) {
              w_p <- w_parcel_map[[spatial_groups[v]]]
              alpha_sp <- if (spatial_strength > 0) {
                min(0.85, spatial_strength / (spatial_strength + 0.02)) *
                  (if (spatial_strength >= 1.0) 1.0 else concordance[idx_e]^2)
              } else if (!keep[idx_e] && sum(w_p[-fallback_kernel]) > 0) {
                0.65
              } else {
                0.0
              }
              w_blend <- (1.0 - alpha_sp) * w_v_own + alpha_sp * w_p
              w_blend <- round(w_blend / 0.02) * 0.02
              if (sum(w_blend) <= 0) {
                w_blend <- rep(0, ncol(kernels))
                w_blend[fallback_kernel] <- 1.0
              } else {
                w_blend <- w_blend / sum(w_blend)
              }
              w_fold[, v] <- w_blend
            } else {
              w_fold[, v] <- w_v_own
            }
            chosen[v] <- which.max(w_fold[, v])
            if (reason[v] == "fallback_prespecified" && sum(w_fold[-fallback_kernel, v]) > 0.05) {
              reason[v] <- if (!is.null(spatial_groups) && spatial_strength > 0) "voxel_soft_spatial" else "voxel_soft"
            } else if (reason[v] == "fallback_parcel" && sum(w_fold[-fallback_kernel, v]) <= 0.05 &&
                       chosen[v] == fallback_kernel && !is.na(anchor[v]) && anchor[v] != fallback_kernel) {
              reason[v] <- "fallback_prespecified"
            }
          }
        }
      }
      if (pooling == "parcel" && selection_mode == "soft_average") {
        if (!is.null(gate) && length(eligible)) {
          incompat_v <- eligible[colSums(gate$compatible[, eligible, drop = FALSE]) == 0L]
          if (length(incompat_v)) {
            w_fold[, incompat_v] <- 0
            w_fold[fallback_kernel, incompat_v] <- 1.0
            reason[incompat_v] <- "fallback_prespecified"
          }
        }
        chosen[eligible] <- max.col(t(w_fold[, eligible, drop = FALSE]), ties.method = "first")
        soft_dep <- eligible[colSums(w_fold[-fallback_kernel, eligible, drop = FALSE]) > 0.05 &
          reason[eligible] == "fallback_prespecified"]
        if (length(soft_dep)) reason[soft_dep] <- "pooled_parcel"
      }
      if (selection_mode == "soft_average") {
        chosen <- max.col(t(w_fold), ties.method = "first")
      }
      selected_score[eligible] <- scores[cbind(chosen[eligible], eligible)]
    }
    if (selection_mode == "hard") {
      w_fold[] <- 0
      w_fold[cbind(chosen, seq_along(voxel_ids))] <- 1.0
    }
    k_mat <- kernels %*% w_fold
    k_divs <- if (kernel_normalize == "positive_lobe") {
      colSums(pmax(k_mat, 0))
    } else if (kernel_normalize == "peak") {
      apply(k_mat, 2L, function(col) {
        pk <- max(col)
        if (pk <= .Machine$double.eps) max(abs(col)) else pk
      })
    } else {
      colSums(k_mat)
    }
    bad_sum <- !is.finite(k_divs) | abs(k_divs) < 1e-6
    if (any(bad_sum)) {
      k_mat[, bad_sum] <- kernels[, fallback_kernel]
      k_divs[bad_sum] <- 1.0
    }
    k_mat <- sweep(k_mat, 2L, k_divs, "/")
    colnames(k_mat) <- voxel_ids
    kernel_weights[[target]] <- w_fold
    selected_kernels[[target]] <- k_mat
    selections[[target]] <- data.frame(target_run = names(bold)[target], voxel = voxel_ids,
      parcel = if (is.null(spatial_groups)) NA_character_ else spatial_groups,
      kernel_index = chosen, kernel = candidate_ids[chosen], status = reason,
      n_training_runs = counts, n_parcel_voxels = parcel_n, parcel_kernel_index = anchor,
      best_kernel_index = best, best_adjusted_partial_r2 = best_score,
      selected_adjusted_partial_r2 = selected_score, score_margin = margin,
      voxel_evidence = ifelse(counts < min_training_runs, "insufficient_training",
        ifelse(best_score < min_r2, "low_task_support", "supported")))
    if (!is.null(gate)) {
      sel <- selections[[target]]
      idx <- cbind(sel$kernel_index, seq_along(voxel_ids))
      sel$selected_task_p <- gate$task_p[idx]
      sel$selected_shape_p <- gate$shape_p[idx]
      sel$n_compatible_shapes <- colSums(gate$compatible)
      sel$n_supported_shapes <- colSums(gate$supported)
      sel$stability <- ifelse(counts < 2L, "not_assessable", ifelse(
        sel$n_compatible_shapes == 0L, "no_common_library_shape", "compatible_training_runs"))
      own_support <- gate$supported[idx]
      sel$voxel_evidence <- ifelse(counts < min_training_runs, "insufficient_training",
        ifelse(sel$n_compatible_shapes == 0L, "inconsistent_or_misspecified",
          ifelse(!own_support, "unidentified",
            ifelse(sel$kernel_index == fallback_kernel, "prespecified_shape", "supported_departure"))))
      selections[[target]] <- sel
    }
    if (return_scores) candidate_scores[[target]] <- scores
  }
  # Post-selection target-excluded noise refit:
  # All HRF selections, weights, parcel anchors, and evidence scores are now
  # strictly frozen from training runs (j != r). Refitting each run's noise model
  # using its frozen cross-run selected_kernels[[r]] removes non-canonical HRF
  # task leakage from the REML covariance and GLS scale without allowing any
  # target-run information to enter HRF selection.
  for (r in seq_along(bold)) {
    mask_r <- !is.na(valid_rows[[r]]) & valid_rows[[r]]
    if (sum(mask_r) < 12L) next
    sel_k_r <- selected_kernels[[r]]
    if (is.null(sel_k_r)) next
    diff_fb <- max(abs(sweep(sel_k_r, 1L, kernels[, fallback_kernel], "-")))
    if (!is.finite(diff_fb) || diff_fb <= 1e-6) next
    z_r <- cbind(intercept = 1, trend = seq(-1, 1, length.out = nrow(bold[[r]])), nuisance[[r]])
    fixed_r <- if (noise_model == "ols") 0 else if (!is.null(noise_rho)) noise_rho[r] else NULL
    noise[[r]] <- hrf_estimate_noise(bold[[r]], events[[r]], tr, offsets[r],
      kernels, fallback_kernel, z_r, mask_r, eff_noise_groups, noise_grid, batch_size,
      fixed_r, selected_kernels = sel_k_r,
      spatial_coords = spatial_coords, spatial_fwhm = spatial_fwhm)
  }
  structure(list(kernels = kernels, parameters = parameters, selections = selections,
    kernel_weights = kernel_weights, selected_kernels = selected_kernels,
    parcels = anchors, training_runs = training_runs,
    signal_status = do.call(rbind, event_manifest), designs = do.call(rbind, design_manifest),
    scores = if (return_scores) candidate_scores else NULL,
    noise = noise, run_evidence = if (return_scores && evidence == "calibrated") run_results else NULL,
    run_ids = names(bold),
    voxel_ids = voxel_ids, n_time = vapply(bold, nrow, integer(1)),
    specification = list(events = events, nuisance = nuisance, valid_rows = valid_rows,
      spatial_groups = spatial_groups, noise_groups = eff_noise_groups),
    settings = list(tr = tr, time_offset = offsets, pooling = pooling, amplitude_model = amplitude_model,
      spatial_strength = spatial_strength, fallback_kernel = as.integer(fallback_kernel),
      min_r2 = min_r2, min_training_runs = min_training_runs, min_parcel_voxels = min_parcel_voxels,
      batch_size = as.integer(batch_size), score = "mean_adjusted_partial_r2",
      rank_tolerance = 1e-8, evidence = evidence, evidence_alpha = evidence_alpha,
      canonical_strength = canonical_strength, selection_mode = selection_mode,
      kernel_normalize = kernel_normalize,
      noise_grid = noise_grid,
      noise_pooling = if (!is.null(noise_groups)) "supplied_noise_groups" else
        if (is.null(spatial_groups)) "run" else "supplied_parcel",
      nuisance = "run_specific_intercept_linear_trend_and_supplied", noise_model = noise_model,
      selection = if (is_single_run) "within_run_subblock_pooled" else "leave_target_run_out",
      parcel_fallback = !is.null(spatial_groups),
      spatial_regularization = if (pooling == "parcel") "shared_parcel_hrf" else
        if (spatial_strength > 0) "discrete_parcel_shape_penalty" else "none")),
    class = "crossrun_hrf")
}

#' Deconvolve runs using frozen HRFs learned from other runs
#'
#' @param bold Named time-by-voxel run matrices with IDs/order/dimensions matching
#'   the HRF fit. These should be the same preprocessed BOLD used for estimation.
#' @param hrf_fit Object returned by estimate_crossrun_hrf().
#' @param lambda Nonnegative regularization strength. If \code{NULL} (the
#'   default for modern \code{ar1}/\code{calibrated} fits), selected
#'   automatically from \code{lambda_grid} by gain-calibrated, task-weighted
#'   adaptive prewhitened generalized cross-validation using \code{gcv_rule}.
#' @param scaling Input scaling: zscore (default), noise_sd, demean, none, or psc.
#' @param penalty Quadratic penalty: diff2, diff1, ridge, matern1, or sobolev.
#' @param backend Fixed linear inverse backend: cpp, r, or scalar wiener.
#' @param batch_size Maximum voxel columns per scaling/multiplication block.
#' @param tail_padding Number of appended zero BOLD observations, default zero.
#'   Use nrow(hrf_fit$kernels) to reproduce the current voxelwise wrapper assumption.
#' @param ridge_floor Numerical diagonal ridge for cpp/r.
#' @param operator_cache_mb Maximum inverse-operator cache size in MiB, default
#'   128. Zero disables reuse across runs. A single operator may exceed the cache
#'   budget while it is being applied, but is not retained in the cache.
#' @param nuisance_mode Jointly fit unpenalized intercept, trend, and supplied
#'   nuisance terms, or ignore them (default, preserving prepared input signals).
#'   Joint projection is applied to both data and design, but can amplify weakly
#'   identified low-frequency activity in a flexible continuous inverse. Validate
#'   it for the intended design rather than assuming it improves recovery.
#' @param observation_mask Use fit's valid_rows in the inverse (fit), or all
#'   observations (all). Defaults to fit for AR/calibrated, all for OLS/legacy.
#' @param lambda_grid Candidate regularization values evaluated when
#'   \code{lambda = NULL}.
#' @param gcv_rule Selection rule (\code{"min"} or \code{"1se"}) used when
#'   \code{lambda = NULL}.
#' @param rwls Robust weighted least-squares frame weighting mode for transient
#'   temporal artifacts (e.g., sub-censoring-threshold motion glitches or RF spikes):
#'   \code{"auto"} (default; enables \code{"woodbury"} on modern AR/calibrated fits and
#'   falls back to \code{"none"} on legacy OLS or \code{backend = "wiener"}),
#'   \code{"none"} (or \code{FALSE}; unweighted GLS),
#'   \code{"woodbury"} (or \code{TRUE}; additive-transient Sparse Bayesian Learning
#'   Woodbury covariance \eqn{\Sigma_{\mathrm{AR}} + \mathrm{diag}(\delta^2)}), or
#'   \code{"innovation"} (Diedrichsen-style innovation-space downweighting).
#' @return A list with activity matrices, per-run scaling/status tables, the HRF
#'   fit, frozen inverse settings, standardized noise innovations, and per-run
#'   observation \code{frame_weights}. Activity matrices retain voxel IDs and carry
#'   tr and time_offset attributes for subsequent event alignment.
#' @details One inverse operator is constructed per selected HRF/run length,
#'   cached for consecutive runs with identical masks/nuisance designs, and applied using block
#'   matrix products for cpp/r. Scalar Wiener uses direct block FFT fits.
#'   HRFs are frozen; no voxelwise factorization or R fitting loop is required.
#'   AR/calibrated fits use a masked GLS inverse with the frozen run-specific
#'   noise covariance. A gap of d acquisitions uses rho^d and the corresponding
#'   accumulated innovation variance. Censored rows contribute no observations;
#'   latent activity is still estimated on the complete original grid. Its
#'   observation_mask attribute identifies modeled versus observed samples.
#'   Nonfinite/constant signals on usable observations remain NA. Nonfinite
#'   censored observations do not invalidate a run. Joint nuisance fitting avoids
#'   projecting only the BOLD and inadvertently changing the task model.
#'   Wiener requires OLS, nuisance_mode = "ignore", observation_mask = "all".
#'   Scaling is applied once before optional padding, and its centers/divisors are
#'   returned. Different run z-scores do not preserve absolute between-run units.
#' @export
deconvolve_crossrun_hrf <- function(bold, hrf_fit, lambda = NULL, scaling = "zscore",
    penalty = c("diff2", "diff1", "ridge", "matern1", "sobolev"), backend = c("cpp", "r", "wiener"),
    batch_size = 2048L, tail_padding = 0L, ridge_floor = 1e-8, operator_cache_mb = 128,
    nuisance_mode = NULL, observation_mask = NULL,
    lambda_grid = 10^seq(-2, 1.0, length.out = 16), gcv_rule = c("min", "1se"),
    rwls = "auto") {
  hrf_validate_bold(bold)
  if (!inherits(hrf_fit, "crossrun_hrf")) stop("hrf_fit must be a crossrun_hrf object")
  if (!identical(names(bold), hrf_fit$run_ids) || !identical(colnames(bold[[1L]]), hrf_fit$voxel_ids) ||
      !identical(vapply(bold, nrow, integer(1)), hrf_fit$n_time))
    stop("BOLD run/voxel IDs and dimensions must match the HRF fit")
  modern <- identical(hrf_fit$settings$noise_model, "ar1") || identical(hrf_fit$settings$evidence, "calibrated")
  if (missing(lambda) && !modern) stop("Specify lambda explicitly (or pass lambda = NULL for GCV tuning)")
  checkmate::assert_number(lambda, lower = 0, finite = TRUE, null.ok = TRUE)
  gcv_rule <- match.arg(gcv_rule)
  scaling <- match.arg(scaling, c("zscore", "demean", "none", "psc", "noise_sd"))
  penalty <- match.arg(penalty); backend <- match.arg(backend)
  checkmate::assert_integerish(batch_size, lower = 1, len = 1L)
  checkmate::assert_integerish(tail_padding, lower = 0, len = 1L)
  checkmate::assert_number(ridge_floor, lower = 0, finite = TRUE)
  checkmate::assert_number(operator_cache_mb, lower = 0, finite = TRUE)
  if (is.null(nuisance_mode)) nuisance_mode <- "ignore"
  if (is.null(observation_mask)) observation_mask <- if (modern) "fit" else "all"
  nuisance_mode <- match.arg(nuisance_mode, c("joint", "ignore"))
  observation_mask <- match.arg(observation_mask, c("fit", "all"))
  if (is.logical(rwls)) {
    checkmate::assert_flag(rwls)
    rwls_mode <- if (isTRUE(rwls)) "woodbury" else "none"
  } else if (is.character(rwls) && length(rwls) == 1L && !is.na(rwls)) {
    rwls_mode <- match.arg(rwls, c("auto", "none", "woodbury", "innovation"))
  } else {
    stop("rwls must be a logical flag or one of 'auto', 'none', 'woodbury', 'innovation'")
  }
  if (backend == "wiener" && (identical(hrf_fit$settings$noise_model, "ar1") ||
      nuisance_mode != "ignore" || observation_mask != "all"))
    stop("Wiener does not implement masked, nuisance-adjusted AR GLS; use cpp/r, or an OLS fit with nuisance_mode = 'ignore', observation_mask = 'all'")
  if (backend == "wiener" && rwls_mode %in% c("woodbury", "innovation")) {
    stop("Wiener backend does not implement rWLS frame weighting; use backend = 'cpp' or 'r'")
  }
  if (rwls_mode == "auto" && (!modern || backend == "wiener")) {
    rwls_mode <- "none"
  }
  rwls_space <- if (rwls_mode == "innovation") "innovation" else "woodbury"

  # Per-run noise SD, standardized innovations, and optional global rWLS frame weights
  # pre-pass (computed once before GCV lambda tuning so tuning also benefits from
  # precomputed GLS noise_sd and artifact downweighting without redundant QR passes).
  run_prep <- vector("list", length(bold))
  frame_weights <- precomputed_nsd <- vector("list", length(bold))
  names(frame_weights) <- names(precomputed_nsd) <- names(bold)
  for (r in seq_along(bold)) {
    y <- bold[[r]]; n <- nrow(y); n_vox <- ncol(y)
    mask <- if (observation_mask == "fit") hrf_fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    strat <- hrf_stratum_kernels(hrf_fit, r)
    rho <- strat$rho
    status <- if (length(rows) >= 3L) decon_signal_status(y[rows, , drop = FALSE]) else
      rep("insufficient_observations", n_vox)
    if (identical(hrf_fit$settings$noise_model, "ar1") && !is.null(hrf_fit$noise)) {
      unknown <- grepl("^insufficient_noise", hrf_fit$noise[[r]]$method)
      status[status == "ok" & unknown] <- "noise_unidentified"
    }
    ok_cols <- which(status == "ok")
    store_innov <- length(rows) >= 3L && n_vox <= 16384L
    z_innov_r <- if (store_innov) {
      matrix(0.0, length(rows), n_vox, dimnames = list(NULL, colnames(y)))
    } else NULL
    coh_nsd_r <- rep(NA_real_, n_vox)

    need_rwls <- (rwls_mode != "none") && length(rows) >= 16L && length(ok_cols) >= 1L
    rwls_cols <- integer(0L)
    rwls_map <- rep(NA_integer_, n_vox)
    z_raw_r <- r_obs_r <- NULL
    if (need_rwls) {
      if (length(ok_cols) <= 4096L) {
        rwls_cols <- ok_cols
      } else {
        r2_r <- if (!is.null(hrf_fit$selections[[r]]$r2)) {
          as.numeric(hrf_fit$selections[[r]]$r2[ok_cols])
        } else {
          rep(1.0, length(ok_cols))
        }
        r2_r[!is.finite(r2_r)] <- -Inf
        thr_r2 <- sort(r2_r, decreasing = TRUE)[4096L]
        rwls_cols <- ok_cols[which(r2_r >= thr_r2)]
        if (length(rwls_cols) > 4096L) rwls_cols <- rwls_cols[seq_len(4096L)]
      }
      rwls_map[rwls_cols] <- seq_along(rwls_cols)
      z_raw_r <- matrix(NA_real_, length(rows), length(rwls_cols))
      r_obs_r <- matrix(NA_real_, length(rows), length(rwls_cols))
    }

    strata <- strat$keys
    task_env_accum <- rep(0.0, length(rows))
    n_task_strata <- 0L
    diff1_cols <- function(mat) {
      if (ncol(mat) == 0L || nrow(mat) < 2L) return(mat)
      rbind(mat[2:nrow(mat), , drop = FALSE], 0.0) -
        rbind(0.0, mat[1:(nrow(mat) - 1L), , drop = FALSE])
    }
    for (key in unique(strata[status == "ok"])) {
      columns <- which(strata == key & status == "ok")
      ar <- rho[columns[1L]]
      ar_op <- if (!is.null(strat$get_ar)) strat$get_ar(columns[1L]) else ar
      k_vec <- as.numeric(strat$get_kernel(columns[1L]))
      x_gls_stratum <- NULL
      x_trial_innov <- NULL
      x_trial_raw <- NULL
      if (!is.null(hrf_fit$specification$events[[r]]) && length(rows) >= 20L) {
        mask_decon_bool <- seq_len(n) %in% rows
        nuis_r <- if (!is.null(hrf_fit$specification$nuisance)) hrf_fit$specification$nuisance[[r]] else matrix(numeric(), n, 0L)
        z_r_base <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n), nuis_r)
        off_r <- if (length(hrf_fit$settings$time_offset) >= r) hrf_fit$settings$time_offset[r] else 0.0
        md_str <- hrf_build_noise_mean_design(n, hrf_fit$specification$events[[r]],
          hrf_fit$settings$tr, off_r, k_vec, z_r_base, mask_decon_bool)
        if (md_str$df_nominal >= 16L) {
          x_gls_stratum <- md_str$design[rows, , drop = FALSE]
        }
        if (!is.null(md_str$trial_design) && ncol(md_str$trial_design) > 0L) {
          td_full <- md_str$trial_design
          td_d1 <- diff1_cols(td_full)
          td_d2 <- diff1_cols(td_d1)
          task_env_accum <- task_env_accum +
            rowSums(td_full[rows, , drop = FALSE]^2) +
            0.5 * rowSums(td_d1[rows, , drop = FALSE]^2) +
            0.25 * rowSums(td_d2[rows, , drop = FALSE]^2)
          n_task_strata <- n_task_strata + 1L
          if (!isTRUE(md_str$used_trial_design)) {
            x_trial_innov <- cbind(td_full, td_d1)[rows, , drop = FALSE]
          }
          x_trial_raw <- if (isTRUE(md_str$used_trial_design)) {
            cbind(td_d1, td_d2)[rows, , drop = FALSE]
          } else {
            x_trial_innov
          }
        }
      }
      for (indices in split(columns, ceiling(seq_along(columns) / batch_size))) {
        y_blk <- y[rows, indices, drop = FALSE]
        idx_rw <- if (need_rwls) which(!is.na(rwls_map[indices])) else integer(0L)
        do_raw_blk <- length(idx_rw) > 0L
        gls_res <- hrf_estimate_gls_noise_sd(
          y_blk, x_gls_stratum, rows, ar_op,
          return_innovations = store_innov && !do_raw_blk,
          return_raw_innovations = do_raw_blk,
          trial_rows = if (do_raw_blk) x_trial_raw else x_trial_innov
        )
        coh_nsd <- if (is.list(gls_res)) gls_res$sd else gls_res
        coh_nsd_r[indices] <- coh_nsd
        if (store_innov && !do_raw_blk && is.list(gls_res) && !is.null(gls_res$z_unit)) {
          z_innov_r[, indices] <- gls_res$z_unit
        }
        if (do_raw_blk && is.list(gls_res) && !is.null(gls_res$z_raw)) {
          dst <- rwls_map[indices[idx_rw]]
          z_raw_r[, dst] <- gls_res$z_raw[, idx_rw, drop = FALSE]
          r_obs_r[, dst] <- gls_res$r_obs[, idx_rw, drop = FALSE]
        }
      }
    }

    w_r <- rep(1.0, length(rows))
    if (need_rwls && !is.null(z_raw_r) && ncol(z_raw_r) >= 1L) {
      med_or0 <- function(vec) {
        if (is.null(vec)) return(0.0)
        val <- stats::median(vec[rwls_cols], na.rm = TRUE)
        if (is.finite(val)) as.numeric(val) else 0.0
      }
      rho_med <- med_or0(rho)
      p1_med <- med_or0(strat$psi1)
      p2_med <- med_or0(strat$psi2)
      p3_med <- med_or0(strat$psi3)
      p4_med <- med_or0(strat$psi4)
      ar_rep_r <- if (abs(p3_med) > 1e-10 || abs(p4_med) > 1e-10) {
        c(rho_med, p1_med, p2_med, p3_med, p4_med)
      } else if (abs(p1_med) > 1e-10 || abs(p2_med) > 1e-10) {
        c(rho_med, p1_med, p2_med)
      } else {
        rho_med
      }
      w_r <- hrf_estimate_rwls_weights(
        z_raw = z_raw_r,
        r_obs = r_obs_r,
        rows = rows,
        ar_vec = ar_rep_r,
        mode = rwls_mode,
        task_env = if (n_task_strata > 0L) (task_env_accum / n_task_strata) else NULL
      )
    }

    # Refresh `coh_nsd_r` and `z_innov_r` in the exact rWLS-whitened GLS space
    # (and without `x_trial_raw` derivative projection when `used_trial_design == TRUE`)
    # so downstream parcel/cross-parcel error covariances (`M-MOMS` shared_error) match
    # the exact post-rWLS residual geometry without unweighted spike-projection smearing.
    if (need_rwls) {
      w_active_r <- any(abs(w_r - 1.0) > 1e-6)
      w_obs_arg <- if (w_active_r) w_r else NULL
      for (key in unique(strata[status == "ok"])) {
        columns <- which(strata == key & status == "ok")
        if (!w_active_r && !any(!is.na(rwls_map[columns]))) next
        ar <- rho[columns[1L]]
        ar_op <- if (!is.null(strat$get_ar)) strat$get_ar(columns[1L]) else ar
        k_vec <- as.numeric(strat$get_kernel(columns[1L]))
        x_gls_stratum <- NULL
        x_trial_innov <- NULL
        if (!is.null(hrf_fit$specification$events[[r]]) && length(rows) >= 20L) {
          mask_decon_bool <- seq_len(n) %in% rows
          nuis_r <- if (!is.null(hrf_fit$specification$nuisance)) hrf_fit$specification$nuisance[[r]] else matrix(numeric(), n, 0L)
          z_r_base <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n), nuis_r)
          off_r <- if (length(hrf_fit$settings$time_offset) >= r) hrf_fit$settings$time_offset[r] else 0.0
          md_str <- hrf_build_noise_mean_design(n, hrf_fit$specification$events[[r]],
            hrf_fit$settings$tr, off_r, k_vec, z_r_base, mask_decon_bool)
          if (md_str$df_nominal >= 16L) {
            x_gls_stratum <- md_str$design[rows, , drop = FALSE]
          }
          if (!isTRUE(md_str$used_trial_design) && !is.null(md_str$trial_design) && ncol(md_str$trial_design) > 0L) {
            td_full <- md_str$trial_design
            td_d1 <- diff1_cols(td_full)
            x_trial_innov <- cbind(td_full, td_d1)[rows, , drop = FALSE]
          }
        }
        for (indices in split(columns, ceiling(seq_along(columns) / batch_size))) {
          y_blk <- y[rows, indices, drop = FALSE]
          gls_clean <- hrf_estimate_gls_noise_sd(
            y_blk, x_gls_stratum, rows, ar_op,
            return_innovations = store_innov,
            return_raw_innovations = FALSE,
            trial_rows = x_trial_innov,
            obs_weights = w_obs_arg
          )
          if (!w_active_r) {
            coh_nsd_r[indices] <- if (is.list(gls_clean)) gls_clean$sd else gls_clean
          }
          if (store_innov && is.list(gls_clean) && !is.null(gls_clean$z_unit)) {
            z_innov_r[, indices] <- gls_clean$z_unit
          }
        }
      }
      if (store_innov && w_active_r && !is.null(z_innov_r)) {
        attr(z_innov_r, "rwls_whitened") <- TRUE
      }
    }
    frame_weights[[r]] <- w_r
    precomputed_nsd[[r]] <- coh_nsd_r
    run_prep[[r]] <- list(
      rows = rows,
      strat = strat,
      rho = rho,
      status = status,
      coh_nsd = coh_nsd_r,
      z_innov = z_innov_r,
      w_r = w_r
    )
  }

  # Automatic gain-calibrated, task-weighted prewhitened GCV tuning when lambda = NULL:
  lambda_tuning <- NULL
  if (is.null(lambda)) {
    checkmate::assert_numeric(lambda_grid, lower = 0, any.missing = FALSE, min.len = 1L, finite = TRUE)
    lambda_grid <- sort(unique(as.numeric(lambda_grid)))
    tuned <- hrf_tune_decon_lambda(
      bold = bold,
      hrf_fit = hrf_fit,
      lambda_grid = lambda_grid,
      scaling = scaling,
      penalty = penalty,
      ridge_floor = ridge_floor,
      observation_mask = observation_mask,
      nuisance_mode = nuisance_mode,
      tail_padding = tail_padding,
      backend = backend,
      gcv_rule = gcv_rule,
      frame_weights = frame_weights,
      precomputed_nsd = precomputed_nsd
    )
    lambda <- tuned$lambda
    lambda_tuning <- tuned$tuning
  }

  settings <- list(lambda = lambda, lambda_tuning = lambda_tuning, scaling = scaling, penalty = penalty,
    backend = backend, tail_padding = as.integer(tail_padding), ridge_floor = ridge_floor,
    batch_size = as.integer(batch_size), operator_cache_mb = operator_cache_mb,
    hrf_selection = "frozen_crossrun", nuisance_mode = nuisance_mode,
    observation_mask = observation_mask, noise_model = hrf_fit$settings$noise_model,
    rwls = rwls_mode, rwls_space = rwls_space)
  activity <- scaling_tables <- noise_innovations <- vector("list", length(bold))
  names(activity) <- names(scaling_tables) <- names(noise_innovations) <- names(bold)
  operators <- new.env(parent = emptyenv())
  cached_design <- NULL
  cache_keys <- character(); cache_bytes <- 0
  for (r in seq_along(bold)) {
    y <- bold[[r]]; n <- nrow(y)
    prep_r <- run_prep[[r]]
    rows <- prep_r$rows
    strat <- prep_r$strat
    rho <- prep_r$rho
    status <- prep_r$status
    coh_nsd_r <- prep_r$coh_nsd
    z_innov_r <- prep_r$z_innov
    w_r <- prep_r$w_r
    w_active <- length(w_r) == length(rows) && any(abs(w_r - 1.0) > 1e-8)
    w_op <- if (w_active) w_r else NULL

    z <- if (nuisance_mode == "joint") cbind(1, seq(-1, 1, length.out = n),
      hrf_fit$specification$nuisance[[r]]) else matrix(numeric(), n, 0L)
    design <- list(n = n, rows = rows, nuisance = z, obs_weights = w_op)
    if (!identical(design, cached_design)) {
      rm(list = ls(operators), envir = operators)
      cache_keys <- character(); cache_bytes <- 0
      cached_design <- design
    }
    out <- matrix(NA_real_, n, ncol(y), dimnames = dimnames(y))
    centers <- divisors <- noise_sds <- rep(NA_real_, ncol(y))
    # Group voxels into strata that share the exact same (kernel / weight vector, rho) pair.
    # Because the regularized GLS inverse operator M = (H' W' Q_Z W H + lambda P)^-1 H' W' Q_Z W
    # depends only on (n, kernel, valid_rows, rho, nuisance, lambda, penalty, obs_weights),
    # all voxels in a stratum share a single n x m_obs inverse operator matrix,
    # turning voxelwise deconvolution into fast BLAS matrix multiplication (M %*% Y_block)
    # and caching operators across runs when observation masks, nuisance, and weights match.
    strata <- strat$keys
    for (key in unique(strata[status == "ok"])) {
      columns <- which(strata == key & status == "ok")
      ar <- rho[columns[1L]]
      ar_op <- if (!is.null(strat$get_ar)) strat$get_ar(columns[1L]) else ar
      ar_psi1 <- if (!is.null(strat$psi1)) strat$psi1[columns[1L]] else 0.0
      ar_psi2 <- if (!is.null(strat$psi2)) strat$psi2[columns[1L]] else 0.0
      ar_psi3 <- if (!is.null(strat$psi3)) strat$psi3[columns[1L]] else 0.0
      ar_psi4 <- if (!is.null(strat$psi4)) strat$psi4[columns[1L]] else 0.0
      k_vec <- as.numeric(strat$get_kernel(columns[1L]))
      if (backend != "wiener") {
        if (exists(key, envir = operators, inherits = FALSE)) operator <- get(key, envir = operators) else {
          operator <- hrf_inverse_operator(n, k_vec, rows, ar_op, z, settings, obs_weights = w_op)
          bytes <- as.numeric(utils::object.size(operator))
          if (bytes <= operator_cache_mb * 1024^2) {
            while (length(cache_keys) && cache_bytes + bytes > operator_cache_mb * 1024^2) {
              cache_bytes <- cache_bytes - as.numeric(utils::object.size(get(cache_keys[1L], envir = operators)))
              rm(list = cache_keys[1L], envir = operators); cache_keys <- cache_keys[-1L]
            }
            assign(key, operator, envir = operators)
            cache_keys <- c(cache_keys, key); cache_bytes <- cache_bytes + bytes
          }
        }
      }
      for (indices in split(columns, ceiling(seq_along(columns) / batch_size))) {
        y_blk <- y[rows, indices, drop = FALSE]
        coh_nsd <- coh_nsd_r[indices]
        scaled <- scale_decon_input(
          y_blk,
          scaling = scaling,
          noise_rho = ar,
          noise_psi1 = ar_psi1,
          noise_psi2 = ar_psi2,
          noise_psi3 = ar_psi3,
          noise_psi4 = ar_psi4,
          noise_sd_override = if (identical(scaling, "noise_sd")) coh_nsd else NULL,
          rows = rows,
          obs_weights = w_op
        )
        out[, indices] <- if (backend == "wiener") {
          fit_input <- rbind(scaled$y, matrix(0, tail_padding, length(indices)))
          deconvolve_reglin(fit_input, k_vec, lambda = lambda,
            scaling = "none", backend = "wiener",
            kernel_normalize = if (!is.null(hrf_fit$settings$kernel_normalize)) hrf_fit$settings$kernel_normalize else "signed_sum")[seq_len(n), , drop = FALSE]
        } else operator %*% scaled$y
        centers[indices] <- scaled$center; divisors[indices] <- scaled$scale
        noise_sds[indices] <- if (identical(scaling, "noise_sd")) scaled$scale else coh_nsd
      }
    }
    attr(out, "tr") <- hrf_fit$settings$tr
    attr(out, "time_offset") <- hrf_fit$settings$time_offset[r]
    attr(out, "observation_mask") <- seq_len(n) %in% rows
    activity[[r]] <- out
    noise_innovations[[r]] <- z_innov_r
    scaling_tables[[r]] <- data.frame(run = names(bold)[r], voxel = colnames(y),
      status = status, center = centers, scale = divisors, noise_sd = noise_sds,
      scaling = scaling, noise_rho = rho, n_observations = length(rows))
  }
  list(
    activity = activity,
    scaling = scaling_tables,
    hrf_fit = hrf_fit,
    settings = settings,
    noise_innovations = noise_innovations,
    frame_weights = frame_weights
  )
}
