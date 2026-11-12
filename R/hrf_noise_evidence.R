# Helpers for the cross-run interface. Covariance parameters may use their own
# run's observations; shape selection and its evidence use OTHER runs only.

# Convert cascaded AR parameters c(rho1, psi_1, ..., psi_q) representing
# (1 - rho1 L)(1 - sum_{j=1}^q psi_j L^j) = 1 - sum_{k=1}^p phi_k L^k
# into stationary causal AR(p) coefficients, unit-marginal ACF, and exact
# innovation-to-marginal variance ratio sigma_u^2 / sigma_e^2.
hrf_ar_cascade_info <- function(ar_vec, max_lag = 16L) {
  max_lag <- max(1L, as.integer(max_lag))
  if (!length(ar_vec)) {
    return(list(phi = numeric(0), acf = c(1.0, rep(0.0, max_lag)), sig_u2 = 1.0))
  }
  rho1 <- if (is.finite(ar_vec[1L])) min(max(ar_vec[1L], -0.95), 0.95) else 0.0
  psi <- if (length(ar_vec) >= 2L) as.numeric(ar_vec[-1L]) else numeric(0)
  psi[!is.finite(psi)] <- 0.0
  while (length(psi) > 0L && abs(psi[length(psi)]) <= 1e-10) {
    psi <- psi[-length(psi)]
  }
  if (!length(psi)) {
    if (abs(rho1) <= 1e-10) {
      return(list(phi = numeric(0), acf = c(1.0, rep(0.0, max_lag)), sig_u2 = 1.0))
    }
    acf_v <- rho1^(0:max_lag)
    return(list(phi = rho1, acf = acf_v, sig_u2 = max(1e-6, 1.0 - rho1 * rho1)))
  }
  # Enforce strict causal stationarity of the residual psi polynomial
  roots_psi <- tryCatch(Mod(polyroot(c(1.0, -psi))), error = function(e) rep(0.0, length(psi)))
  if (any(roots_psi <= 1.02)) {
    rad <- 0.97 / max(1e-6, max(1.0 / pmax(roots_psi, 1e-6)))
    psi <- psi * (rad^seq_along(psi))
  }
  q <- length(psi)
  phi <- numeric(q + 1L)
  phi[1L] <- rho1 + psi[1L]
  if (q >= 2L) {
    for (k in 2L:q) {
      phi[k] <- psi[k] - rho1 * psi[k - 1L]
    }
  }
  phi[q + 1L] <- -rho1 * psi[q]
  while (length(phi) > 0L && abs(phi[length(phi)]) <= 1e-12) {
    phi <- phi[-length(phi)]
  }
  p <- length(phi)
  if (p == 0L) {
    return(list(phi = numeric(0), acf = c(1.0, rep(0.0, max_lag)), sig_u2 = 1.0))
  }
  roots_phi <- tryCatch(Mod(polyroot(c(1.0, -phi))), error = function(e) rep(0.0, p))
  if (any(roots_phi <= 1.015)) {
    rad <- 0.975 / max(1e-6, max(1.0 / pmax(roots_phi, 1e-6)))
    phi <- phi * (rad^seq_len(p))
  }
  lag_need <- max(max_lag, p)
  acf_v <- as.numeric(stats::ARMAacf(ar = phi, lag.max = lag_need))
  sig_u2 <- max(1e-6, 1.0 - sum(phi * acf_v[2L:(p + 1L)]))
  list(phi = phi, acf = acf_v[seq_len(max_lag + 1L)], sig_u2 = sig_u2)
}

# Build the nuisance + free-trial + condition-derivative + DF-adaptive DCT (< 1/128 Hz)
# mean design matrix for Step 1 noise estimation and GLS residual scale estimation.
hrf_build_noise_mean_design <- function(n_obs, events, tr, offset, kernel_col, z, mask) {
  n_mask <- sum(mask)
  rank_fn <- function(a) {
    if (ncol(a) == 0L || n_mask == 0L) return(0L)
    a <- a[mask, , drop = FALSE]
    norm <- sqrt(colSums(a^2)); norm[norm == 0] <- 1
    s <- svd(sweep(a, 2L, norm, "/"), nu = 0, nv = 0)$d
    if (length(s) && s[1L] > 0) sum(s > 1e-8 * s[1L]) else 0L
  }
  if (is.null(events) || nrow(events) == 0L) {
    return(list(
      design = z,
      df_nominal = n_mask - rank_fn(z),
      trial_design = matrix(0, n_obs, 0L),
      used_trial_design = TRUE
    ))
  }
  durations <- if (is.null(events$duration)) rep(0, nrow(events)) else events$duration
  k_mat <- matrix(as.numeric(kernel_col), ncol = 1L)
  x <- hrf_event_design_cpp(n_obs, tr, offset, events$onset, durations,
    seq_len(nrow(events)), nrow(events), k_mat)[[1L]]
  cond_levels <- unique(events$condition)
  cond_ids <- match(events$condition, cond_levels)
  x_cond <- if (length(cond_levels) > 0L) {
    hrf_event_design_cpp(n_obs, tr, offset, events$onset, durations,
      cond_ids, length(cond_levels), k_mat)[[1L]]
  } else {
    matrix(0, n_obs, 0L)
  }
  derivative <- function(a) {
    if (ncol(a) == 0L || nrow(a) < 2L) return(a)
    rbind(a[2:nrow(a), , drop = FALSE], 0) - rbind(0, a[1:(nrow(a) - 1L), , drop = FALSE])
  }
  z_dct <- z
  if (n_obs * tr >= 96 && n_mask >= 40L) {
    cap_cos <- max(6L, min(24L, as.integer(floor(0.15 * n_mask))))
    n_cos <- max(2L, min(cap_cos, as.integer(floor(2 * n_obs * tr / 128.0))))
    t_idx <- seq(0, n_obs - 1L)
    dct_cols <- vapply(seq_len(n_cos), function(k) {
      cos(pi * (2 * t_idx + 1) * k / (2 * n_obs))
    }, numeric(n_obs))
    z_dct <- cbind(z, dct_cols)
  }
  aug_trial_cond <- cbind(z_dct, x, derivative(x_cond), derivative(derivative(x_cond)))
  aug_trial <- cbind(z_dct, x)
  aug_cond <- cbind(z_dct, x_cond, derivative(x_cond), derivative(derivative(x_cond)))
  min_df_full <- max(36L, as.integer(0.68 * n_mask))
  used_trial <- FALSE
  mean_design <- if (n_mask - rank_fn(aug_trial_cond) >= min_df_full) {
    used_trial <- TRUE
    aug_trial_cond
  } else if (n_mask - rank_fn(aug_trial) >= min_df_full) {
    used_trial <- TRUE
    aug_trial
  } else if (n_mask - rank_fn(cbind(z, x)) >= 20L) {
    used_trial <- TRUE
    cbind(z, x)
  } else if (n_mask - rank_fn(aug_cond) >= 20L) {
    aug_cond
  } else if (n_mask - rank_fn(cbind(z, x_cond)) >= 12L) {
    cbind(z, x_cond)
  } else {
    z
  }
  list(
    design = mean_design,
    df_nominal = n_mask - rank_fn(mean_design),
    trial_design = x,
    used_trial_design = used_trial
  )
}

# Coherent, spike-robust leverage-standardized GLS/REML residual scale estimator
# (Woolrich et al. 2001 Eq. 4 / Seber 1977, robustified against un-scrubbed motion spikes).
# Under N(0, sigma^2 Sigma_cor), unit-marginal whitened residuals e = (I - Q Q') W_unit y
# satisfy Var(e_t) = sigma^2 (1 - h_t), where h_t = ||Q[t, ]||_2^2 is the GLS leverage.
# Leverage-standardized innovations z_t = e_t / sqrt(1 - h_t) have exact variance sigma^2;
# Winsorizing at c = 3.0 * s_mad with exact Gaussian consistency beta(3.0) = E[min(Z^2, 9)]
# retains >99% REML efficiency on clean Gaussian data while bounding motion spikes.
hrf_estimate_gls_noise_sd <- function(y_rows, x_rows, rows, ar_vec, fallback_d2_sd = NULL,
                                      return_innovations = FALSE, trial_rows = NULL) {
  y_rows <- as.matrix(y_rows)
  m <- length(rows)
  p_vox <- ncol(y_rows)
  make_d2_unit <- function(mat) {
    z_fb <- matrix(0.0, m, p_vox)
    if (m < 3L || p_vox == 0L) return(z_fb)
    ok_c <- which(colSums(!is.finite(mat)) == 0L)
    if (!length(ok_c)) return(z_fb)
    tri_idx <- which((rows[3:m] - rows[1:(m - 2L)]) == 2L)
    if (length(tri_idx) >= 2L) {
      d2_c <- mat[tri_idx + 2L, ok_c, drop = FALSE] -
        2.0 * mat[tri_idx + 1L, ok_c, drop = FALSE] +
        mat[tri_idx, ok_c, drop = FALSE]
      out_rows <- tri_idx + 2L
    } else {
      d2_c <- mat[3:m, ok_c, drop = FALSE] - 2.0 * mat[2:(m - 1L), ok_c, drop = FALSE] + mat[1:(m - 2L), ok_c, drop = FALSE]
      out_rows <- 3:m
    }
    d2_c <- sweep(d2_c, 2L, colMeans(d2_c), "-")
    nrm_c <- sqrt(colSums(d2_c^2))
    good_c <- is.finite(nrm_c) & nrm_c > .Machine$double.eps
    if (any(good_c)) {
      z_fb[out_rows, ok_c[good_c]] <- sweep(d2_c[, good_c, drop = FALSE], 2L, nrm_c[good_c], "/")
    }
    z_fb
  }
  if (m < 4L || p_vox == 0L) {
    sd_na <- rep(NA_real_, p_vox)
    return(if (isTRUE(return_innovations)) list(sd = sd_na, z_unit = matrix(0.0, m, p_vox)) else sd_na)
  }
  if (is.null(fallback_d2_sd)) {
    r1 <- if (length(ar_vec) >= 1L) ar_vec[1L] else 0.0
    p1 <- if (length(ar_vec) >= 2L) ar_vec[2L] else 0.0
    p2 <- if (length(ar_vec) >= 3L) ar_vec[3L] else 0.0
    p3 <- if (length(ar_vec) >= 4L) ar_vec[4L] else 0.0
    p4 <- if (length(ar_vec) >= 5L) ar_vec[5L] else 0.0
    fallback_d2_sd <- estimate_decon_noise_sd(
      y_rows, noise_rho = r1, noise_psi1 = p1, noise_psi2 = p2, noise_psi3 = p3, noise_psi4 = p4,
      rows = rows
    )
  }
  if (is.null(x_rows) || m < 20L) {
    sd_fb <- as.numeric(fallback_d2_sd)
    return(if (isTRUE(return_innovations)) list(sd = sd_fb, z_unit = make_d2_unit(y_rows)) else sd_fb)
  }

  x_rows <- as.matrix(x_rows)
  ok_cols <- which(colSums(!is.finite(y_rows)) == 0L)
  want_trial_proj <- isTRUE(return_innovations) && !is.null(trial_rows) && ncol(trial_rows) > 0L
  w_bundle <- hrf_apply_whitening(
    list(
      wx = x_rows,
      wy = if (length(ok_cols)) y_rows[, ok_cols, drop = FALSE] else NULL,
      wt = if (want_trial_proj) as.matrix(trial_rows) else NULL
    ),
    rows = rows, ar_vec = ar_vec, unit_marginal = TRUE
  )
  wx <- w_bundle$wx
  nrm <- sqrt(colSums(wx^2)); nrm[nrm == 0] <- 1.0
  sv <- svd(sweep(wx, 2L, nrm, "/"), nv = 0)
  rk <- if (length(sv$d) && sv$d[1L] > 0) sum(sv$d > 1e-8 * sv$d[1L]) else 0L
  df_res <- m - rk
  if (df_res < 16L) {
    sd_fb <- as.numeric(fallback_d2_sd)
    return(if (isTRUE(return_innovations)) list(sd = sd_fb, z_unit = make_d2_unit(y_rows)) else sd_fb)
  }

  q <- if (rk > 0L) sv$u[, seq_len(rk), drop = FALSE] else matrix(0.0, m, 0L)
  out_sd <- rep(NA_real_, p_vox)
  z_unit <- if (isTRUE(return_innovations)) matrix(0.0, m, p_vox) else NULL
  if (!length(ok_cols)) {
    return(if (isTRUE(return_innovations)) list(sd = out_sd, z_unit = z_unit) else out_sd)
  }

  wy <- w_bundle$wy
  e_mat <- if (rk > 0L) wy - q %*% crossprod(q, wy) else wy
  lev <- if (rk > 0L) rowSums(q^2) else rep(0.0, m)
  keep_t <- (1.0 - lev) > 0.10
  use_lev <- sum(keep_t) >= 16L
  inv_sd_t <- if (use_lev) 1.0 / sqrt(1.0 - lev[keep_t]) else NULL
  c_win <- 3.0
  beta_win <- (2 * stats::pnorm(c_win) - 1) - 2 * c_win * stats::dnorm(c_win) +
    2 * (c_win^2) * (1 - stats::pnorm(c_win))

  # For spatial/cross-parcel innovation covariance (`return_innovations = TRUE`), if `x_rows`
  # fell back to a condition-mean design (dense rapid ITI), project out the top passband
  # singular modes of `trial_rows` so trial-by-trial neural fluctuations never enter `z_unit`.
  e_innov <- e_mat
  keep_innov <- keep_t
  inv_sd_innov <- inv_sd_t
  if (want_trial_proj && (df_res > 24L)) {
    wt <- w_bundle$wt
    if (rk > 0L) wt <- wt - q %*% crossprod(q, wt)
    nrm_t <- sqrt(colSums(wt^2))
    pos_t <- nrm_t > 1e-8
    if (any(pos_t)) {
      sv_t <- svd(sweep(wt[, pos_t, drop = FALSE], 2L, nrm_t[pos_t], "/"), nv = 0)
      rk_max <- max(0L, df_res - 24L)
      rk_t <- if (length(sv_t$d) && sv_t$d[1L] > 0) min(rk_max, sum(sv_t$d > 0.06 * sv_t$d[1L])) else 0L
      if (rk_t > 0L) {
        qt <- sv_t$u[, seq_len(rk_t), drop = FALSE]
        e_innov <- e_innov - qt %*% crossprod(qt, e_innov)
        lev_tot <- lev + rowSums(qt^2)
        keep_innov <- (1.0 - lev_tot) > 0.10
        if (sum(keep_innov) >= 16L) {
          inv_sd_innov <- 1.0 / sqrt(1.0 - lev_tot[keep_innov])
        } else {
          keep_innov <- keep_t
        }
      }
    }
  }
  use_lev_innov <- sum(keep_innov) >= 16L

  for (idx_k in seq_along(ok_cols)) {
    j <- ok_cols[idx_k]
    ej <- e_mat[, idx_k]
    s_gls <- sqrt(sum(ej^2) / df_res)
    if (!is.finite(s_gls) || s_gls <= .Machine$double.eps) {
      out_sd[j] <- fallback_d2_sd[j]
      next
    }
    s_est <- s_gls
    if (use_lev) {
      zj <- ej[keep_t] * inv_sd_t
      s0 <- stats::mad(zj, center = 0)
      if (is.finite(s0) && s0 > .Machine$double.eps) {
        u <- pmin(abs(zj), c_win * s0)
        s_est <- sqrt(mean(u^2) / beta_win)
      }
    }
    if (isTRUE(return_innovations)) {
      ej_in <- e_innov[, idx_k]
      if (use_lev_innov) {
        zi <- ej_in[keep_innov] * inv_sd_innov
        si0 <- stats::mad(zi, center = 0)
        ui <- if (is.finite(si0) && si0 > .Machine$double.eps) {
          sign(zi) * pmin(abs(zi), c_win * si0)
        } else zi
        ui <- ui - mean(ui)
        nrm_u <- sqrt(sum(ui^2))
        if (is.finite(nrm_u) && nrm_u > .Machine$double.eps) {
          z_unit[keep_innov, j] <- ui / nrm_u
        }
      } else {
        ui <- ej_in - mean(ej_in)
        nrm_u <- sqrt(sum(ui^2))
        if (is.finite(nrm_u) && nrm_u > .Machine$double.eps) {
          z_unit[, j] <- ui / nrm_u
        }
      }
    }
    # Guard against severe omitted-event misspecification inflating GLS residuals far above Delta^2 y
    d2_j <- fallback_d2_sd[j]
    if (is.finite(d2_j) && d2_j > .Machine$double.eps && s_est > 1.35 * d2_j) {
      s_est <- 1.25 * d2_j
    }
    out_sd[j] <- if (is.finite(s_est) && s_est > .Machine$double.eps) s_est else d2_j
  }
  if (isTRUE(return_innovations)) list(sd = out_sd, z_unit = z_unit) else out_sd
}

hrf_estimate_noise <- function(y, events, tr, offset, kernels, fallback, z, mask,
    groups, rho_grid, batch_size, fixed_rho = NULL, selected_kernels = NULL,
    spatial_coords = NULL, spatial_fwhm = NULL) {
  ids <- colnames(y)
  n_vox <- ncol(y)
  labels <- if (is.null(groups)) rep("run", n_vox) else unname(as.character(groups))
  n_obs <- nrow(y)
  rows <- which(mask)
  n_mask <- length(rows)

  if (!is.null(fixed_rho)) {
    r_fix <- unname(fixed_rho)
    nsd_fix <- rep(NA_real_, n_vox)
    if (n_mask >= 4L) {
      if (is.null(selected_kernels)) {
        md_info <- hrf_build_noise_mean_design(n_obs, events, tr, offset, kernels[, fallback], z, mask)
        x_r <- if (md_info$df_nominal >= 16L) md_info$design[rows, , drop = FALSE] else NULL
        nsd_fix <- hrf_estimate_gls_noise_sd(y[rows, , drop = FALSE], x_r, rows, r_fix)
      } else {
        k_sigs <- apply(round(selected_kernels, 4L), 2L, paste, collapse = ",")
        for (sg in unique(k_sigs)) {
          idx_s <- which(k_sigs == sg)
          md_info <- hrf_build_noise_mean_design(n_obs, events, tr, offset, selected_kernels[, idx_s[1L]], z, mask)
          x_r <- if (md_info$df_nominal >= 16L) md_info$design[rows, , drop = FALSE] else NULL
          nsd_fix[idx_s] <- hrf_estimate_gls_noise_sd(y[rows, idx_s, drop = FALSE], x_r, rows, r_fix)
        }
      }
    }
    return(data.frame(voxel = ids, group = labels,
      rho = rep(r_fix, n_vox), psi1 = 0, psi2 = 0, psi3 = 0, psi4 = 0,
      noise_sd = nsd_fix, method = "supplied", n_voxels = NA_integer_))
  }
  if (n_mask < 12L) return(data.frame(voxel = ids, group = labels,
    rho = 0, psi1 = 0, psi2 = 0, psi3 = 0, psi4 = 0, noise_sd = NA_real_,
    method = "insufficient_noise_df", n_voxels = 0L))

  # Partition voxels by frozen HRF kernel stratum (either fallback in Pass 1 or
  # cross-run selected_kernels in Pass 2) so non-canonical HRFs do not leak task
  # signal into noise REML profiles or GLS residual variance.
  if (is.null(selected_kernels)) {
    k_sigs <- rep("fb", n_vox)
    k_list <- list(fb = kernels[, fallback])
  } else {
    k_sigs <- apply(round(selected_kernels, 4L), 2L, paste, collapse = ",")
    u_sigs <- unique(k_sigs)
    k_list <- lapply(u_sigs, function(sg) selected_kernels[, which(k_sigs == sg)[1L]])
    names(k_list) <- u_sigs
  }
  md_list <- lapply(k_list, function(k_col) {
    hrf_build_noise_mean_design(n_obs, events, tr, offset, k_col, z, mask)
  })
  min_df_any <- min(vapply(md_list, `[[`, integer(1), "df_nominal"))
  if (min_df_any < 10L) return(data.frame(voxel = ids, group = labels,
    rho = 0, psi1 = 0, psi2 = 0, psi3 = 0, psi4 = 0, noise_sd = NA_real_,
    method = "insufficient_noise_df", n_voxels = 0L))

  des_list <- lapply(md_list, `[[`, "design")
  profiles <- matrix(NA_real_, length(rho_grid), n_vox)
  for (sg in names(des_list)) {
    idx_s <- which(k_sigs == sg)
    profiles[, idx_s] <- hrf_noise_profile_cpp(
      y[, idx_s, drop = FALSE], des_list[[sg]], rows, rho_grid, as.integer(batch_size)
    )
  }

  out <- data.frame(voxel = ids, group = labels, rho = 0, psi1 = 0, psi2 = 0,
    psi3 = 0, psi4 = 0, noise_sd = NA_real_,
    method = "insufficient_noise_voxels", n_voxels = 0L)

  q_use <- if (n_mask >= 96L && tr < 0.95) 4L else 2L
  lag_idx <- lapply(seq_len(q_use), function(L) {
    if (n_mask > L) which((rows[(1L + L):n_mask] - rows[1L:(n_mask - L)]) == L) else integer(0)
  })

  interp_rho <- function(prof_vec) {
    best <- which.min(prof_vec)
    rho_hat <- rho_grid[best]
    if (best > 1L && best < length(rho_grid) && length(rho_grid) >= 5L) {
      x1 <- rho_grid[best - 1L]; x2 <- rho_grid[best]; x3 <- rho_grid[best + 1L]
      y1 <- prof_vec[best - 1L]; y2 <- prof_vec[best]; y3 <- prof_vec[best + 1L]
      d1 <- x2 - x1; d2 <- x3 - x2
      if (d1 > 0 && d2 > 0) {
        g1 <- (y2 - y1) / d1; g2 <- (y3 - y2) / d2
        if (is.finite(g1) && is.finite(g2) && (g2 - g1) > 1e-12) {
          rho_interp <- x2 - (d1 * g2 + d2 * g1) / (2 * (g2 - g1))
          rho_hat <- round(min(max(rho_interp, x1), x3), 4L)
        }
      }
    }
    list(rho = rho_hat, best = best)
  }

  estimate_psi_weighted <- function(vox_idx, vox_weights, rho_val) {
    psi_zero <- rep(0.0, 4L)
    if (n_mask < 48L || length(lag_idx[[1L]]) < 24L || length(lag_idx[[2L]]) < 24L) {
      return(psi_zero)
    }
    q_eff <- if (q_use == 4L && length(lag_idx[[3L]]) >= 24L && length(lag_idx[[4L]]) >= 24L) 4L else 2L
    r_lags <- numeric(q_eff)
    w_tot <- 0.0
    df_acc <- 0.0
    u_chunks <- list()
    w_chunks <- list()
    vox_sigs <- k_sigs[vox_idx]
    u_sigs_psi <- unique(vox_sigs)
    in_psi_list <- vector("list", 2L * length(u_sigs_psi))
    for (si in seq_along(u_sigs_psi)) {
      sg <- u_sigs_psi[si]
      sub_v <- vox_idx[vox_sigs == sg]
      in_psi_list[[2L * si - 1L]] <- des_list[[sg]][rows, , drop = FALSE]
      in_psi_list[[2L * si]] <- y[rows, sub_v, drop = FALSE]
    }
    w_psi_list <- hrf_apply_whitening(in_psi_list, rows, rho_val, unit_marginal = FALSE)
    for (si in seq_along(u_sigs_psi)) {
      sg <- u_sigs_psi[si]
      sub_pos <- which(vox_sigs == sg)
      sub_w <- vox_weights[sub_pos]
      wz <- w_psi_list[[2L * si - 1L]]
      norm_z <- sqrt(colSums(wz^2)); norm_z[norm_z == 0] <- 1
      sv <- svd(sweep(wz, 2L, norm_z, "/"), nv = 0)
      rk <- if (length(sv$d) && sv$d[1L] > 0) sum(sv$d > 1e-8 * sv$d[1L]) else 0L
      df_res <- n_mask - rk
      if (df_res < 20L) next
      q_des <- if (rk > 0L) sv$u[, seq_len(rk), drop = FALSE] else matrix(0, n_mask, 0L)
      wy <- w_psi_list[[2L * si]]
      res1 <- if (rk > 0L) wy - q_des %*% crossprod(q_des, wy) else wy
      v_norms <- sqrt(colSums(res1^2))
      ok_v <- is.finite(v_norms) & v_norms > 1e-10
      if (!any(ok_v)) next
      res1_ok <- res1[, ok_v, drop = FALSE]
      res1_n <- sweep(res1_ok, 2L, v_norms[ok_v], "/")
      sw_ok <- sub_w[ok_v]
      sw_sum <- sum(sw_ok)
      if (sw_sum <= 0) next
      m_diag <- if (rk > 0L) pmax(0.0, 1.0 - rowSums(q_des^2)) else rep(1.0, n_mask)
      for (L in seq_len(q_eff)) {
        li <- lag_idx[[L]]
        # Normalize each lag-L product by the local energy of the participating gap-free
        # index pairs (li, li + L) so scattered motion/censor gaps never attenuate emp_L.
        d_exp_L <- sqrt(max(sum(m_diag[li + L]) * sum(m_diag[li]), 1e-12))
        bL <- if (rk > 0L) (-sum(q_des[li + L, , drop = FALSE] * q_des[li, , drop = FALSE])) / d_exp_L else 0.0
        r_hi <- res1_ok[li + L, , drop = FALSE]
        r_lo <- res1_ok[li, , drop = FALSE]
        den_L <- sqrt(pmax(colSums(r_hi^2) * colSums(r_lo^2), 1e-24))
        emp_L <- colSums(r_hi * r_lo) / den_L
        r_lags[L] <- r_lags[L] + sum(sw_ok * (emp_L - bL))
      }
      u_chunks[[length(u_chunks) + 1L]] <- res1_n
      w_chunks[[length(w_chunks) + 1L]] <- sw_ok
      w_tot <- w_tot + sw_sum
      df_acc <- df_acc + sw_sum * df_res
    }
    if (w_tot <= 0) return(psi_zero)
    r_lags <- pmin(pmax(r_lags / w_tot, -0.80), 0.80)
    df_eff <- max(20.0, df_acc / w_tot)
    pair_frac_12 <- min(length(lag_idx[[1L]]), length(lag_idx[[2L]])) / max(1.0, n_mask - 2.0)
    df_eff_12 <- max(16.0, df_eff * min(1.0, pair_frac_12))

    # Exact Isserlis fourth-moment variance factor sum_{i,j} w_i w_j (u_i' u_j)^2
    # accounting for both parcel size and cross-voxel spatial noise correlation.
    u_all <- if (length(u_chunks) == 1L) u_chunks[[1L]] else do.call(cbind, u_chunks)
    w_norm_u <- unlist(w_chunks, use.names = FALSE) / w_tot
    u_w <- sweep(u_all, 2L, sqrt(pmax(w_norm_u, 0)), "*")
    sp_var_factor <- if (ncol(u_w) <= nrow(u_w)) {
      g_mat <- crossprod(u_w)
      sum(g_mat * g_mat)
    } else {
      m_mat <- tcrossprod(u_w)
      sum(m_mat * m_mat)
    }
    sigma_k2_12 <- max(1e-6, sp_var_factor / df_eff_12)
    tau12 <- max(0.060, sqrt(5.9915 * sigma_k2_12))

    if (q_eff == 2L) {
      r1 <- r_lags[1L]; r2 <- r_lags[2L]
      k1 <- r1
      k2 <- min(max((r2 - r1 * r1) / max(1.0 - r1 * r1, 0.10), -0.80), 0.80)
      k_norm2 <- k1 * k1 + k2 * k2
      if (!is.finite(k_norm2) || k_norm2 <= tau12 * tau12) return(psi_zero)
      shrink_eb <- max(0.0, 1.0 - (tau12 * tau12) / k_norm2)
      k1_s <- min(max(shrink_eb * k1, -0.70), 0.70)
      k2_s <- min(max(shrink_eb * k2, -0.70), 0.70)
      psi2 <- round(k2_s, 4L)
      psi1 <- round(k1_s * (1.0 - psi2), 4L)
      return(c(psi1, psi2, 0.0, 0.0))
    }

    # Step-up Levinson-Durbin up to q_eff = 4 with blockwise Empirical-Bayes PACF shrinkage
    pacf_k <- numeric(q_eff)
    phi_ld <- numeric(0)
    pred_v <- 1.0
    for (m_ord in seq_len(q_eff)) {
      num <- if (m_ord == 1L) r_lags[1L] else r_lags[m_ord] - sum(phi_ld * rev(r_lags[1L:(m_ord - 1L)]))
      km <- min(max(num / max(pred_v, 0.12), -0.75), 0.75)
      pacf_k[m_ord] <- km
      phi_ld <- if (m_ord == 1L) km else c(phi_ld - km * rev(phi_ld), km)
      pred_v <- pred_v * (1.0 - km * km)
    }
    pair_frac_34 <- min(length(lag_idx[[3L]]), length(lag_idx[[4L]])) / max(1.0, n_mask - 4.0)
    df_eff_34 <- max(16.0, df_eff * min(1.0, pair_frac_34))
    sigma_k2_34 <- max(1e-6, sp_var_factor / df_eff_34)
    tau34 <- max(0.090, sqrt(9.2103 * sigma_k2_34))
    n12 <- sum(pacf_k[1:2]^2)
    s12 <- if (is.finite(n12) && n12 > tau12^2) max(0.0, 1.0 - (tau12^2) / n12) else 0.0
    pacf_k[1:2] <- pmin(pmax(s12 * pacf_k[1:2], -0.70), 0.70)
    n34 <- sum(pacf_k[3:4]^2)
    s34 <- if (is.finite(n34) && n34 > tau34^2) max(0.0, 1.0 - (tau34^2) / n34) else 0.0
    pacf_k[3:4] <- pmin(pmax(s34 * pacf_k[3:4], -0.55), 0.55)
    if (all(abs(pacf_k) <= 1e-10)) return(psi_zero)
    if (all(abs(pacf_k[3:4]) <= 1e-10)) {
      psi2 <- round(pacf_k[2L], 4L)
      psi1 <- round(pacf_k[1L] * (1.0 - psi2), 4L)
      return(c(psi1, psi2, 0.0, 0.0))
    }
    psi_cur <- numeric(0)
    for (m_ord in seq_len(4L)) {
      km <- pacf_k[m_ord]
      psi_cur <- if (m_ord == 1L) km else c(psi_cur - km * rev(psi_cur), km)
    }
    round(psi_cur, 4L)
  }

  # Damped 2nd-pass alternating REML: pre-whiten demeaned y and design by AR(q)(psi_0),
  # re-profile rho, and refresh psi so short-lag physiological/MA structure does not pull rho.
  refine_ar_fit <- function(vox_idx, vox_weights) {
    w_norm <- vox_weights / sum(vox_weights)
    prof_0 <- as.numeric(profiles[, vox_idx, drop = FALSE] %*% w_norm)
    ip0 <- interp_rho(prof_0)
    rho_hat <- ip0$rho
    best_idx <- ip0$best
    psi_hat <- estimate_psi_weighted(vox_idx, w_norm, rho_hat)
    if (n_mask >= 64L && any(abs(psi_hat) > 1e-4)) {
      prof_2 <- matrix(NA_real_, length(rho_grid), length(vox_idx))
      vox_sigs <- k_sigs[vox_idx]
      u_sigs_ref <- unique(vox_sigs)
      in_ref_list <- vector("list", 2L * length(u_sigs_ref))
      for (si in seq_along(u_sigs_ref)) {
        sg <- u_sigs_ref[si]
        sub_v <- vox_idx[vox_sigs == sg]
        y_rows_sub <- y[rows, sub_v, drop = FALSE]
        in_ref_list[[2L * si - 1L]] <- sweep(y_rows_sub, 2L, colMeans(y_rows_sub), "-")
        in_ref_list[[2L * si]] <- des_list[[sg]][rows, , drop = FALSE]
      }
      w_ref_list <- hrf_apply_whitening(in_ref_list, rows, c(0.0, psi_hat), unit_marginal = TRUE)
      for (si in seq_along(u_sigs_ref)) {
        sg <- u_sigs_ref[si]
        sub_pos <- which(vox_sigs == sg)
        sub_v <- vox_idx[sub_pos]
        y_f <- y[, sub_v, drop = FALSE]
        d_f <- des_list[[sg]]
        y_f[rows, ] <- w_ref_list[[2L * si - 1L]]
        d_f[rows, ] <- w_ref_list[[2L * si]]
        prof_2[, sub_pos] <- hrf_noise_profile_cpp(
          y_f, d_f, rows, rho_grid, as.integer(batch_size)
        )
      }
      ok2 <- colSums(is.finite(prof_2)) == length(rho_grid)
      if (any(ok2)) {
        sw2 <- w_norm[ok2] / sum(w_norm[ok2])
        ip1 <- interp_rho(as.numeric(prof_2[, ok2, drop = FALSE] %*% sw2))
        rho_hat <- round(0.75 * ip1$rho + 0.25 * rho_hat, 4L)
        best_idx <- ip1$best
        psi_hat <- estimate_psi_weighted(vox_idx, w_norm, rho_hat)
      }
    }
    list(rho = rho_hat, psi = psi_hat, boundary = best_idx %in% c(1L, length(rho_grid)))
  }

  spatial_w <- if (!is.null(spatial_coords)) {
    hrf_compute_spatial_weights(spatial_coords, labels, spatial_fwhm)
  } else NULL

  # Deterministic 1D quantization into at most K clusters (never touches .Random.seed)
  cluster_1d_det <- function(x, max_k = 3L) {
    ux <- sort(unique(round(x, 2L)))
    k <- min(as.integer(max_k), length(ux))
    if (k <= 1L) return(rep(1L, length(x)))
    probs <- seq(0.15, 0.85, length.out = k)
    ctrs <- unique(as.numeric(stats::quantile(x, probs = probs, names = FALSE)))
    if (length(ctrs) <= 1L) return(rep(1L, length(x)))
    for (iter in 1:10) {
      d_mat <- abs(outer(x, ctrs, "-"))
      cl <- max.col(-d_mat, ties.method = "first")
      new_ctrs <- vapply(seq_along(ctrs), function(ci) {
        m_i <- which(cl == ci)
        if (length(m_i)) mean(x[m_i]) else ctrs[ci]
      }, numeric(1))
      if (max(abs(new_ctrs - ctrs)) < 1e-6) break
      ctrs <- new_ctrs
    }
    as.integer(factor(cl))
  }

  for (g in unique(labels)) {
    index <- which(labels == g)
    usable <- index[colSums(is.finite(profiles[, index, drop = FALSE])) == length(rho_grid)]
    if (!length(usable)) next

    # Pool likelihood PROFILES across concordant voxels, not BOLD signals.
    # Woolrich (2001) SUSAN principle: when a parcel mixes distinct noise regimes
    # (e.g., GM core vs CSF/pial-vein or WM partial-volume border voxels) whose
    # 10-90% profile-argmin spread exceeds homogeneous sampling noise, weight each
    # voxel's neighborhood by bilateral similarity + spatial affinity and quantize
    # into <= 3 strata so operator caching remains compact.
    use_susan_g <- FALSE
    if (length(usable) >= 10L && n_mask >= 36L) {
      ind_rhos <- vapply(usable, function(v) interp_rho(profiles[, v])$rho, numeric(1))
      df_est <- max(20.0, mean(vapply(md_list, `[[`, integer(1), "df_nominal")))
      r_med <- stats::median(ind_rhos)
      exp_sd <- sqrt(max(0.05, 1.0 - r_med * r_med) / df_est)
      q90_10 <- as.numeric(diff(stats::quantile(ind_rhos, probs = c(0.10, 0.90), names = FALSE)))
      if (is.finite(q90_10) && q90_10 > max(0.22, 3.8 * exp_sd)) {
        use_susan_g <- TRUE
      }
    }

    if (!use_susan_g) {
      fit_g <- refine_ar_fit(usable, rep(1.0, length(usable)))
      out$rho[index] <- fit_g$rho
      out$psi1[index] <- fit_g$psi[1L]
      out$psi2[index] <- fit_g$psi[2L]
      out$psi3[index] <- fit_g$psi[3L]
      out$psi4[index] <- fit_g$psi[4L]
      out$method[index] <- if (fit_g$boundary) "pooled_reml_grid_boundary" else "pooled_reml"
      out$n_voxels[index] <- length(usable)
    } else {
      # Restrict SUSAN's auxiliary second-difference scale feature to gap-free
      # contiguous triplets on the acquisition clock so drift/task jumps across
      # censored motion segments never contaminate `raw_lsd`.
      tri_idx <- if (n_mask >= 3L) which((rows[3:n_mask] - rows[1:(n_mask - 2L)]) == 2L) else integer(0)
      if (length(tri_idx) >= 8L) {
        d2_u <- y[rows[tri_idx + 2L], usable, drop = FALSE] -
          2 * y[rows[tri_idx + 1L], usable, drop = FALSE] +
          y[rows[tri_idx], usable, drop = FALSE]
      } else {
        pair_idx <- if (n_mask >= 2L) which((rows[2:n_mask] - rows[1:(n_mask - 1L)]) == 1L) else integer(0)
        if (length(pair_idx) >= 6L) {
          d2_u <- y[rows[pair_idx + 1L], usable, drop = FALSE] -
            y[rows[pair_idx], usable, drop = FALSE]
        } else {
          d2_u <- y[rows[3:n_mask], usable, drop = FALSE] -
            2 * y[rows[2:(n_mask - 1L)], usable, drop = FALSE] +
            y[rows[1:(n_mask - 2L)], usable, drop = FALSE]
        }
      }
      raw_lsd <- log(pmax(apply(d2_u, 2L, stats::mad, center = 0), 1e-12))
      lsd_scale <- max(stats::sd(raw_lsd), 0.25)
      rho_scale <- max(0.12, 0.50 * stats::sd(ind_rhos))
      n_u <- length(usable)
      aff_mat <- matrix(0.0, n_u, n_u)
      smooth_rhos <- numeric(n_u)
      for (ii in seq_len(n_u)) {
        d_rho <- (ind_rhos - ind_rhos[ii]) / rho_scale
        d_lsd <- (raw_lsd - raw_lsd[ii]) / lsd_scale
        w_sim <- exp(-0.5 * (d_rho * d_rho + 0.5 * d_lsd * d_lsd))
        if (!is.null(spatial_w)) {
          w_sim <- w_sim * (0.25 + 0.75 * spatial_w[usable[ii], usable])
        }
        if (sum(w_sim) <= 1e-12) w_sim[ii] <- 1.0
        aff_mat[ii, ] <- w_sim / sum(w_sim)
        smooth_rhos[ii] <- interp_rho(as.numeric(profiles[, usable, drop = FALSE] %*% aff_mat[ii, ]))$rho
      }
      cl <- cluster_1d_det(smooth_rhos, max_k = 3L)
      susan_pfx <- if (length(unique(cl)) > 1L) "susan_" else ""
      for (cid in unique(cl)) {
        m_pos <- which(cl == cid)
        w_cl <- colMeans(aff_mat[m_pos, , drop = FALSE])
        fit_c <- refine_ar_fit(usable, w_cl)
        v_target <- usable[m_pos]
        out$rho[v_target] <- fit_c$rho
        out$psi1[v_target] <- fit_c$psi[1L]
        out$psi2[v_target] <- fit_c$psi[2L]
        out$psi3[v_target] <- fit_c$psi[3L]
        out$psi4[v_target] <- fit_c$psi[4L]
        out$method[v_target] <- paste0(susan_pfx, if (fit_c$boundary) "pooled_reml_grid_boundary" else "pooled_reml")
        out$n_voxels[v_target] <- length(m_pos)
      }
      missing_v <- setdiff(index, usable)
      if (length(missing_v)) {
        maj_cid <- as.integer(names(which.max(table(cl))))
        ref_v <- usable[which(cl == maj_cid)[1L]]
        out$rho[missing_v] <- out$rho[ref_v]
        out$psi1[missing_v] <- out$psi1[ref_v]
        out$psi2[missing_v] <- out$psi2[ref_v]
        out$psi3[missing_v] <- out$psi3[ref_v]
        out$psi4[missing_v] <- out$psi4[ref_v]
        out$method[missing_v] <- out$method[ref_v]
        out$n_voxels[missing_v] <- out$n_voxels[ref_v]
      }
    }

    # Populate coherent spike-robust GLS residual noise_sd per (HRF, AR) stratum
    g_keys <- paste(k_sigs[index], out$rho[index], out$psi1[index], out$psi2[index],
                    out$psi3[index], out$psi4[index], sep = ":")
    for (gk in unique(g_keys)) {
      sub_v <- index[g_keys == gk]
      v1 <- sub_v[1L]
      ar_v <- c(out$rho[v1], out$psi1[v1], out$psi2[v1], out$psi3[v1], out$psi4[v1])
      out$noise_sd[sub_v] <- hrf_estimate_gls_noise_sd(
        y[rows, sub_v, drop = FALSE],
        des_list[[k_sigs[v1]]][rows, , drop = FALSE],
        rows,
        ar_v
      )
    }
  }
  out
}

# Estimate the effective number of independent candidate tests across the HRF
# library. For compact libraries (K <= 9, such as "standard"), exact Bonferroni
# multiplicity (K for task, K - 1 for non-fallback shape departure) is used.
# For extended or dense libraries (K > 9) where candidate HRF shapes are
# strongly positively correlated, standard Bonferroni is overly conservative;
# we compute the Galwey (2009) spectral effective number of tests
# M_eff = (sum_i sqrt(lambda_i))^2 / sum_i lambda_i from the eigenvalues of the
# absolute cosine-similarity matrix |H_norm' H_norm|, floored at 9.
hrf_effective_candidates <- function(kernels) {
  if (is.null(kernels) || !is.matrix(kernels) || ncol(kernels) <= 9L) {
    k <- if (is.matrix(kernels)) ncol(kernels) else 1L
    return(list(task = k, shape = max(1L, k - 1L)))
  }
  k <- ncol(kernels)
  norms <- sqrt(colSums(kernels^2))
  norms[norms <= 0] <- 1
  h_norm <- sweep(kernels, 2L, norms, "/")
  r_mat <- abs(crossprod(h_norm))
  eig <- eigen(r_mat, symmetric = TRUE, only.values = TRUE)$values
  eig <- pmax(eig, 0)
  # Galwey (2009) spectral effective number of independent tests for dense/extended libraries (K > 9).
  m_eff <- if (sum(eig) > 0) (sum(sqrt(eig))^2) / sum(eig) else k
  m_eff <- max(9, min(k, m_eff))
  list(task = m_eff, shape = max(1, m_eff - 1))
}

# Cauchy combination test (Liu & Xie, 2020) for pooling p-values across
# spatially correlated voxels within a parcel. Because the standard Cauchy tail
# is invariant to arbitrary dependency structures, this provides a valid
# parcel-level omnibus test without assuming spatial independence across voxels.
# Conservative upper-tail point masses (e.g., p = 1 from one-sided/clamped tests)
# are bounded at 0.90 so they do not generate -Inf tangents that mask true signals.
hrf_cauchy_p <- function(p, upper_bound = 0.90) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-15), upper_bound)
  t_stat <- mean(tan((0.5 - p) * pi))
  0.5 - atan(t_stat) / pi
}

# Combine per-run calibrated GLS F-test p-values across training runs in a
# leave-target-run-out fold using Fisher's method (-2 * sum log p ~ chi^2_{2R}),
# apply library multiplicity correction (hrf_effective_candidates), and evaluate
# the three-gate calibrated decision rule:
#   1. compatible: no single training run rejects candidate k in favor of a
#      better-scoring alternative at alpha / (R * (K - 1)).
#   2. task_supported: Fisher-combined task F-test p <= alpha AND mean adjusted
#      partial R^2 >= min_r2 AND compatible across runs.
#   3. depart: Fisher-combined pairwise encompassing F-test vs fallback_kernel
#      has p <= alpha AND mean score strictly exceeds fallback_kernel (the
#      fallback kernel itself is exempt from the departure requirement).
hrf_fold_evidence <- function(results, train, counts, scores, fallback, alpha, min_r2,
                              kernels = NULL, subblock_mode = FALSE) {
  k <- nrow(scores); v <- ncol(scores)
  eff <- if (is.null(kernels)) list(task = k, shape = max(1L, k - 1L)) else hrf_effective_candidates(kernels)
  combine <- function(field, multiplicity) {
    logs <- matrix(0, k, v)
    cauchy_sum <- if (isTRUE(subblock_mode)) matrix(0, k, v) else NULL
    for (r in train) {
      p <- results[[r]][[field]]
      if (is.null(p)) next
      ok_p <- is.finite(p)
      p[!ok_p] <- 1
      logs <- logs + log(pmax(p, .Machine$double.xmin))
      if (isTRUE(subblock_mode)) {
        p_c <- pmin(pmax(p, 1e-15), 0.90)
        t_c <- tan((0.5 - p_c) * pi)
        t_c[!ok_p] <- 0
        cauchy_sum <- cauchy_sum + t_c
      }
    }
    raw <- matrix(stats::pchisq(as.vector(-2 * logs),
      df = rep(2 * pmax(counts, 1L), each = k), lower.tail = FALSE), k, v)
    if (isTRUE(subblock_mode)) {
      # Guard against residual within-run sub-block dependence (shared AR fit and
      # boundary carryover) using Brown-scaled Fisher + Cauchy tail bound.
      c_scale <- ifelse(counts >= 2L, 1.15, 1.0)
      raw_brown <- matrix(stats::pchisq(as.vector(sweep(-2 * logs, 2L, c_scale, "/")),
        df = rep(2 * pmax(counts, 1L) / c_scale, each = k), lower.tail = FALSE), k, v)
      raw_cauchy <- 0.5 - atan(sweep(cauchy_sum, 2L, pmax(counts, 1L), "/")) / pi
      raw <- pmax(raw, pmin(raw_brown, raw_cauchy))
    }
    raw[, counts == 0] <- NA_real_
    out <- pmin(raw * multiplicity, 1)
    out[, counts == 0] <- NA_real_
    list(adj = out, raw = raw)
  }
  task_res <- combine("task_p", eff$task)
  shape_res <- combine("fallback_p", eff$shape)
  task <- task_res$adj
  shape <- shape_res$adj
  compatible <- matrix(TRUE, k, v)
  # Simultaneous candidate compatibility within each voxel over all its usable
  # training runs. No averaging can hide one run's strong rejection of a shape.
  # Cap competitor multiplicity by 2 * eff$shape so dense/oversampled HRF grids
  # do not dilute incompatibility sensitivity through near-duplicate bases.
  m_incompat <- min(max(1L, k - 1L), max(eff$shape, as.integer(ceiling(2.0 * eff$shape))))
  cutoff <- alpha / (pmax(counts, 1L) * m_incompat)
  for (r in train) {
    p <- results[[r]]$incompatible_p
    if (is.null(p)) next
    p[!is.finite(p)] <- 1
    compatible <- compatible & sweep(p, 2L, cutoff, ">")
  }
  task_supported <- compatible & task <= alpha & scores >= min_r2
  task_supported[is.na(task_supported)] <- FALSE
  depart <- shape <= alpha & sweep(scores, 2L, scores[fallback, ], ">")
  depart[is.na(depart)] <- FALSE
  depart[fallback, ] <- TRUE
  list(task_p = task, shape_p = shape,
    task_p_raw = task_res$raw, shape_p_raw = shape_res$raw,
    eff_task = eff$task, eff_shape = eff$shape,
    compatible = compatible, task_supported = task_supported,
    depart = depart, supported = task_supported & depart)
}

# Pool calibrated evidence and compute Empirical-Bayes departure confidence omega
# across spatially correlated voxels within a parcel (or local neighborhood).
hrf_parcel_pooled_evidence <- function(gate, scores, eligible, fallback_kernel, alpha,
    min_r2, min_parcel_voxels, h_norm, kernel_cos, dist_fb, canonical_strength) {
  scores_el <- scores[, eligible, drop = FALSE]
  K <- nrow(scores_el)
  V <- ncol(scores_el)
  pooled <- if (V >= 1L) rowMeans(scores_el, na.rm = TRUE) else rep(NA_real_, K)
  if (V < min_parcel_voxels || all(!is.finite(pooled)) || max(pooled, na.rm = TRUE) < min_r2) {
    return(list(
      admissible = if (is.null(gate)) rep(TRUE, K) else (seq_len(K) == fallback_kernel),
      comp_ok = is.finite(pooled) & (pooled >= min_r2),
      omega = if (is.null(gate) && any(is.finite(pooled)) && max(pooled, na.rm = TRUE) >= min_r2) 1.0 else 0.0
    ))
  }
  if (is.null(gate)) {
    comp_ok <- is.finite(pooled) & (pooled >= min_r2)
    comp_ok[fallback_kernel] <- TRUE
    return(list(admissible = rep(TRUE, K), comp_ok = comp_ok, omega = 1.0))
  }

  # Estimate spatial correlation among parcel voxels from candidate score profiles
  cor_mat <- tryCatch(stats::cor(scores_el, use = "pairwise.complete.obs"),
    error = function(e) matrix(0.5, V, V))
  rho_bar <- if (V > 1L) mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE) else 1.0
  if (!is.finite(rho_bar)) rho_bar <- 0.5
  rho_bar <- max(0.25, min(0.95, rho_bar))
  denom_z <- sqrt(V + V * (V - 1) * rho_bar)

  eff_t <- if (!is.null(gate$eff_task)) gate$eff_task else K
  eff_s <- if (!is.null(gate$eff_shape)) gate$eff_shape else max(1L, K - 1L)
  p_task_raw <- if (!is.null(gate$task_p_raw)) gate$task_p_raw[, eligible, drop = FALSE] else
    (gate$task_p[, eligible, drop = FALSE] / eff_t)
  p_shape_raw <- if (!is.null(gate$shape_p_raw)) gate$shape_p_raw[, eligible, drop = FALSE] else
    (gate$shape_p[, eligible, drop = FALSE] / eff_s)
  p_t_clamped <- matrix(pmin(pmax(p_task_raw, 1e-15), 0.99), nrow = K, ncol = V)
  p_s_clamped <- matrix(pmin(pmax(p_shape_raw, 1e-15), 0.99), nrow = K, ncol = V)

  # Spatially-correlated Stouffer/Brown z-test combined with bounded Cauchy test
  z_task <- rowSums(stats::qnorm(1 - p_t_clamped), na.rm = TRUE) / denom_z
  z_shape <- rowSums(stats::qnorm(1 - p_s_clamped), na.rm = TRUE) / denom_z
  p_task_stouffer <- pmin(1.0, stats::pnorm(z_task, lower.tail = FALSE) * eff_t)
  p_shape_stouffer <- pmin(1.0, stats::pnorm(z_shape, lower.tail = FALSE) * eff_s)
  p_task_cauchy <- pmin(1.0, apply(p_t_clamped, 1L, hrf_cauchy_p) * eff_t)
  p_shape_cauchy <- pmin(1.0, apply(p_s_clamped, 1L, hrf_cauchy_p) * eff_s)

  p_task_pool <- pmin(p_task_stouffer, p_task_cauchy, na.rm = TRUE)
  p_shape_pool <- pmin(p_shape_stouffer, p_shape_cauchy, na.rm = TRUE)
  p_shape_pool[fallback_kernel] <- 0.0

  comp_cnt <- rowSums(gate$compatible[, eligible, drop = FALSE])
  task_cnt <- rowSums(gate$task_supported[, eligible, drop = FALSE])
  beat_fb <- rowSums(sweep(scores_el, 2L, scores_el[fallback_kernel, ], ">="), na.rm = TRUE)

  adm_strict <- rowSums(gate$supported[, eligible, drop = FALSE]) >= min_parcel_voxels
  comp_ok <- (comp_cnt >= ceiling(0.6 * V)) &
    ((task_cnt >= ceiling(0.5 * min_parcel_voxels)) | (!is.na(p_task_pool) & p_task_pool <= alpha)) &
    is.finite(pooled) & (pooled >= min_r2)

  adm_pooled <- if (K > 1L && any(gate$supported[-fallback_kernel, eligible, drop = FALSE])) {
    (comp_cnt >= min_parcel_voxels) &
      (task_cnt >= ceiling(0.5 * min_parcel_voxels)) &
      (beat_fb >= min_parcel_voxels) &
      !is.na(p_task_pool) & (p_task_pool <= alpha) &
      !is.na(p_shape_pool) & (p_shape_pool <= alpha)
  } else {
    rep(FALSE, K)
  }
  admissible <- adm_strict | adm_pooled

  # Continuous Empirical-Bayes departure weight omega in [0, 1] for soft_average:
  cand_non_fb <- setdiff(which(comp_ok), fallback_kernel)
  if (!length(cand_non_fb) || K <= 1L) {
    omega <- 0.0
  } else {
    obj_non_fb <- pooled[cand_non_fb] - canonical_strength * dist_fb[cand_non_fb]
    best_c <- cand_non_fb[which.max(obj_non_fb)]
    tau_best <- h_norm[, best_c] - kernel_cos[fallback_kernel, best_c] * h_norm[, fallback_kernel]
    tau_best <- tau_best / max(sqrt(sum(tau_best^2)), 1e-8)
    opp_cands <- vapply(cand_non_fb, function(k) {
      tk <- h_norm[, k] - kernel_cos[fallback_kernel, k] * h_norm[, fallback_kernel]
      tk <- tk / max(sqrt(sum(tk^2)), 1e-8)
      sum(tau_best * tk)
    }, numeric(1))
    spread_all <- max(pooled[comp_ok], na.rm = TRUE) - min(pooled[comp_ok], na.rm = TRUE)
    max_d_fb <- max(dist_fb[comp_ok], 1e-6)
    expected_drop <- 0.35 * spread_all * (dist_fb[best_c] / max_d_fb)
    diffs <- (scores_el[best_c, ] - canonical_strength * dist_fb[best_c]) -
      (scores_el[fallback_kernel, ] - expected_drop)
    if (any(opp_cands < -0.5)) {
      opp_c <- cand_non_fb[which.min(opp_cands)]
      diffs_asym <- scores_el[best_c, ] - scores_el[opp_c, ]
      if (mean(diffs_asym, na.rm = TRUE) <= canonical_strength * dist_fb[best_c]) {
        diffs[] <- -1
      }
    }
    if (!any(is.finite(diffs)) || mean(diffs, na.rm = TRUE) <= 0) {
      omega <- 0.0
    } else {
      sd_d <- stats::sd(diffs, na.rm = TRUE)
      eff_v <- V / (1 + (V - 1) * rho_bar * 0.65)
      t_diff <- if (is.finite(sd_d) && sd_d > 1e-8) mean(diffs, na.rm = TRUE) / (sd_d / sqrt(eff_v)) else 0.0
      z_eff <- max(z_shape[best_c], 0.75 * t_diff, na.rm = TRUE)
      omega <- max(0.0, min(1.0, (z_eff - 1.35) / (2.5 - 1.35)))
    }
  }
  comp_ok[fallback_kernel] <- TRUE
  list(admissible = admissible, comp_ok = comp_ok, omega = omega,
    p_task_pool = p_task_pool, p_shape_pool = p_shape_pool)
}

# Spread-normalized local Boltzmann evidence weights blended with Empirical-Bayes
# fallback prior weight (1 - omega).
hrf_soft_boltzmann_weights <- function(scores_vec, kernel_cos, dist_fb,
    fallback_kernel, comp_ok, omega, canonical_strength = 0.01,
    top_k = NULL, temp = 0.20, quantize_step = 0.02) {
  K <- length(scores_vec)
  w_fb <- rep(0.0, K)
  w_fb[fallback_kernel] <- 1.0
  if (K <= 1L || !is.finite(omega) || omega <= 1e-4) return(w_fb)
  comp_ok <- !is.na(comp_ok) & comp_ok & is.finite(scores_vec)
  if (!any(comp_ok)) return(w_fb)
  if (is.null(top_k)) {
    top_k <- max(4L, min(10L, as.integer(ceiling(0.25 * sum(comp_ok)))))
  }

  obj <- scores_vec - canonical_strength * dist_fb
  obj[!comp_ok] <- -Inf
  if (!any(is.finite(obj))) return(w_fb)

  k_top <- which.max(obj)
  ord <- order(obj, decreasing = TRUE)
  ord <- ord[is.finite(obj[ord]) & (kernel_cos[k_top, ord] >= 0.85)]
  if (!length(ord)) return(w_fb)
  top_idx <- utils::head(ord, min(as.integer(top_k), length(ord)))

  sc_ok <- scores_vec[comp_ok]
  spread <- max(sc_ok, na.rm = TRUE) - stats::median(sc_ok, na.rm = TRUE)
  spread <- max(spread, 1e-4)
  delta <- (obj[top_idx] - max(obj[top_idx])) / (temp * spread)
  w_sub <- exp(pmax(-20, delta))
  w_data <- rep(0.0, K)
  w_data[top_idx] <- w_sub / sum(w_sub)

  w_fin <- omega * w_data + (1.0 - omega) * w_fb
  if (quantize_step > 0) {
    w_q <- round(w_fin / quantize_step) * quantize_step
    if (sum(w_q) > 0) w_fin <- w_q / sum(w_q)
  }
  w_fin / sum(w_fin)
}

# Build stratum keys and per-voxel HRF kernel accessor for deconvolve_crossrun_hrf()
# and decon_epoch_regression(). In "hard" mode, preserves exact integer kernel_index
# strata; in "soft_average" mode, groups voxels that share identical composite weight
# vectors and returns their normalized composite HRF kernel.
hrf_stratum_kernels <- function(hrf_fit, r) {
  selected <- hrf_fit$selections[[r]]$kernel_index
  noise_r <- if (is.null(hrf_fit$noise)) NULL else hrf_fit$noise[[r]]
  rho <- if (is.null(noise_r)) rep(0, length(selected)) else noise_r$rho
  psi1 <- if (!is.null(noise_r) && "psi1" %in% names(noise_r)) noise_r$psi1 else rep(0, length(selected))
  psi2 <- if (!is.null(noise_r) && "psi2" %in% names(noise_r)) noise_r$psi2 else rep(0, length(selected))
  psi3 <- if (!is.null(noise_r) && "psi3" %in% names(noise_r)) noise_r$psi3 else rep(0, length(selected))
  psi4 <- if (!is.null(noise_r) && "psi4" %in% names(noise_r)) noise_r$psi4 else rep(0, length(selected))
  has_p34 <- abs(psi3) > 1e-10 | abs(psi4) > 1e-10
  has_p12 <- abs(psi1) > 1e-10 | abs(psi2) > 1e-10 | has_p34
  psi_sfx <- ifelse(
    has_p34,
    paste0(":", sprintf("%.6g", psi1), ":", sprintf("%.6g", psi2),
           ":", sprintf("%.6g", psi3), ":", sprintf("%.6g", psi4)),
    ifelse(
      has_p12,
      paste0(":", sprintf("%.6g", psi1), ":", sprintf("%.6g", psi2)),
      ""
    )
  )
  get_ar <- function(voxel_idx) {
    idx <- voxel_idx[1L]
    if (has_p34[idx]) {
      c(rho[idx], psi1[idx], psi2[idx], psi3[idx], psi4[idx])
    } else if (has_p12[idx]) {
      c(rho[idx], psi1[idx], psi2[idx])
    } else {
      rho[idx]
    }
  }
  is_soft <- identical(hrf_fit$settings$selection_mode, "soft_average") &&
    !is.null(hrf_fit$selected_kernels) && !is.null(hrf_fit$selected_kernels[[r]])
  if (!is_soft) {
    keys <- paste0(paste(selected, sprintf("%.17g", rho), sep = ":"), psi_sfx)
    get_kernel <- function(voxel_idx) {
      hrf_fit$kernels[, selected[voxel_idx], drop = FALSE]
    }
    return(list(
      keys = keys, rho = rho, psi1 = psi1, psi2 = psi2, psi3 = psi3, psi4 = psi4,
      selected = selected, is_soft = FALSE, get_kernel = get_kernel, get_ar = get_ar
    ))
  }
  w_mat <- if (!is.null(hrf_fit$kernel_weights) && !is.null(hrf_fit$kernel_weights[[r]])) {
    hrf_fit$kernel_weights[[r]]
  } else {
    hrf_fit$selected_kernels[[r]]
  }
  w_sig <- apply(round(w_mat, 4L), 2L, paste, collapse = ",")
  keys <- paste0(paste(w_sig, sprintf("%.17g", rho), sep = ":"), psi_sfx)
  sel_k_mat <- hrf_fit$selected_kernels[[r]]
  get_kernel <- function(voxel_idx) {
    sel_k_mat[, voxel_idx, drop = FALSE]
  }
  list(
    keys = keys, rho = rho, psi1 = psi1, psi2 = psi2, psi3 = psi3, psi4 = psi4,
    selected = selected, is_soft = TRUE, get_kernel = get_kernel, get_ar = get_ar
  )
}

# Construct the exact gap-aware AR(1) or cascaded AR(1)xAR(q) Cholesky whitening
# matrix W on observed acquisition indices `rows`.
# - When `unit_marginal = FALSE` (default), W %*% Sigma_cor %*% t(W) = sig_u2 * I_m
#   (reducing identically to (1 - rho^2) * I_m under pure AR(1) when psi = 0).
# - When `unit_marginal = TRUE`, W %*% Sigma_cor %*% t(W) = I_m for exact GLS/REML
#   variance estimation.
# On contiguous grids (all(diff(rows) == 1)), builds the exact O(m * p) stationary
# AR(p) Cholesky factor with exact p x p Toeplitz initial block; across arbitrary
# scrub/censor gaps with higher-order psi != 0, builds the exact subset Toeplitz
# Cholesky factor so zero off-diagonal correlation leaks across censor boundaries.
hrf_whitening_matrix <- function(rows, rho, unit_marginal = FALSE) {
  m <- length(rows)
  w <- diag(m)
  if (!length(rho) || m <= 1L) return(w)
  has_psi <- length(rho) >= 2L && any(is.finite(rho[-1L]) & abs(rho[-1L]) > 1e-10)
  if (!has_psi) {
    rho1 <- rho[1L]
    if (!is.finite(rho1) || abs(rho1) <= 1e-12) return(w)
    rho1 <- min(max(rho1, -0.95), 0.95)
    s_u <- sqrt(max(1.0 - rho1 * rho1, 1e-8))
    w[1L, 1L] <- s_u
    for (i in seq.int(2L, m)) {
      a <- rho1^(rows[i] - rows[i - 1L])
      multiplier <- s_u / sqrt(max(1.0 - a * a, 1e-12))
      w[i, i] <- multiplier
      w[i, i - 1L] <- -a * multiplier
    }
    return(if (isTRUE(unit_marginal)) w / s_u else w)
  }

  is_contig <- all(diff(rows) == 1L)
  lag_span <- if (is_contig) max(8L, length(rho) + 2L) else max(8L, max(rows) - min(rows))
  info <- hrf_ar_cascade_info(rho, max_lag = lag_span)
  phi <- info$phi
  p <- length(phi)
  if (p == 0L) return(w)
  sig_u <- sqrt(info$sig_u2)

  if (is_contig && m > p) {
    w_unit <- matrix(0.0, m, m)
    Gamma_p <- matrix(info$acf[abs(outer(seq_len(p), seq_len(p), "-")) + 1L], p, p)
    L_p <- tryCatch(t(chol(Gamma_p)), error = function(e) NULL)
    if (!is.null(L_p)) {
      w_unit[seq_len(p), seq_len(p)] <- forwardsolve(L_p, diag(p))
      inv_su <- 1.0 / sig_u
      for (i in seq.int(p + 1L, m)) {
        w_unit[i, i] <- inv_su
        for (k in seq_len(p)) {
          w_unit[i, i - k] <- -phi[k] * inv_su
        }
      }
      return(if (isTRUE(unit_marginal)) w_unit else w_unit * sig_u)
    }
  }

  idx <- abs(outer(rows, rows, "-")) + 1L
  Sigma_sub <- matrix(info$acf[idx], m, m)
  R_chol <- tryCatch(chol(Sigma_sub), error = function(e) {
    diag(Sigma_sub) <- diag(Sigma_sub) + 1e-6
    chol(Sigma_sub)
  })
  w_unit <- t(backsolve(R_chol, diag(m)))
  if (isTRUE(unit_marginal)) w_unit else w_unit * sig_u
}

# Apply exact AR(1) or cascaded AR(p) Cholesky whitening to a matrix (or list of
# matrices) on observed acquisition indices `rows` via O(m * p) C++ ArWhitener,
# matching `hrf_whitening_matrix(rows, ar_vec, unit_marginal) %*% X` without
# allocating a dense m x m matrix in R.
hrf_apply_whitening <- function(target, rows, ar_vec, unit_marginal = FALSE) {
  rows <- as.integer(rows)
  m <- length(rows)
  as_mat_target <- function(obj) {
    if (is.list(obj)) {
      lapply(obj, function(el) if (is.null(el)) NULL else as.matrix(el))
    } else {
      as.matrix(obj)
    }
  }
  if (m <= 1L || !length(ar_vec)) return(as_mat_target(target))
  has_psi <- length(ar_vec) >= 2L && any(is.finite(ar_vec[-1L]) & abs(ar_vec[-1L]) > 1e-10)
  if (!has_psi) {
    rho1 <- ar_vec[1L]
    if (!is.finite(rho1) || abs(rho1) <= 1e-12) return(as_mat_target(target))
    rho1 <- min(max(rho1, -0.95), 0.95)
    return(hrf_whiten_apply_cpp(
      target = as_mat_target(target), rows = rows, rho = rho1,
      unit_marginal = isTRUE(unit_marginal)
    ))
  }
  is_contig <- all(diff(rows) == 1L)
  lag_span <- if (is_contig) max(8L, length(ar_vec) + 2L) else max(8L, max(rows) - min(rows))
  info <- hrf_ar_cascade_info(ar_vec, max_lag = lag_span)
  if (!length(info$phi)) return(as_mat_target(target))
  hrf_whiten_apply_cpp(
    target = as_mat_target(target), rows = rows, rho = 0.0,
    ar_phi = info$phi, ar_acf = info$acf, sig_u2 = info$sig_u2,
    unit_marginal = isTRUE(unit_marginal)
  )
}

hrf_score_run_ar <- function(y_sub, design, z, rows, batch_size, ar_vec,
    evidence = FALSE, fallback_index = 1L) {
  rho1 <- if (length(ar_vec)) ar_vec[1L] else 0.0
  has_psi <- length(ar_vec) >= 2L && any(is.finite(ar_vec[-1L]) & abs(ar_vec[-1L]) > 1e-10)
  if (!has_psi || length(rows) < 3L) {
    return(hrf_score_run_cpp(y_sub, design, z, rows, as.integer(batch_size),
      rho = rho1, evidence = evidence, fallback_index = as.integer(fallback_index)))
  }
  is_contig <- all(diff(rows) == 1L)
  lag_span <- if (is_contig) max(8L, length(ar_vec) + 2L) else max(8L, max(rows) - min(rows))
  info <- hrf_ar_cascade_info(ar_vec, max_lag = lag_span)
  if (!length(info$phi)) {
    return(hrf_score_run_cpp(y_sub, design, z, rows, as.integer(batch_size),
      rho = 0.0, evidence = evidence, fallback_index = as.integer(fallback_index)))
  }
  hrf_score_run_cpp(y_sub, design, z, rows, as.integer(batch_size),
    rho = 0.0, evidence = evidence, fallback_index = as.integer(fallback_index),
    ar_phi = info$phi, ar_acf = info$acf, sig_u2 = info$sig_u2)
}

hrf_inverse_operator <- function(n, kernel, rows, rho, nuisance, settings) {
  has_psi <- length(rho) >= 2L && any(is.finite(rho[-1L]) & abs(rho[-1L]) > 1e-10)
  if (settings$backend == "cpp") {
    if (!has_psi) {
      return(hrf_inverse_operator_cpp(as.integer(n), kernel,
        rows, rho[1L], nuisance, settings$lambda, settings$penalty, settings$ridge_floor,
        as.integer(settings$tail_padding)))
    }
    all_rows <- c(rows, if (settings$tail_padding) n + seq_len(settings$tail_padding))
    is_contig <- all(diff(all_rows) == 1L)
    lag_span <- if (is_contig) max(8L, length(rho) + 2L) else max(8L, max(all_rows) - min(all_rows))
    info <- hrf_ar_cascade_info(rho, max_lag = lag_span)
    if (!length(info$phi)) {
      return(hrf_inverse_operator_cpp(as.integer(n), kernel,
        rows, 0.0, nuisance, settings$lambda, settings$penalty, settings$ridge_floor,
        as.integer(settings$tail_padding)))
    }
    return(hrf_inverse_operator_cpp(as.integer(n), kernel,
      rows, 0.0, nuisance, settings$lambda, settings$penalty, settings$ridge_floor,
      as.integer(settings$tail_padding),
      ar_phi = info$phi, ar_acf = info$acf, sig_u2 = info$sig_u2))
  }
  # Exact masked, nuisance-adjusted GLS reference (used for backend = "r").
  measured <- length(rows)
  all_rows <- c(rows, if (settings$tail_padding) n + seq_len(settings$tail_padding))
  h <- decon_convolution_matrix(n + settings$tail_padding, kernel)[all_rows, , drop = FALSE]
  w <- hrf_whitening_matrix(all_rows, rho)
  z <- rbind(nuisance[rows, , drop = FALSE], matrix(0, settings$tail_padding, ncol(nuisance)))
  z <- w %*% z; h <- w %*% h
  rhs <- w[, seq_len(measured), drop = FALSE]
  if (ncol(z)) {
    norm <- sqrt(colSums(z^2)); norm[norm == 0] <- 1
    fit <- svd(sweep(z, 2L, norm, "/"), nv = 0)
    rank <- sum(fit$d > 1e-8 * fit$d[1L])
    if (rank) {
      q <- fit$u[, seq_len(rank), drop = FALSE]
      h <- h - q %*% crossprod(q, h); rhs <- rhs - q %*% crossprod(q, rhs)
    }
  }
  lhs <- crossprod(h) + settings$lambda * decon_penalty_crossprod(ncol(h), settings$penalty)
  diag(lhs) <- diag(lhs) + settings$ridge_floor
  solve(lhs, crossprod(h, rhs))[length(kernel) - 1L + seq_len(n), , drop = FALSE]
}

# Gain-calibrated, task-weighted adaptive prewhitened GCV selection for the
# deconvolution penalty lambda in deconvolve_crossrun_hrf().
#
# Standard observation-space GCV minimizes prediction error ||W H (s_hat - s)||_2^2
# in convolved BOLD space, which severely over-smooths event-related MEDuSA
# deconvolution for three reasons:
#   1. Unfiltered null/low-SNR voxels have strictly decreasing E[GCV(lambda)],
#      pulling a global lambda toward the upper grid boundary.
#   2. Because the passband singular values of W H are small (~0.15-0.30),
#      shrinking s_hat -> 0 reduces observation-space noise faster than signal
#      loss whenever per-TR event power is smaller than per-TR noise variance.
#   3. In Stage 3 (decon_epoch_regression), physical scale is restored by
#      dividing by the diagonal resolution gain g_lvl(lambda) = b_0' h_0.
#      Consequently, the gain-restored variance of both a single-lag readout
#      V_lvl(lambda) = ||b_0||_2^2 / g_lvl(lambda)^2 and a 1-TR temporal contrast
#      V_diff(lambda) = ||b_0 - b_1||_2^2 / g_diff(lambda)^2 are U-shaped in
#      lambda: beyond a finite variance minimum lambda_var*(TR, rho, penalty),
#      increasing lambda increases BOTH smoothing bias AND gain-restored variance.
#   4. Discrete derivative penalties scale with TR^4 (diff2/sobolev) and AR(1)
#      whitening scales signal energy by (1 - rho)^2, so the optimal penalty
#      shifts systematically across acquisition TRs and noise autocorrelation.
#
# This helper pools task-supported voxels/parcels from hrf_fit$selections,
# evaluates the prewhitened GCV trace alongside the exact interior gain-restored
# level/contrast variances and relative modulation transfer function C_rel(lambda)
# on a compact local window, and selects the lambda that minimizes the combined
# gain-calibrated prewhitened GCV criterion gcv_adaptive.
hrf_tune_decon_lambda <- function(bold, hrf_fit, lambda_grid, scaling = "zscore",
    penalty = "sobolev", ridge_floor = 1e-8, observation_mask = "fit",
    backend = "cpp", gcv_rule = c("min", "1se"), nuisance_mode = "ignore", tail_padding = 0L) {
  gcv_rule <- match.arg(gcv_rule)
  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  if (length(lambda_grid) == 1L) {
    return(list(
      lambda = lambda_grid[1L],
      tuning = data.frame(
        kernel = 1L, lambda = lambda_grid[1L], df = NA_real_, rss = NA_real_,
        gcv = NA_real_, gcv_se = NA_real_, var_level = NA_real_,
        var_contrast = NA_real_, mtf_rel = 1.0, gcv_adaptive = 1.0, snr_task = NA_real_
      )
    ))
  }

  n_runs <- length(bold)
  n_vox <- ncol(bold[[1L]])
  r2_mat <- matrix(0.0, nrow = n_vox, ncol = n_runs)
  tp_mat <- matrix(1.0, nrow = n_vox, ncol = n_runs)
  has_r2 <- FALSE
  if (!is.null(hrf_fit$selections)) {
    for (r in seq_len(min(n_runs, length(hrf_fit$selections)))) {
      sel_r <- hrf_fit$selections[[r]]
      if (is.data.frame(sel_r) && nrow(sel_r) == n_vox) {
        if ("best_adjusted_partial_r2" %in% names(sel_r)) {
          v_r2 <- as.numeric(sel_r$best_adjusted_partial_r2)
          v_r2[!is.finite(v_r2)] <- 0.0
          r2_mat[, r] <- pmax(0.0, v_r2)
          has_r2 <- TRUE
        }
        if ("selected_task_p" %in% names(sel_r)) {
          v_tp <- as.numeric(sel_r$selected_task_p)
          v_tp[!is.finite(v_tp)] <- 1.0
          tp_mat[, r] <- pmin(1.0, pmax(0.0, v_tp))
        }
      }
    }
  }
  r2_v <- rowMeans(r2_mat, na.rm = TRUE)
  r2_v[!is.finite(r2_v)] <- 0.0
  tp_v <- apply(tp_mat, 1L, min, na.rm = TRUE)
  tp_v[!is.finite(tp_v)] <- 1.0

  tr <- if (!is.null(hrf_fit$settings$tr) && is.finite(hrf_fit$settings$tr) && hrf_fit$settings$tr > 0) {
    as.numeric(hrf_fit$settings$tr)
  } else {
    1.0
  }

  run_tuning_list <- list()
  run_weights <- numeric(0)
  wiener_lambdas <- numeric(0)
  wiener_tunings <- list()

  for (r_tune in seq_along(bold)) {
    y_t <- bold[[r_tune]]
    n_t <- nrow(y_t)
    mask_t <- if (observation_mask == "fit" && !is.null(hrf_fit$specification$valid_rows)) {
      hrf_fit$specification$valid_rows[[r_tune]]
    } else {
      rep(TRUE, n_t)
    }
    rows_t <- which(!is.na(mask_t) & mask_t)
    if (length(rows_t) < 3L) next
    st_t <- decon_signal_status(y_t[rows_t, , drop = FALSE])
    ok_cols <- which(st_t == "ok")
    if (!length(ok_cols)) next

    r2_r_vec <- if (has_r2 && ncol(r2_mat) >= r_tune) r2_mat[, r_tune] else r2_v
    tp_r_vec <- if (ncol(tp_mat) >= r_tune) tp_mat[, r_tune] else tp_v
    r2_ok <- r2_r_vec[ok_cols]
    tp_ok <- tp_r_vec[ok_cols]
    pos_r2 <- r2_ok[r2_ok > 0]
    if (length(pos_r2)) {
      r2_cut <- max(0.005, as.numeric(stats::quantile(pos_r2, probs = 0.35, na.rm = TRUE)))
      active_cols <- ok_cols[r2_ok >= r2_cut | tp_ok <= 0.15]
    } else {
      active_cols <- integer(0)
    }
    if (length(active_cols) < min(8L, length(ok_cols))) {
      ord_ok <- order(r2_ok, -tp_ok, decreasing = TRUE)
      active_cols <- ok_cols[utils::head(ord_ok, min(64L, length(ok_cols)))]
    } else if (length(active_cols) > 128L) {
      ord_act <- order(r2_r_vec[active_cols], decreasing = TRUE)
      active_cols <- active_cols[utils::head(ord_act, 128L)]
    }

    w_act <- pmax(r2_r_vec[active_cols], 1e-3)
    w_act <- w_act / sum(w_act)
    sel_t <- hrf_fit$selections[[r_tune]]$kernel_index[active_cols]
    k_tab <- tapply(w_act, sel_t, sum)
    k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
    if (!length(k_mode) || is.na(k_mode)) k_mode <- 1L

    norm_mode_tune <- if (!is.null(hrf_fit$settings$kernel_normalize)) {
      hrf_fit$settings$kernel_normalize
    } else {
      "signed_sum"
    }
    k_tune <- if (identical(hrf_fit$settings$selection_mode, "soft_average") &&
        !is.null(hrf_fit$selected_kernels) && !is.null(hrf_fit$selected_kernels[[r_tune]])) {
      k_avg <- as.numeric(hrf_fit$selected_kernels[[r_tune]][, active_cols, drop = FALSE] %*% w_act)
      tryCatch(
        normalize_decon_kernel(k_avg, mode = norm_mode_tune),
        error = function(e) as.numeric(hrf_fit$kernels[, k_mode])
      )
    } else {
      as.numeric(hrf_fit$kernels[, k_mode])
    }
    attr(k_tune, "normalize") <- norm_mode_tune

    noise_rt <- if (is.null(hrf_fit$noise)) NULL else hrf_fit$noise[[r_tune]]
    rho_t <- if (is.null(noise_rt)) {
      0.0
    } else {
      stats::median(noise_rt$rho[active_cols], na.rm = TRUE)
    }
    if (!is.finite(rho_t)) rho_t <- 0.0
    psi1_t <- if (!is.null(noise_rt) && "psi1" %in% names(noise_rt)) {
      stats::median(noise_rt$psi1[active_cols], na.rm = TRUE)
    } else {
      0.0
    }
    psi2_t <- if (!is.null(noise_rt) && "psi2" %in% names(noise_rt)) {
      stats::median(noise_rt$psi2[active_cols], na.rm = TRUE)
    } else {
      0.0
    }
    psi3_t <- if (!is.null(noise_rt) && "psi3" %in% names(noise_rt)) {
      stats::median(noise_rt$psi3[active_cols], na.rm = TRUE)
    } else {
      0.0
    }
    psi4_t <- if (!is.null(noise_rt) && "psi4" %in% names(noise_rt)) {
      stats::median(noise_rt$psi4[active_cols], na.rm = TRUE)
    } else {
      0.0
    }
    if (!is.finite(psi1_t)) psi1_t <- 0.0
    if (!is.finite(psi2_t)) psi2_t <- 0.0
    if (!is.finite(psi3_t)) psi3_t <- 0.0
    if (!is.finite(psi4_t)) psi4_t <- 0.0
    ar_t <- if (abs(psi3_t) > 1e-10 || abs(psi4_t) > 1e-10) {
      c(rho_t, psi1_t, psi2_t, psi3_t, psi4_t)
    } else if (abs(psi1_t) > 1e-10 || abs(psi2_t) > 1e-10) {
      c(rho_t, psi1_t, psi2_t)
    } else {
      rho_t
    }

    nsd_tune <- NULL
    if (identical(scaling, "noise_sd") && !is.null(hrf_fit$specification$events[[r_tune]])) {
      mask_tune_bool <- seq_len(n_t) %in% rows_t
      nuis_t <- if (!is.null(hrf_fit$specification$nuisance)) hrf_fit$specification$nuisance[[r_tune]] else matrix(numeric(), n_t, 0L)
      z_t_base <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n_t), nuis_t)
      off_t <- if (length(hrf_fit$settings$time_offset) >= r_tune) hrf_fit$settings$time_offset[r_tune] else 0.0
      md_t <- hrf_build_noise_mean_design(n_t, hrf_fit$specification$events[[r_tune]], tr, off_t, k_tune, z_t_base, mask_tune_bool)
      if (md_t$df_nominal >= 16L) {
        nsd_tune <- hrf_estimate_gls_noise_sd(
          y_t[rows_t, active_cols, drop = FALSE],
          md_t$design[rows_t, , drop = FALSE],
          rows_t,
          ar_t
        )
      }
    }

    sc_vox <- scale_decon_input(
      y_t[rows_t, active_cols, drop = FALSE],
      scaling = scaling,
      noise_rho = rho_t,
      noise_psi1 = psi1_t,
      noise_psi2 = psi2_t,
      noise_psi3 = psi3_t,
      noise_psi4 = psi4_t,
      noise_sd_override = nsd_tune,
      rows = rows_t
    )$y
    if (identical(nuisance_mode, "joint")) {
      nuis_j <- if (!is.null(hrf_fit$specification$nuisance)) hrf_fit$specification$nuisance[[r_tune]] else matrix(numeric(), n_t, 0L)
      z_j <- cbind(1, seq(-1, 1, length.out = n_t), nuis_j)[rows_t, , drop = FALSE]
      nrm_j <- sqrt(colSums(z_j^2)); nrm_j[nrm_j == 0] <- 1
      sv_j <- svd(sweep(z_j, 2L, nrm_j, "/"), nv = 0)
      rk_j <- if (length(sv_j$d) && sv_j$d[1L] > 0) sum(sv_j$d > 1e-8 * sv_j$d[1L]) else 0L
      if (rk_j > 0L && length(rows_t) - rk_j >= 8L) {
        qj <- sv_j$u[, seq_len(rk_j), drop = FALSE]
        sc_vox <- sc_vox - qj %*% crossprod(qj, sc_vox)
      }
    }
    y_tune <- sc_vox
    groups_all <- hrf_fit$specification$spatial_groups
    if (!is.null(groups_all) && length(groups_all) == n_vox) {
      g_act <- as.character(groups_all[active_cols])
      u_groups <- unique(g_act[!is.na(g_act)])
      parcel_cols <- list()
      for (g in u_groups) {
        idx_g <- which(g_act == g)
        if (length(idx_g) >= 2L) {
          wg <- sqrt(pmax(r2_r_vec[active_cols[idx_g]], 1e-3))
          wg <- wg / sum(wg)
          sub_g <- sc_vox[, idx_g, drop = FALSE]
          # Align voxel signs deterministically via leading SVD direction so mixed-polarity
          # task voxels inside a parcel do not cancel out during GCV pooling.
          sgn_g <- tryCatch({
            vg <- svd(sub_g, nu = 0, nv = 1)$v[, 1L]
            if (vg[which.max(wg)] < 0) vg <- -vg
            ifelse(vg >= 0, 1.0, -1.0)
          }, error = function(e) rep(1.0, length(idx_g)))
          p_series <- as.numeric(sub_g %*% (wg * sgn_g))
          if (all(is.finite(p_series)) && stats::sd(p_series) > 1e-8) {
            parcel_cols[[g]] <- p_series
          }
        }
      }
      if (length(parcel_cols)) {
        p_mat <- do.call(cbind, parcel_cols)
        p_scaled <- scale_decon_input(
          p_mat,
          scaling = "noise_sd",
          noise_rho = rho_t,
          noise_psi1 = psi1_t,
          noise_psi2 = psi2_t,
          noise_psi3 = psi3_t,
          noise_psi4 = psi4_t,
          rows = rows_t
        )$y
        y_tune <- cbind(p_scaled, sc_vox[, seq_len(min(32L, ncol(sc_vox))), drop = FALSE])
      }
    }
    if (tail_padding > 0L) {
      y_tune <- rbind(y_tune, matrix(0.0, as.integer(tail_padding), ncol(y_tune)))
    }

    tune_res <- deconvolve_reglin(
      BOLDobs = y_tune,
      kernel = k_tune,
      lambda = NULL,
      lambda_grid = lambda_grid,
      penalty = penalty,
      tune_by = "global",
      normalize = FALSE,
      demean = FALSE,
      ridge_floor = ridge_floor,
      return_diagnostics = TRUE,
      prewhiten_gcv = rho_t != 0,
      prewhiten_rho = rho_t,
      backend = if (backend == "wiener") "wiener" else "cpp",
      gcv_rule = gcv_rule,
      scaling = "none"
    )
    tun <- tune_res$tuning
    if (backend == "wiener" || !is.data.frame(tun) || nrow(tun) != length(lambda_grid)) {
      lam_out <- as.numeric(tune_res$lambda[1L])
      if (!is.finite(lam_out)) lam_out <- lambda_grid[1L]
      wiener_lambdas <- c(wiener_lambdas, lam_out)
      wiener_tunings[[length(wiener_tunings) + 1L]] <- tun
      run_weights <- c(run_weights, as.numeric(length(rows_t)))
      next
    }

    # Evaluate exact interior gain-restored variance and modulation transfer on a
    # compact window using a fast 2-column normal-equation solve.
    l_k <- length(k_tune)
    n_win <- min(length(rows_t), max(64L, 2L * l_k + 16L))
    if (n_win < l_k + 8L) n_win <- max(length(rows_t), l_k + 8L)
    t0 <- max(3L, min(n_win - 2L, as.integer(floor(n_win * 0.45))))
    t1 <- min(n_win, t0 + 1L)

    h_full <- decon_convolution_matrix(n_win, k_tune)[seq_len(n_win), , drop = FALSE]
    wh <- hrf_apply_whitening(h_full, seq_len(n_win), ar_t, unit_marginal = FALSE)
    hth <- crossprod(wh)
    p_pen <- decon_penalty_crossprod(ncol(wh), penalty)
    idx_lat <- (l_k - 1L) + seq_len(n_win)
    idx0 <- idx_lat[t0]
    idx1 <- idx_lat[t1]
    h0 <- wh[, idx0]
    h1 <- wh[, idx1]
    h_diff <- h0 - h1

    e_rhs <- matrix(0.0, nrow = ncol(wh), ncol = 2L)
    e_rhs[idx0, 1L] <- 1.0
    e_rhs[idx0, 2L] <- 1.0
    e_rhs[idx1, 2L] <- -1.0

    n_grid <- length(lambda_grid)
    v_lvl <- v_diff <- c_mtf <- g_rel <- numeric(n_grid)
    rf <- max(ridge_floor, 1e-8)
    d_isi <- max(2L, as.integer(round(2.5 / tr)))
    far_mask <- abs(seq_len(n_win) - t0) >= d_isi
    for (i in seq_len(n_grid)) {
      lhs <- hth + lambda_grid[i] * p_pen
      diag(lhs) <- diag(lhs) + rf
      q_mat <- tryCatch(solve(lhs, e_rhs), error = function(e) MASS::ginv(lhs) %*% e_rhs)
      b0 <- as.numeric(wh %*% q_mat[, 1L])
      bdiff <- as.numeric(wh %*% q_mat[, 2L])
      r_row <- as.numeric(crossprod(q_mat[, 1L], hth[, idx_lat, drop = FALSE]))

      g_lvl <- max(sum(b0 * h0), 1e-4)
      g_diff <- max(sum(bdiff * h_diff), 1e-4)
      leak_isi <- if (any(far_mask)) {
        min(0.85, sqrt(sum(r_row[far_mask]^2)) / g_lvl)
      } else {
        0.0
      }
      g_eff <- g_lvl * max(0.15, 1.0 - 0.65 * (leak_isi^2))
      v_lvl[i] <- sum(b0^2) / (g_eff^2)
      v_diff[i] <- sum(bdiff^2) / (g_diff^2)
      g_rel[i] <- g_lvl

      conc_0 <- g_lvl / sqrt(max(sum(r_row^2), 1e-12))
      step_rel <- g_diff / (2.0 * g_lvl)
      c_mtf[i] <- 0.5 * conc_0 + 0.5 * step_rel
    }

    g_rel <- pmin(1.0, pmax(0.0, g_rel / max(g_rel[1L], 1e-6)))
    c_rel <- pmin(1.0, pmax(0.0, 0.65 * (c_mtf / max(c_mtf[1L], 1e-6)) + 0.35 * g_rel))
    p_exp <- min(1.0, tr^0.65)
    v_norm <- 0.50 * (v_lvl / min(v_lvl))^p_exp + 0.50 * (v_diff / min(v_diff))^p_exp

    gcv_vec <- as.numeric(tun$gcv)
    gcv_min <- min(gcv_vec, na.rm = TRUE)
    gcv_end <- gcv_vec[n_grid]
    snr_gcv <- if (is.finite(gcv_min) && gcv_min > 0 && is.finite(gcv_end)) {
      max(0.0, (gcv_end - gcv_min) / gcv_min)
    } else {
      0.0
    }
    r2_top <- if (has_r2 && length(r2_ok)) {
      mean(utils::head(sort(r2_ok, decreasing = TRUE), min(12L, length(r2_ok))), na.rm = TRUE)
    } else {
      0.0
    }
    if (!is.finite(r2_top)) r2_top <- 0.0
    snr_r2 <- (8.0 / max(tr, 0.25)) * r2_top / max(1.0 - r2_top, 0.10)
    snr_task <- min(6.0, max(snr_r2, if (has_r2) 0.35 * snr_gcv else 2.0 * snr_gcv))

    idx_gmin <- which.min(gcv_vec)
    gcv_excess <- rep(0.0, n_grid)
    if (length(idx_gmin) == 1L && is.finite(gcv_min) && gcv_min > 0) {
      s_hi <- seq_len(n_grid) >= idx_gmin
      gcv_excess[s_hi] <- pmax(0.0, (gcv_vec[s_hi] - gcv_min) / gcv_min)
    }

    w_bias <- 0.25 + 0.35 * snr_task
    gcv_adaptive <- v_norm + w_bias * (1.0 - c_rel)^2 + 0.30 * (gcv_excess^2)

    tun$var_level <- v_lvl
    tun$var_contrast <- v_diff
    tun$mtf_rel <- c_rel
    tun$gcv_adaptive <- gcv_adaptive
    tun$snr_task <- rep(snr_task, n_grid)
    run_tuning_list[[length(run_tuning_list) + 1L]] <- tun
    run_weights <- c(run_weights, as.numeric(length(rows_t)))
  }

  if (length(wiener_lambdas) > 0L && !length(run_tuning_list)) {
    w_norm <- run_weights / sum(run_weights)
    ord_w <- order(w_norm, wiener_lambdas)
    lam_pool <- exp(sum(w_norm[ord_w] * log(pmax(wiener_lambdas[ord_w], 1e-12))))
    best_w <- which.min(abs(log(pmax(lambda_grid, 1e-12)) - log(max(lam_pool, 1e-12))))
    return(list(lambda = lambda_grid[best_w], tuning = wiener_tunings[[ord_w[1L]]]))
  }

  if (length(run_tuning_list) > 0L) {
    n_grid <- length(lambda_grid)
    w_norm <- run_weights / sum(run_weights)
    # Sort run contributions canonically so floating-point summation is strictly
    # permutation-invariant across arbitrary run orderings.
    sig_sort <- vapply(run_tuning_list, function(df_r) as.numeric(df_r$gcv_adaptive[1L]), numeric(1L))
    ord_r <- order(w_norm, sig_sort)
    tun_pool <- run_tuning_list[[ord_r[1L]]]
    avg_cols <- intersect(
      c("df", "rss", "gcv", "gcv_se", "var_level", "var_contrast", "mtf_rel", "gcv_adaptive", "snr_task"),
      names(tun_pool)
    )
    for (cn in avg_cols) {
      acc_c <- numeric(n_grid)
      for (k in ord_r) {
        acc_c <- acc_c + w_norm[k] * as.numeric(run_tuning_list[[k]][[cn]])
      }
      tun_pool[[cn]] <- acc_c
    }
    gcv_adaptive_pool <- tun_pool$gcv_adaptive
    best_idx <- which.min(gcv_adaptive_pool)
    if (!length(best_idx) || !is.finite(best_idx)) best_idx <- 1L
    if (gcv_rule == "1se" && best_idx < n_grid) {
      var_min_idx <- max(which.min(tun_pool$var_level), which.min(tun_pool$var_contrast))
      if (best_idx < var_min_idx) {
        best_idx <- min(var_min_idx, best_idx + 1L)
      }
    }
    return(list(lambda = lambda_grid[best_idx], tuning = tun_pool))
  }

  list(
    lambda = lambda_grid[1L],
    tuning = data.frame(
      kernel = 1L, lambda = lambda_grid, df = NA_real_, rss = NA_real_,
      gcv = NA_real_, gcv_se = NA_real_, var_level = NA_real_,
      var_contrast = NA_real_, mtf_rel = NA_real_, gcv_adaptive = NA_real_, snr_task = NA_real_
    )
  )
}
