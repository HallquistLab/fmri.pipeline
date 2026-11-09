# Helpers for the cross-run interface. Covariance parameters may use their own
# run's observations; shape selection and its evidence use OTHER runs only.
hrf_estimate_noise <- function(y, events, tr, offset, kernels, fallback, z, mask,
    groups, rho_grid, batch_size, fixed_rho = NULL) {
  ids <- colnames(y)
  labels <- if (is.null(groups)) rep("run", ncol(y)) else unname(groups)
  if (!is.null(fixed_rho)) return(data.frame(voxel = ids, group = labels,
    rho = rep(unname(fixed_rho), length(ids)), method = "supplied", n_voxels = NA_integer_))
  if (sum(mask) < 12L) return(data.frame(voxel = ids, group = labels,
    rho = 0, method = "insufficient_noise_df", n_voxels = 0L))
  # Free event amplitudes keep trial variation out of the covariance estimate.
  # Derivatives allow modest shape mismatch, when residual df permit them.
  durations <- if (is.null(events$duration)) rep(0, nrow(events)) else events$duration
  x <- hrf_event_design_cpp(nrow(y), tr, offset, events$onset, durations,
    seq_len(nrow(events)), nrow(events), kernels[, fallback, drop = FALSE])[[1L]]
  derivative <- function(a) rbind(a[2:nrow(a), , drop = FALSE], 0) - rbind(0, a[1:(nrow(a) - 1L), , drop = FALSE])
  augmented <- cbind(z, x, derivative(x), derivative(derivative(x)))
  # A rank check is preferable to a column-count check for repeated events.
  rank <- function(a) {
    a <- a[mask, , drop = FALSE]
    norm <- sqrt(colSums(a^2)); norm[norm == 0] <- 1
    s <- svd(sweep(a, 2L, norm, "/"), nu = 0, nv = 0)$d
    if (length(s) && s[1L] > 0) sum(s > 1e-8 * s[1L]) else 0L
  }
  mean_design <- if (sum(mask) - rank(augmented) >= 20L) augmented else cbind(z, x)
  if (sum(mask) - rank(mean_design) < 10L) return(data.frame(voxel = ids, group = labels,
    rho = 0, method = "insufficient_noise_df", n_voxels = 0L))
  profiles <- hrf_noise_profile_cpp(y, mean_design, which(mask), rho_grid, as.integer(batch_size))
  out <- data.frame(voxel = ids, group = labels, rho = 0,
    method = "insufficient_noise_voxels", n_voxels = 0L)
  for (g in unique(labels)) {
    index <- which(labels == g)
    usable <- index[colSums(is.finite(profiles[, index, drop = FALSE])) == length(rho_grid)]
    if (!length(usable)) next
    # Pool likelihood PROFILES, not BOLD signals. This regularizes covariance
    # without canceling signed neural responses or assuming independent voxels
    # for significance testing. A group with no data does not borrow other groups.
    best <- which.min(rowMeans(profiles[, usable, drop = FALSE]))
    out$rho[index] <- rho_grid[best]
    out$method[index] <- if (best %in% c(1L, length(rho_grid))) "pooled_reml_grid_boundary" else "pooled_reml"
    out$n_voxels[index] <- length(usable)
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
hrf_fold_evidence <- function(results, train, counts, scores, fallback, alpha, min_r2, kernels = NULL) {
  k <- nrow(scores); v <- ncol(scores)
  eff <- if (is.null(kernels)) list(task = k, shape = max(1L, k - 1L)) else hrf_effective_candidates(kernels)
  combine <- function(field, multiplicity) {
    logs <- matrix(0, k, v)
    for (r in train) {
      p <- results[[r]][[field]]
      if (is.null(p)) next
      p[!is.finite(p)] <- 1
      logs <- logs + log(pmax(p, .Machine$double.xmin))
    }
    raw <- matrix(stats::pchisq(as.vector(-2 * logs),
      df = rep(2 * pmax(counts, 1L), each = k), lower.tail = FALSE), k, v)
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
  cutoff <- alpha / (pmax(counts, 1L) * max(1L, k - 1L))
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
# vectors and returns their unit-sum composite HRF kernel.
hrf_stratum_kernels <- function(hrf_fit, r) {
  selected <- hrf_fit$selections[[r]]$kernel_index
  rho <- if (is.null(hrf_fit$noise)) rep(0, length(selected)) else hrf_fit$noise[[r]]$rho
  is_soft <- identical(hrf_fit$settings$selection_mode, "soft_average") &&
    !is.null(hrf_fit$selected_kernels) && !is.null(hrf_fit$selected_kernels[[r]])
  if (!is_soft) {
    keys <- paste(selected, sprintf("%.17g", rho), sep = ":")
    get_kernel <- function(voxel_idx) {
      hrf_fit$kernels[, selected[voxel_idx], drop = FALSE]
    }
    return(list(keys = keys, rho = rho, selected = selected, is_soft = FALSE, get_kernel = get_kernel))
  }
  w_mat <- if (!is.null(hrf_fit$kernel_weights) && !is.null(hrf_fit$kernel_weights[[r]])) {
    hrf_fit$kernel_weights[[r]]
  } else {
    hrf_fit$selected_kernels[[r]]
  }
  w_sig <- apply(round(w_mat, 4L), 2L, paste, collapse = ",")
  keys <- paste(w_sig, sprintf("%.17g", rho), sep = ":")
  sel_k_mat <- hrf_fit$selected_kernels[[r]]
  get_kernel <- function(voxel_idx) {
    sel_k_mat[, voxel_idx, drop = FALSE]
  }
  list(keys = keys, rho = rho, selected = selected, is_soft = TRUE, get_kernel = get_kernel)
}

# Construct the exact gap-aware AR(1) Cholesky whitening matrix W on observed
# acquisition indices `rows`. When censoring creates a gap of d = rows[i] - rows[i-1]
# TRs between consecutive usable observations, the lag-d correlation is rho^d
# and the accumulated innovation variance is (1 - rho^(2d)) / (1 - rho^2),
# ensuring whitening remains exact across scrub/censor gaps.
hrf_whitening_matrix <- function(rows, rho) {
  m <- length(rows)
  w <- diag(m)
  if (rho == 0) return(w)
  w[1L, 1L] <- sqrt(1 - rho^2)
  for (i in seq.int(2L, m)) {
    a <- rho^(rows[i] - rows[i - 1L])
    multiplier <- sqrt((1 - rho^2) / (1 - a^2))
    w[i, i] <- multiplier; w[i, i - 1L] <- -a * multiplier
  }
  w
}

hrf_inverse_operator <- function(n, kernel, rows, rho, nuisance, settings) {
  if (settings$backend == "cpp") return(hrf_inverse_operator_cpp(as.integer(n), kernel,
    rows, rho, nuisance, settings$lambda, settings$penalty, settings$ridge_floor,
    as.integer(settings$tail_padding)))
  # Independent R reference for the exact same masked, nuisance-adjusted GLS.
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
    backend = "cpp", gcv_rule = c("min", "1se")) {
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

    r2_ok <- r2_v[ok_cols]
    tp_ok <- tp_v[ok_cols]
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
      ord_act <- order(r2_v[active_cols], decreasing = TRUE)
      active_cols <- active_cols[ utils::head(ord_act, 128L) ]
    }

    w_act <- pmax(r2_v[active_cols], 1e-3)
    w_act <- w_act / sum(w_act)
    sel_t <- hrf_fit$selections[[r_tune]]$kernel_index[active_cols]
    k_tab <- tapply(w_act, sel_t, sum)
    k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
    if (!length(k_mode) || is.na(k_mode)) k_mode <- 1L

    k_tune <- if (identical(hrf_fit$settings$selection_mode, "soft_average") &&
        !is.null(hrf_fit$selected_kernels) && !is.null(hrf_fit$selected_kernels[[r_tune]])) {
      k_avg <- as.numeric(hrf_fit$selected_kernels[[r_tune]][, active_cols, drop = FALSE] %*% w_act)
      if (all(is.finite(k_avg)) && abs(sum(k_avg)) > 1e-6) k_avg / sum(k_avg) else hrf_fit$kernels[, k_mode]
    } else {
      hrf_fit$kernels[, k_mode]
    }

    rho_t <- if (is.null(hrf_fit$noise) || is.null(hrf_fit$noise[[r_tune]])) {
      0.0
    } else {
      stats::median(hrf_fit$noise[[r_tune]]$rho[active_cols], na.rm = TRUE)
    }
    if (!is.finite(rho_t)) rho_t <- 0.0

    sc_vox <- scale_decon_input(y_t[rows_t, active_cols, drop = FALSE], scaling = scaling, noise_rho = rho_t)$y
    y_tune <- sc_vox
    groups_all <- hrf_fit$specification$spatial_groups
    if (!is.null(groups_all) && length(groups_all) == n_vox) {
      g_act <- as.character(groups_all[active_cols])
      u_groups <- unique(g_act[!is.na(g_act)])
      parcel_cols <- list()
      for (g in u_groups) {
        idx_g <- which(g_act == g)
        if (length(idx_g) >= 2L) {
          wg <- sqrt(pmax(r2_v[active_cols[idx_g]], 1e-3))
          wg <- wg / sum(wg)
          p_series <- as.numeric(sc_vox[, idx_g, drop = FALSE] %*% wg)
          if (all(is.finite(p_series)) && stats::sd(p_series) > 1e-8) {
            parcel_cols[[g]] <- p_series
          }
        }
      }
      if (length(parcel_cols)) {
        p_mat <- do.call(cbind, parcel_cols)
        p_scaled <- scale_decon_input(p_mat, scaling = "noise_sd", noise_rho = rho_t)$y
        y_tune <- cbind(p_scaled, sc_vox[, seq_len(min(32L, ncol(sc_vox))), drop = FALSE])
      }
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
      return(list(lambda = lam_out, tuning = tun))
    }

    # Evaluate exact interior gain-restored variance and modulation transfer on a
    # compact window using a fast 2-column normal-equation solve.
    l_k <- length(k_tune)
    n_win <- min(length(rows_t), max(64L, 2L * l_k + 16L))
    if (n_win < l_k + 8L) n_win <- max(length(rows_t), l_k + 8L)
    t0 <- max(3L, min(n_win - 2L, as.integer(floor(n_win * 0.45))))
    t1 <- min(n_win, t0 + 1L)

    h_full <- decon_convolution_matrix(n_win, k_tune)[seq_len(n_win), , drop = FALSE]
    w_win <- hrf_whitening_matrix(seq_len(n_win), rho_t)
    wh <- w_win %*% h_full
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

    best_idx <- which.min(gcv_adaptive)
    if (!length(best_idx) || !is.finite(best_idx)) best_idx <- 1L
    if (gcv_rule == "1se" && best_idx < n_grid) {
      var_min_idx <- max(which.min(v_lvl), which.min(v_diff))
      if (best_idx < var_min_idx) {
        best_idx <- min(var_min_idx, best_idx + 1L)
      }
    }

    tun$var_level <- v_lvl
    tun$var_contrast <- v_diff
    tun$mtf_rel <- c_rel
    tun$gcv_adaptive <- gcv_adaptive
    tun$snr_task <- rep(snr_task, n_grid)
    return(list(lambda = lambda_grid[best_idx], tuning = tun))
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
