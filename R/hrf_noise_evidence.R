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
hrf_cauchy_p <- function(p) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
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
    out <- matrix(stats::pchisq(as.vector(-2 * logs),
      df = rep(2 * pmax(counts, 1L), each = k), lower.tail = FALSE), k, v)
    out <- pmin(out * multiplicity, 1)
    out[, counts == 0] <- NA_real_
    out
  }
  task <- combine("task_p", eff$task)
  shape <- combine("fallback_p", eff$shape)
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
  list(task_p = task, shape_p = shape, compatible = compatible,
    task_supported = task_supported, depart = depart,
    supported = task_supported & depart)
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
