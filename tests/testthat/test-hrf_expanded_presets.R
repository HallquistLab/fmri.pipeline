test_that("hrf_library supports manifold, enet, Inf undershoot, and positive_lobe normalization", {
  set.seed(101)
  seed_before <- .Random.seed

  man_1 <- hrf_library(1.0, preset = "manifold")
  man_2 <- hrf_library(2.0, preset = "manifold")
  enet_45 <- hrf_library(1.0, preset = "enet", n_enet = 45L)
  enet_25 <- hrf_library(1.0, preset = "enet", n_enet = 25L)

  # Deterministic and RNG-free
  expect_identical(.Random.seed, seed_before)
  expect_identical(enet_45, hrf_library(1.0, preset = "enet", n_enet = 45L))

  expect_equal(ncol(man_1), 44L)
  expect_equal(ncol(enet_45), 45L)
  expect_equal(ncol(enet_25), 25L)

  # Canonical (6, 1, 6) is column 1
  expect_equal(attr(man_1, "parameters")[1L, ],
    data.frame(peak_delay = 6, dispersion = 1, undershoot_ratio = 6))
  expect_equal(attr(enet_45, "parameters")[1L, ],
    data.frame(peak_delay = 6, dispersion = 1, undershoot_ratio = 6))

  # TR-aware dispersion floor prevents sub-Nyquist narrower-than-TR aliasing at slow TR
  expect_gte(min(attr(man_2, "parameters")$dispersion), 0.80)
  expect_equal(min(attr(man_1, "parameters")$dispersion), 0.60)

  # Auto normalization uses positive_lobe for medusa/manifold/enet and signed_sum for fixed-ratio-6 presets
  expect_equal(attr(man_1, "normalize"), "positive_lobe")
  expect_equal(attr(enet_45, "normalize"), "positive_lobe")
  expect_identical(hrf_library(1.0, preset = "medusa"), man_1)
  expect_equal(attr(hrf_library(1.0, preset = "canonical_27"), "normalize"), "signed_sum")
  expect_equal(ncol(hrf_library(1.0, preset = "canonical_27")), 27L)
  expect_equal(colSums(pmax(man_1, 0)), rep(1, ncol(man_1)), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(colSums(pmax(enet_45, 0)), rep(1, ncol(enet_45)), tolerance = 1e-12, ignore_attr = TRUE)

  # Inf undershoot produces pure single-gamma nonnegative HRF matching large finite limit
  cust_inf <- hrf_library(0.8, peak_delays = c(3.0, 6.0), dispersions = 0.8,
    undershoot_ratios = c(6, Inf), normalize = "positive_lobe")
  pars_inf <- attr(cust_inf, "parameters")
  inf_cols <- which(is.infinite(pars_inf$undershoot_ratio))
  expect_length(inf_cols, 2L)
  expect_true(all(cust_inf[, inf_cols] >= 0))

  cust_huge <- hrf_library(0.8, peak_delays = 3.0, dispersions = 0.8,
    undershoot_ratios = 1e12, normalize = "positive_lobe")
  idx_3_inf <- which(pars_inf$peak_delay == 3.0 & is.infinite(pars_inf$undershoot_ratio))
  expect_equal(cust_inf[, idx_3_inf], cust_huge[, 1L], tolerance = 1e-10)

  # Peak normalization
  cust_pk <- hrf_library(0.8, peak_delays = c(3.5, 6.0), dispersions = 1.0,
    undershoot_ratios = c(4, Inf), normalize = "peak")
  expect_equal(apply(cust_pk, 2L, max), rep(1, ncol(cust_pk)), tolerance = 1e-12, ignore_attr = TRUE)
})

test_that("C++ cascaded AR(p) whitener and inverse operator match dense R reference across gaps and penalties", {
  n <- 96L
  kernel <- as.numeric(hrf_library(1.0, peak_delays = 4.5, dispersions = 0.8,
    undershoot_ratios = Inf, normalize = "positive_lobe")[, 1L])
  rows_contig <- seq_len(n)
  rows_gap <- setdiff(seq_len(n), c(1L, 12:16, 48L, 80:83))
  z_none <- matrix(numeric(0), n, 0L)
  z_joint <- cbind(1, seq(-1, 1, length.out = n), sin(seq_len(n) / 8), cos(seq_len(n) / 11))

  ar_cases <- list(
    ar3 = c(0.45, 0.22, -0.14),
    ar5 = c(0.35, 0.18, -0.11, 0.07, -0.04)
  )

  for (ar_vec in ar_cases) {
    for (rows in list(rows_contig, rows_gap)) {
      for (tail_padding in c(0L, 4L)) {
        for (penalty in c("diff2", "sobolev", "ridge")) {
          s_cpp <- list(backend = "cpp", lambda = 0.05, penalty = penalty,
            ridge_floor = 1e-8, tail_padding = tail_padding)
          s_r <- modifyList(s_cpp, list(backend = "r"))

          op_cpp_ig <- fmri.pipeline:::hrf_inverse_operator(n, kernel, rows, ar_vec, z_none, s_cpp)
          op_r_ig <- fmri.pipeline:::hrf_inverse_operator(n, kernel, rows, ar_vec, z_none, s_r)
          expect_equal(op_cpp_ig, op_r_ig, tolerance = 1e-9)

          op_cpp_jt <- fmri.pipeline:::hrf_inverse_operator(n, kernel, rows, ar_vec, z_joint, s_cpp)
          op_r_jt <- fmri.pipeline:::hrf_inverse_operator(n, kernel, rows, ar_vec, z_joint, s_r)
          expect_equal(op_cpp_jt, op_r_jt, tolerance = 1e-7)
          expect_lt(max(abs(op_cpp_jt %*% z_joint[rows, ])), 5e-7)
        }
      }

      # Also verify Step-1 scoring equivalence against explicit dense whitening matrix
      set.seed(202)
      lib_small <- hrf_library(1.0, peak_delays = c(3.5, 5.0, 6.0), dispersions = c(0.8, 1.0),
        undershoot_ratios = c(6, Inf), normalize = "positive_lobe")
      ev <- data.frame(onset = seq(10, 82, by = 9), condition = rep(c("A", "B"), length.out = 9))
      des <- fmri.pipeline:::hrf_event_design_cpp(n, 1.0, 0.0, ev$onset, rep(0, nrow(ev)),
        match(ev$condition, c("A", "B")), 2L, lib_small)
      y <- matrix(rnorm(n * 8L), nrow = n, ncol = 8L)
      y[, 1:4] <- y[, 1:4] + des[[2L]][, 1L] * 1.5
      z_run <- cbind(1, seq(-1, 1, length.out = n))

      got_cpp <- fmri.pipeline:::hrf_score_run_ar(y, des, z_run, rows, 64L,
        ar_vec = ar_vec, evidence = TRUE, fallback_index = 1L)

      w_mat <- fmri.pipeline:::hrf_whitening_matrix(rows, ar_vec)
      m <- length(rows)
      w_int <- as.numeric(w_mat %*% rep(1.0, m))
      y_ref <- y
      y_c <- sweep(y[rows, , drop = FALSE], 2L, colMeans(y[rows, , drop = FALSE]), "-")
      wy <- w_mat %*% y_c
      wy <- wy - outer(w_int / mean(w_int), colMeans(wy))
      y_ref[rows, ] <- wy
      z_ref <- z_run
      z_ref[rows, ] <- w_mat %*% z_run[rows, , drop = FALSE]
      des_ref <- lapply(des, function(d) {
        d[rows, ] <- w_mat %*% d[rows, , drop = FALSE]
        d
      })
      got_ref <- fmri.pipeline:::hrf_score_run_cpp(y_ref, des_ref, z_ref, rows, 64L,
        rho = 0.0, evidence = TRUE, fallback_index = 1L)

      expect_equal(got_cpp$scores, got_ref$scores, tolerance = 1e-10)
      expect_equal(got_cpp$task_p, got_ref$task_p, tolerance = 1e-10)
      expect_equal(got_cpp$fallback_p, got_ref$fallback_p, tolerance = 1e-10)
      expect_equal(got_cpp$incompatible_p, got_ref$incompatible_p, tolerance = 1e-10)
    }
  }
})

test_that("manifold and enet presets recover fast narrow zero-undershoot HRFs and preserve cross-region amplitude parity", {
  set.seed(404)
  tr <- 1.0
  n_time <- 180L
  runs <- c("run1", "run2", "run3")
  groups <- stats::setNames(rep(c("subcort_fast", "cort_deep", "null_parcel"), each = 8L),
    paste0("v", seq_len(24L)))

  true_fast <- hrf_library(tr, peak_delays = 3.0, dispersions = 0.80,
    undershoot_ratios = Inf, normalize = "positive_lobe")[, 1L, drop = FALSE]
  true_deep <- hrf_library(tr, peak_delays = 4.5, dispersions = 0.80,
    undershoot_ratios = 3.0, normalize = "positive_lobe")[, 1L, drop = FALSE]

  bold <- events <- cov_list <- vector("list", length(runs))
  names(bold) <- names(events) <- names(cov_list) <- runs
  for (r in seq_along(runs)) {
    onsets <- seq(12, 156, by = 8)
    mod <- as.numeric(scale(rnorm(length(onsets))))
    ev <- data.frame(onset = onsets, duration = 0, condition = "stim")
    cov_list[[r]] <- data.frame(onset = onsets, mod = mod)
    d_fast <- fmri.pipeline:::hrf_event_design_cpp(n_time, tr, 0, onsets, rep(0, length(onsets)),
      seq_along(onsets), length(onsets), true_fast)[[1L]]
    d_deep <- fmri.pipeline:::hrf_event_design_cpp(n_time, tr, 0, onsets, rep(0, length(onsets)),
      seq_along(onsets), length(onsets), true_deep)[[1L]]
    amps <- 1.0 + 0.45 * mod
    s_fast <- as.numeric(d_fast %*% amps)
    s_deep <- as.numeric(d_deep %*% amps)

    y <- matrix(0, nrow = n_time, ncol = 24L, dimnames = list(NULL, names(groups)))
    for (v in seq_len(24L)) {
      eps <- as.numeric(stats::filter(rnorm(n_time + 40L, sd = 0.08), 0.35,
        method = "recursive"))[-seq_len(40L)]
      sig <- if (v <= 8L) s_fast else if (v <= 16L) s_deep else 0
      y[, v] <- 100 + sig + eps
    }
    bold[[r]] <- y
    events[[r]] <- ev
  }

  lib_man <- hrf_library(tr, preset = "manifold")
  fit_man <- estimate_crossrun_hrf(
    bold = bold, events = events, tr = tr, kernel_grid = lib_man,
    spatial_groups = groups, pooling = "parcel", selection_mode = "soft_average"
  )

  expect_equal(fit_man$settings$kernel_normalize, "positive_lobe")
  expect_equal(colSums(pmax(fit_man$selected_kernels$run1, 0)),
    rep(1, 24L), tolerance = 1e-6, ignore_attr = TRUE)

  # Subcortical fast parcel recovers fast peak (~3s) with high waveform cosine (> 0.98)
  k_est_fast <- fit_man$selected_kernels$run1[, 1L]
  cos_fast <- sum(k_est_fast * true_fast[, 1L]) /
    sqrt(sum(k_est_fast^2) * sum(true_fast[, 1L]^2))
  expect_gt(cos_fast, 0.98)

  # Null parcel stays at canonical fallback
  expect_true(all(fit_man$selections$run1$status[17:24] == "fallback_prespecified"))

  # Stage-3 gain-restored modulator slope parity: positive_lobe prevents the ~33-50%
  # undershoot-ratio scale distortion seen under signed_sum normalization
  dec_man <- deconvolve_crossrun_hrf(bold, fit_man, lambda = 0.05, scaling = "demean")
  ep_pos <- decon_epoch_regression(
    deconvolution = dec_man, covariates = cov_list, events = events,
    target_condition = "stim", lags = c(0, 1), parcel_groups = groups,
    aggregate_by_parcel = TRUE, restore_scale = TRUE, demix_other_events = TRUE
  )
  cf_pos <- ep_pos$coefficients
  b_fast_pos <- cf_pos$beta_corrected[cf_pos$voxel == "subcort_fast" & cf_pos$term == "mod" & cf_pos$lag == 0]
  b_deep_pos <- cf_pos$beta_corrected[cf_pos$voxel == "cort_deep" & cf_pos$term == "mod" & cf_pos$lag == 0]
  rel_gap_pos <- abs(b_fast_pos - b_deep_pos) / max(abs(b_fast_pos), abs(b_deep_pos))
  expect_lt(rel_gap_pos, 0.15)

  lib_man_signed <- hrf_library(tr, preset = "manifold", normalize = "signed_sum")
  fit_signed <- estimate_crossrun_hrf(
    bold = bold, events = events, tr = tr, kernel_grid = lib_man_signed,
    spatial_groups = groups, pooling = "parcel", selection_mode = "soft_average"
  )
  dec_signed <- deconvolve_crossrun_hrf(bold, fit_signed, lambda = 0.05, scaling = "demean")
  ep_signed <- decon_epoch_regression(
    deconvolution = dec_signed, covariates = cov_list, events = events,
    target_condition = "stim", lags = c(0, 1), parcel_groups = groups,
    aggregate_by_parcel = TRUE, restore_scale = TRUE, demix_other_events = TRUE
  )
  cf_sig <- ep_signed$coefficients
  b_fast_sig <- cf_sig$beta_corrected[cf_sig$voxel == "subcort_fast" & cf_sig$term == "mod" & cf_sig$lag == 0]
  b_deep_sig <- cf_sig$beta_corrected[cf_sig$voxel == "cort_deep" & cf_sig$term == "mod" & cf_sig$lag == 0]
  rel_gap_sig <- abs(b_fast_sig - b_deep_sig) / max(abs(b_fast_sig), abs(b_deep_sig))
  expect_gt(rel_gap_sig, rel_gap_pos + 0.12)
})

test_that("hrf_library rejects undershoot_ratios <= 1 to prevent zero-sum or inverted kernels", {
  expect_error(
    hrf_library(1.0, undershoot_ratios = 1.0),
    "strictly greater than 1"
  )
  expect_error(
    hrf_library(1.0, undershoot_ratios = 0.75),
    "strictly greater than 1"
  )
  expect_error(
    hrf_library(1.0, undershoot_ratios = c(6, 1.0)),
    "strictly greater than 1"
  )
})

test_that("GCV lambda = NULL tuning and diagnose_crossrun_recovery preserve kernel_normalize", {
  set.seed(20261112)
  tr <- 1.0
  n_time <- 120L
  onsets <- seq(10, 106, by = 8)
  ev <- list(
    run1 = data.frame(onset = onsets, duration = 0, condition = "stim"),
    run2 = data.frame(onset = onsets + 1, duration = 0, condition = "stim")
  )

  # Single HRF with pronounced undershoot (ratio = 3) so positive_lobe and
  # signed_sum differ substantially in norm (~1.45x amplitude, ~2.1x H'H scale).
  lib_sig <- hrf_library(tr, peak_delays = 6, dispersions = 1, undershoot_ratios = 3,
    normalize = "signed_sum")
  lib_pos <- hrf_library(tr, peak_delays = 6, dispersions = 1, undershoot_ratios = 3,
    normalize = "positive_lobe")
  scale_c <- as.numeric(lib_pos[which.max(lib_pos[, 1L]), 1L] / lib_sig[which.max(lib_sig[, 1L]), 1L])
  expect_lt(scale_c, 0.80)

  groups <- setNames(rep("p1", 8L), paste0("v", seq_len(8L)))
  bold <- lapply(names(ev), function(rn) {
    d <- fmri.pipeline:::hrf_event_design_cpp(n_time, tr, 0, ev[[rn]]$onset,
      rep(0, nrow(ev[[rn]])), rep(1L, nrow(ev[[rn]])), 1L, lib_pos)[[1L]]
    y <- matrix(0, n_time, 8L, dimnames = list(NULL, names(groups)))
    for (v in seq_len(8L)) {
      y[, v] <- 100 + 1.5 * d[, 1L] + rnorm(n_time, sd = 0.25)
    }
    y
  })
  names(bold) <- names(ev)

  base_grid <- 10^seq(-1, 2, length.out = 9L)
  for (smode in c("hard", "soft_average")) {
    fit_sig <- estimate_crossrun_hrf(
      bold = bold, events = ev, tr = tr, kernel_grid = lib_sig,
      spatial_groups = groups, pooling = "parcel", selection_mode = smode,
      noise_model = "ols"
    )
    fit_pos <- estimate_crossrun_hrf(
      bold = bold, events = ev, tr = tr, kernel_grid = lib_pos,
      spatial_groups = groups, pooling = "parcel", selection_mode = smode,
      noise_model = "ols"
    )

    dec_sig <- deconvolve_crossrun_hrf(
      bold, fit_sig, lambda = NULL, lambda_grid = base_grid,
      ridge_floor = 0, scaling = "demean"
    )
    dec_pos <- deconvolve_crossrun_hrf(
      bold, fit_pos, lambda = NULL, lambda_grid = (scale_c^2) * base_grid,
      ridge_floor = 0, scaling = "demean"
    )

    # Because k_tune inside GCV honors kernel_normalize, scaling H by scale_c and
    # lambda_grid by scale_c^2 yields identical effective degrees of freedom & GCV!
    expect_equal(dec_pos$tuning$df, dec_sig$tuning$df, tolerance = 1e-5)
    expect_equal(dec_pos$tuning$gcv, dec_sig$tuning$gcv, tolerance = 1e-5)
    expect_equal(dec_pos$settings$lambda, (scale_c^2) * dec_sig$settings$lambda, tolerance = 1e-8)
  }

  # diagnose_crossrun_recovery matches fitted normalization even when bare
  # subsetted matrix (stripped of attributes by `[`) is supplied to source_kernel_grid
  fit_pos_hard <- estimate_crossrun_hrf(
    bold = bold, events = ev, tr = tr, kernel_grid = lib_pos,
    spatial_groups = groups, pooling = "parcel", selection_mode = "hard"
  )
  dec_pos_hard <- deconvolve_crossrun_hrf(bold, fit_pos_hard, lambda = 0.1, scaling = "demean")
  diag_default <- diagnose_crossrun_recovery(dec_pos_hard)
  bare_k1 <- fit_pos_hard$kernels[, 1L, drop = FALSE]
  expect_null(attr(bare_k1, "normalize"))
  diag_bare <- diagnose_crossrun_recovery(dec_pos_hard, source_kernel_grid = bare_k1)
  expect_equal(diag_bare$summary$mean_gain, diag_default$summary$mean_gain, tolerance = 1e-10)
})

test_that("cross-undershoot scale inflation is prevented broadly across custom libraries, unattributed grids, reglin, and impulse diagnostics", {
  tr <- 1.0

  # 1. hrf_library(normalize = "auto") switches to positive_lobe for custom multi-undershoot,
  #    Inf undershoot, or non-canonical undershoot ratios
  cust_multi <- hrf_library(tr, peak_delays = c(4, 6), dispersions = 1.0, undershoot_ratios = c(3, 12))
  cust_inf <- hrf_library(tr, peak_delays = 3.5, dispersions = 0.8, undershoot_ratios = Inf)
  cust_deep <- hrf_library(tr, peak_delays = 5.0, dispersions = 1.0, undershoot_ratios = 3)
  expect_equal(attr(cust_multi, "normalize"), "positive_lobe")
  expect_equal(attr(cust_inf, "normalize"), "positive_lobe")
  expect_equal(attr(cust_deep, "normalize"), "positive_lobe")
  expect_equal(colSums(pmax(cust_multi, 0)), rep(1, ncol(cust_multi)), tolerance = 1e-12, ignore_attr = TRUE)

  # 2. Unattributed multi-undershoot grid (e.g., cbind stripping attributes) automatically
  #    promotes to positive_lobe in deconvolve_reglin() and estimate_crossrun_hrf()
  h_can <- as.numeric(hrf_library(tr, peak_delays = 6, dispersions = 1, undershoot_ratios = 6)[, 1L])
  h_inf <- as.numeric(cust_inf[, 1L])
  h_deep <- as.numeric(cust_deep[, 1L])
  bare_multi <- cbind(can = h_can, subcort = h_inf, deep = h_deep)
  expect_null(attr(bare_multi, "normalize"))

  set.seed(404)
  n_time <- 90L
  onsets <- seq(8, 76, by = 8)
  ev <- list(
    run1 = data.frame(onset = onsets, duration = 0, condition = "stim"),
    run2 = data.frame(onset = onsets + 1, duration = 0, condition = "stim")
  )
  bold <- lapply(ev, function(e) {
    matrix(100 + rnorm(n_time * 4L, sd = 0.3), n_time, 4L,
      dimnames = list(NULL, paste0("v", 1:4)))
  })

  fit_bare <- estimate_crossrun_hrf(bold, ev, tr, bare_multi, noise_model = "ols")
  expect_equal(fit_bare$settings$kernel_normalize, "positive_lobe")
  expect_equal(colSums(pmax(fit_bare$kernels, 0)), rep(1, 3L), tolerance = 1e-12, ignore_attr = TRUE)

  reg_bare <- deconvolve_reglin(bold$run1, kernel = h_can, kernel_grid = bare_multi,
    lambda = 0.1, return_diagnostics = TRUE)
  expect_equal(reg_bare$kernel_normalize, "positive_lobe")
  expect_equal(vapply(reg_bare$kernels, function(k) sum(pmax(k, 0)), numeric(1L)),
    rep(1, 3L), tolerance = 1e-12)

  # Wiener backend inside deconvolve_crossrun_hrf also honors kernel_normalize
  dec_wiener <- deconvolve_crossrun_hrf(bold, fit_bare, lambda = 0.1,
    backend = "wiener", scaling = "demean",
    nuisance_mode = "ignore", observation_mask = "all")
  ref_w <- deconvolve_reglin(bold$run1[, 1L] - mean(bold$run1[, 1L]),
    fit_bare$kernels[, fit_bare$selections$run1$kernel_index[1L]],
    lambda = 0.1, scaling = "none", backend = "wiener",
    kernel_normalize = "positive_lobe")[, 1L]
  expect_equal(as.numeric(dec_wiener$activity$run1[, 1L]), as.numeric(ref_w), tolerance = 1e-10)

  # 3. deconvolution_impulse_diagnostics() auto-switches to positive_lobe when
  #    kernel and forward_kernel differ in undershoot depth
  ev_imp <- data.frame(trial = c("t1", "t2"), onset = c(12, 36))
  imp_diag <- deconvolution_impulse_diagnostics(
    event_data = ev_imp, kernel = h_can, forward_kernel = h_inf,
    n_time = 70L, tr = tr, lambda = 0.05
  )
  expect_equal(imp_diag$settings$kernel_normalize, "positive_lobe")
  expect_equal(sum(pmax(imp_diag$settings$kernel, 0)), 1.0, tolerance = 1e-12)
  expect_equal(sum(pmax(imp_diag$settings$forward_kernel, 0)), 1.0, tolerance = 1e-12)
})
