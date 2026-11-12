test_that("P1-1: hrf_build_noise_mean_design handles empty events and degrades monotonically in short/dense runs", {
  tr <- 1.0
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  k_vec <- lib[, 1L]

  n_obs <- 60L
  z_base <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n_obs))
  mask_all <- rep(TRUE, n_obs)

  # 1. NULL and 0-row events return full structured list
  res_null <- hrf_build_noise_mean_design(n_obs, NULL, tr, 0, k_vec, z_base, mask_all)
  ev_empty <- data.frame(onset = numeric(0), duration = numeric(0), condition = character(0))
  res_empty <- hrf_build_noise_mean_design(n_obs, ev_empty, tr, 0, k_vec, z_base, mask_all)

  for (res in list(res_null, res_empty)) {
    expect_true(is.list(res))
    expect_true(all(c("design", "df_nominal", "trial_design", "used_trial_design") %in% names(res)))
    expect_equal(nrow(res$design), n_obs)
    expect_equal(res$df_nominal, n_obs - 2L)
    expect_equal(ncol(res$trial_design), 0L)
    expect_true(isTRUE(res$used_trial_design))
  }

  # 2. Monotonic fallback ladder as trial density increases relative to masked frames
  ev_20 <- data.frame(onset = seq(2, by = 2.5, length.out = 20L), duration = 0, condition = "fb")
  res_trial <- hrf_build_noise_mean_design(n_obs, ev_20, tr, 0, k_vec, z_base, mask_all)
  expect_true(res_trial$used_trial_design)
  expect_gte(res_trial$df_nominal, 20L)

  # Dense trials (44 trials in 60 frames): trial-wise design has < 20 DF -> falls back cleanly to condition design
  ev_44 <- data.frame(onset = seq(1, by = 1.25, length.out = 44L), duration = 0, condition = "fb")
  res_cond <- hrf_build_noise_mean_design(n_obs, ev_44, tr, 0, k_vec, z_base, mask_all)
  expect_false(res_cond$used_trial_design)
  expect_gte(res_cond$df_nominal, 20L)
  expect_equal(ncol(res_cond$trial_design), 44L)

  # Heavily scrubbed run (18 valid frames) with 4 conditions -> falls back to z + x_cond (12 <= df < 20) or z
  mask_short <- c(rep(TRUE, 18L), rep(FALSE, n_obs - 18L))
  ev_short <- data.frame(
    onset = c(1, 4, 7, 10, 13, 15),
    duration = 0,
    condition = c("c1", "c2", "c3", "c4", "c1", "c2")
  )
  res_short <- hrf_build_noise_mean_design(n_obs, ev_short, tr, 0, k_vec, z_base, mask_short)
  expect_gte(res_short$df_nominal, 12L)

  # Ultra-short mask (10 valid frames) -> falls back safely to z_base
  mask_tiny <- c(rep(TRUE, 10L), rep(FALSE, n_obs - 10L))
  res_tiny <- hrf_build_noise_mean_design(n_obs, ev_short, tr, 0, k_vec, z_base, mask_tiny)
  expect_equal(ncol(res_tiny$design), ncol(z_base))
  expect_equal(res_tiny$df_nominal, 8L)
})

test_that("P1-2: soft_average with pooling = 'voxel' synchronizes chosen, status, kernel_weights, and parcel fallback", {
  set.seed(20260502L)
  tr <- 1.0
  n_time <- 160L
  n_trials <- 22L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_early <- which(abs(pars$peak_delay - 4.5) < 1e-6)[1L]
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6)[1L]
  idx_late <- which(abs(pars$peak_delay - 7.5) < 1e-6)[1L]

  # Parcel P1 has 5 strong early-HRF voxels, 1 pure-noise low-R2 voxel (v6),
  # and 1 cross-run incompatible voxel (v7: early in run1, late in run2-3).
  # Parcel P2 has 4 pure-noise voxels (unsupported parcel).
  groups <- stats::setNames(c(rep("P1", 7L), rep("P2", 4L)), paste0("v", 1:11))
  bold <- events <- list()
  for (r in 1:3) {
    t_ev <- seq(6, by = 6.5, length.out = n_trials)
    amps <- 1.0 + 0.45 * as.numeric(scale(stats::rnorm(n_trials)))
    des_e <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
      seq_len(n_trials), n_trials, lib[, idx_early, drop = FALSE])[[1L]]
    des_l <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials),
      seq_len(n_trials), n_trials, lib[, idx_late, drop = FALSE])[[1L]]
    ymat <- matrix(0.0, n_time, 11L, dimnames = list(NULL, names(groups)))
    for (v in 1:5) ymat[, v] <- 100 + as.numeric(des_e %*% amps) + stats::rnorm(n_time, sd = 0.12)
    ymat[, 6L] <- 100 + stats::rnorm(n_time, sd = 0.50) # low-R2 voxel in supported parcel P1
    ymat[, 7L] <- 100 + as.numeric((if (r == 1L) des_e else des_l) %*% amps) + stats::rnorm(n_time, sd = 0.03)
    for (v in 8:11) ymat[, v] <- 100 + stats::rnorm(n_time, sd = 0.50) # unsupported parcel P2
    bold[[paste0("run", r)]] <- ymat
    events[[paste0("run", r)]] <- data.frame(onset = t_ev, duration = 0, condition = "fb")
  }
  # Test run 3 (trained on runs 1 and 2, where v7 is early vs late -> incompatible!)
  fit_vox <- estimate_crossrun_hrf(
    bold, events, tr, lib,
    spatial_groups = groups, pooling = "voxel", min_parcel_voxels = 3L,
    fallback_kernel = idx_can, selection_mode = "soft_average", noise_model = "ols"
  )
  sel3 <- fit_vox$selections$run3
  w3 <- fit_vox$kernel_weights$run3
  k3 <- fit_vox$selected_kernels$run3

  # Global synchronization invariant: kernel_index == max-weight index for every voxel
  expect_equal(sel3$kernel_index, max.col(t(w3), ties.method = "first"))
  expect_equal(unname(colSums(w3)), rep(1.0, 11L), tolerance = 1e-10)
  expect_equal(unname(colSums(k3)), rep(1.0, 11L), tolerance = 1e-10)

  # Low-R2 voxel v6 inside reliable parcel P1 retains parcel soft-average anchor (idx_early), NOT wiped to idx_can
  expect_equal(sel3$kernel_index[6L], idx_early)
  expect_equal(sel3$status[6L], "fallback_parcel")
  expect_gt(w3[idx_early, 6L], 0.50)

  # Incompatible voxel v7 falls back cleanly to fallback_kernel in BOTH kernel_weights and selections
  expect_equal(sel3$kernel_index[7L], idx_can)
  expect_equal(sel3$status[7L], "fallback_prespecified")
  expect_equal(as.numeric(w3[, 7L]), as.numeric(seq_len(ncol(lib)) == idx_can))

  # Unsupported parcel P2 voxels (v8:v11) have synchronized fallback_prespecified status and weights
  expect_equal(sel3$kernel_index[8:11], rep(idx_can, 4L))
  expect_equal(sel3$status[8:11], rep("fallback_prespecified", 4L))
})

test_that("P1-3: per-parcel focal-condition screening recovers zero-mean RPE HRFs and single-run subblocks remain calibrated", {
  set.seed(20260503L)
  tr <- 1.0
  n_time <- 190L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_early <- which(abs(pars$peak_delay - 4.5) < 1e-6)[1L]
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6)[1L]

  # 16 visual voxels respond to 'cue' (canonical HRF, trial-varying amplitude)
  # 4 striatal voxels have ZERO mean and pure trial-varying RPEs on closely-coupled 'feedback' (+2.2s) with EARLY HRF (4.5s)
  groups <- stats::setNames(c(rep("R_vis", 16L), rep("R_str", 4L)), paste0("v", 1:20))
  bold <- events <- list()
  for (r in 1:2) {
    t_cue <- seq(7, by = 7.2, length.out = n_trials)
    t_fb <- t_cue + 2.2
    cue_amp <- 1.2 + 0.70 * as.numeric(scale(stats::rnorm(n_trials)))
    rpe_amp <- 0.85 * as.numeric(scale(stats::rnorm(n_trials))) # exact zero mean!
    des_cue_can <- hrf_event_design_cpp(n_time, tr, 0, t_cue, rep(0, n_trials),
      seq_len(n_trials), n_trials, lib[, idx_can, drop = FALSE])[[1L]]
    des_fb_early <- hrf_event_design_cpp(n_time, tr, 0, t_fb, rep(0, n_trials),
      seq_len(n_trials), n_trials, lib[, idx_early, drop = FALSE])[[1L]]

    ymat <- matrix(0.0, n_time, 20L, dimnames = list(NULL, names(groups)))
    for (v in 1:16) ymat[, v] <- 100 + as.numeric(des_cue_can %*% cue_amp) + stats::rnorm(n_time, sd = 0.12)
    for (v in 17:20) ymat[, v] <- 100 + as.numeric(des_fb_early %*% rpe_amp) + stats::rnorm(n_time, sd = 0.10)
    bold[[paste0("run", r)]] <- ymat
    events[[paste0("run", r)]] <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue"),
      data.frame(onset = t_fb, duration = 0, condition = "feedback")
    )
  }

  fit_grp <- estimate_crossrun_hrf(
    bold, events, tr, lib,
    spatial_groups = groups, pooling = "parcel", min_parcel_voxels = 3L,
    fallback_kernel = idx_can, amplitude_model = "trial", noise_model = "ols"
  )
  expect_equal(fit_grp$selections$run1$kernel_index[1:16], rep(idx_can, 16L))
  expect_equal(fit_grp$selections$run1$kernel_index[17:20], rep(idx_early, 4L))

  # Ungrouped voxel mode also recovers zero-mean RPE voxels
  fit_ung <- estimate_crossrun_hrf(
    bold, events, tr, lib,
    pooling = "voxel", fallback_kernel = idx_can,
    amplitude_model = "trial", noise_model = "ols"
  )
  expect_equal(fit_ung$selections$run1$kernel_index[17:20], rep(idx_early, 4L))

  # Single-run sub-block calibration under AR(1) null noise (rho = 0.50) + boundary-straddling events
  t_single <- seq(5, 173, by = 6.0)
  ev_1 <- list(run1 = data.frame(onset = t_single, duration = 0, condition = "fb"))
  W_1 <- hrf_whitening_matrix(seq_len(180L), 0.50)
  n_null_vox <- 80L
  y_null_1 <- matrix(0.0, 180L, n_null_vox, dimnames = list(NULL, paste0("nv", seq_len(n_null_vox))))
  for (v in seq_len(n_null_vox)) {
    y_null_1[, v] <- 100 + sqrt(1 - 0.50^2) * as.numeric(solve(W_1, stats::rnorm(180L)))
  }
  fit_null_1 <- estimate_crossrun_hrf(
    list(run1 = y_null_1), ev_1, tr, lib,
    fallback_kernel = idx_can, amplitude_model = "trial", noise_model = "ar1"
  )
  fp_rate_1 <- mean(fit_null_1$selections$run1$kernel_index != idx_can)
  expect_lte(fp_rate_1, 0.08)
})

test_that("P1-4: build_other_demix_operator handles NA conditions, edge-truncated events, and non-impulse trial-varying confounds", {
  set.seed(20260504L)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 22L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]

  # Cue precedes choice by 1.8s; cue has strong trial-varying amplitude correlated r = 0.75 with x,
  # while choice has ZERO true modulation (true_beta = 0) and sustained RT duration!
  # Also include an edge cue at t = 179.6s (truncated at run boundary -> !all(ok_all)) and an NA-condition distractor.
  bold_dur <- ev_dur <- x_list <- list()
  for (r in 1:2) {
    t_cue <- seq(6, by = 7.2, length.out = n_trials)
    t_choice <- t_cue + 1.8
    x <- as.numeric(scale(stats::rnorm(n_trials)))
    z_ind <- as.numeric(scale(stats::rnorm(n_trials)))
    amp_cue <- 1.5 + 1.10 * (0.75 * x + sqrt(1 - 0.75^2) * z_ind)
    rt <- runif(n_trials, 0.8, 1.8)
    amp_choice <- rep(1.0, n_trials) # TRUE choice modulation = 0!

    des_cue <- hrf_event_design_cpp(n_time, tr, 0, t_cue, rep(0, n_trials),
      seq_len(n_trials), n_trials, hrf_k)[[1L]]
    des_ch <- hrf_event_design_cpp(n_time, tr, 0, t_choice, rt,
      seq_len(n_trials), n_trials, hrf_k)[[1L]]
    sig <- as.numeric(des_cue %*% amp_cue + des_ch %*% amp_choice)
    ymat <- replicate(6L, 100 + sig + stats::rnorm(n_time, sd = 0.12))
    colnames(ymat) <- paste0("v", 1:6)
    bold_dur[[paste0("run", r)]] <- ymat
    ev_dur[[paste0("run", r)]] <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue", stringsAsFactors = FALSE),
      data.frame(onset = t_choice, duration = rt, condition = "choice", stringsAsFactors = FALSE),
      data.frame(onset = c(3.0, 179.8), duration = c(0, 0), condition = c(NA_character_, "cue"), stringsAsFactors = FALSE)
    )
    x_list[[paste0("run", r)]] <- x
  }

  # Fit HRF on clean non-NA copy (estimate_crossrun_hrf validates event table), then test Stage 3 with NA & edge events
  ev_clean <- lapply(ev_dur, function(df) {
    df$condition[is.na(df$condition)] <- "distractor"
    df
  })
  fit <- estimate_crossrun_hrf(bold_dur, ev_clean, tr, lib, noise_model = "ols")
  dec <- deconvolve_crossrun_hrf(bold_dur, fit, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")

  # Condition-mean other-event demixing (demix_ridge = 0.05) leaves large cue bleed into choice;
  # joint trial-varying demixing (demix_ridge = 0.02) works cleanly with NA condition + edge event + duration/window models!
  reg_mean <- decon_epoch_regression(
    dec, x_list, events = ev_dur, target_condition = "choice",
    lags = 0, source_model = "duration", demix_ridge = 0.05, restore_scale = TRUE
  )
  reg_joint_dur <- decon_epoch_regression(
    dec, x_list, events = ev_dur, target_condition = "choice",
    lags = 0, source_model = "duration", demix_ridge = 0.02, restore_scale = TRUE
  )
  reg_joint_win <- decon_epoch_regression(
    dec, x_list, events = ev_dur, target_condition = "choice",
    lags = 0, source_model = "window", demix_ridge = 0.02, restore_scale = TRUE
  )
  reg_joint_imp <- decon_epoch_regression(
    dec, x_list, events = ev_dur, target_condition = "choice",
    lags = 0, source_model = "impulse", demix_ridge = 0.02, restore_scale = TRUE
  )

  expect_true(all(is.finite(reg_joint_dur$coefficients$beta_corrected)))
  expect_true(all(is.finite(reg_joint_win$coefficients$beta_corrected)))
  expect_true(all(is.finite(reg_joint_imp$coefficients$beta_corrected)))
  expect_lt(
    abs(mean(reg_joint_dur$coefficients$beta_corrected)),
    0.45 * abs(mean(reg_mean$coefficients$beta_corrected))
  )
})

test_that("P1-5: duration source_model does not inject +alpha0 baseline into pre-onset or post-offset lags", {
  set.seed(20260505L)
  tr <- 1.0
  n_time <- 220L
  n_trials <- 18L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  alpha0 <- 2.2

  bold <- events <- x_list <- list()
  for (r in 1:2) {
    t_ev <- seq(14, by = 10.5, length.out = n_trials)
    x <- as.numeric(scale(stats::rnorm(n_trials)))
    rt <- pmin(2.4, pmax(0.5, 1.4 - 0.45 * x + stats::rnorm(n_trials, sd = 0.25)))
    des <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rt, seq_len(n_trials), n_trials, hrf_k)[[1L]]
    sig <- as.numeric(des %*% rep(alpha0, n_trials))
    ymat <- replicate(4L, 100 + sig + stats::rnorm(n_time, sd = 0.08))
    colnames(ymat) <- paste0("v", 1:4)
    bold[[paste0("run", r)]] <- ymat
    events[[paste0("run", r)]] <- data.frame(onset = t_ev, duration = rt, condition = "choice")
    x_list[[paste0("run", r)]] <- x
  }

  fit <- estimate_crossrun_hrf(bold, events, tr, lib, noise_model = "ols")
  dec <- deconvolve_crossrun_hrf(bold, fit, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")
  reg <- decon_epoch_regression(
    dec, x_list, events = events, target_condition = "choice",
    lags = c(-4, 0, 5), source_model = "duration", restore_scale = TRUE, return_panel = TRUE
  )

  p_df <- reg$panel
  m_pre <- mean(p_df$decon_demixed[p_df$evt_time == -4], na.rm = TRUE)
  m_on <- mean(p_df$decon_demixed[p_df$evt_time == 0], na.rm = TRUE)
  m_post <- mean(p_df$decon_demixed[p_df$evt_time == 5], na.rm = TRUE)

  # Pre-onset (-4s) and post-offset (+5s) mean activity remain near 0 rather than +alpha0 (2.2),
  # while lag 0 recovers +alpha0 and RT-confounded slope remains near 0 across all lags!
  expect_lt(abs(m_pre), 0.35)
  expect_lt(abs(m_post), 0.45)
  expect_gt(m_on, 1.5)
  expect_lt(max(abs(reg$summary$mean_beta_corrected)), 0.08)

  # Zero-offset BOLD (ym - 100) with scaling = "none" also shows zero baseline bleed at pre/post lags
  bold_zero_dc <- lapply(bold, function(ym) ym - 100)
  dec_none <- deconvolve_crossrun_hrf(bold_zero_dc, fit, lambda = 0.8, scaling = "none", penalty = "sobolev")
  reg_none <- decon_epoch_regression(
    dec_none, x_list, events = events, target_condition = "choice",
    lags = c(-4, 0, 5), source_model = "duration", restore_scale = FALSE, return_panel = TRUE
  )
  expect_lt(abs(mean(reg_none$panel$decon_demixed[reg_none$panel$evt_time == -4], na.rm = TRUE)), 0.35)
  expect_lt(abs(mean(reg_none$panel$decon_demixed[reg_none$panel$evt_time == 5], na.rm = TRUE)), 0.45)
})

test_that("P1-6: short-run and NA-covariate subsetting preserve exact centered error covariance and unbiased REPAIR EB weights", {
  set.seed(20260506L)
  tr <- 1.0
  n_time <- 90L
  n_trials <- 6L # short run (m_r = 6) where m/(m-1) = 1.20 would cause 20% over-subtraction
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]

  n_subj <- 16L
  panels <- vector("list", n_subj)
  tau2_sub <- numeric(n_subj)
  for (s in seq_len(n_subj)) {
    bold <- events <- covs <- list()
    for (r in 1:2) {
      t_ev <- seq(10, by = 11.0, length.out = n_trials)
      x <- as.numeric(scale(stats::rnorm(n_trials)))
      if (r == 1L) x[3L] <- NA_real_ # induce m_r = 5 < m_ok = 6 subsetting!
      ymat <- matrix(100 + stats::rnorm(n_time * 2L, sd = 0.40), nrow = n_time, ncol = 2L,
        dimnames = list(NULL, c("R_seed", "R_targ")))
      bold[[paste0("run", r)]] <- ymat
      events[[paste0("run", r)]] <- data.frame(onset = t_ev, duration = 0, condition = "ev")
      covs[[paste0("run", r)]] <- data.frame(x = x)
    }
    fit <- estimate_crossrun_hrf(bold, events, tr, lib, fallback_kernel = 1L, noise_model = "ols")
    dec <- deconvolve_crossrun_hrf(bold, fit, lambda = 1.0, scaling = "noise_sd", penalty = "sobolev")
    reg <- decon_epoch_regression(
      dec, covs, events = events, target_condition = "ev", lags = 0,
      whiten_panel = TRUE, return_panel = TRUE, restore_scale = TRUE
    )
    tau2_sub[s] <- mean(reg$coefficients$tau2_bio, na.rm = TRUE)
    p_s <- reg$panel
    p_s$id <- sprintf("sub%02d", s)
    p_s$parcel <- p_s$voxel
    panels[[s]] <- p_s
  }
  # Under pure noise (true tau2_bio = 0), median tau2_bio across short runs is near 0
  expect_lt(stats::median(tau2_sub), 0.05)

  # Whitened SEs match unwhitened SEs closely on well-separated trials even when m_r = 5 < m_ok = 6
  p1 <- panels[[1L]]
  ok_r1 <- p1$run == "run1" & is.finite(p1$se_whitened) & is.finite(p1$se_demixed)
  expect_equal(sum(ok_r1), 10L) # 5 non-NA trials x 2 parcels
  expect_lt(max(abs(p1$se_whitened[ok_r1] / p1$se_demixed[ok_r1] - 1.0)), 0.15)

  # REPAIR EB predictor mode under pure noise remains unbiased around 0 on short runs
  full_panel <- do.call(rbind, panels)
  rep_eb <- decon_connectivity_repair(
    full_panel, seed_region = "R_seed", target_region = "R_targ",
    seed_lag = 0, target_lag = 0, covariates = "x", method = "eb_predictor"
  )
  expect_true(is.finite(rep_eb$group$estimate[1L]))
  expect_lt(abs(rep_eb$group$estimate[1L]), 0.25)
})

test_that("P1-7: deconvolve_reglin updates noise_sd after Pass-1 prewhiten_rho, is batch-invariant, and supports signal 1-SE GCV", {
  set.seed(20260507L)
  n_time <- 240L
  hrf <- as.numeric(spm_hrf(1.0)$hrf)
  hrf <- hrf / sum(hrf)
  conv_m <- decon_convolution_matrix(n_time, hrf)

  # 1. Pass-1 prewhiten_rho updates scaling = "noise_sd" divisor
  s_true <- numeric(n_time)
  s_true[seq(12, 220, by = 9)] <- 1.4
  W_ar <- hrf_whitening_matrix(seq_len(n_time), 0.55)
  Y <- matrix(0.0, n_time, 6L)
  for (j in 1:6) {
    eps <- 0.35 * sqrt(1 - 0.55^2) * as.numeric(solve(W_ar, stats::rnorm(n_time)))
    Y[, j] <- 100 + as.numeric(conv_m[seq_len(n_time), seq_len(n_time) + length(hrf) - 1L] %*% s_true) + eps
  }
  fit_auto <- deconvolve_reglin(
    Y, kernel = hrf, lambda_grid = c(0.01, 0.1, 1.0, 5.0), penalty = "sobolev",
    scaling = "noise_sd", prewhiten_gcv = TRUE, prewhiten_rho = NULL,
    backend = "cpp", return_diagnostics = TRUE
  )
  fit_explicit <- deconvolve_reglin(
    Y, kernel = hrf, lambda_grid = c(0.01, 0.1, 1.0, 5.0), penalty = "sobolev",
    scaling = "noise_sd", prewhiten_gcv = TRUE, prewhiten_rho = fit_auto$ar1_rho,
    backend = "cpp", return_diagnostics = TRUE
  )
  expect_equal(as.numeric(fit_auto$scaled_scale), as.numeric(fit_explicit$scaled_scale), tolerance = 1e-10)
  expect_equal(fit_auto$activity, fit_explicit$activity, ignore_attr = TRUE, tolerance = 1e-8)

  # 2. Strict column-batch invariance of estimate_decon_noise_sd at n >= 200 with supplied noise_rho
  Y_het <- cbind(
    Y[, 1:2],
    100 + 2.5 * sin(seq_len(n_time) / 18) + matrix(stats::rnorm(n_time * 4L, sd = 0.25), n_time, 4L)
  )
  rhos_het <- c(0.55, 0.55, 0.70, 0.70, 0.65, 0.40)
  sd_all <- estimate_decon_noise_sd(Y_het, noise_rho = rhos_het)
  sd_each <- vapply(seq_len(ncol(Y_het)), function(j) {
    estimate_decon_noise_sd(Y_het[, j, drop = FALSE], noise_rho = rhos_het[j])
  }, numeric(1L))
  expect_equal(sd_all, sd_each, tolerance = 1e-13)

  # 3. Per-signal 1-SE GCV rule (tune_by = "signal", gcv_rule = "1se") matches across R and C++ and >= min rule
  lam_grid <- 10^seq(-3, 1, length.out = 9L)
  min_cpp <- deconvolve_reglin(
    Y, kernel = hrf, lambda_grid = lam_grid, penalty = "sobolev",
    tune_by = "signal", gcv_rule = "min", backend = "cpp", return_diagnostics = TRUE
  )
  onese_cpp <- deconvolve_reglin(
    Y, kernel = hrf, lambda_grid = lam_grid, penalty = "sobolev",
    tune_by = "signal", gcv_rule = "1se", backend = "cpp", return_diagnostics = TRUE
  )
  onese_r <- deconvolve_reglin(
    Y, kernel = hrf, lambda_grid = lam_grid, penalty = "sobolev",
    tune_by = "signal", gcv_rule = "1se", backend = "r", return_diagnostics = TRUE
  )
  onese_w <- deconvolve_reglin(
    Y, kernel = hrf, lambda_grid = lam_grid,
    tune_by = "signal", gcv_rule = "1se", backend = "wiener", return_diagnostics = TRUE
  )
  min_w <- deconvolve_reglin(
    Y, kernel = hrf, lambda_grid = lam_grid,
    tune_by = "signal", gcv_rule = "min", backend = "wiener", return_diagnostics = TRUE
  )

  expect_equal(as.numeric(onese_cpp$lambda), as.numeric(onese_r$lambda), tolerance = 1e-12)
  expect_equal(onese_cpp$activity, onese_r$activity, ignore_attr = TRUE, tolerance = 1e-8)
  expect_true(all(as.numeric(onese_cpp$lambda) >= as.numeric(min_cpp$lambda)))
  expect_true(any(as.numeric(onese_cpp$lambda) > as.numeric(min_cpp$lambda)))
  expect_true(all(as.numeric(onese_w$lambda) >= as.numeric(min_w$lambda)))
})

test_that("P1-8: parse_decon_identity ignores directory words, voxelwise_deconvolution deduplicates metadata, and mixed_by normalizes complete-case weights", {
  # 1. parse_decon_identity handles /subjects/, /sessions/, /runs/, /runtime/, flattened filenames, and directory fallbacks
  p_bids <- "/Users/lab/subjects/sub-01/sessions/ses-02/runs/sub-01_ses-02_run-03_bold.nii.gz"
  expect_equal(parse_decon_identity(p_bids), list(id = "01", run = "03", session = "02"))

  p_flat <- "decon_.Users.lab.subjects.sub-01.sessions.ses-02.runs.sub-01_ses-02_run-03_bold.nii.gz_schaefer_nii_atlas_deconvolved.csv.gz"
  expect_equal(parse_decon_identity(p_flat), list(id = "01", run = "03", session = "02"))

  p_dir_only <- "/scratch/runtime/sub-104/ses-baseline/run-2/bold_preproc.nii.gz"
  expect_equal(parse_decon_identity(p_dir_only), list(id = "104", run = "2", session = "baseline"))

  p_parent_pilot <- "/scratch/run2_pilot/sub-09_run-01_bold.nii.gz"
  expect_equal(parse_decon_identity(p_parent_pilot), list(id = "09", run = "01"))

  expect_equal(parse_decon_identity("subA_run1_deconvolved.csv.gz"), list(id = "A", run = "1"))
  expect_error(parse_decon_identity("sub-A_sub-B_run1.csv"), "Ambiguous id")

  # 2. voxelwise_deconvolution escapes .nii(.gz)?$ properly and deduplicates metadata on re-run
  tmp_dir <- tempfile("p1_vox_decon_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  arr4d <- array(100 + stats::rnorm(2 * 2 * 1 * 20, sd = 2), dim = c(2, 2, 1, 20))
  atlas_arr <- array(c(1L, 1L, 0L, 0L), dim = c(2, 2, 1))
  nii_file <- file.path(tmp_dir, "sub-01_run-01_bold.nii.gz")
  # Name contains '_nii_' inside stem to test regex escaping!
  atlas_file <- file.path(tmp_dir, "schaefer_nii_mask.nii.gz")
  RNifti::writeNifti(arr4d, nii_file)
  RNifti::writeNifti(atlas_arr, atlas_file)

  out_dir <- file.path(tmp_dir, "out")
  meta_df <- data.frame(id = "01", run = "01", group = "ctrl", stringsAsFactors = FALSE)
  decon_set <- list(kernel = c(0.1, 0.5, 0.3, 0.1), lambda = 0.5)

  voxelwise_deconvolution(
    niftis = nii_file, add_metadata = meta_df, out_dir = out_dir,
    TR = 1.0, atlas_files = atlas_file, algorithm = "reglin",
    decon_settings = decon_set, nprocs = 1L
  )
  # Re-run on same out_dir: should skip existing deconvolved CSV and NOT duplicate metadata row
  voxelwise_deconvolution(
    niftis = nii_file, add_metadata = meta_df, out_dir = out_dir,
    TR = 1.0, atlas_files = atlas_file, algorithm = "reglin",
    decon_settings = decon_set, nprocs = 1L
  )

  expected_subdir <- file.path(out_dir, "schaefer_nii_mask")
  expect_true(dir.exists(expected_subdir))
  meta_path <- file.path(expected_subdir, "schaefer_nii_mask_metadata.csv")
  expect_true(file.exists(meta_path))
  meta_read <- readr::read_csv(meta_path, show_col_types = FALSE)
  expect_equal(nrow(meta_read), 1L)

  # 3. mixed_by normalizes weights to mean 1.0 on complete.cases rows when NAs exist in outcome/covariates
  set.seed(20260508L)
  df_mb <- expand.grid(trial = 1:20, id = paste0("s", 1:8), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df_mb$evt_time <- 0
  df_mb$x <- stats::rnorm(nrow(df_mb))
  df_mb$y <- 0.5 * df_mb$x + rep(stats::rnorm(8, sd = 0.3), each = 20) + stats::rnorm(nrow(df_mb), sd = 0.4)
  # Extreme weights on rows that have NA in x or y would previously skew mean(w_vec) away from 1.0!
  df_mb$wt <- runif(nrow(df_mb), 0.5, 1.5)
  na_rows <- 1:25
  df_mb$x[1:12] <- NA_real_
  df_mb$y[13:25] <- NA_real_
  df_mb$wt[na_rows] <- 500.0

  res_mb <- mixed_by(
    data = df_mb,
    outcomes = "y",
    rhs_model_formulae = ~ x + (1 | id),
    split_on = "evt_time",
    weights = "wt",
    padjust_by = NULL,
    return_models = TRUE,
    calculate = c("parameter_estimates_reml")
  )
  fit_lmer <- res_mb$models$model[[1L]]
  while (is.list(fit_lmer)) fit_lmer <- fit_lmer[[1L]]
  used_wts <- stats::weights(fit_lmer)
  expect_equal(mean(fit_lmer@resp$weights), 1.0, tolerance = 1e-10)
  expect_equal(mean(used_wts, na.rm = TRUE), 1.0, tolerance = 1e-10)

  # Matches fitting directly on complete rows where weights have mean 1.0
  df_comp <- df_mb[stats::complete.cases(df_mb[, c("y", "x", "id")]), ]
  df_comp$wt <- df_comp$wt / mean(df_comp$wt)
  ref_lmer <- lmerTest::lmer(y ~ x + (1 | id), data = df_comp, weights = df_comp$wt)
  expect_equal(as.numeric(lme4::VarCorr(fit_lmer)$id), as.numeric(lme4::VarCorr(ref_lmer)$id), tolerance = 1e-8)
  expect_equal(stats::sigma(fit_lmer), stats::sigma(ref_lmer), tolerance = 1e-8)
})
