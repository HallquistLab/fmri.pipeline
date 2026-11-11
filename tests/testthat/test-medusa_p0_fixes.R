test_that("P0-1: post-demixing covariate operator M_post eliminates autocorrelated-covariate bias and calibrates SEs in 1-stage and 2-stage workflows", {
  set.seed(20260401L)
  tr <- 1.0
  n_time <- 190L
  n_trials <- 42L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]

  # Rapid event-related design (ISI 2.3-3.6s) with strongly autocorrelated covariate (rho_x = 0.70)
  n_sched <- 16L
  n_rep_per_sched <- 12L
  beta_true <- 0.80

  est_act <- se_act <- df_act <- numeric(0L)
  est_null <- se_null <- df_null <- numeric(0L)
  est_2stage <- numeric(0L)

  for (sc in seq_len(n_sched)) {
    bold <- events <- covs <- list()
    for (r in 1:2) {
      isis <- stats::runif(n_trials - 1L, 2.3, 3.6)
      onsets <- cumsum(c(8.0, isis))
      x_ar <- as.numeric(stats::arima.sim(model = list(ar = 0.70), n = n_trials))
      x_ar <- as.numeric(scale(x_ar))
      G_r <- hrf_event_design_cpp(
        n_time, tr, 0, onsets, rep(0, n_trials), seq_len(n_trials), n_trials, hrf_k
      )[[1L]]

      n_vox <- 2L * n_rep_per_sched
      Y_r <- matrix(0.0, n_time, n_vox, dimnames = list(NULL, paste0("v", seq_len(n_vox))))
      for (vi in seq_len(n_vox)) {
        eps_v <- as.numeric(stats::arima.sim(model = list(ar = 0.35), n = n_time, sd = 0.30))
        b_v <- if (vi <= n_rep_per_sched) beta_true else 0.0
        amp_v <- 0.50 + b_v * x_ar + stats::rnorm(n_trials, 0, 0.20)
        Y_r[, vi] <- 100.0 + as.numeric(G_r %*% amp_v) + eps_v
      }
      rn <- paste0("run", r)
      bold[[rn]] <- Y_r
      events[[rn]] <- data.frame(onset = onsets, duration = 0, condition = "cue")
      covs[[rn]] <- data.frame(x = x_ar)
    }

    fit_h <- estimate_crossrun_hrf(bold, events, tr, lib, noise_model = "ar1")

    dec <- deconvolve_crossrun_hrf(bold, fit_h, lambda = 0.7, scaling = "noise_sd", penalty = "sobolev")
    reg1 <- decon_epoch_regression(
      dec, covs, events = events, target_condition = "cue",
      lags = c(0, 1), restore_scale = TRUE, whiten_panel = TRUE, return_panel = sc <= 4L
    )
    cf0 <- reg1$coefficients[abs(reg1$coefficients$lag) < 1e-6, ]
    act_idx <- seq_len(n_rep_per_sched)
    nul_idx <- (n_rep_per_sched + 1L):(2L * n_rep_per_sched)

    est_act <- c(est_act, cf0$beta_gls[act_idx])
    se_act <- c(se_act, cf0$se_gls[act_idx])
    df_act <- c(df_act, cf0$df[act_idx])

    est_null <- c(est_null, cf0$beta_gls[nul_idx])
    se_null <- c(se_null, cf0$se_gls[nul_idx])
    df_null <- c(df_null, cf0$df[nul_idx])

    # Check 2-stage workflow (covariates = NULL in Stage 3, merged onto panel for Stage 4)
    if (sc <= 4L) {
      reg0 <- decon_epoch_regression(
        dec, covariates = NULL, events = events, target_condition = "cue",
        lags = c(0, 1), restore_scale = TRUE, whiten_panel = TRUE, return_panel = TRUE
      )
      x_map <- rbind(
        data.frame(run = "run1", trial = seq_len(n_trials), x = covs$run1$x),
        data.frame(run = "run2", trial = seq_len(n_trials), x = covs$run2$x)
      )
      p0 <- reg0$panel
      m_ord <- match(paste(p0$run, p0$trial), paste(x_map$run, x_map$trial))
      p0$x <- x_map$x[m_ord]
      reg0$panel <- p0
      ml2 <- medusa_multilag_infer(reg0, covariates = "x", n_perm = 200L, n_boot_latency = 0L, seed = 100L + sc)
      cf2 <- ml2$coefficients[abs(lag) < 1e-6 & voxel %in% paste0("v", act_idx)]
      est_2stage <- c(est_2stage, cf2$estimate)
    }
  }

  # 1. Unbiased recovery under rho_x = 0.70 (eliminates ~20% pre-fix attenuation)
  expect_lt(abs(mean(est_act) - beta_true), 0.045)
  expect_lt(abs(mean(est_2stage) - beta_true), 0.055)

  # 2. Nominal 95% confidence interval coverage on active slopes
  t_crit_a <- stats::qt(0.975, df = df_act)
  cov_95 <- mean(abs(est_act - beta_true) <= t_crit_a * se_act)
  expect_gte(cov_95, 0.89)
  expect_lte(cov_95, 0.99)

  # 3. Calibrated Type I error (~5%) on null slopes under autocorrelated covariate
  p_null <- 2.0 * stats::pt(abs(est_null / se_null), df = df_null, lower.tail = FALSE)
  type1_rate <- mean(p_null < 0.05)
  expect_gte(type1_rate, 0.015)
  expect_lte(type1_rate, 0.090)
})

test_that("P0-2: parcel pre-aggregation uses task-free innovation covariance, preserves negative spatial correlation, and mixes heterogeneous strata", {
  set.seed(20260402L)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 30L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)

  # 1. Task-free spatial variance: strong task curvature must NOT inflate parcel se_demixed
  groups <- stats::setNames(
    c(rep("P_strong", 5L), rep("P_null", 5L), rep("P_anticor", 2L)),
    paste0("v", 1:12)
  )
  bold <- events <- covs <- list()
  L_pos <- chol(0.45 * matrix(1, 5, 5) + 0.55 * diag(5))
  L_neg <- chol(matrix(c(1.0, -0.65, -0.65, 1.0), 2, 2))

  for (r in 1:2) {
    onsets <- cumsum(c(7.0, stats::runif(n_trials - 1L, 4.2, 6.5)))
    xv <- as.numeric(scale(stats::rnorm(n_trials)))
    G_can <- hrf_event_design_cpp(
      n_time, tr, 0, onsets, rep(0, n_trials), seq_len(n_trials), n_trials, lib[, 2L, drop = FALSE]
    )[[1L]]
    sig_strong <- as.numeric(G_can %*% (2.8 + 3.2 * xv))

    eps_strong <- matrix(stats::rnorm(n_time * 5L, sd = 0.25), n_time, 5L) %*% L_pos
    eps_null <- matrix(stats::rnorm(n_time * 5L, sd = 0.25), n_time, 5L) %*% L_pos
    eps_anti <- matrix(stats::rnorm(n_time * 2L, sd = 0.25), n_time, 2L) %*% L_neg

    Y_r <- cbind(
      sweep(eps_strong, 1L, 100 + sig_strong, "+"),
      100 + eps_null,
      100 + eps_anti
    )
    colnames(Y_r) <- names(groups)
    rn <- paste0("run", r)
    bold[[rn]] <- Y_r
    events[[rn]] <- data.frame(onset = onsets, duration = 0, condition = "stim")
    covs[[rn]] <- data.frame(x = xv)
  }

  fit_p <- estimate_crossrun_hrf(
    bold, events, tr, lib, spatial_groups = groups,
    pooling = "parcel", fallback_kernel = 2L, noise_model = "ar1",
    evidence = "calibrated", return_scores = TRUE
  )
  dec_p <- deconvolve_crossrun_hrf(bold, fit_p, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
  expect_true(is.list(dec_p$noise_innovations))

  reg_p <- decon_epoch_regression(
    dec_p, covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE
  )
  reg_v <- decon_epoch_regression(
    dec_p, covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = FALSE
  )

  se_strong <- reg_p$coefficients$mean_se_trial[reg_p$coefficients$parcel == "P_strong"]
  se_null <- reg_p$coefficients$mean_se_trial[reg_p$coefficients$parcel == "P_null"]
  expect_lt(abs(se_strong - se_null) / se_null, 0.16)

  # 2. Negative spatial correlation reduces parcel SE strictly below independent-voxel 1/sqrt(2) benchmark
  se_anti_parcel <- reg_p$coefficients$mean_se_trial[reg_p$coefficients$parcel == "P_anticor"]
  se_anti_vox <- reg_v$coefficients$mean_se_trial[reg_v$coefficients$voxel %in% c("v11", "v12")]
  se_indep_bench <- sqrt(sum((0.5 * se_anti_vox)^2))
  expect_lt(se_anti_parcel, 0.85 * se_indep_bench)

  # 3. restore_scale = FALSE precision_ev weights do not double-divide by s2_v
  reg_noscale <- decon_epoch_regression(
    dec_p, covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = FALSE, aggregate_by_parcel = TRUE
  )
  expect_true(all(vapply(reg_noscale$parcel_weights, function(w) abs(sum(w) - 1.0) < 1e-10, logical(1L))))

  # 4. Multi-stratum parcel operator (heterogeneous fast + slow HRF voxels aggregated into one parcel)
  bold_het <- list()
  het_groups <- stats::setNames(rep("P_het", 6L), paste0("h", 1:6))
  for (r in 1:2) {
    onsets <- events[[r]]$onset
    xv <- covs[[r]]$x
    G_fast <- hrf_event_design_cpp(
      n_time, tr, 0, onsets, rep(0, n_trials), seq_len(n_trials), n_trials, lib[, 1L, drop = FALSE]
    )[[1L]]
    G_slow <- hrf_event_design_cpp(
      n_time, tr, 0, onsets, rep(0, n_trials), seq_len(n_trials), n_trials, lib[, 3L, drop = FALSE]
    )[[1L]]
    amps <- 0.8 + 1.0 * xv
    Yh <- cbind(
      replicate(3L, 100 + as.numeric(G_fast %*% amps) + stats::rnorm(n_time, sd = 0.12)),
      replicate(3L, 100 + as.numeric(G_slow %*% amps) + stats::rnorm(n_time, sd = 0.12))
    )
    colnames(Yh) <- names(het_groups)
    bold_het[[paste0("run", r)]] <- Yh
  }
  fit_v <- estimate_crossrun_hrf(
    bold_het, events, tr, lib, pooling = "voxel", noise_model = "ar1", evidence = "calibrated"
  )
  dec_v <- deconvolve_crossrun_hrf(bold_het, fit_v, lambda = 0.5, scaling = "noise_sd", penalty = "sobolev")
  reg_het <- decon_epoch_regression(
    dec_v, covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = het_groups
  )
  expect_lt(abs(reg_het$coefficients$beta_gls[1L] - 1.0), 0.10)
})

test_that("P0-3: medusa_multilag_infer partials correlated covariates, estimates Level-2 traits, retains partial-lag subjects, and calibrates small-N omnibus & bands", {
  set.seed(20260403L)
  lags <- seq(-1, 6, by = 1.0)
  L <- length(lags)
  ops_ref <- medusa_reference_lag_operators(lags = lags, tr = 1.0)
  L_err <- chol(ops_ref$R_op + 1e-6 * diag(L))

  # 1. Joint multiple regression on correlated within-run covariates (cor(x1, x2) = 0.65, beta2 == 0)
  #    + Level-2 between-subject covariate (trait_age, gamma_age > 0; trait_null, gamma == 0)
  #    + Available-case retention when half of subjects have boundary lag NA
  S_n <- 20L
  n_tr_run <- 24L
  w_x1_true <- 1.10 * exp(-0.5 * ((lags - 2.0) / 1.1)^2)
  w_age_true <- 0.55 * as.numeric(lags >= 1.0 & lags <= 4.0)

  trait_age_s <- as.numeric(scale(stats::rnorm(S_n)))
  trait_null_s <- as.numeric(scale(stats::rnorm(S_n)))

  panel_parts <- vector("list", S_n)
  for (si in seq_len(S_n)) {
    sid <- sprintf("sub_%02d", si)
    b_x1_s <- w_x1_true + stats::rnorm(1L, 0, 0.15)
    b_int_s <- 0.40 + trait_age_s[si] * w_age_true + stats::rnorm(L, 0, 0.18)
    run_parts <- vector("list", 2L)
    for (r in 1:2) {
      z1 <- stats::rnorm(n_tr_run)
      x1 <- as.numeric(scale(z1))
      x2 <- as.numeric(scale(0.65 * z1 + sqrt(1 - 0.65^2) * stats::rnorm(n_tr_run)))
      E_tr <- matrix(stats::rnorm(n_tr_run * L, sd = 0.45), n_tr_run, L) %*% L_err
      Y_tr <- matrix(b_int_s, n_tr_run, L, byrow = TRUE) + outer(x1, b_x1_s) + E_tr
      # Half of subjects are missing the final lag (+6s) due to boundary gating
      if (si <= (S_n / 2L)) {
        Y_tr[, L] <- NA_real_
      }
      dt_r <- data.table::CJ(trial = seq_len(n_tr_run), lag_idx = seq_len(L))
      dt_r[, `:=`(
        id = sid,
        run = paste0("run", r),
        parcel = "R1",
        evt_time = lags[lag_idx],
        decon_demixed = Y_tr[cbind(trial, lag_idx)],
        precision_weight = 1.0,
        precision_weight_cohort = 1.0,
        x1 = x1[trial],
        x2 = x2[trial],
        trait_age = trait_age_s[si],
        trait_null = trait_null_s[si]
      )]
      dt_r[, lag_idx := NULL]
      run_parts[[r]] <- dt_r
    }
    panel_parts[[si]] <- data.table::rbindlist(run_parts)
  }
  panel_multi <- data.table::rbindlist(panel_parts)

  ml <- medusa_multilag_infer(
    panel_multi,
    lag_resolution = diag(L),
    lag_error_cor = ops_ref$R_op,
    trajectory_mode = "pspline",
    n_perm = 1500L,
    n_boot_latency = 50L,
    seed = 4031L
  )
  cf <- ml$coefficients

  # All 4 covariates (within-run x1, x2 AND between-subject trait_age, trait_null) are estimated
  expect_true(all(c("x1", "x2", "trait_age", "trait_null") %in% unique(cf$term)))
  # All 20 subjects retained despite 10 subjects having boundary lag NA
  expect_equal(length(unique(ml$subjects[term == "x1"]$id)), S_n)

  # Active x1 and trait_age strongly detected; correlated null x2 and trait_null controlled
  expect_lt(min(cf[term == "x1"]$p_tfce), 0.005)
  expect_gt(min(cf[term == "x2"]$p_maxT), 0.05)
  expect_gt(cf[term == "x2"]$omnibus_p[1L], 0.05)
  expect_lt(max(abs(cf[term == "x2"]$estimate)), 0.16)

  expect_lt(min(cf[term == "trait_age"]$p_tfce), 0.02)
  expect_lt(abs(mean(cf[term == "trait_age" & lag >= 1 & lag <= 4]$estimate) - 0.55), 0.15)
  expect_gt(cf[term == "trait_null"]$omnibus_p[1L], 0.05)

  # 2. Monte Carlo calibration at N = 15, L = 8: pivotal omnibus_p Type I error (~5%),
  #    simultaneous multilag band coverage (~95%), and unbiased unsmoothed contrast CIs
  n_mc <- 160L
  S_mc <- 15L
  omni_null_p <- numeric(n_mc)
  maxT_null_p <- numeric(n_mc)
  band_cov <- logical(n_mc)
  con_est <- numeric(n_mc)
  con_cov <- logical(n_mc)

  w_sharp <- 0.90 * exp(-0.5 * ((lags - 2.0) / 0.75)^2)
  c_late_early <- ifelse(lags %in% c(1, 2, 3), 1 / 3, ifelse(lags %in% c(-1, 0), -1 / 2, 0))
  true_con <- sum(c_late_early * w_sharp)

  for (mc in seq_len(n_mc)) {
    E_s <- matrix(stats::rnorm(S_mc * L, sd = 0.40), S_mc, L) %*% L_err
    dt_null <- data.table::CJ(si = seq_len(S_mc), li = seq_len(L))
    dt_null[, `:=`(
      id = sprintf("s%02d", si),
      parcel = "P0",
      term = "x",
      lag = lags[li],
      estimate = E_s[cbind(si, li)],
      std.error = 0.35
    )]
    res_0 <- medusa_multilag_infer(
      dt_null, lag_resolution = diag(L), lag_error_cor = ops_ref$R_op,
      trajectory_mode = "pspline", n_perm = 600L, n_boot_latency = 0L
    )
    omni_null_p[mc] <- res_0$coefficients$omnibus_p[1L]
    maxT_null_p[mc] <- min(res_0$coefficients$p_maxT)

    dt_sig <- data.table::copy(dt_null)
    dt_sig[, estimate := estimate + w_sharp[li]]
    res_1 <- medusa_multilag_infer(
      dt_sig, lag_resolution = diag(L), lag_error_cor = ops_ref$R_op,
      trajectory_mode = "pspline", contrasts = list(le = c_late_early),
      contrast_source = "unsmoothed", n_perm = 600L, n_boot_latency = 0L
    )
    cf1 <- res_1$coefficients[order(lag)]
    band_cov[mc] <- all(cf1$conf.low.multilag.sim <= w_sharp & w_sharp <= cf1$conf.high.multilag.sim)
    cn1 <- res_1$contrasts[contrast == "le"]
    con_est[mc] <- cn1$estimate
    con_cov[mc] <- (cn1$conf.low <= true_con) && (true_con <= cn1$conf.high)
  }

  expect_gte(mean(omni_null_p < 0.05), 0.015)
  expect_lte(mean(omni_null_p < 0.05), 0.095)
  expect_gte(mean(maxT_null_p < 0.05), 0.015)
  expect_lte(mean(maxT_null_p < 0.05), 0.095)
  expect_gte(mean(band_cov), 0.88)
  expect_lt(abs(mean(con_est) - true_con), 0.035)
  expect_gte(mean(con_cov), 0.89)
})

test_that("P0-4: decon_connectivity_repair eliminates shared cross-parcel & cross-lag error bias (even post-rbindlist) and calibrates unequal-weight HKSJ", {
  set.seed(20260404L)

  # 1. End-to-end BOLD simulation with strong shared physiological/global noise (r_shared = 0.60)
  #    and ZERO true neural connectivity or PPI, combined via data.table::rbindlist (stripping attributes)
  tr <- 1.0
  n_time <- 160L
  n_trials <- 26L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ"), each = 3L), paste0("v", 1:6))

  panels <- vector("list", 18L)
  for (s in seq_len(18L)) {
    bold_s <- ev_s <- x_s <- list()
    for (r in 1:2) {
      t_ev <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.3, 6.8)))
      xv <- as.numeric(scale(stats::rnorm(n_trials)))
      # Strictly independent latent neural variations (true connectivity = 0, true PPI = 0)
      eta_s <- stats::rnorm(n_trials, 0, 0.55)
      eta_t <- stats::rnorm(n_trials, 0, 0.55)
      des <- hrf_event_design_cpp(
        n_time, tr, 0, t_ev, rep(0, n_trials), seq_len(n_trials), n_trials, hrf_k
      )[[1L]]
      # Shared cross-parcel observation noise e_global
      e_glob <- stats::rnorm(n_time, sd = 0.30)
      ymat <- cbind(
        replicate(3L, 100 + as.numeric(des %*% (0.7 + 0.3 * xv + eta_s)) + e_glob + stats::rnorm(n_time, sd = 0.20)),
        replicate(3L, 100 + as.numeric(des %*% (0.7 + 0.3 * xv + eta_t)) + e_glob + stats::rnorm(n_time, sd = 0.20))
      )
      colnames(ymat) <- names(groups)
      rn <- paste0("run", r)
      bold_s[[rn]] <- ymat
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0, condition = "fb")
      x_s[[rn]] <- data.frame(x = xv)
    }
    fh <- estimate_crossrun_hrf(bold_s, ev_s, tr, lib, spatial_groups = groups,
      pooling = "parcel", noise_model = "ar1", evidence = "calibrated")
    dc <- deconvolve_crossrun_hrf(bold_s, fh, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
    rg <- decon_epoch_regression(dc, x_s, events = ev_s, target_condition = "fb",
      lags = c(0, 1), restore_scale = TRUE, aggregate_by_parcel = TRUE, whiten_panel = TRUE, return_panel = TRUE)
    ps <- rg$panel
    ps$id <- sprintf("sub%02d", s)
    panels[[s]] <- ps
  }
  panel_rb <- data.table::rbindlist(panels)
  expect_null(attr(panel_rb, "unit_error_cor"))
  expect_true("shared_noise_loading" %in% names(panel_rb))

  fit_zero <- decon_connectivity_repair(
    panel_rb, "R_seed", "R_targ", test_ppi = TRUE, method = "mmoms", shared_error = "zero"
  )
  fit_auto <- decon_connectivity_repair(
    panel_rb, "R_seed", "R_targ", test_ppi = TRUE, method = "mmoms",
    shared_error = "auto", return_augmented_panel = TRUE
  )

  # Neglecting shared_error ("zero") produces severe positive connectivity bias; "auto" removes it
  expect_gt(fit_zero$group[term == "seed_c"]$estimate, 0.20)
  expect_lt(abs(fit_auto$group[term == "seed_c"]$estimate), 0.14)
  expect_gt(fit_auto$group[term == "seed_c"]$p.value, 0.05)
  expect_gt(fit_auto$group[term == "seed_c:x_c"]$p.value, 0.05)
  expect_gt(mean(fit_auto$subjects$r_shared_error), 0.35)

  # Moment-preserving REPAIR predictors reproduce subject M-MOMS slopes in WLS even when err_xy != 0
  p_s1 <- fit_auto$panel_repair[id == "sub01"]
  lm_s1 <- stats::lm(targ_c ~ seed_c_repair + x_c_repair + `seed_c:x_c_repair`,
    data = p_s1, weights = repair_weight)
  expect_equal(unname(stats::coef(lm_s1)[-1L]), fit_auto$subjects[id == "sub01"]$estimate, tolerance = 1e-5)

  # Granger lagged connectivity (seed_lag = 0 -> target_lag = 1 controlling target at lag 0)
  fit_granger <- decon_connectivity_repair(
    panel_rb, "R_seed", "R_targ", seed_lag = 0, target_lag = 1,
    control_target_seed_lag = TRUE, method = "mmoms", shared_error = "auto"
  )
  expect_gt(fit_granger$group[term == "seed_c"]$p.value, 0.05)

  # 2. Monte Carlo calibration under shared error + unequal subject trial counts (HKSJ leverage correction)
  n_mc <- 100L
  p_mc_conn <- numeric(n_mc)
  p_mc_ppi <- numeric(n_mc)
  n_tr_vec <- c(90L, 18L, 18L, 18L, 18L, 18L, 18L, 18L)
  S_mc <- length(n_tr_vec)

  for (mc in seq_len(n_mc)) {
    dt_parts <- vector("list", S_mc)
    for (si in seq_len(S_mc)) {
      nt <- n_tr_vec[si]
      xv <- as.numeric(scale(stats::rnorm(nt)))
      eta_s <- stats::rnorm(nt, 0, 0.55)
      eta_t <- stats::rnorm(nt, 0, 0.55)
      e_sh <- stats::rnorm(nt, 0, 0.35)
      e_s <- e_sh + stats::rnorm(nt, 0, 0.25)
      e_t <- e_sh + stats::rnorm(nt, 0, 0.25)
      se_tot <- sqrt(0.35^2 + 0.25^2)
      r_sh <- 0.35^2 / (0.35^2 + 0.25^2)
      sid <- sprintf("s%02d", si)
      dt_parts[[si]] <- rbind(
        data.table::data.table(
          id = sid, run = "r1", trial = seq_len(nt), evt_time = 0, parcel = "S",
          decon_demixed = eta_s + e_s, se_demixed = se_tot, precision_weight = 1.0,
          shared_noise_loading = sqrt(r_sh), x = xv
        ),
        data.table::data.table(
          id = sid, run = "r1", trial = seq_len(nt), evt_time = 0, parcel = "T",
          decon_demixed = eta_t + e_t, se_demixed = se_tot, precision_weight = 1.0,
          shared_noise_loading = sqrt(r_sh), x = xv
        )
      )
    }
    res_mc <- decon_connectivity_repair(
      data.table::rbindlist(dt_parts), "S", "T", test_ppi = TRUE, method = "mmoms", shared_error = "auto"
    )
    p_mc_conn[mc] <- res_mc$group[term == "seed_c"]$p.value
    p_mc_ppi[mc] <- res_mc$group[term == "seed_c:x_c"]$p.value
  }
  expect_lte(mean(p_mc_conn < 0.05), 0.095)
  expect_lte(mean(p_mc_ppi < 0.05), 0.095)
})

test_that("P0-5: hrf_tune_decon_lambda pools all runs invariantly to run order, aligns mixed-polarity parcels, and forwards nuisance_mode", {
  set.seed(20260405L)
  tr <- 1.0
  n_time <- 140L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep("P_mix", 6L), paste0("v", 1:6))

  # Run 1 has pure noise (zero task SNR), while Runs 2 and 3 have strong opposite-polarity parcel signals (+1 on v1:v3, -1 on v4:v6)
  ev_base <- data.frame(onset = seq(10, 118, by = 9.0), duration = 0, condition = "go")
  G_go <- hrf_event_design_cpp(
    n_time, tr, 0, ev_base$onset, rep(0, nrow(ev_base)), seq_len(nrow(ev_base)), nrow(ev_base), hrf_k
  )[[1L]]

  Y_r1 <- matrix(100 + stats::rnorm(n_time * 6L, sd = 0.25), n_time, 6L, dimnames = list(NULL, names(groups)))
  sig2 <- as.numeric(G_go %*% stats::rnorm(nrow(ev_base), 1.5, 0.5))
  sig3 <- as.numeric(G_go %*% stats::rnorm(nrow(ev_base), 1.5, 0.5))
  pol <- c(1, 1, 1, -1, -1, -1)
  Y_r2 <- sweep(matrix(stats::rnorm(n_time * 6L, sd = 0.20), n_time, 6L), 2L, pol, FUN = function(e, s) e) +
    tcrossprod(sig2, pol) + 100
  Y_r3 <- tcrossprod(sig3, pol) + matrix(stats::rnorm(n_time * 6L, sd = 0.20), n_time, 6L) + 100
  colnames(Y_r2) <- colnames(Y_r3) <- names(groups)

  bold_123 <- list(run1 = Y_r1, run2 = Y_r2, run3 = Y_r3)
  events_123 <- list(run1 = ev_base, run2 = ev_base, run3 = ev_base)

  bold_321 <- list(run3 = Y_r3, run2 = Y_r2, run1 = Y_r1)
  events_321 <- list(run3 = ev_base, run2 = ev_base, run1 = ev_base)

  fit_123 <- estimate_crossrun_hrf(bold_123, events_123, tr, lib, spatial_groups = groups, pooling = "parcel")
  fit_321 <- estimate_crossrun_hrf(bold_321, events_321, tr, lib, spatial_groups = groups, pooling = "parcel")

  dec_123 <- deconvolve_crossrun_hrf(bold_123, fit_123, lambda = NULL, scaling = "noise_sd", nuisance_mode = "joint")
  dec_321 <- deconvolve_crossrun_hrf(bold_321, fit_321, lambda = NULL, scaling = "noise_sd", nuisance_mode = "joint")

  # Strictly invariant to run permutation (Run 1 no longer short-circuits lambda tuning)
  expect_equal(dec_123$settings$lambda, dec_321$settings$lambda, tolerance = 1e-12)

  # Mixed-polarity parcel alignment matches all-positive parcel lambda (preventing cancellation to max lambda = 8)
  Y_r2_pos <- tcrossprod(sig2, rep(1, 6L)) + (Y_r2 - 100 - tcrossprod(sig2, pol)) + 100
  Y_r3_pos <- tcrossprod(sig3, rep(1, 6L)) + (Y_r3 - 100 - tcrossprod(sig3, pol)) + 100
  colnames(Y_r2_pos) <- colnames(Y_r3_pos) <- names(groups)
  bold_pos <- list(run1 = Y_r1, run2 = Y_r2_pos, run3 = Y_r3_pos)
  fit_pos <- estimate_crossrun_hrf(bold_pos, events_123, tr, lib, spatial_groups = groups, pooling = "parcel")
  dec_pos <- deconvolve_crossrun_hrf(bold_pos, fit_pos, lambda = NULL, scaling = "noise_sd", nuisance_mode = "joint")
  expect_lt(dec_123$settings$lambda, 6.0)
  expect_equal(dec_123$settings$lambda, dec_pos$settings$lambda, tolerance = 1e-6)
})

test_that("P0-6: causal BOLD observability gating masks unobservable boundary/censored trial lags without dropping short-RT interior trials", {
  set.seed(20260406L)
  tr <- 1.0
  n_time <- 120L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]

  # Include short-RT interior trials (duration = 0.30s) plus a final trial at t = 117.5s
  onsets <- c(seq(10, 100, by = 9.0), 117.5)
  n_trials <- length(onsets)
  durations <- rep(0.30, n_trials)
  xv <- as.numeric(scale(seq_len(n_trials)))

  G_0 <- hrf_event_design_cpp(
    n_time, tr, 0, onsets, durations, seq_len(n_trials), n_trials, hrf_k
  )[[1L]]
  Y <- cbind(
    v1 = 100 + as.numeric(G_0 %*% (1.0 + 1.5 * xv)) + stats::rnorm(n_time, sd = 0.15),
    v2 = 100 + stats::rnorm(n_time, sd = 0.15)
  )
  bold <- list(run1 = Y)
  events <- list(run1 = data.frame(onset = onsets, duration = durations, condition = "resp"))
  # Censor the causal HRF window (t = 46..57s, rows 47..58) of trial 5 (onset = 46s)
  valid_rows <- list(run1 = !(seq_len(n_time) %in% 47:58))

  fit <- estimate_crossrun_hrf(bold, events, tr, lib, valid_rows = valid_rows, noise_model = "ar1")
  dec <- deconvolve_crossrun_hrf(bold, fit, lambda = 0.6, scaling = "noise_sd")
  reg <- decon_epoch_regression(
    dec, list(run1 = data.frame(x = xv)), events = events, target_condition = "resp",
    lags = c(-12, 0, 8), source_model = "duration", duration_scale = "rate",
    restore_scale = TRUE, return_panel = TRUE
  )
  p <- reg$panel[voxel == "v1"]

  # 1. Final trial at t = 117.5s at lag +8s (125.5s > 120s) has 0 observable BOLD energy -> gated to NA
  expect_true(is.na(p[trial == n_trials & evt_time == 8]$decon_demixed))
  # 2. First trial at t = 10s at lag -12s (-2s) -> pre-run, check gated if outside causal window or trial 5 censored
  expect_true(is.na(p[trial == 5L & evt_time == 0]$decon_demixed))
  # 3. Uncensored short-RT (0.30s) interior trials at lag 0 are 100% retained
  interior_ok <- setdiff(2:(n_trials - 2L), 5L)
  expect_true(all(is.finite(p[trial %in% interior_ok & evt_time == 0]$decon_demixed)))
})

test_that("P0-7 & P1-1: gap-aware estimate_psi_weighted eliminates censoring attenuation, propagates psi3..psi4, and handles 0-event runs", {
  set.seed(20260407L)

  # 1. Gap-aware estimate_psi_weighted inside hrf_estimate_noise under 25% random censoring at fast TR (0.72s)
  n_obs <- 420L
  n_vox <- 120L
  tr_fast <- 0.72
  ar_phi <- c(0.35, -0.08, 0.04, 0.28)
  Z <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n_obs))
  U <- matrix(0.0, n_obs, n_vox, dimnames = list(NULL, paste0("v", seq_len(n_vox))))
  for (v in seq_len(n_vox)) {
    U[, v] <- 100.0 + as.numeric(stats::arima.sim(model = list(ar = ar_phi), n = n_obs, sd = 0.25))
  }

  lib_fast <- hrf_library(tr_fast, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  ev_0 <- data.frame(onset = numeric(0), duration = numeric(0), condition = character(0))
  rho_grid <- seq(-0.20, 0.85, by = 0.05)

  mask_full <- rep(TRUE, n_obs)
  censor_idx <- sort(sample(2:(n_obs - 1L), size = round(0.25 * n_obs), replace = FALSE))
  mask_cens <- !(seq_len(n_obs) %in% censor_idx)

  nf_full <- hrf_estimate_noise(
    U, ev_0, tr_fast, 0.0, lib_fast, 1L, Z, mask_full, NULL, rho_grid, 256L
  )
  nf_cens <- hrf_estimate_noise(
    U, ev_0, tr_fast, 0.0, lib_fast, 1L, Z, mask_cens, NULL, rho_grid, 256L
  )

  # Under 25% random censoring (~76% of 5-TR lag-4 windows broken!), gap-local normalization prevents attenuation
  expect_gt(nf_full$psi4[1L], 0.12)
  expect_lt(abs(nf_cens$psi4[1L] - nf_full$psi4[1L]), 0.06)
  expect_lt(abs(nf_cens$rho[1L] - nf_full$rho[1L]), 0.06)

  # 2. psi3 and psi4 propagate into estimate_decon_noise_sd and scale_decon_input
  sd_ar2_only <- estimate_decon_noise_sd(
    U[, 1:20], noise_psi1 = rep(0.35, 20), noise_psi2 = rep(-0.08, 20)
  )
  sd_ar4_full <- estimate_decon_noise_sd(
    U[, 1:20], noise_psi1 = rep(0.35, 20), noise_psi2 = rep(-0.08, 20),
    noise_psi3 = rep(0.04, 20), noise_psi4 = rep(0.28, 20)
  )
  expect_true(all(is.finite(sd_ar4_full)))
  expect_gt(mean(abs(sd_ar4_full - sd_ar2_only)), 1e-3)

  # 3. P1-1: 0-event run does not crash estimate_crossrun_hrf or deconvolve_crossrun_hrf
  lib <- hrf_library(1.0, peak_delays = c(5.0, 6.5), dispersions = 1.0, undershoot_ratios = 6)
  ev1 <- data.frame(onset = c(10, 28, 46, 64), duration = 0, condition = "cue")
  ev_empty <- data.frame(onset = numeric(0), duration = numeric(0), condition = character(0))
  G1 <- hrf_event_design_cpp(90L, 1.0, 0, ev1$onset, rep(0, 4L), 1:4, 4L, lib[, 1L, drop = FALSE])[[1L]]
  bold_0ev <- list(
    run1 = cbind(v1 = 100 + as.numeric(G1 %*% c(1, 1.2, 0.9, 1.1)) + stats::rnorm(90, sd = 0.1),
                 v2 = 100 + stats::rnorm(90, sd = 0.1)),
    run2 = cbind(v1 = 100 + stats::rnorm(90, sd = 0.1),
                 v2 = 100 + stats::rnorm(90, sd = 0.1))
  )
  fit_0ev <- estimate_crossrun_hrf(bold_0ev, list(run1 = ev1, run2 = ev_empty), tr = 1.0, lib)
  dec_0ev <- deconvolve_crossrun_hrf(bold_0ev, fit_0ev, lambda = NULL, scaling = "noise_sd")
  expect_true(is.list(dec_0ev$activity))
  expect_equal(length(dec_0ev$activity), 2L)
})

test_that("P0-8: compute_zca_whitening_matrix is smooth across spectrum, leaves centered space invariant to mean mode, and exports cohort precision weights", {
  set.seed(20260408L)

  # 1. Smoothness across eigenvalue spectrum (no cliff at 0.35 * mean_diag) & run-mean null-mode invariance
  m <- 20L
  Q_m <- diag(m) - matrix(1 / m, m, m)
  V <- qr.Q(qr(Q_m %*% matrix(stats::rnorm(m * (m - 1L)), m, m - 1L)))[, seq_len(m - 1L)]
  base_eigs <- seq(0.50, 2.0, length.out = m - 2L)

  eval_at <- function(lam_small) {
    S_cent <- V %*% (c(lam_small, base_eigs) * t(V))
    Phi <- compute_zca_whitening_matrix(S_cent, eig_floor = 0.08, offdiag_scale = 0.90)
    Phi
  }
  Phi_below <- eval_at(0.349)
  Phi_above <- eval_at(0.351)
  expect_lt(max(abs(Phi_above - Phi_below)), 0.01)

  # Adding an arbitrary multiple of 1*1'/m does not change centered whitening Q_m %*% Phi %*% Q_m when < mean_diag
  S0 <- V %*% (c(0.40, base_eigs) * t(V))
  S_plus <- S0 + 0.25 * matrix(1 / m, m, m)
  expect_equal(
    Q_m %*% compute_zca_whitening_matrix(S0) %*% Q_m,
    Q_m %*% compute_zca_whitening_matrix(S_plus) %*% Q_m,
    tolerance = 1e-12
  )

  # 2. Cohort precision weights reflect absolute noise precision across low-noise vs high-noise units
  tr <- 1.0
  n_time <- 130L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  onsets <- seq(10, 110, by = 6.5)
  nt <- length(onsets)
  xv <- as.numeric(scale(stats::rnorm(nt)))
  G <- hrf_event_design_cpp(n_time, tr, 0, onsets, rep(0, nt), seq_len(nt), nt, lib[, 1L, drop = FALSE])[[1L]]
  Y <- cbind(
    v_low_noise = 100 + as.numeric(G %*% (0.6 + 0.5 * xv)) + stats::rnorm(n_time, sd = 0.12),
    v_high_noise = 100 + as.numeric(G %*% (0.6 + 0.5 * xv)) + stats::rnorm(n_time, sd = 0.65)
  )
  ev <- list(run1 = data.frame(onset = onsets, duration = 0, condition = "c"))
  fit <- estimate_crossrun_hrf(list(run1 = Y), ev, tr, lib, noise_model = "ar1")
  dec <- deconvolve_crossrun_hrf(list(run1 = Y), fit, lambda = 0.6, scaling = "noise_sd")
  reg <- decon_epoch_regression(
    dec, list(run1 = data.frame(x = xv)), events = ev, target_condition = "c",
    lags = 0, restore_scale = TRUE, whiten_panel = TRUE, return_panel = TRUE
  )
  p <- reg$panel
  expect_equal(mean(p[voxel == "v_low_noise"]$precision_weight), 1.0, tolerance = 1e-10)
  expect_equal(mean(p[voxel == "v_high_noise"]$precision_weight), 1.0, tolerance = 1e-10)
  expect_gt(
    mean(p[voxel == "v_low_noise"]$precision_weight_cohort),
    3.0 * mean(p[voxel == "v_high_noise"]$precision_weight_cohort)
  )
})
