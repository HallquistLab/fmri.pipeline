library(testthat)
library(dplyr)

data.table::setDTthreads(1L)

# ==============================================================================
# 1. Unified Contrast Engine & YAML Parity (L1 / L2 / L3)
# ==============================================================================

test_that("specify_contrasts and get_contrasts_from_spec support ~0 + A:B, interaction_diffs, conditioned '|', slope_diffs, and custom contrasts", {
  set.seed(101)
  df <- expand.grid(
    A = c("a1", "a2"),
    B = c("b1", "b2"),
    rep = 1:4,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  df$id <- paste0("s", seq_len(nrow(df)))
  df$session <- 1L
  df$age <- rnorm(nrow(df), mean = 25, sd = 4)

  # 1. Cell-means parameterization (~ 0 + A:B) preserves cat_vars and generates interaction_diffs & conditioned contrasts
  m_cell <- fmri.pipeline:::create_new_hi_model(
    data = df,
    level = 3L,
    cur_model_names = NULL,
    spec_list = list(
      name = "cell_ab",
      model_formula = "~ 0 + A:B",
      diagonal = FALSE,
      cell_means = TRUE,
      cond_means = "A",
      pairwise_diffs = c("A", "A | B"),
      interaction_diffs = "A:B",
      custom = list(corner_diff = c("Aa1:Bb1" = 1, "Aa2:Bb2" = -1))
    ),
    lg = lgr::get_logger("test")
  )

  expect_equal(sort(m_cell$contrast_spec$cat_vars), c("A", "B"))
  cmat <- m_cell$contrasts
  expect_true(is.matrix(cmat) && nrow(cmat) > 0L)
  # Interaction diff on cell means: (a1b1 - a2b1) - (a1b2 - a2b2) = [1, -1, -1, 1]
  ix_rows <- grep("^A_x_B\\.", rownames(cmat), value = TRUE)
  expect_length(ix_rows, 1L)
  expect_equal(
    as.numeric(cmat[ix_rows, c("Aa1:Bb1", "Aa2:Bb1", "Aa1:Bb2", "Aa2:Bb2")]),
    c(1, -1, -1, 1),
    tolerance = 1e-7
  )
  # Conditioned simple effects "A | B" present
  expect_true(any(grepl("^B\\.b1_A\\.", rownames(cmat))))
  expect_true(any(grepl("^B\\.b2_A\\.", rownames(cmat))))
  # Custom contrast preserved
  expect_true("corner_diff" %in% rownames(cmat))
  expect_equal(as.numeric(cmat["corner_diff", c("Aa1:Bb1", "Aa2:Bb2")]), c(1, -1))

  # 2. Continuous x categorical interaction supports slope_diffs and survives respecify_l3_model
  m_slope <- fmri.pipeline:::create_new_hi_model(
    data = df,
    level = 3L,
    cur_model_names = NULL,
    spec_list = list(
      name = "slope_mod",
      model_formula = "~ A * age",
      diagonal = FALSE,
      simple_slopes = list(age = "A"),
      slope_diffs = list(age = "A"),
      custom = list(age_main = c("age" = 1))
    ),
    lg = lgr::get_logger("test")
  )
  expect_true(any(grepl("^age\\.A\\.", rownames(m_slope$contrasts))))
  expect_true("age_main" %in% rownames(m_slope$contrasts))

  # Survives L3 refit on subset of subjects
  refit_slope <- fmri.pipeline:::respecify_l3_model(m_slope, df[1:12, c("id", "session"), drop = FALSE])
  expect_true("age_main" %in% rownames(refit_slope$contrasts))
  expect_true(any(grepl("^age\\.A\\.", rownames(refit_slope$contrasts))))
})

test_that("factor_coding = 'sum' survives relevel() and L2/L3 refits, and special-character levels do not collide", {
  df <- expand.grid(
    run_number = 1:6,
    session = 1L,
    id = c("s1", "s2"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  df$cue <- rep(c("CS+", "CS-", "CS0"), length.out = nrow(df))
  df$valence <- rep(c("high", "low"), length.out = nrow(df))

  m_l2 <- fmri.pipeline:::create_new_hi_model(
    data = df,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_sum",
      model_formula = "~ cue * valence",
      diagonal = TRUE,
      cond_means = "cue",
      pairwise_diffs = "cue",
      factor_coding = "sum",
      reference_level = c(valence = "low")
    ),
    lg = lgr::get_logger("test")
  )

  # Sum coding has entries in {-1, 0, 1} rather than {0, 1}
  expect_true(any(m_l2$model_matrix < -0.5))
  # CS+, CS-, CS0 produce 3 distinct condition means + 3 distinct pairwise differences = 6 non-zero rows!
  cue_rows <- grep("^cue\\.", rownames(m_l2$contrasts), value = TRUE)
  expect_length(cue_rows, 6L)
  expect_true(all(rowSums(abs(m_l2$contrasts[cue_rows, , drop = FALSE])) > 0))

  # Survives L2 subject refit with sum coding and all contrasts intact
  refit_l2 <- fmri.pipeline:::respecify_l2_models_by_subject(m_l2, df, split_on = c("id", "session"))
  s1_mm <- refit_l2$by_subject$model_matrix[[1L]]
  expect_true(any(s1_mm < -0.5))
  expect_equal(nrow(refit_l2$by_subject$contrasts[[1L]]), nrow(m_l2$contrasts))
})

test_that("YAML round-trip preserves custom contrasts, interaction_diffs, slope_diffs, wi_contrasts simple_slopes, factor_coding, and empty [] vectors", {
  skip_if_not_installed("yaml")

  gpa <- create_mock_gpa(n_subjects = 4, n_runs = 2, n_trials = 8, output_directory = tempdir())
  gpa$trial_data$cue <- rep(c("A", "B"), length.out = nrow(gpa$trial_data))
  gpa$trial_data$rew <- rep(c("hi", "lo"), each = 2, length.out = nrow(gpa$trial_data))
  gpa$trial_data$rt_mod <- seq_len(nrow(gpa$trial_data)) / 10
  gpa$subject_data$group <- rep(c("pat", "ctl"), length.out = nrow(gpa$subject_data))
  gpa$subject_data$age <- c(21, 25, 31, 36)

  spec <- list(
    onsets = list("onset"),
    durations = list("duration"),
    events = list(stim = list(onset = "onset", duration = "duration")),
    signals = list(
      stim_wi = list(
        event = "stim",
        normalization = "none",
        value_fixed = 1,
        wi_factors = list("cue", "rew"),
        wi_formula = "~ cue * rew",
        wi_expansion = "auto"
      ),
      stim_slope = list(
        event = "stim",
        normalization = "none",
        value_fixed = 1,
        wi_factors = list("cue", "rt_mod"),
        wi_formula = "~ cue * rt_mod",
        wi_expansion = "formula"
      )
    ),
    l1_models = list(
      m1 = list(
        signals = list("stim_wi", "stim_slope"),
        contrasts = list(
          diagonal = FALSE,
          wi_contrasts = list(
            stim_wi = list(
              diagonal = FALSE,
              cell_means = TRUE,
              cond_means = list(), # YAML [] -> list() should coerce cleanly without crashing
              pairwise_diffs = list("cue"),
              interaction_diffs = list("cue:rew")
            ),
            stim_slope = list(
              diagonal = FALSE,
              simple_slopes = list(rt_mod = list("cue")),
              slope_diffs = list(rt_mod = list("cue"))
            )
          )
        )
      )
    ),
    l3_models = list(
      g1 = list(
        model_formula = "~ group * age",
        factor_coding = "sum",
        contrasts = list(
          diagonal = FALSE,
          cond_means = list(),
          pairwise_diffs = list("group"),
          simple_slopes = list(age = list("group")),
          slope_diffs = list(age = list("group")),
          custom = list(age_overall = list(age = 1))
        )
      )
    )
  )

  yaml_in <- tempfile(fileext = ".yaml")
  writeLines(yaml::as.yaml(spec), yaml_in)

  gpa <- build_l1_models(gpa, from_spec_file = yaml_in)
  gpa <- build_l3_models(gpa, from_spec_file = yaml_in)

  out_yaml <- tempfile(fileext = ".yaml")
  export_glm_config(gpa, file = out_yaml)
  parsed <- yaml::read_yaml(out_yaml)

  expect_equal(parsed$l1_models$m1$contrasts$wi_contrasts$stim_wi$interaction_diffs, "cue:rew")
  expect_equal(parsed$l1_models$m1$contrasts$wi_contrasts$stim_slope$simple_slopes$rt_mod, "cue")
  expect_equal(parsed$l1_models$m1$contrasts$wi_contrasts$stim_slope$slope_diffs$rt_mod, "cue")
  expect_equal(parsed$l3_models$g1$factor_coding, "sum")
  expect_true("age_overall" %in% names(parsed$l3_models$g1$contrasts$custom))

  # Re-import from exported YAML and verify identical contrast matrices
  gpa2 <- create_mock_gpa(n_subjects = 4, n_runs = 2, n_trials = 8, output_directory = tempdir())
  gpa2$trial_data <- gpa$trial_data
  gpa2$subject_data <- gpa$subject_data
  gpa2 <- build_l1_models(gpa2, from_spec_file = out_yaml)
  gpa2 <- build_l3_models(gpa2, from_spec_file = out_yaml)

  expect_equal(gpa2$l1_models$models$m1$contrasts, gpa$l1_models$models$m1$contrasts, tolerance = 1e-7)
  expect_equal(gpa2$l3_models$models$g1$contrasts, gpa$l3_models$models$g1$contrasts, tolerance = 1e-7)
})

# ==============================================================================
# 2. Level 1 Factorial Specification (wi_factors, derivatives, empty EVs)
# ==============================================================================

test_that("wi_expansion = 'auto' uses cell means for saturated ~A*B (matching 'formula' COPEs) and preserves unsaturated ~A+B", {
  set.seed(202)
  td <- expand.grid(
    cue = c("c1", "c2"),
    rew = c("r1", "r2"),
    rep = 1:3,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  td$id <- "sub1"
  td$session <- 1L
  td$run_number <- 1L
  td$trial <- seq_len(nrow(td))
  td$onset <- seq(2, by = 4, length.out = nrow(td))
  td$duration <- 1.0
  td$event <- "stim"
  td$value <- 1.0

  wi_cspec <- list(
    diagonal = FALSE,
    cell_means = TRUE,
    cond_means = "cue",
    pairwise_diffs = c("cue", "rew", "cue | rew"),
    interaction_diffs = "cue:rew"
  )

  make_sig_and_contrasts <- function(exp_mode) {
    sig <- list(
      name = "stim",
      event = "stim",
      value_type = "unit",
      value_fixed = 1,
      value = td[, c("id", "session", "run_number", "trial", "value", "cue", "rew")],
      normalization = "none",
      beta_series = FALSE,
      wi_factors = c("cue", "rew"),
      wi_formula = "~ cue * rew",
      wi_expansion = exp_mode
    )
    sig <- fmri.pipeline:::fit_wi_model(sig)
    wi_obj <- structure(
      list(
        level = 1L, lmfit = sig$wi_model, signal_name = "stim",
        wi_factors = sig$wi_factors, wi_formula = sig$wi_formula,
        wi_expansion = sig$wi_expansion, wi_expansion_resolved = sig$wi_expansion_resolved
      ),
      class = c("list", "l1_wi_spec")
    )
    wi_con <- fmri.pipeline:::specify_contrasts(wi_obj, spec_list = wi_cspec)
    list(sig = sig, wi_con = wi_con)
  }

  res_auto <- make_sig_and_contrasts("auto")
  res_form <- make_sig_and_contrasts("formula")

  expect_equal(res_auto$sig$wi_expansion_resolved, "cell_means")
  expect_equal(res_form$sig$wi_expansion_resolved, "formula")
  expect_true(grepl("^~\\s*0\\s*\\+", res_auto$sig$wi_ev_formula))

  # Verify exact mathematical COPE equivalence on arbitrary BOLD response y = X_cell %*% true_cell_betas + noise
  exp_auto <- fmri.pipeline:::expand_signal(res_auto$sig)
  exp_form <- fmri.pipeline:::expand_signal(res_form$sig)

  rd <- data.frame(run_number = 1L, run_volumes = 60L, drop_volumes = 0L, tr = 1.0)
  ev_df <- td[, c("event", "run_number", "trial", "onset", "duration")]

  bdm_auto <- build_design_matrix(events = ev_df, signals = exp_auto, tr = 1.0, run_data = rd, write_timing_files = NULL, plot = FALSE)
  bdm_form <- build_design_matrix(events = ev_df, signals = exp_form, tr = 1.0, run_data = rd, write_timing_files = NULL, plot = FALSE)

  X_auto <- cbind(1, as.matrix(bdm_auto$design_convolved[[1L]]))
  X_form <- cbind(1, as.matrix(bdm_form$design_convolved[[1L]]))
  y <- rnorm(nrow(X_auto)) + X_auto[, -1, drop = FALSE] %*% c(2.0, -1.0, 0.5, 3.0)

  b_auto <- stats::coef(stats::lm.fit(X_auto, y))[-1L]
  b_form <- stats::coef(stats::lm.fit(X_form, y))[-1L]

  copes_auto <- as.numeric(res_auto$wi_con$contrasts %*% b_auto)
  copes_form <- as.numeric(res_form$wi_con$contrasts %*% b_form)
  names(copes_auto) <- rownames(res_auto$wi_con$contrasts)
  names(copes_form) <- rownames(res_form$wi_con$contrasts)

  expect_equal(copes_auto, copes_form[names(copes_auto)], tolerance = 1e-6)

  # Unsaturated additive 3x3 (~ A + B) stays in 'formula' mode (5 EVs instead of 9 cell EVs)
  td3 <- expand.grid(A = c("a1", "a2", "a3"), B = c("b1", "b2", "b3"), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  td3$id <- "sub1"; td3$session <- 1L; td3$run_number <- 1L; td3$trial <- seq_len(nrow(td3)); td3$value <- 1.0
  sig_add <- list(
    name = "stim", event = "stim", value_type = "unit", value = td3, normalization = "none", beta_series = FALSE,
    wi_factors = c("A", "B"), wi_formula = "~ A + B", wi_expansion = "auto"
  )
  sig_add <- fmri.pipeline:::fit_wi_model(sig_add)
  expect_equal(sig_add$wi_expansion_resolved, "formula")
  expect_equal(length(stats::coef(sig_add$wi_model)), 5L)
})

test_that("expand_signal handles NA rows in wi_factors without misaligning trial onsets and supports <2 observed subject levels", {
  td_all <- data.frame(
    id = c(rep("s1", 4), rep("s2", 3)),
    session = 1L,
    run_number = 1L,
    trial = c(1:4, 1:3),
    value = 1.0,
    cue = c("A", NA, "B", "A", "A", "A", "A"),
    stringsAsFactors = FALSE
  )

  sig <- list(
    name = "cue_sig",
    event = "stim",
    value_type = "unit",
    value = td_all,
    normalization = "none",
    beta_series = FALSE,
    wi_factors = "cue",
    wi_formula = "~ cue",
    wi_expansion = "auto"
  )
  sig <- fmri.pipeline:::fit_wi_model(sig)

  # s1 has NA on trial 2: trial 3 ("B") must map to cueB with trial = 3, NOT trial 2!
  sig_s1 <- sig
  sig_s1$value <- td_all[td_all$id == "s1", , drop = FALSE]
  expect_warning(exp_s1 <- fmri.pipeline:::expand_signal(sig_s1), "missing \\(NA\\)")
  s1_names <- vapply(exp_s1, `[[`, character(1), "name")
  b_val <- exp_s1[[grep("cueB", s1_names)]]$value
  expect_equal(b_val$trial, 3L)
  expect_true(!2L %in% exp_s1[[grep("cueA", s1_names)]]$value$trial)

  # s2 observes ONLY level "A" (< 2 levels): expand_signal must not crash, and cueB has 0 trials
  sig_s2 <- sig
  sig_s2$value <- td_all[td_all$id == "s2", , drop = FALSE]
  exp_s2 <- fmri.pipeline:::expand_signal(sig_s2)
  s2_names <- vapply(exp_s2, `[[`, character(1), "name")
  expect_equal(nrow(exp_s2[[grep("cueA", s2_names)]]$value), 3L)
  expect_equal(nrow(exp_s2[[grep("cueB", s2_names)]]$value), 0L)
})

test_that("fsl_l1_model emits shape=10 for empty EVs, zeros non-estimable contrasts, writes _cope_validity.tsv, and guards weights='cells'", {
  tmp_dir <- tempfile("l1_empty_cell_")
  dir.create(tmp_dir, recursive = TRUE)

  # 2 runs: Run 1 has all 4 cells of A x B; Run 2 is missing cell (a2, b2)
  r1 <- expand.grid(A = c("a1", "a2"), B = c("b1", "b2"), rep = 1:2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  r1$run_number <- 1L
  r2 <- subset(r1, !(A == "a2" & B == "b2"))
  r2$run_number <- 2L
  td <- rbind(r1, r2)
  td$id <- "sub1"
  td$session <- 1L
  td$trial <- ave(td$run_number, td$run_number, FUN = seq_along)
  td$onset <- td$trial * 4
  td$duration <- 1.0

  # Verify weights = "cells" warning and marginal masking when fitted directly on incomplete data
  m_incomp <- structure(list(
    level = 2L,
    contrast_spec = list(
      diagonal = FALSE,
      cell_means = TRUE,
      cond_means = "A",
      pairwise_diffs = c("A", "A | B"),
      weights = "cells"
    )
  ), class = c("hi_model_spec", "list"))
  m_incomp <- fmri.pipeline:::mobj_fit_lm(m_incomp, model_formula = "~ A * B", data = cbind(td[td$run_number == 2L, ], dummy = 1), id_cols = c("id", "session", "run_number"))
  expect_warning(
    m_incomp <- fmri.pipeline:::get_contrasts_from_spec(m_incomp, m_incomp$lmfit, allow_interactive_dedupe = FALSE),
    "zero observations"
  )
  # Unconditioned A.a2 and A.a1-a2 touch empty cell (a2, b2) and are dropped,
  # while A.a1 (both b1 and b2 observed), cell mean A.a2_B.b1, and conditioned B.b1_A.a1-a2 survive!
  expect_false(any(c("A.a2", "A.a1-a2") %in% rownames(m_incomp$contrasts)))
  expect_true("A.a1" %in% rownames(m_incomp$contrasts))
  expect_true("A.a2_B.b1" %in% rownames(m_incomp$contrasts))
  expect_true(any(grepl("^B\\.b1_A\\.", rownames(m_incomp$contrasts))))

  # Now test full L1 setup via build_design_matrix + fsl_l1_model in BOTH 'auto' and 'formula' modes
  for (exp_mode in c("auto", "formula")) {
    mode_dir <- file.path(tmp_dir, exp_mode)
    dir.create(mode_dir, recursive = TRUE)
    nifti_files <- file.path(mode_dir, paste0("bold_run", 1:2, ".nii.gz"))
    for (nf in nifti_files) {
      RNifti::writeNifti(array(0, dim = c(2, 2, 2, 60)), file = nf)
    }

    gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, output_directory = mode_dir)
    gpa$trial_data <- td
    gpa$subject_data$mr_dir <- mode_dir
    gpa$run_data <- data.frame(
      id = "sub1",
      session = 1L,
      run_number = 1:2,
      mr_dir = mode_dir,
      run_nifti = nifti_files,
      run_nifti_present = TRUE,
      run_volumes = 60L,
      drop_volumes = 0L,
      tr = 1.0,
      nvoxels = 1000L,
      exclude_run = FALSE,
      exclude_subject = FALSE,
      stringsAsFactors = FALSE
    )
    gpa$use_preconvolve <- TRUE
    gpa$glm_software <- "fsl"
    gpa$lgr_threshold <- "warn"
    gpa$log_json <- FALSE
    gpa$log_txt <- FALSE
    gpa$parallel$l1_setup_cores <- 1L
    gpa$output_locations$feat_l1_directory <- file.path(mode_dir, "feat_l1", "sub-{id}", "ses-{session}", "{l1_model}")
    gpa$output_locations$feat_l1_fsf <- "FEAT_LVL1_run{run_number}.fsf"
    gpa$glm_settings <- list(fsl = list(
      force_l1_creation = TRUE,
      smooth_fwhm = 5,
      hp_filter = 100,
      prewhiten = TRUE,
      z_thresh = 2.3,
      p_thresh = 0.05
    ))
    gpa$confound_settings <- list(l1_confound_regressors = NULL)

    spec <- list(
      onsets = list("onset"),
      durations = list("duration"),
      events = list(stim = list(onset = "onset", duration = "duration")),
      signals = list(
        task = list(
          event = "stim",
          normalization = "none",
          value_fixed = 1,
          wi_factors = list("A", "B"),
          wi_formula = "~ A * B",
          wi_expansion = exp_mode
        )
      ),
      l1_models = list(
        m1 = list(
          signals = list("task"),
          contrasts = list(
            diagonal = FALSE,
            wi_contrasts = list(
              task = list(
                diagonal = FALSE,
                cell_means = TRUE,
                pairwise_diffs = list("A"),
                interaction_diffs = list("A:B")
              )
            )
          )
        )
      )
    )
    spec_f <- tempfile(fileext = ".yaml")
    writeLines(yaml::as.yaml(spec), spec_f)
    gpa <- build_l1_models(gpa, from_spec_file = spec_f)

    sig_sub <- gpa$l1_models$signals$task
    exp_sigs <- fmri.pipeline:::expand_signal(sig_sub)
    bdm <- build_design_matrix(
      events = gpa$l1_models$events$stim$data,
      signals = exp_sigs,
      tr = 1.0,
      run_data = gpa$run_data,
      write_timing_files = "convolved",
      output_directory = file.path(mode_dir, "timing_files"),
      plot = FALSE
    )
    bdm$run_nifti <- nifti_files

    l1_df <- fsl_l1_model(
      id = "sub1",
      session = 1L,
      d_obj = bdm,
      gpa = gpa,
      model_name = "m1",
      run_nifti = nifti_files,
      nvoxels = c(1000L, 1000L)
    )

    expect_equal(nrow(l1_df), 2L)
    fsf_r2 <- readLines(l1_df$feat_fsf[2L])
    # Run 2 has an empty/aliased EV -> must emit shape 10 (FSL empty EV)
    expect_true(any(grepl("^set fmri\\(shape[0-9]+\\) 10$", fsf_r2)), info = exp_mode)

    # Sidecar _cope_validity.tsv written and marks A.a1_B.b1 estimable in run 2 while A.a2_B.b2 & interaction are invalid
    sidecar_r2 <- file.path(dirname(l1_df$feat_fsf[2L]), paste0(sub("\\.fsf$", "", basename(l1_df$feat_fsf[2L])), "_cope_validity.tsv"))
    expect_true(file.exists(sidecar_r2), info = exp_mode)
    sc_df <- utils::read.delim(sidecar_r2, stringsAsFactors = FALSE)
    expect_true(sc_df$estimable[grep("A\\.a1_B\\.b1$", sc_df$l1_cope_name)], info = exp_mode)
    expect_false(sc_df$estimable[grep("A\\.a2_B\\.b2$", sc_df$l1_cope_name)], info = exp_mode)
    expect_false(sc_df$estimable[grep("A_x_B\\.", sc_df$l1_cope_name)], info = exp_mode)

    # Create mock FEAT completion for run 1 & run 2 and verify build_l1_cope_validity honors sidecar
    for (rr in 1:2) {
      dir.create(file.path(l1_df$feat_dir[rr], "stats"), recursive = TRUE)
      file.create(file.path(l1_df$feat_dir[rr], ".feat_complete"))
      for (cc in seq_len(nrow(sc_df))) {
        file.create(file.path(l1_df$feat_dir[rr], "stats", paste0("cope", cc, ".nii.gz")))
        file.create(file.path(l1_df$feat_dir[rr], "stats", paste0("varcope", cc, ".nii.gz")))
      }
    }
    val_df <- fmri.pipeline:::build_l1_cope_validity(l1_df, gpa, lg = lgr::get_logger("test"))
    r2_val <- val_df[val_df$run_number == 2L, ]
    expect_true(r2_val$valid_for_higher_level[grep("A\\.a1_B\\.b1$", r2_val$l1_cope_name)], info = exp_mode)
    expect_false(r2_val$valid_for_higher_level[grep("A\\.a2_B\\.b2$", r2_val$l1_cope_name)], info = exp_mode)
  }
})

test_that("formula-mode parametric modulator pre-centers per run, zeros 1-event/zero-SD modulators, orthogonalizes add_deriv, and computes ':' VIFs", {
  # Run 1 has 4 events with varying RT; Run 2 has only 1 event (zero variance after centering)
  td <- data.frame(
    id = "sub1",
    session = 1L,
    event = "stim",
    run_number = c(1L, 1L, 1L, 1L, 2L),
    trial = c(1:4, 1L),
    onset = c(4, 12, 20, 28, 6),
    duration = 1.0,
    cond = c("A", "A", "B", "B", "A"),
    rew = c("r1", "r2", "r1", "r2", "r1"),
    rt = c(0.4, 0.8, 0.5, 0.9, 0.7),
    stringsAsFactors = FALSE
  )

  # Signal 1: categorical A:B cell-means with ':' in regressor names + add_deriv = TRUE
  sig_evt <- list(
    name = "task_evt",
    event = "stim",
    value_type = "unit",
    value = data.frame(id = td$id, session = td$session, run_number = td$run_number, trial = td$trial, value = 1.0, cond = td$cond, rew = td$rew, stringsAsFactors = FALSE),
    normalization = "none",
    beta_series = FALSE,
    add_deriv = TRUE,
    wi_factors = c("cond", "rew"),
    wi_formula = "~ cond * rew",
    wi_expansion = "auto"
  )
  sig_evt <- fmri.pipeline:::fit_wi_model(sig_evt)

  # Signal 2: formula-mode parametric modulator (~ cond with sum coding) to test per-run pre-centering & 1-event zeroing
  sig_pm <- list(
    name = "task_rt",
    event = "stim",
    value_type = "parametric",
    parametric_modulator = "rt",
    value = data.frame(id = td$id, session = td$session, run_number = td$run_number, trial = td$trial, value = td$rt, cond = td$cond, stringsAsFactors = FALSE),
    normalization = "none",
    beta_series = FALSE,
    wi_factors = "cond",
    wi_formula = "~ cond",
    wi_expansion = "formula",
    factor_coding = "sum"
  )
  sig_pm <- fmri.pipeline:::fit_wi_model(sig_pm)

  expanded <- c(fmri.pipeline:::expand_signal(sig_evt), fmri.pipeline:::expand_signal(sig_pm))

  rd <- data.frame(run_number = 1:2, run_volumes = c(45L, 45L), drop_volumes = c(0L, 0L), tr = 1.0)
  bdm <- build_design_matrix(
    events = td[, c("event", "run_number", "trial", "onset", "duration")],
    signals = expanded,
    tr = 1.0,
    run_data = rd,
    center_values = TRUE,
    write_timing_files = NULL,
    plot = FALSE
  )

  # 1. Derivative regressor in Run 1 is orthogonalized temporal derivative (NOT identical copy with cor == 1)
  r1_mat <- bdm$design_convolved[[1L]]
  base_col <- "task_evt.condA:rewr1"
  deriv_col <- "d_task_evt.condA:rewr1"
  expect_true(all(c(base_col, deriv_col) %in% names(r1_mat)))
  expect_gt(stats::sd(r1_mat[[deriv_col]]), 0)
  expect_lt(abs(stats::cor(r1_mat[[base_col]], r1_mat[[deriv_col]])), 1e-5)

  # 2. ':' interaction column names compute valid finite VIFs in Run 1
  expect_true(all(is.finite(bdm$collinearity$run_number1$vif)))

  # 3. Run 2 has only 1 event -> centered parametric modulators must be all zeros (marked empty)
  r2_mat <- bdm$design_convolved[[2L]]
  pm_cols <- grep("^task_rt\\.", names(r2_mat), value = TRUE)
  expect_true(length(pm_cols) >= 1L)
  for (pc in pm_cols) {
    expect_equal(max(abs(r2_mat[[pc]])), 0)
    expect_true(pc %in% bdm$runs_empty_evs[[2L]])
  }
})

# ==============================================================================
# 3. Level 2 Factorial Specification, Passthrough & Metadata Alignment
# ==============================================================================

test_that("respecify_l2_models_by_subject handles <2 surviving levels, missing reference level, num2fac, and zero-SD zscore", {
  # Balanced 2x3x2 parent dataset so all cond x block_num cells exist in parent
  rd <- expand.grid(
    cond = c("ref", "alt1", "alt2"),
    block_num = c(1, 2),
    id = c("s1", "s2", "s3"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  rd$session <- 1L
  rd$run_number <- ave(rd$session, rd$id, FUN = seq_along)
  rd$cov_z <- ifelse(rd$id == "s2", 5, rd$run_number * 2 + ifelse(rd$block_num == 2, 1, 0))

  m_l2 <- fmri.pipeline:::create_new_hi_model(
    data = rd,
    level = 2L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l2_fac",
      model_formula = "~ cond + block_num + cov_z",
      diagonal = TRUE,
      cond_means = "cond",
      pairwise_diffs = "cond",
      num2fac = "block_num",
      covariate_transform = c(cov_z = "zscore"),
      reference_level = c(cond = "ref")
    ),
    lg = lgr::get_logger("test")
  )

  # Drop runs so s3 has ONLY "alt1" and "alt2" (missing reference level "ref"!) and only block_num == 1 (<2 levels!)
  retained <- rd[!(rd$id == "s3" & (rd$cond == "ref" | rd$block_num == 2)), , drop = FALSE]
  # Give s3 2 replicates of alt1 and alt2 with block_num == 1 so it has 4 runs and df > 0
  s3_extra <- retained[retained$id == "s3", , drop = FALSE]
  s3_extra$run_number <- s3_extra$run_number + 10L
  s3_extra$cov_z <- s3_extra$cov_z + 3
  retained <- rbind(retained, s3_extra)

  expect_warning(
    refit <- fmri.pipeline:::respecify_l2_models_by_subject(m_l2, retained, split_on = c("id", "session")),
    "zero observations"
  )

  # s2: zero-SD cov_z is zeroed and pruned without NaN, while estimable contrasts survive
  s2_dmat <- refit$by_subject[id == "s2"]$model_matrix[[1L]]
  expect_false(anyNA(s2_dmat))
  expect_false("cov_z" %in% colnames(s2_dmat))
  expect_gt(nrow(refit$by_subject[id == "s2"]$contrasts[[1L]]), 0L)

  # s3: missing reference level 'ref' and single level of block_num does not crash;
  # EV_Intercept is non-estimable and pruned!
  s3_cmat <- refit$by_subject[id == "s3"]$contrasts[[1L]]
  expect_false("EV_Intercept" %in% rownames(s3_cmat))
  expect_equal(qr(refit$by_subject[id == "s3"]$model_matrix[[1L]])$rank, ncol(refit$by_subject[id == "s3"]$model_matrix[[1L]]))
})

test_that("single-run L2 passthrough allows pooling contrast, drops slope contrast, synthesizes masked tdof_t1, and rejects difference-only L2", {
  skip_if_not_installed("RNifti")

  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, include_l1_models = TRUE, include_l2_models = FALSE)
  gpa$run_data$exclude_run <- FALSE
  gpa$run_data$exclude_subject <- FALSE
  gpa$run_data$run_centered <- c(-0.5, 0.5)
  gpa$run_data$cond <- c("A", "B")
  gpa$lgr_threshold <- "warn"
  gpa$log_json <- FALSE
  gpa$log_txt <- FALSE
  gpa$parallel$l2_setup_cores <- 1L
  gpa$glm_settings <- list(fsl = list(force_l2_creation = TRUE))
  gpa$output_locations$setup_l2_log_json <- tempfile(fileext = ".json")
  gpa$output_locations$setup_l2_log_txt <- tempfile(fileext = ".txt")
  gpa$output_locations$setup_l2_l1_cope_validity_tsv <- tempfile(fileext = ".tsv")
  gpa$output_locations$feat_l2_directory <- file.path(gpa$output_directory, "feat_l2", "sub-{id}", "ses-{session}", "L2m-{l2_model}")

  contrasts <- matrix(1, nrow = 1, ncol = 1, dimnames = list("copeA", "regA"))
  gpa$l1_models$models$model1$contrasts <- contrasts
  gpa$l1_models$n_contrasts <- c(model1 = 1L)
  class(gpa$l1_models) <- c("l1_model_set", "list")
  gpa$l1_cope_names <- list(model1 = "copeA")

  # Create real tiny NIfTI mask + cope + stats/dof in run1.feat; mark run 2 empty
  root <- tempfile("l2_pt_tdof_")
  dir.create(root, recursive = TRUE)
  feat_dirs <- file.path(root, paste0("run", 1:2, ".feat"))
  feat_fsfs <- file.path(root, paste0("run", 1:2, ".fsf"))
  mask_arr <- array(c(1, 1, 0, 1, 0, 0, 1, 1), dim = c(2, 2, 2))
  cope_arr <- array(2.5, dim = c(2, 2, 2))
  for (ii in 1:2) {
    dir.create(file.path(feat_dirs[ii], "stats"), recursive = TRUE)
    file.create(file.path(feat_dirs[ii], ".feat_complete"))
    RNifti::writeNifti(mask_arr, file = file.path(feat_dirs[ii], "mask.nii.gz"))
    RNifti::writeNifti(cope_arr, file = file.path(feat_dirs[ii], "stats", "cope1.nii.gz"))
    RNifti::writeNifti(cope_arr, file = file.path(feat_dirs[ii], "stats", "varcope1.nii.gz"))
    writeLines("94", file.path(feat_dirs[ii], "stats", "dof"))
    writeLines(
      if (ii == 1L) c('set fmri(shape1) 3', 'set fmri(evtitle1) "regA"') else c('set fmri(shape1) 10', 'set fmri(evtitle1) "regA"'),
      feat_fsfs[ii]
    )
  }
  gpa$l1_model_setup <- structure(list(fsl = data.frame(
    id = "sub1", session = 1L, run_number = 1:2, l1_model = "model1",
    feat_fsf = feat_fsfs, feat_dir = feat_dirs, feat_complete = TRUE, feat_failed = FALSE,
    stringsAsFactors = FALSE
  )), class = c("l1_setup", "list"))

  # Model 1: continuous covariate with EV_Intercept (pooling) + EV_run_centered (slope)
  m_pool_slope <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data, level = 2L, cur_model_names = NULL,
    spec_list = list(name = "l2_pool_slope", model_formula = "~ run_centered", diagonal = TRUE, covariate_transform = c(run_centered = "mean")),
    lg = lgr::get_logger("test")
  )
  # Model 2: categorical difference-only model (B - A) -> must NOT passthrough on 1 surviving run!
  m_diff_only <- fmri.pipeline:::create_new_hi_model(
    data = gpa$run_data, level = 2L, cur_model_names = NULL,
    spec_list = list(name = "l2_diff_only", model_formula = "~ cond", diagonal = FALSE, pairwise_diffs = "cond"),
    lg = lgr::get_logger("test")
  )
  gpa$l2_models <- structure(list(models = list(l2_pool_slope = m_pool_slope, l2_diff_only = m_diff_only)), class = c("hi_model_set", "list"))

  res <- setup_l2_models(gpa, l1_model_names = "model1", l2_model_names = c("l2_pool_slope", "l2_diff_only"), backend = "fsl")
  l2_df <- res$l2_model_setup$fsl

  expect_equal(nrow(l2_df), 1L)
  expect_equal(l2_df$l2_model, "l2_pool_slope")
  expect_equal(l2_df$cope_list[[1L]]$l2_cope_name, "EV_Intercept")

  # Masked tdof_t1.nii.gz synthesized with dof = 94 inside mask and 0 outside mask!
  tdof_file <- file.path(feat_dirs[1L], "stats", "tdof_t1.nii.gz")
  expect_true(file.exists(tdof_file))
  tdof_img <- as.array(RNifti::readNifti(tdof_file))
  expect_equal(as.numeric(tdof_img[mask_arr == 1]), rep(94, sum(mask_arr == 1)))
  expect_equal(as.numeric(tdof_img[mask_arr == 0]), rep(0, sum(mask_arr == 0)))
})

test_that("lookup_feat_l3_metadata deduplicates shifted l2_cope_number across subjects", {
  gpa <- list(
    subject_data = data.frame(id = c("s1", "s2"), session = c(1L, 1L), exclude_subject = FALSE, stringsAsFactors = FALSE),
    l1_cope_names = list(m1 = "copeA"),
    l2_model_setup = structure(list(fsl = tibble::tibble(
      id = c("s1", "s2"),
      session = c(1L, 1L),
      l1_model = c("m1", "m1"),
      l1_cope_number = c(1L, 1L),
      l1_cope_name = c("copeA", "copeA"),
      l2_model = c("l2m", "l2m"),
      l2_input_mode = c("cope_files", "cope_files"),
      cope_list = list(
        data.frame(id = "s1", session = 1L, l2_cope_number = 1:2, l2_cope_name = c("EV_Intercept", "shared_con"), stringsAsFactors = FALSE),
        data.frame(id = "s2", session = 1L, l2_cope_number = 1L, l2_cope_name = "shared_con", stringsAsFactors = FALSE)
      )
    )), class = c("l2_setup", "list")),
    l3_models = list(models = list(l3m = list(
      contrasts = matrix(1, nrow = 1, ncol = 1, dimnames = list("group_mean", "(Intercept)"))
    )))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  meta <- fmri.pipeline:::lookup_feat_l3_metadata(
    gpa,
    data.frame(
      l1_model = "m1",
      l1_cope_name = "copeA",
      l2_model = "l2m",
      l2_cope_name = "shared_con",
      l3_model = "l3m",
      stringsAsFactors = FALSE
    )
  )
  expect_equal(sum(meta$l2_cope_name == "shared_con"), 1L)
})

# ==============================================================================
# 4. Level 3 Pure FSL & Hybrid FSL + AFNI (3dLMEr & 3dClusterize)
# ==============================================================================

test_that("respecify_l3_model handles pooled_sessions_subject_ev mixed 2x2 ANOVA and num2fac='session' cleanly", {
  sdat <- expand.grid(
    id = c("s1", "s2", "s3", "s4"),
    session = 1:2,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  sdat$group <- ifelse(sdat$id %in% c("s1", "s2"), "pat", "ctl")
  sdat$session_label <- ifelse(sdat$session == 1L, "pre", "post")

  m_mixed <- fmri.pipeline:::create_new_hi_model(
    data = sdat,
    level = 3L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l3_mixed",
      model_formula = "~ group * session_label",
      l3_input_mode = "pooled_sessions_subject_ev",
      diagonal = FALSE,
      pairwise_diffs = c("group", "session_label"),
      interaction_diffs = "group:session_label"
    ),
    lg = lgr::get_logger("test")
  )

  expect_warning(
    refit <- fmri.pipeline:::respecify_l3_model(m_mixed, sdat[, c("id", "session")]),
    "pure between-subject effects are aliased"
  )

  # Full-rank model_matrix_noalias used for FSL FEAT
  dmat <- refit$model_matrix_noalias %||% refit$model_matrix
  expect_equal(qr(dmat)$rank, ncol(dmat))
  # Subject indicator EV_id* diagonals suppressed AND aliased group main-effect contrast pruned!
  expect_false(any(grepl("^EV_id", rownames(refit$contrasts))))
  expect_false(any(grepl("^group\\.", rownames(refit$contrasts))))
  # Within-subject session_label main effect AND group x session_label interaction survive!
  expect_true(any(grepl("^session_label\\.", rownames(refit$contrasts))))
  expect_true(any(grepl("^group_x_session_label\\.", rownames(refit$contrasts))))

  # Also verify num2fac = "session" joins cleanly against integer new_data$session
  m_sess_fac <- fmri.pipeline:::create_new_hi_model(
    data = sdat,
    level = 3L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l3_sess_fac",
      model_formula = "~ session",
      num2fac = "session",
      diagonal = FALSE,
      pairwise_diffs = "session"
    ),
    lg = lgr::get_logger("test")
  )
  refit_sess <- fmri.pipeline:::respecify_l3_model(m_sess_fac, data.frame(id = sdat$id, session = as.integer(sdat$session)))
  expect_equal(nrow(refit_sess$model_matrix), nrow(sdat))
  expect_true(any(grepl("^session\\.", rownames(refit_sess$contrasts))))
})

test_that("emmeans_to_3dlmer_glt and validate_3dlmer_glt_table support full factorial & continuous slope specs", {
  dt <- expand.grid(
    Subj = paste0("s", 1:6),
    session_label = c("pre", "post"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dt$id <- dt$Subj
  dt$session <- ifelse(dt$session_label == "pre", 1L, 2L)
  dt$group <- ifelse(dt$Subj %in% c("s1", "s2", "s3"), "pat", "ctl")
  dt$age <- seq_len(nrow(dt)) + 20
  dt$InputFile <- paste0("cope", seq_len(nrow(dt)), ".nii.gz")

  mobj <- fmri.pipeline:::create_new_hi_model(
    data = dt,
    level = 3L,
    cur_model_names = NULL,
    spec_list = list(
      name = "lmer_full",
      model_formula = "~ group * session_label + group * age",
      l3_input_mode = "3dlmer",
      lmer_formula = "group * session_label + group * age + (1 | Subj)",
      diagonal = FALSE,
      cond_means = c("group", "session_label | group"),
      cell_means = TRUE,
      overall_response = TRUE,
      pairwise_diffs = c("group", "session_label | group"),
      interaction_diffs = "group:session_label",
      simple_slopes = list(age = "group"),
      slope_diffs = list(age = "group"),
      custom = list(age_overall = c("age" = 1)),
      delete = "overall",
      lmer_glt_codes = list(manual_pat_slope = "group : 1*pat age :")
    ),
    lg = lgr::get_logger("test")
  )

  glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, dt, qVars = "age", raw_glt_codes = mobj$lmer_glt_codes, context = "unit test")
  expect_true("age_overall" %in% glts$label_afni)
  expect_true("manual_pat_slope" %in% glts$label_afni)
  expect_false("overall" %in% glts$label_afni) # deleted via spec$delete
  expect_true(any(grepl("group : 1\\*.* -1\\*.* session_label : 1\\*.* -1\\*", glts$code)))
  expect_true(any(grepl("group : 1\\*pat age :", glts$code)))
  expect_true(any(grepl("group : 1\\*ctl -1\\*pat age :|group : 1\\*pat -1\\*ctl age :", glts$code)))

  # All generated + manual GLTs pass validate_3dlmer_glt_table!
  validated <- fmri.pipeline:::validate_3dlmer_glt_table(glts, dt, qVars = "age", context = "unit test")
  expect_equal(nrow(validated), nrow(glts))
})

test_that("afni_3dlmer_setup emits qVarCenters, numeric session qVars, random_effects_vars, NA filtering, and single-run L1->3dLMEr harvesting", {
  skip_if_not_installed("RNifti")
  tmp_dir <- tempfile("lmer_single_run_")
  dir.create(tmp_dir, recursive = TRUE)

  # Single-run multi-session longitudinal dataset (4 subjects x 2 sessions, 1 subject has NA age on session 2, 1 run has invalid L1 cope)
  sdat <- expand.grid(id = paste0("s", 1:4), session = 1:2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  sdat$exclude_subject <- FALSE
  sdat$age <- c(20, 22, 24, 26, 20, 22, 24, NA_real_)
  sdat$site <- rep(c("siteA", "siteB"), length.out = nrow(sdat))

  rdat <- sdat[, c("id", "session")]
  rdat$run_number <- 1L
  rdat$exclude_run <- FALSE
  rdat$exclude_subject <- FALSE

  mask_arr <- array(1, dim = c(2, 2, 2))
  cope_arr <- array(1.5, dim = c(2, 2, 2))
  feat_dirs <- character(nrow(rdat))
  feat_fsfs <- character(nrow(rdat))
  for (ii in seq_len(nrow(rdat))) {
    fdir <- file.path(tmp_dir, sprintf("%s_ses%d.feat", rdat$id[ii], rdat$session[ii]))
    fsf <- file.path(tmp_dir, sprintf("%s_ses%d.fsf", rdat$id[ii], rdat$session[ii]))
    dir.create(file.path(fdir, "reg_standard", "stats"), recursive = TRUE)
    dir.create(file.path(fdir, "stats"), recursive = TRUE)
    file.create(file.path(fdir, ".feat_complete"))
    RNifti::writeNifti(mask_arr, file = file.path(fdir, "reg_standard", "mask.nii.gz"))
    RNifti::writeNifti(cope_arr, file = file.path(fdir, "reg_standard", "stats", "cope1.nii.gz"))
    RNifti::writeNifti(cope_arr, file = file.path(fdir, "reg_standard", "stats", "varcope1.nii.gz"))
    RNifti::writeNifti(cope_arr, file = file.path(fdir, "stats", "cope1.nii.gz"))
    RNifti::writeNifti(cope_arr, file = file.path(fdir, "stats", "varcope1.nii.gz"))
    # Mark s3 session 2 as empty EV so build_l1_cope_validity drops it!
    is_empty <- rdat$id[ii] == "s3" && rdat$session[ii] == 2L
    writeLines(c(sprintf("set fmri(shape1) %d", if (is_empty) 10L else 3L), 'set fmri(evtitle1) "regA"'), fsf)
    feat_dirs[ii] <- fdir
    feat_fsfs[ii] <- fsf
  }

  l3_mod <- fmri.pipeline:::create_new_hi_model(
    data = sdat[!is.na(sdat$age), ],
    level = 3L,
    cur_model_names = NULL,
    spec_list = list(
      name = "l3_long",
      model_formula = "~ session + age",
      l3_input_mode = "3dlmer",
      lmer_formula = "session + age + (session | Subj) + (1 | site)",
      covariate_transform = c(age = "mean"),
      lmer_glt_codes = list(sess_slope = "session :", age_slope = "age :")
    ),
    lg = lgr::get_logger("test")
  )

  gpa <- list(
    multi_run = FALSE,
    output_directory = tmp_dir,
    output_locations = list(
      afni_3dlmer_directory = file.path("{gpa$output_directory}", "afni_3dlmer", "L3m-{l3_model}", "L1m-{l1_model}", "l1c-{l1_cope_name}"),
      feat_l3_directory = file.path("{gpa$output_directory}", "feat_l3")
    ),
    parallel = list(afni = list(l3_lmer_njobs = 2L), fsl = list(l3_setup_cores = 1L)),
    lgr_threshold = "warn",
    log_json = FALSE,
    log_txt = FALSE,
    subject_data = sdat,
    run_data = rdat,
    l1_cope_names = list(m1 = "copeA"),
    l1_models = structure(list(
      models = list(m1 = list(contrasts = matrix(1, 1, 1, dimnames = list("copeA", "regA")))),
      n_contrasts = c(m1 = 1L)
    ), class = c("l1_model_set", "list")),
    l1_model_setup = structure(list(fsl = data.frame(
      id = rdat$id, session = rdat$session, run_number = 1L, l1_model = "m1",
      feat_fsf = feat_fsfs, feat_dir = feat_dirs, feat_complete = TRUE, feat_failed = FALSE,
      stringsAsFactors = FALSE
    )), class = c("l1_setup", "list")),
    l2_models = NULL,
    l3_models = structure(list(models = list(l3_long = l3_mod)), class = c("hi_model_set", "list"))
  )
  class(gpa) <- c("glm_pipeline_arguments", "list")

  res <- fmri.pipeline:::afni_3dlmer_setup(
    gpa = gpa,
    backend = NULL,
    lg = lgr::get_logger("test"),
    l1_model_names = "m1",
    l2_model_names = NULL,
    l3_model_names = "l3_long",
    l2_l3_pairs = NULL,
    subj_df = sdat[, c("id", "session")],
    requires_l2 = FALSE
  )

  expect_equal(nrow(res$data), 1L)
  script_path <- res$data$afni_script[1L]
  dt_written <- utils::read.delim(file.path(dirname(script_path), "dataTable.txt"), stringsAsFactors = FALSE)
  # 8 total subject-sessions minus 1 empty L1 EV (s3 ses 2) minus 1 NA age (s4 ses 2) = 6 valid rows!
  expect_equal(nrow(dt_written), 6L)
  expect_true("site" %in% names(dt_written))
  expect_false(anyNA(dt_written$age))
  expect_true(all(grepl("reg_standard/stats/cope1\\.nii\\.gz$", dt_written$InputFile)))

  sh_lines <- paste(readLines(script_path), collapse = "\n")
  # Both numeric session and age in -qVars, NOT site! And -qVarCenters emitted!
  expect_match(sh_lines, "-qVars 'age,session'|-qVars 'session,age'")
  expect_match(sh_lines, "-qVarCenters")
  expect_true(file.exists(file.path(dirname(script_path), "intersection_mask.nii.gz")))
})

test_that("afni_3dclusterize emits -1sided <SIDE> <THRESH> in AFNI-required argument order", {
  tmp_nii <- tempfile(fileext = ".nii.gz")
  file.create(tmp_nii)
  obj <- afni_3dclusterize$new(
    threshold_file = tmp_nii,
    onesided = TRUE,
    one_thresh = 3.09,
    one_tail = "RIGHT_TAIL",
    clust_nvox = 20L
  )
  call_str <- obj$get_call()
  expect_match(call_str, "-1sided RIGHT_TAIL 3\\.09")
})
