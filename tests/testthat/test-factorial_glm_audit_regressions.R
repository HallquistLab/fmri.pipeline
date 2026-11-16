test_that("unbalanced formula-coded L1 means recover known cell amplitudes", {
  values <- data.frame(
    run_number = 1L, trial = 1:5, value = 1,
    A = c("a", "a", "a", "b", "a")
  )
  events <- data.frame(
    event = "stim", run_number = 1L, trial = 1:5,
    onset = c(4, 20, 36, 52, 104), duration = 1
  )
  runs <- data.frame(run_number = 1L, run_volumes = 100L, drop_volumes = 0L, tr = 1)
  # Generate the response from two independently specified unit event trains.
  cell_signals <- lapply(c("a", "b"), function(lv) {
    list(name = lv, event = "stim", value = values[values$A == lv, 1:3],
         normalization = "none", beta_series = FALSE)
  })
  cell_bdm <- build_design_matrix(events, cell_signals, tr = 1, run_data = runs, plot = FALSE)
  y <- 0.5 + as.matrix(cell_bdm$design_convolved[[1L]]) %*% c(5, -1)

  for (coding in c("sum", "helmert", "poly")) {
    sig <- fmri.pipeline:::fit_wi_model(list(
      name = "task", event = "stim", value = values, value_type = "unit",
      normalization = "none", beta_series = FALSE,
      wi_factors = "A", wi_formula = "~ A", wi_expansion = "formula", factor_coding = coding
    ))
    wi <- structure(list(level = 1L, lmfit = sig$wi_model, signal_name = "task"),
                    class = c("l1_wi_spec", "list"))
    wi <- fmri.pipeline:::specify_contrasts(wi, spec_list = list(diagonal = FALSE, cond_means = "A"))
    for (within_run in c(TRUE, FALSE)) {
      bdm <- build_design_matrix(
        events, list(sig), tr = 1, run_data = runs, plot = FALSE,
        center_values = TRUE, convolve_wi_run = within_run
      )
      beta <- stats::coef(stats::lm.fit(cbind(1, as.matrix(bdm$design_convolved[[1L]])), y))[-1L]
      expect_equal(as.numeric(wi$contrasts %*% beta), c(5, -1), tolerance = 1e-7,
                   info = paste(coding, within_run))
    }
  }
})

test_that("canonical factor padding never enters emmeans covariate reductions", {
  for (transform in c("none", "min", "max", "mean", "zscore")) {
    x <- c(10, 20, 30, 40, 50, 60)
    x <- switch(transform, min = x - min(x), max = x - max(x),
                mean = x - mean(x), zscore = as.numeric(scale(x)), x)
    df <- data.frame(A = factor(rep("a", 6L), levels = c("a", "b")), x = x)
    df$dummy <- 4 + 2 * x + c(-1, 1, 0, 0, 1, -1)
    result <- fmri.pipeline:::fit_lm_with_canonical_levels(
      dummy ~ A * x, df, canonical_levels = list(A = c("a", "b"))
    )
    fit <- result$lmfit
    rg <- emmeans::ref_grid(fit)
    expect_equal(unique(rg@grid$x), mean(x), tolerance = 1e-10, info = transform)
    expect_equal(rg@grid$.wgt., c(6, 0))
    expect_equal(nrow(stats::model.matrix(fit)[result$obs_idx, , drop = FALSE]), nrow(df))
    expected <- stats::lm(dummy ~ x, df)
    expect_equal(unname(stats::coef(fit)[c("(Intercept)", "x")]), unname(stats::coef(expected)))
    expect_equal(unname(stats::fitted(fit)[result$obs_idx]), unname(stats::fitted(expected)))
    expect_equal(summary(fit)$sigma, summary(expected)$sigma)
    means <- as.data.frame(emmeans::emmeans(fit, ~ A))
    expect_equal(means$emmean[means$A == "a"], mean(df$dummy), tolerance = 1e-7)
    expect_true(is.na(means$emmean[means$A == "b"]))
  }
  # Numeric binary predictors must retain their original reference-grid support.
  df <- data.frame(A = factor(rep("a", 4L), levels = letters[1:5]),
                   x = c(1, 1, 2, 2), dummy = c(1, 3, 4, 3))
  fit <- fmri.pipeline:::fit_lm_with_canonical_levels(dummy ~ A + x, df)$lmfit
  expect_equal(sort(unique(emmeans::ref_grid(fit)@grid$x)), c(1, 2))
})

test_that("FSL L3 writes estimable contrasts on the same basis as the reduced design", {
  df <- data.frame(id = paste0("s", 1:8), session = 1L, A = rep(c("a", "b"), each = 4L))
  parent <- structure(list(name = "m3", level = 3L, model_variables = "A", l3_input_mode = "per_session"),
                      class = c("hi_model_spec", "list"))
  parent <- fmri.pipeline:::mobj_fit_lm(parent, "~ A", df, c("id", "session"))
  parent <- fmri.pipeline:::specify_contrasts(parent, spec_list = list(diagonal = TRUE))
  inputs <- transform(df[1:4, ], l1_model = "m1", l1_cope_name = "signal", l3_model = "m3",
                      cope_file = paste0("cope", 1:4, ".nii.gz"))
  outdir <- tempfile("l3_audit_")
  dir.create(outdir)
  gpa <- structure(list(
    lgr_threshold = "warn", multi_run = FALSE, l3_models = list(models = list(m3 = parent)),
    output_locations = list(feat_l3_fsf = "audit.fsf"), additional = list(feat_l3_args = NULL),
    glm_settings = list(fsl = list(force_l3_creation = TRUE))
  ), class = c("glm_pipeline_arguments", "list"))
  local_mocked_bindings(get_output_directory = function(...) outdir, .package = "fmri.pipeline")
  result <- fmri.pipeline:::fsl_l3_model(inputs, gpa)
  fsf <- readLines(result$feat_fsf)
  expect_true("set fmri(evs_real) 1" %in% fsf)
  expect_true("set fmri(ncon_real) 1" %in% fsf)
  expect_true("set fmri(con_real1.1) 1" %in% fsf)
  expect_false(any(grepl("^set fmri\\(con_real[0-9]+\\.2\\)", fsf)))
  expect_identical(result$cope_list[[1L]]$l3_cope_name, "EV_Intercept")
})

test_that("FSL L3 metadata preserves each analysis and session's COPE mapping", {
  foreach::registerDoSEQ()
  units <- expand.grid(id = paste0("s", 1:4), session = 1:2, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  analyses <- expand.grid(l1_model = c("m1", "m2"), l1_cope_name = c("C1", "C2"),
                          l2_model = c("pool", "difference"), l2_cope_name = c("mean", "slope"),
                          l3_model = c("per", "pooled"), stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  cfg <- dplyr::cross_join(units, analyses) %>%
    dplyr::mutate(l3_cope_number = 1L, l3_cope_name = "parent")
  inputs <- split(cfg, interaction(cfg$l1_model, cfg$l1_cope_name, cfg$l2_model, cfg$l2_cope_name, cfg$l3_model,
                                  ifelse(cfg$l3_model == "per", cfg$session, 0L), drop = TRUE))
  backend <- list(l3_setup = function(info, gpa) {
    result <- info[1L, c("l1_model", "l1_cope_name", "l2_model", "l2_cope_name", "l3_model")]
    result$session <- if (result$l3_model == "pooled") 0L else info$session[1L]
    result$cope_list <- list(if (result$l1_model == "m1" && result$l1_cope_name == "C1" &&
                                result$l2_model == "pool" && result$l2_cope_name == "mean" && result$session != 2L) {
      data.frame(l3_cope_number = 1:2, l3_cope_name = c("x", "y"))
    } else data.frame(l3_cope_number = 1L, l3_cope_name = "y"))
    result
  })
  gpa <- structure(list(multi_run = TRUE), class = c("glm_pipeline_arguments", "list"))
  local_mocked_bindings(
    get_fsl_l3_model_df = function(...) cfg,
    get_feat_l3_inputs = function(...) inputs,
    .package = "fmri.pipeline"
  )
  result <- fmri.pipeline:::setup_l3_backend_fsl(
    gpa, backend, lgr::get_logger("test"), c("m1", "m2"), c("pool", "difference"),
    c("per", "pooled"), NULL, units, requires_l2 = TRUE
  )
  for (ii in seq_len(nrow(result$data))) {
    setup <- result$data[ii, ]
    meta <- result$metadata %>%
      dplyr::filter(.data$l1_model == setup$l1_model, .data$l1_cope_name == setup$l1_cope_name,
                    .data$l2_model == setup$l2_model, .data$l2_cope_name == setup$l2_cope_name,
                    .data$l3_model == setup$l3_model)
    if (setup$session != 0L) meta <- dplyr::filter(meta, .data$session == setup$session)
    actual <- dplyr::distinct(meta, .data$l3_cope_number, .data$l3_cope_name)
    expected <- if (setup$l1_model == "m1" && setup$l1_cope_name == "C1" &&
                    setup$l2_model == "pool" && setup$l2_cope_name == "mean" && setup$session != 2L) {
      data.frame(l3_cope_number = 1:2, l3_cope_name = c("x", "y"))
    } else data.frame(l3_cope_number = 1L, l3_cope_name = "y")
    expect_equal(dplyr::arrange(actual, .data$l3_cope_number), expected)
    expect_false(anyDuplicated(meta[, c("id", "session", "l3_cope_number")]) > 0L)
    if (setup$session == 0L) expect_equal(sort(unique(meta$session)), 1:2)
  }
})

test_that("AFNI overall GLTs preserve supported weights and reject nonseparable cell means", {
  make_model <- function(counts, weights) {
    df <- expand.grid(A = factor(c("a", "b")), B = factor(c("pre", "post")))
    df$dummy <- c(10, 0, 0, 10)
    df <- df[rep(1:4, counts), ]
    df$dummy <- df$dummy + rep(c(-1, 1), length.out = nrow(df))
    list(data = df, model = structure(list(
      lmfit = stats::lm(dummy ~ A * B, df),
      contrast_spec = list(diagonal = FALSE, overall_response = TRUE, weights = weights, cat_vars = c("A", "B"))
    ), class = c("hi_model_spec", "list")))
  }
  for (weights in c("cells", "proportional")) {
    bad <- make_model(c(8, 2, 2, 8), weights)
    expect_error(fmri.pipeline:::emmeans_to_3dlmer_glt(bad$model, bad$data), "joint cell weights are not separable")
  }
  # Evaluate the generated factor-level blocks as a contrast of fitted cell means.
  for (weights in c("equal", "flat", "cells", "proportional")) {
    fit <- make_model(c(8, 4, 4, 2), weights)
    glt <- fmri.pipeline:::emmeans_to_3dlmer_glt(fit$model, fit$data)
    cells <- emmeans::emmeans(fit$model$lmfit, ~ A * B)
    blocks <- fmri.pipeline:::parse_3dlmer_glt_blocks(glt$code, names(fit$data))$blocks
    joint <- rep(1, nrow(cells@grid))
    for (block in blocks) {
      tokens <- strsplit(block$expr, "\\s+")[[1L]]
      ws <- stats::setNames(as.numeric(sub("\\*.*", "", tokens)), sub("^[^*]*\\*", "", tokens))
      joint <- joint * ws[as.character(cells@grid[[block$var]])]
    }
    expected <- emmeans::emmeans(fit$model$lmfit, ~ 1, weights = weights)@linfct
    expect_equal(unname(as.numeric(joint %*% cells@linfct)), unname(as.numeric(expected)), tolerance = 1e-8)
  }
  flat <- make_model(c(8, 2, 2, 8), "flat")
  expect_no_error(fmri.pipeline:::emmeans_to_3dlmer_glt(flat$model, flat$data))
})

test_that("AFNI slope differences match emtrends within every conditioning slice", {
  df <- expand.grid(A = factor(c("a", "b")), B = factor(c("pre", "post")),
                    C = factor(c("left", "right")), age = 1:4)
  df$dummy <- with(df, age * ifelse(A == "a", ifelse(B == "pre", 2, 7), ifelse(B == "pre", -1, 3)))
  df$dummy <- df$dummy + rep(c(-1, 1, 1, -1), each = 8L)
  m <- structure(list(
    lmfit = stats::lm(dummy ~ A * B * C * age, df),
    contrast_spec = list(diagonal = FALSE, slope_diffs = list(age = "A | B:C"),
                         weights = "equal", cat_vars = c("A", "B", "C"))
  ), class = c("hi_model_spec", "list"))
  glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(m, df, qVars = "age")
  expect_equal(nrow(glts), 4L)
  trends <- suppressMessages(emmeans::emtrends(m$lmfit, ~ A * B * C, var = "age"))
  for (ii in seq_len(nrow(glts))) {
    blocks <- fmri.pipeline:::parse_3dlmer_glt_blocks(glts$code[ii], names(df))$blocks
    weights <- rep(1, nrow(trends@grid))
    for (block in blocks) {
      if (block$var == "age") next
      tokens <- strsplit(block$expr, "\\s+")[[1L]]
      ws <- stats::setNames(as.numeric(sub("\\*.*", "", tokens)), sub("^[^*]*\\*", "", tokens))
      w <- ws[as.character(trends@grid[[block$var]])]
      w[is.na(w)] <- 0
      weights <- weights * w
    }
    b <- sub(".*B : 1\\*([^ ]+).*", "\\1", glts$code[ii])
    c_level <- sub(".*C : 1\\*([^ ]+).*", "\\1", glts$code[ii])
    pair <- suppressMessages(emmeans::contrast(trends, "pairwise", by = c("B", "C")))
    expected <- pair@linfct[pair@grid$B == b & pair@grid$C == c_level, , drop = FALSE]
    expect_equal(as.numeric(weights %*% trends@linfct), as.numeric(expected), tolerance = 1e-7)
  }
  # A missing arm must remove only its own conditioned comparison.
  df <- subset(df, !(A == "b" & B == "pre" & C == "left"))
  partial <- fmri.pipeline:::emmeans_to_3dlmer_glt(m, df, qVars = "age")
  expect_equal(nrow(partial), 3L)
  m$contrast_spec$slope_diffs <- list(age = "A:B | C")
  expect_error(fmri.pipeline:::emmeans_to_3dlmer_glt(m, df, qVars = "age"), "one focal factor")
})

test_that("alias-reduced observed L1 means survive validity checks without hiding missing inputs", {
  tmp_dir <- tempfile("l1_alias_validity_")
  dir.create(tmp_dir)
  td <- data.frame(
    id = "sub1", session = 1L, run_number = c(rep(1L, 4L), rep(2L, 2L)),
    trial = c(1:4, 1:2), onset = c(4, 12, 20, 28, 4, 12), duration = 1,
    A = c("a", "b", "a", "b", "b", "b")
  )
  niftis <- file.path(tmp_dir, paste0("run", 1:2, ".nii.gz"))
  for (nf in niftis) RNifti::writeNifti(array(0, c(2, 2, 2, 60)), file = nf)
  gpa <- create_mock_gpa(n_subjects = 1, n_runs = 2, output_directory = tmp_dir)
  gpa$trial_data <- td
  gpa$subject_data$mr_dir <- tmp_dir
  gpa$run_data <- data.frame(
    id = "sub1", session = 1L, run_number = 1:2, mr_dir = tmp_dir, run_nifti = niftis,
    run_nifti_present = TRUE, run_volumes = 60L, drop_volumes = 0L, tr = 1,
    nvoxels = 1000L, exclude_run = FALSE, exclude_subject = FALSE
  )
  gpa$use_preconvolve <- TRUE
  gpa$lgr_threshold <- "warn"
  gpa$log_json <- gpa$log_txt <- FALSE
  gpa$output_locations$feat_l1_directory <- file.path(tmp_dir, "sub-{id}", "ses-{session}", "{l1_model}")
  gpa$glm_settings <- list(fsl = list(force_l1_creation = TRUE, smooth_fwhm = 5, hp_filter = 100,
                                    prewhiten = TRUE, z_thresh = 2.3, p_thresh = 0.05))
  gpa$confound_settings <- list(l1_confound_regressors = NULL)
  spec <- list(
    onsets = list("onset"), durations = list("duration"),
    events = list(stim = list(onset = "onset", duration = "duration")),
    signals = list(task = list(event = "stim", normalization = "none", value_fixed = 1,
                              wi_factors = list("A"), wi_formula = "~ A", wi_expansion = "formula")),
    l1_models = list(m1 = list(signals = list("task"), contrasts = list(
      diagonal = FALSE, wi_contrasts = list(task = list(diagonal = FALSE, cond_means = "A"))
    )))
  )
  spec_file <- file.path(tmp_dir, "spec.yaml")
  yaml::write_yaml(spec, spec_file)
  gpa <- build_l1_models(gpa, from_spec_file = spec_file)
  bdm <- build_design_matrix(
    events = gpa$l1_models$events$stim$data, signals = list(gpa$l1_models$signals$task),
    tr = 1, run_data = gpa$run_data, write_timing_files = "convolved",
    output_directory = file.path(tmp_dir, "timing"), plot = FALSE
  )
  bdm$run_nifti <- niftis
  setup <- fsl_l1_model(id = "sub1", session = 1L, d_obj = bdm, gpa = gpa,
                        model_name = "m1", run_nifti = niftis, nvoxels = c(1000L, 1000L))
  validity <- fmri.pipeline:::build_l1_cope_validity(setup, gpa)
  r2 <- subset(validity, run_number == 2L)
  expect_true(r2$valid_for_higher_level[r2$l1_cope_name == "task.A.b"])
  expect_false(r2$valid_for_higher_level[r2$l1_cope_name == "task.A.a"])

  # Completing a run still requires the real COPE and VARCOPE files.
  setup$feat_complete <- TRUE
  validity <- fmri.pipeline:::build_l1_cope_validity(setup, gpa)
  expect_true(all(validity$validity_state == "invalid_missing_file"))

  setup$feat_complete <- FALSE
  fsf <- readLines(setup$feat_fsf[2L])
  writeLines(sub("^set fmri\\(shape1\\) [0-9]+$", "set fmri(shape1) 10", fsf), setup$feat_fsf[2L])
  validity <- fmri.pipeline:::build_l1_cope_validity(setup, gpa)
  expect_false(any(subset(validity, run_number == 2L)$valid_for_higher_level))
  writeLines(fsf, setup$feat_fsf[2L])

  sidecar_file <- sub("\\.fsf$", "_cope_validity.tsv", setup$feat_fsf[2L])
  sidecar <- utils::read.delim(sidecar_file)
  sidecar$l1_cope_name <- "stale_name"
  utils::write.table(sidecar, sidecar_file, sep = "\t", row.names = FALSE, quote = FALSE)
  validity <- fmri.pipeline:::build_l1_cope_validity(setup, gpa)
  expect_false(any(subset(validity, run_number == 2L)$valid_for_higher_level))
})
