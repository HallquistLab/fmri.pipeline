make_lookup_feat_outputs_gpa <- function(root, l1_feat = NULL, l2_feat = NULL, l3_feat = NULL) {
  if (is.null(l1_feat)) {
    l1_feat <- file.path(root, "feat_l1", "sub-01", "ses-1", "facehouse", "FEAT_LVL1_run1.feat")
  }
  if (is.null(l2_feat)) {
    l2_feat <- file.path(
      root, "feat_l2", "sub-01", "ses-1", "L1m-facehouse",
      "l1c-01_EV_face", "L2m-l2_session.gfeat"
    )
  }
  if (is.null(l3_feat)) {
    l3_feat <- file.path(
      root, "feat_l3", "L3m-l3_group", "L1m-facehouse",
      "L2m-l2_session", "l2c-overall", "FEAT_l1c-EV_face.gfeat"
    )
  }

  l1_setup <- data.frame(
    id = "01",
    session = 1L,
    run_number = 1L,
    run_volumes = 100L,
    run_nifti = file.path(root, "run.nii.gz"),
    l1_model = "facehouse",
    l1_confound_file = NA_character_,
    to_run = FALSE,
    feat_fsf = sub("\\.feat$", ".fsf", l1_feat),
    feat_dir = l1_feat,
    feat_complete = TRUE,
    feat_failed = FALSE,
    feat_dir_exists = dir.exists(l1_feat),
    feat_fsf_exists = FALSE,
    stringsAsFactors = FALSE
  )

  l2_setup <- data.frame(
    id = "01",
    session = 1L,
    l1_model = "facehouse",
    l1_cope_number = 1L,
    l1_cope_name = "EV_face",
    l2_model = "l2_session",
    l2_scope = "id_session",
    l2_input_mode = "cope_files",
    l2_passthrough = FALSE,
    n_l2_copes = 2L,
    passthrough_cope_file = NA_character_,
    feat_fsf = sub("\\.gfeat$", ".fsf", l2_feat),
    feat_dir = l2_feat,
    feat_complete = TRUE,
    feat_failed = FALSE,
    feat_dir_exists = dir.exists(l2_feat),
    feat_fsf_exists = FALSE,
    stringsAsFactors = FALSE
  )
  l2_setup$cope_list <- list(data.frame(
    id = c("01", "01"),
    session = c(1L, 1L),
    l2_cope_number = c(1L, 2L),
    l2_cope_name = c("overall", "stress"),
    stringsAsFactors = FALSE
  ))

  l3_setup <- data.frame(
    l1_model = "facehouse",
    l1_cope_name = "EV_face",
    l2_model = "l2_session",
    l2_cope_name = "overall",
    l3_model = "l3_group",
    session = 0L,
    l3_input_mode = "all_sessions",
    feat_fsf = sub("\\.gfeat$", ".fsf", l3_feat),
    feat_dir = l3_feat,
    feat_complete = TRUE,
    feat_failed = FALSE,
    feat_dir_exists = dir.exists(l3_feat),
    feat_fsf_exists = FALSE,
    stringsAsFactors = FALSE
  )

  l3_metadata <- data.frame(
    id = c("01", "01"),
    session = c(1L, 1L),
    l1_model = c("facehouse", "facehouse"),
    l1_cope_number = c(1L, 1L),
    l1_cope_name = c("EV_face", "EV_face"),
    l2_model = c("l2_session", "l2_session"),
    l2_cope_number = c(1L, 1L),
    l2_cope_name = c("overall", "overall"),
    l3_model = c("l3_group", "l3_group"),
    l3_cope_number = c(1L, 2L),
    l3_cope_name = c("group", "age"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      output_directory = root,
      lgr_threshold = "warn",
      multi_run = TRUE,
      subject_data = data.frame(id = "01", session = 1L, exclude_subject = FALSE),
      run_data = data.frame(id = "01", session = 1L, run_number = 1L),
      l1_cope_names = list(facehouse = c("EV_face", "EV_house")),
      l1_models = list(models = list(facehouse = list(
        contrasts = matrix(
          1,
          nrow = 2L,
          ncol = 1L,
          dimnames = list(c("EV_face", "EV_house"), "Intercept")
        )
      ))),
      l2_models = list(models = list(l2_session = list(l2_scope = "id_session"))),
      l3_models = list(models = list(l3_group = list(
        contrasts = matrix(
          1,
          nrow = 2L,
          ncol = 1L,
          dimnames = list(c("group", "age"), "Intercept")
        )
      ))),
      l1_model_setup = structure(list(fsl = l1_setup), class = c("l1_setup", "list")),
      l2_model_setup = structure(list(fsl = l2_setup), class = c("l2_setup", "list")),
      l3_model_setup = structure(
        list(metadata = l3_metadata, fsl = l3_setup),
        class = c("l3_setup", "list")
      )
    ),
    class = "glm_pipeline_arguments"
  )
}

make_lookup_feat_stats <- function(feat_dir, contrast_names) {
  stats_dir <- file.path(feat_dir, "stats")
  dir.create(stats_dir, recursive = TRUE)
  for (ii in seq_along(contrast_names)) {
    file.create(file.path(stats_dir, paste0("cope", ii, ".nii.gz")))
    file.create(file.path(stats_dir, paste0("varcope", ii, ".nii.gz")))
    file.create(file.path(stats_dir, paste0("zstat", ii, ".nii.gz")))
    file.create(file.path(stats_dir, paste0("tstat", ii, ".nii.gz")))
  }
  writeLines(
    paste0("/ContrastName", seq_along(contrast_names), "\t", contrast_names),
    file.path(feat_dir, "design.con")
  )
  invisible(feat_dir)
}

test_that("lookup_feat_outputs maps L1 FEAT contrast outputs", {
  root <- tempfile("lookup_feat_outputs_")
  l1_feat <- file.path(root, "feat_l1", "sub-01", "ses-1", "facehouse", "FEAT_LVL1_run1.feat")
  dir.create(file.path(l1_feat, "stats"), recursive = TRUE)
  file.create(file.path(l1_feat, "stats", "cope1.nii.gz"))

  gpa <- make_lookup_feat_outputs_gpa(root, l1_feat = l1_feat)

  result <- lookup_feat_outputs(gpa, level = 1L, what = c("cope", "zstat"))
  out <- result$l1

  expect_s3_class(result, "feat_output_lookup")
  expect_named(result, c("l1", "l2", "l3"))
  expect_equal(nrow(out), 4L)
  expect_equal(unique(out$output_level), 1L)
  expect_true(all(c("EV_face", "EV_house") %in% out$l1_cope_name))
  expect_true(any(out$image_exists))
  expect_true(file.path(l1_feat, "stats", "cope1.nii.gz") %in% out$image_file)

  long <- lookup_feat_outputs(
    gpa, level = 1L, what = c("cope", "zstat"), format = "long"
  )
  expect_s3_class(long, "data.frame")
  expect_equal(nrow(long), 4L)
  expect_equal(unique(long$output_model), "facehouse")
  expect_equal(long, as.data.frame(result))
})

test_that("lookup_feat_outputs maps L2 FEAT cope-list outputs", {
  root <- tempfile("lookup_feat_outputs_")
  l2_feat <- file.path(
    root, "feat_l2", "sub-01", "ses-1", "L1m-facehouse",
    "l1c-01_EV_face", "L2m-l2_session.gfeat"
  )
  dir.create(file.path(l2_feat, "cope1.feat", "stats"), recursive = TRUE)
  file.create(file.path(l2_feat, "cope1.feat", "stats", "zstat2.nii.gz"))

  gpa <- make_lookup_feat_outputs_gpa(root, l2_feat = l2_feat)

  out <- lookup_feat_outputs(gpa, level = 2L, what = "zstat")$l2

  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$output_level), 2L)
  expect_equal(sort(out$l2_cope_name), c("overall", "stress"))
  expect_true(file.path(l2_feat, "cope1.feat", "stats", "zstat2.nii.gz") %in% out$image_file)
  expect_true(out$image_exists[out$l2_cope_name == "stress"])
  expect_false("cope_list" %in% names(out))
  expect_equal(unique(out$l2_input_mode), "cope_files")
  expect_true(all(!out$is_passthrough))

  debug_out <- lookup_feat_outputs(
    gpa, level = 2L, what = "zstat", include_internal = TRUE
  )$l2
  expect_true("l2_input_mode" %in% names(debug_out))
  expect_false("cope_list" %in% names(debug_out))
})

test_that("lookup_feat_outputs maps L3 FEAT metadata outputs", {
  root <- tempfile("lookup_feat_outputs_")
  l3_feat <- file.path(
    root, "feat_l3", "L3m-l3_group", "L1m-facehouse",
    "L2m-l2_session", "l2c-overall", "FEAT_l1c-EV_face.gfeat"
  )
  dir.create(file.path(l3_feat, "cope1.feat", "stats"), recursive = TRUE)
  file.create(file.path(l3_feat, "cope1.feat", "stats", "zstat2.nii.gz"))

  gpa <- make_lookup_feat_outputs_gpa(root, l3_feat = l3_feat)

  out <- lookup_feat_outputs(gpa, level = 3L, what = "zstat")$l3

  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$output_level), 3L)
  expect_equal(sort(out$l3_cope_name), c("age", "group"))
  expect_true(file.path(l3_feat, "cope1.feat", "stats", "zstat2.nii.gz") %in% out$image_file)
  expect_true(out$image_exists[out$l3_cope_name == "age"])
})

test_that("lookup_feat_outputs falls back to scheduler caches", {
  root <- tempfile("lookup_feat_outputs_")
  l2_feat <- file.path(
    root, "feat_l2", "sub-01", "ses-1", "L1m-facehouse",
    "l1c-01_EV_face", "L2m-l2_session.gfeat"
  )
  dir.create(file.path(l2_feat, "cope1.feat", "stats"), recursive = TRUE)
  file.create(file.path(l2_feat, "cope1.feat", "stats", "zstat1.nii.gz"))

  cached_gpa <- make_lookup_feat_outputs_gpa(root, l2_feat = l2_feat)
  gpa <- cached_gpa
  gpa$l1_model_setup <- NULL
  gpa$l2_model_setup <- NULL
  gpa$l3_model_setup <- NULL

  cache_dir <- file.path(root, "scheduler_scripts", "batch_abc")
  dir.create(cache_dir, recursive = TRUE)
  cache_file <- file.path(cache_dir, "run_pipeline_cache_fsl.RData")
  gpa_to_save <- cached_gpa
  gpa <- gpa_to_save
  save(gpa, file = cache_file)
  gpa <- cached_gpa
  gpa$l1_model_setup <- NULL
  gpa$l2_model_setup <- NULL
  gpa$l3_model_setup <- NULL

  out <- lookup_feat_outputs(
    gpa, level = 2L, what = "zstat", cache_dir = dirname(cache_dir)
  )$l2

  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$lookup_source), "cache")
  expect_true(any(grepl("run_pipeline_cache_fsl\\.RData$", out$cache_file)))
})

test_that("lookup_feat_outputs can crawl L1 FEAT folders without setup tables", {
  root <- tempfile("lookup_feat_outputs_")
  l1_feat <- file.path(root, "feat_l1", "sub-01", "ses-1", "facehouse", "FEAT_LVL1_run1.feat")
  make_lookup_feat_stats(l1_feat, contrast_names = c("EV_face", "EV_house"))

  gpa <- make_lookup_feat_outputs_gpa(root, l1_feat = l1_feat)
  gpa$l1_model_setup <- NULL
  gpa$l2_model_setup <- NULL
  gpa$l3_model_setup <- NULL

  out <- lookup_feat_outputs(
    gpa, level = 1L, what = "zstat", source = "filesystem"
  )$l1

  expect_equal(nrow(out), 2L)
  expect_equal(unique(out$lookup_source), "filesystem")
  expect_equal(sort(out$l1_cope_name), c("EV_face", "EV_house"))
  expect_true(all(out$image_exists))
})

test_that("filesystem lookup scans requested statistics independently", {
  root <- tempfile("lookup_feat_outputs_")
  l1_feat <- file.path(root, "feat_l1", "sub-01", "ses-1", "facehouse", "FEAT_LVL1_run1.feat")
  stats_dir <- file.path(l1_feat, "stats")
  dir.create(stats_dir, recursive = TRUE)
  file.create(file.path(stats_dir, "cope1.nii.gz"))
  file.create(file.path(stats_dir, "zstat2.nii.gz"))
  writeLines(
    c(
      "/ContrastName1\tface-house interaction",
      "/ContrastName2\tcontrol",
      "/NumContrasts\t2"
    ),
    file.path(l1_feat, "design.con")
  )

  gpa <- make_lookup_feat_outputs_gpa(root, l1_feat = l1_feat)
  gpa$l1_model_setup <- NULL
  gpa$l2_model_setup <- NULL
  gpa$l3_model_setup <- NULL

  out <- lookup_feat_outputs(
    gpa, level = 1L, what = "cope",
    source = "filesystem", include_missing = TRUE
  )$l1

  expect_equal(nrow(out), 2L)
  expect_equal(out$l1_cope_name, c("face-house interaction", "control"))
  expect_equal(out$image_exists, c(TRUE, FALSE))
  expect_equal(
    out$image_file,
    file.path(stats_dir, paste0("cope", 1:2, ".nii.gz"))
  )

  existing <- lookup_feat_outputs(
    gpa, level = 1L, what = "cope",
    source = "filesystem", include_missing = FALSE
  )$l1
  expect_equal(nrow(existing), 1L)
  expect_equal(existing$l1_cope_name, "face-house interaction")

  zstat <- lookup_feat_outputs(
    gpa, level = 1L, what = "zstat",
    source = "filesystem", include_missing = FALSE
  )$l1
  expect_equal(nrow(zstat), 1L)
  expect_equal(zstat$l1_cope_name, "control")
  expect_equal(zstat$image_file, file.path(stats_dir, "zstat2.nii.gz"))
})

test_that("L2 passthroughs are explicit logical L2 rows", {
  root <- tempfile("lookup_feat_outputs_")
  l1_feat <- file.path(
    root, "feat_l1", "sub-01", "ses-1", "facehouse",
    "FEAT_LVL1_run1.feat"
  )
  passthrough_file <- file.path(l1_feat, "stats", "cope1.nii.gz")
  dir.create(dirname(passthrough_file), recursive = TRUE)
  file.create(passthrough_file)

  gpa <- make_lookup_feat_outputs_gpa(root, l1_feat = l1_feat)
  pass_setup <- data.frame(
    id = "01",
    session = 0L,
    l1_model = "facehouse",
    l1_cope_number = 1L,
    l1_cope_name = "EV_face",
    l2_model = "l2_session",
    l2_scope = "id",
    l2_input_mode = "l1_cope_file_passthrough",
    l2_passthrough = TRUE,
    n_l2_copes = 1L,
    n_input_files = 1L,
    passthrough_cope_file = passthrough_file,
    feat_fsf = NA_character_,
    feat_dir = l1_feat,
    feat_complete = TRUE,
    feat_failed = FALSE,
    feat_dir_exists = TRUE,
    feat_fsf_exists = FALSE,
    stringsAsFactors = FALSE
  )
  pass_setup$cope_list <- list(data.frame(
    id = "01",
    session = 0L,
    l2_cope_number = 1L,
    l2_cope_name = "EV_face",
    stringsAsFactors = FALSE
  ))
  gpa$l2_model_setup <- structure(
    list(fsl = pass_setup),
    class = c("l2_setup", "list")
  )

  result <- lookup_feat_outputs(
    gpa, level = 2L,
    what = c("cope", "varcope", "zstat", "tstat"),
    source = "setup"
  )
  out <- result$l2

  expect_equal(nrow(out), 1L)
  expect_equal(out$output_level, 2L)
  expect_equal(out$materialized_level, 1L)
  expect_true(out$is_passthrough)
  expect_equal(out$l2_input_mode, "l1_cope_file_passthrough")
  expect_equal(out$statistic, "cope")
  expect_equal(out$image_cope_number, 1L)
  expect_equal(out$image_file, passthrough_file)
  expect_equal(out$image_feat_dir, l1_feat)
  expect_true(is.na(out$session))
  expect_true(is.na(out$analysis_fsf))
  expect_true(is.na(out$analysis_dir))
  expect_equal(out$analysis_status, "passthrough")
})

test_that("lookup results have stable empty level schemas", {
  root <- tempfile("lookup_feat_outputs_")
  gpa <- make_lookup_feat_outputs_gpa(root)
  gpa$l1_model_setup <- NULL
  gpa$l2_model_setup <- NULL
  gpa$l3_model_setup <- NULL

  result <- lookup_feat_outputs(
    gpa, level = c(1L, 2L), what = "cope", source = "setup"
  )

  expect_s3_class(result, "feat_output_lookup")
  expect_named(result, c("l1", "l2", "l3"))
  expect_equal(nrow(result$l1), 0L)
  expect_equal(nrow(result$l2), 0L)
  expect_equal(nrow(result$l3), 0L)
  expect_true(all(c("output_level", "l1_model", "image_file") %in% names(result$l1)))
  expect_true(all(c("l2_input_mode", "is_passthrough") %in% names(result$l2)))
  expect_true(all(c("l3_model", "l3_input_mode") %in% names(result$l3)))
  expect_equal(nrow(as.data.frame(result)), 0L)
})

test_that("higher-level estimator fields reach default lookup from refreshed setup and disk", {
  root <- withr::local_tempdir()
  location <- file.path(root, "feat_l3", "L3m-l3_group", "L1m-facehouse",
    "L2m-l2_session", "l2c-overall")
  fit <- make_valid_feat(location, "FEAT_l1c-EV_face", level = 3L)
  syntax <- c(fit$syntax, "set fmri(mixed_yn) 1", "set fmri(robust_yn) 1")
  feat_save_setup(fit$fsf, syntax, character(), 3L)
  writeLines("flameo --ld=stats0000 --runmode=flame12", file.path(fit$fit, ".flame"))
  writeLines(c("ld_dir\tfallback_status", "stats0000\t0"), file.path(fit$fit, "flame_runner_fallbacks.tsv"))
  writeLines("Successful retry", file.path(fit$out, ".feat_auto_retry_warning"))
  feat_validate_outputs(fit$fsf)
  writeLines(rep("2026-01-01 10:00:00", 2), file.path(fit$out, ".feat_complete"))
  gpa <- make_lookup_feat_outputs_gpa(root, l3_feat = fit$out)
  for (source in c("setup", "filesystem")) {
    rows <- lookup_feat_outputs(gpa, level = 3L, what = "cope", source = source,
      refresh_status = source == "setup")$l3
    expect_gt(nrow(rows), 0L)
    expect_identical(unique(rows$feat_requested_estimator), "flame12")
    expect_identical(unique(rows$feat_actual_estimator), "flame1")
    expect_true(all(rows$feat_requested_outlier_deweighting))
    expect_false(any(rows$feat_actual_outlier_deweighting))
    expect_true(all(rows$feat_auto_retried))
    expect_identical(unique(rows$feat_flame1_fallback_count), 1L)
  }
})
