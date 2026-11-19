# Build a complete backend spec by merging test-specific overrides onto defaults
make_test_backend_specs <- function(...) {
  defaults <- fmri.pipeline:::default_glm_backend_specs()
  overrides <- list(...)
  for (nm in names(overrides)) {
    if (nm %in% names(defaults)) {
      defaults[[nm]] <- modifyList(defaults[[nm]], overrides[[nm]])
    } else {
      defaults[[nm]] <- overrides[[nm]]
    }
  }
  defaults
}

make_mock_logger <- function() {
  list(
    set_threshold = function(...) NULL,
    info = function(...) NULL,
    warn = function(...) NULL,
    debug = function(...) NULL,
    error = function(...) NULL
  )
}

make_mock_job_new <- function() {
  function(...) {
    obj <- list(...)
    obj$depends_on_parents <- NULL
    obj$wait_for_children <- FALSE
    obj$input_rdata_file <- NULL
    obj$output_rdata_file <- NULL
    obj$copy <- function(...) {
      args <- list(...)
      new <- modifyList(obj, args)
      new$copy <- obj$copy
      new
    }
    obj
  }
}

make_mock_sequence_new <- function(store_env) {
  function(joblist, sequence_id) {
    store_env$joblist <- joblist
    store_env$sequence_id <- sequence_id
    list(submit = function(...) NULL)
  }
}

flatten_jobs <- function(joblist) {
  out <- list()
  for (item in joblist) {
    if (is.null(item)) next
    if (is.list(item) && !is.null(item$job_name)) {
      out <- c(out, list(item))
    } else if (is.list(item)) {
      out <- c(out, flatten_jobs(item))
    }
  }
  out
}

make_model_list <- function(l1 = NULL, l2 = NULL, l3 = NULL) {
  list(
    l1_model_names = l1,
    l2_model_names = l2,
    l3_model_names = l3,
    l1_model_name = l1,
    l2_model_name = l2,
    l3_model_name = l3
  )
}

test_that("L1 none + L3 requested waits for split_backend_caches", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- FALSE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))
  gpa$glm_backend_specs <- make_test_backend_specs()

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = NULL, l3_model_names = "l3_model1"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  expect_true(length(jobs) > 0L, info = "No jobs were captured; check backend resolution and mock binding setup.")
  split_job <- jobs[vapply(jobs, function(j) j$job_name == "split_backend_caches", logical(1))]
  l3_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l3_fsl", logical(1))]

  expect_null(result)
  expect_length(split_job, 1L)
  expect_length(l3_job, 1L)
  expect_equal(split_job[[1]]$depends_on_parents, "finalize_configuration")
  expect_equal(l3_job[[1]]$depends_on_parents, "split_backend_caches")
})

test_that("L1 none + L2 requested waits for split_backend_caches", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- TRUE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))
  gpa$glm_backend_specs <- make_test_backend_specs()

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "l2_model1", l3_model_names = "none"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  expect_true(length(jobs) > 0L, info = "No jobs were captured; check backend resolution and mock binding setup.")
  split_job <- jobs[vapply(jobs, function(j) j$job_name == "split_backend_caches", logical(1))]
  l2_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l2_fsl", logical(1))]

  expect_null(result)
  expect_length(split_job, 1L)
  expect_length(l2_job, 1L)
  expect_equal(l2_job[[1]]$depends_on_parents, "split_backend_caches")
})

test_that("FSL L3 waits for backend-specific FSL L2 producer job", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- TRUE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))
  gpa$glm_backend_specs <- make_test_backend_specs()

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = "none",
      l2_model_names = "l2_model1",
      l3_model_names = "l3_model1"
    ),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  l2_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l2_fsl", logical(1))]
  l3_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l3_fsl", logical(1))]

  expect_null(result)
  expect_length(l2_job, 1L)
  expect_length(l3_job, 1L)
  expect_equal(l3_job[[1]]$depends_on_parents, "setup_run_l2_fsl")
})

test_that("L2 requested but cannot run completes without L2 job", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- FALSE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))
  gpa$glm_backend_specs <- make_test_backend_specs()

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "l2_model1", l3_model_names = "none"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  l2_jobs <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l2_fsl", logical(1))]

  expect_null(result)
  expect_length(l2_jobs, 0L)
})

test_that("AFNI L3 multi-run path schedules FSL prerequisites and cache sync", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "afni"
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))

  l1_name <- names(gpa$l1_models$models)[1L]
  l2_name <- names(gpa$l2_models$models)[1L]
  l3_name <- names(gpa$l3_models$models)[1L]

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = l1_name,
      l2_model_names = l2_name,
      l3_model_names = l3_name
    ),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  run_l1_fsl <- jobs[vapply(jobs, function(j) j$job_name == "run_l1_fsl", logical(1))]
  l2_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l2_fsl", logical(1))]
  sync_job <- jobs[vapply(jobs, function(j) j$job_name == "sync_l2_backend_caches", logical(1))]
  l3_afni <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l3_afni", logical(1))]

  expect_null(result)
  expect_length(run_l1_fsl, 1L)
  expect_length(l2_job, 1L)
  expect_length(sync_job, 1L)
  expect_length(l3_afni, 1L)
  expect_equal(l2_job[[1]]$depends_on_parents, "run_l1_fsl")
  expect_equal(sync_job[[1]]$depends_on_parents, "setup_run_l2_fsl")
  expect_equal(l3_afni[[1]]$depends_on_parents, "sync_l2_backend_caches")
  expect_match(l2_job[[1]]$input_rdata_file, "run_pipeline_cache_fsl\\.RData$")
  expect_match(l3_afni[[1]]$input_rdata_file, "run_pipeline_cache_afni\\.RData$")
})

test_that("SPM multi-run L3 path uses L1 concat outputs without standalone L2", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "spm"
  gpa$multi_run <- TRUE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    spm = list(l1_spm_alljobs_time = "0:10:00")
  ))

  l1_name <- names(gpa$l1_models$models)[1L]
  l3_name <- names(gpa$l3_models$models)[1L]

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = l1_name,
      l2_model_names = NULL,
      l3_model_names = l3_name
    ),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  run_l1_spm <- jobs[vapply(jobs, function(j) j$job_name == "run_l1_spm", logical(1))]
  l2_job <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l2", logical(1))]
  l3_spm <- jobs[vapply(jobs, function(j) j$job_name == "setup_run_l3_spm", logical(1))]

  expect_null(result)
  expect_length(run_l1_spm, 1L)
  expect_length(l2_job, 0L)
  expect_length(l3_spm, 1L)
  expect_equal(l3_spm[[1]]$depends_on_parents, "run_l1_spm")
})

test_that("run_glm_pipeline honors per-level backend selection", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- c("fsl", "spm", "afni")
  gpa$level_backends <- list(
    l1 = "fsl",
    l2 = "fsl",
    l3 = "afni"
  )
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))

  l1_name <- names(gpa$l1_models$models)[1L]
  l2_name <- names(gpa$l2_models$models)[1L]
  l3_name <- names(gpa$l3_models$models)[1L]

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = l1_name,
      l2_model_names = l2_name,
      l3_model_names = l3_name
    ),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  job_names <- vapply(jobs, `[[`, character(1), "job_name")

  expect_null(result)
  expect_true("run_l1_fsl" %in% job_names)
  expect_true("setup_run_l3_afni" %in% job_names)
  expect_false("run_l1_spm" %in% job_names)
  expect_false("setup_run_l3_spm" %in% job_names)
  expect_false("setup_run_l3_fsl" %in% job_names)
})

test_that("run_glm_pipeline honors model-specific execution backend overrides", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa$l3_models$models$l3_model1$producer_backend <- "fsl"
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))

  l1_name <- names(gpa$l1_models$models)[1L]
  l2_name <- names(gpa$l2_models$models)[1L]
  l3_name <- names(gpa$l3_models$models)[1L]

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = l1_name,
      l2_model_names = l2_name,
      l3_model_names = l3_name
    ),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  job_names <- vapply(jobs, `[[`, character(1), "job_name")

  expect_null(result)
  expect_true("run_l1_fsl" %in% job_names)
  expect_true("setup_run_l2_fsl" %in% job_names)
  expect_true("setup_run_l3_afni" %in% job_names)
  expect_false("setup_run_l3_fsl" %in% job_names)
})

test_that("run_glm_pipeline emits backend preflight report for resolved mixed backend plan", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa$l3_models$models$l3_model1$producer_backend <- "fsl"
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))

  l1_name <- names(gpa$l1_models$models)[1L]
  l2_name <- names(gpa$l2_models$models)[1L]
  l3_name <- names(gpa$l3_models$models)[1L]

  store <- new.env(parent = emptyenv())
  report_store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = l1_name,
      l2_model_names = l2_name,
      l3_model_names = l3_name
    ),
    log_backend_preflight_report = function(report, lg) {
      report_store$report <- report
      invisible(report)
    },
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  expect_null(result)
  expect_true(is.data.frame(report_store$report))
  expect_equal(nrow(report_store$report), 3L)
  expect_identical(report_store$report$level, c(1L, 2L, 3L))
  expect_identical(report_store$report$execution_backend, c("fsl", "fsl", "afni"))
  expect_identical(report_store$report$producer_backend[3L], "fsl")
  expect_identical(report_store$report$producer_level[3L], 2L)
})

test_that("No L3 requested completes without L3 job", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- TRUE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))
  gpa$glm_backend_specs <- make_test_backend_specs()

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "l2_model1", l3_model_names = "none"),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  l3_jobs <- jobs[vapply(jobs, function(j) grepl("^setup_run_l3_", j$job_name), logical(1))]

  expect_null(result)
  expect_length(l3_jobs, 0L)
})

test_that("reassign_unrunnable_models reassigns AFNI-assigned L1 model to FSL when AFNI L3 needs FSL copes", {
  # AFNI can't run L1. When L3 is AFNI (which requires FSL L2 copes),
  # the pipeline should reassign L1 from AFNI to FSL as a fallback.
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- c("fsl", "afni")
  # Explicitly assign L1 to AFNI (which can't run L1) and L3 to AFNI
  gpa$level_backends <- list(l1 = "afni", l2 = "afni", l3 = "afni")
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))

  l1_name <- names(gpa$l1_models$models)[1L]
  l2_name <- names(gpa$l2_models$models)[1L]
  l3_name <- names(gpa$l3_models$models)[1L]

  store <- new.env(parent = emptyenv())
  result <- testthat::with_mocked_bindings(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = l1_name,
      l2_model_names = l2_name,
      l3_model_names = l3_name
    ),
    R_batch_job = list(new = make_mock_job_new()),
    R_batch_sequence = list(new = make_mock_sequence_new(store)),
    .package = "fmri.pipeline"
  )

  jobs <- flatten_jobs(store$joblist)
  job_names <- vapply(jobs, `[[`, character(1), "job_name")

  expect_null(result)
  # L1 should be reassigned to FSL (the fallback from L3 dependency resolution)
  expect_true("run_l1_fsl" %in% job_names)
  expect_false("run_l1_afni" %in% job_names)
  # L2 should also be reassigned to FSL
  expect_true("setup_run_l2_fsl" %in% job_names)
  # L3 remains AFNI
  expect_true("setup_run_l3_afni" %in% job_names)
  expect_false("setup_run_l3_fsl" %in% job_names)
})

test_that("run_glm_pipeline rejects explicit FSL execution for 3dlmer models", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- c("fsl", "afni")
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$l3_models$models$l3_model1$execution_backend <- "fsl"
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))

  expect_error(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = "model1",
      l2_model_names = "l2_model1",
      l3_model_names = "l3_model1"
    ),
    "must execute only with backend 'afni'"
  )
})

test_that("run_glm_pipeline rejects AFNI execution for non-3dlmer models", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- c("fsl", "afni")
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "per_session"
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00",
    l1_setup_time = "0:10:00",
    l1_setup_memgb = "1G",
    l2_setup_cores = 1L,
    l2_setup_run_time = "0:10:00",
    l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")
  ))

  expect_error(
    fmri.pipeline:::run_glm_pipeline(
      gpa,
      l1_model_names = "model1",
      l2_model_names = "l2_model1",
      l3_model_names = "l3_model1"
    ),
    "AFNI L3 execution currently supports only l3_input_mode='3dlmer'"
  )
})

test_that("configured inference waits for collected L3 outputs and supports opt-out", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- FALSE
  gpa$output_locations$sqlite_db <- NULL
  gpa$parallel <- modifyList(gpa$parallel, list(
    finalize_time = "0:10:00", l1_setup_time = "0:10:00", l1_setup_memgb = "1G",
    l2_setup_cores = 1L, l2_setup_run_time = "0:10:00", l3_setup_run_time = "0:10:00",
    fsl = list(l1_feat_alljobs_time = "0:10:00")))
  gpa$glm_backend_specs <- make_test_backend_specs()
  gpa$l3_models$models$l3_model1$inference <- list(method = "ptfce")
  gpa$batch_run <- list(fwe_specs = list(old_correction = list()))
  store <- new.env(parent = emptyenv())
  invoke <- function(enabled, l3 = "l3_model1") {
    testthat::with_mocked_bindings(
      run_glm_pipeline(gpa, l1_model_names = "none", l2_model_names = "none", l3_model_names = l3, run_fwe = enabled),
      R_batch_job = list(new = make_mock_job_new()),
      R_batch_sequence = list(new = make_mock_sequence_new(store)),
      .package = "fmri.pipeline")
    flatten_jobs(store$joblist)
  }
  jobs <- invoke(TRUE)
  correction <- Filter(function(x) x$job_name == "run_fwe", jobs)
  expect_length(correction, 1L)
  expect_identical(correction[[1L]]$depends_on_parents, "cleanup_glm")
  expect_true(correction[[1L]]$wait_for_children)
  expect_false(correction[[1L]]$all_children_successful) # persist status before reporting a child failure
  expect_match(correction[[1L]]$post_children_r_code[1L], "require_complete = TRUE", fixed = TRUE)
  expect_match(correction[[1L]]$post_children_r_code[2L], "!isTRUE(success)", fixed = TRUE)
  cache <- file.path(jobs[[1L]]$batch_directory, "run_pipeline_cache.RData")
  env <- new.env()
  load(cache, env)
  expect_named(env$gpa$batch_run$fwe_specs, "l3_model1_ptfce")
  expect_length(Filter(function(x) x$job_name == "run_fwe", invoke(FALSE)), 0L)
  expect_length(Filter(function(x) x$job_name == "run_fwe", invoke(TRUE, "none")), 0L)
})
