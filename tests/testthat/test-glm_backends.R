# Tests for GLM backend registry helpers

test_that("get_glm_backends returns resolved backends from gpa", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("fsl", "spm")
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  backends <- fmri.pipeline:::get_glm_backends(gpa)

  expect_type(backends, "list")
  expect_true(all(c("fsl", "spm") %in% names(backends)))
  expect_true(is.function(backends$fsl$l1_setup))
  expect_true(is.function(backends$spm$l1_setup))
})

test_that("producer requirement helpers capture AFNI and SPM L1 producers", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("afni", "spm", "fsl")
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  backends <- fmri.pipeline:::get_glm_backends(gpa)

  afni_req <- fmri.pipeline:::resolve_l3_producer_requirements(
    backends = backends,
    l3_backend = backends$afni,
    multi_run = TRUE
  )
  spm_req <- fmri.pipeline:::resolve_l3_producer_requirements(
    backends = backends,
    l3_backend = backends$spm,
    multi_run = TRUE
  )

  expect_equal(nrow(afni_req), 1L)
  expect_identical(afni_req$producer_backend, "afni")
  expect_identical(afni_req$producer_level, 1L)
  expect_equal(nrow(spm_req), 1L)
  expect_identical(spm_req$producer_backend, "spm")
  expect_identical(spm_req$producer_level, 1L)
  expect_identical(
    fmri.pipeline:::get_requirement_producer_backends(afni_req),
    "afni"
  )
  expect_identical(
    fmri.pipeline:::get_requirement_producer_backends(afni_req, producer_level = 1L),
    "afni"
  )
})

test_that("backend capability metadata describes per-level execution strategy", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("fsl", "spm", "afni")
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  backends <- fmri.pipeline:::get_glm_backends(gpa)

  expect_true(fmri.pipeline:::backend_runs_level(backends$fsl, 1L))
  expect_true(fmri.pipeline:::backend_runs_level(backends$fsl, 2L))
  expect_true(fmri.pipeline:::backend_runs_level(backends$fsl, 3L))

  expect_true(fmri.pipeline:::backend_runs_level(backends$spm, 1L))
  expect_false(fmri.pipeline:::backend_runs_level(backends$spm, 2L))
  expect_true(fmri.pipeline:::backend_runs_level(backends$spm, 3L))

  expect_true(fmri.pipeline:::backend_runs_level(backends$afni, 1L))
  expect_false(fmri.pipeline:::backend_runs_level(backends$afni, 2L))
  expect_true(fmri.pipeline:::backend_runs_level(backends$afni, 3L))

  expect_identical(fmri.pipeline:::backend_multi_run_strategy(backends$fsl), "explicit_l2")
  expect_identical(fmri.pipeline:::backend_multi_run_strategy(backends$spm), "concat_in_l1")
  expect_identical(fmri.pipeline:::backend_multi_run_strategy(backends$afni), "concat_in_l1")
})

test_that("artifact helpers expose current producers and AFNI 3dLMEr inputs", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("fsl", "spm", "afni")
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  backends <- fmri.pipeline:::get_glm_backends(gpa)

  expect_true(fmri.pipeline:::backend_produces_artifact(backends$fsl, "subject_session_contrasts"))
  expect_true(fmri.pipeline:::backend_produces_artifact(backends$spm, "subject_session_contrasts"))
  expect_true(fmri.pipeline:::backend_produces_artifact(backends$afni, "subject_session_contrasts"))

  expect_identical(
    fmri.pipeline:::backend_l3_required_artifacts(backends$afni),
    "subject_session_contrasts"
  )
  expect_identical(
    fmri.pipeline:::backend_l3_input_provider_backends(backends$afni),
    "afni"
  )
  expect_identical(
    fmri.pipeline:::get_backends_producing_artifact(backends, "subject_session_contrasts"),
    c("fsl", "spm", "afni")
  )
  expect_identical(
    fmri.pipeline:::resolve_l3_input_producers(backends, backends$afni),
    "afni"
  )
})

test_that("artifact helpers resolve producing level from backend strategy", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("fsl", "spm", "afni")
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)
  backends <- fmri.pipeline:::get_glm_backends(gpa)

  expect_identical(
    fmri.pipeline:::backend_artifact_production_level(
      backends$fsl,
      artifact = "subject_session_contrasts",
      multi_run = TRUE
    ),
    2L
  )
  expect_identical(
    fmri.pipeline:::backend_artifact_production_level(
      backends$fsl,
      artifact = "subject_session_contrasts",
      multi_run = FALSE
    ),
    1L
  )
  expect_identical(
    fmri.pipeline:::backend_artifact_production_level(
      backends$spm,
      artifact = "subject_session_contrasts",
      multi_run = TRUE
    ),
    1L
  )
})

test_that("model-level L3 requirements resolve producer stage from execution and producer backends", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa$l3_models$models$l3_model1$producer_backend <- "fsl"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  req <- fmri.pipeline:::resolve_model_l3_requirements(
    gpa,
    l3_model_names = "l3_model1"
  )

  expect_equal(nrow(req), 1L)
  expect_identical(req$model_name, "l3_model1")
  expect_identical(req$execution_backend, "afni")
  expect_identical(req$producer_backend, "fsl")
  expect_identical(req$artifact, "subject_session_contrasts")
  expect_identical(req$producer_level, 2L)
})

test_that("L3 backend fields round-trip through YAML model building", {
  gpa <- create_mock_gpa(n_subjects = 2, n_runs = 1, include_l2_models = TRUE)
  gpa$glm_software <- c("fsl", "afni")
  gpa$level_backends <- list(l1 = "fsl", l2 = "fsl", l3 = c("fsl", "afni"))
  gpa$subject_data <- expand.grid(
    id = paste0("sub", 1:2),
    session = 1:3,
    stringsAsFactors = FALSE
  )
  gpa$subject_data$session_label <- rep(c("baseline", "month1", "month2"), times = 2)
  gpa$l2_models$models <- list(
    l2_id_session = list(
      name = "l2_id_session",
      level = 2L,
      l2_scope = "id_session"
    )
  )
  class(gpa$l2_models) <- c("list", "hi_model_set")

  spec <- list(
    l3_models = list(
      l3_pooled_subject_ev = list(
        level = 3L,
        l3_input_mode = "pooled_sessions_subject_ev",
        execution_backend = "fsl",
        producer_backend = "fsl",
        model_formula = "~ session_label",
        reference_level = list(session_label = "baseline"),
        contrasts = list(diagonal = TRUE)
      ),
      l3_3dlmer_session = list(
        level = 3L,
        l3_input_mode = "3dlmer",
        execution_backend = "afni",
        producer_backend = "fsl",
        model_formula = "~ session_label",
        random_effects = "(1 | Subj)",
        reference_level = list(session_label = "baseline"),
        contrasts = list(diagonal = TRUE)
      )
    )
  )
  spec_file <- tempfile(fileext = ".yaml")
  yaml::write_yaml(spec, spec_file)

  gpa <- build_l3_models(gpa, from_spec_file = spec_file)

  exec_map <- fmri.pipeline:::get_effective_model_backends(
    gpa,
    level = 3L,
    model_names = names(spec$l3_models)
  )
  producer_map <- fmri.pipeline:::get_effective_model_backends(
    gpa,
    level = 3L,
    model_names = names(spec$l3_models),
    type = "producer"
  )

  expect_identical(exec_map$l3_pooled_subject_ev, "fsl")
  expect_identical(producer_map$l3_pooled_subject_ev, "fsl")
  expect_identical(exec_map$l3_3dlmer_session, "afni")
  expect_identical(producer_map$l3_3dlmer_session, "fsl")
})

test_that("backend preflight report summarizes execution and producer resolution", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa$l3_models$models$l3_model1$producer_backend <- "fsl"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  l1_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 1L, model_names = "model1")
  l2_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 2L, model_names = "l2_model1")
  l3_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1")
  req <- fmri.pipeline:::resolve_model_l3_requirements(gpa, l3_model_names = "l3_model1")

  report <- fmri.pipeline:::build_backend_preflight_report(
    gpa = gpa,
    l1_model_backend_map = l1_map,
    l2_model_backend_map = l2_map,
    l3_model_backend_map = l3_map,
    l3_requirement_df = req
  )

  expect_equal(nrow(report), 3L)
  expect_identical(report$level, c(1L, 2L, 3L))
  expect_identical(report$execution_backend, c("fsl", "fsl", "afni"))
  expect_true(is.na(report$producer_backend[1L]))
  expect_true(is.na(report$producer_backend[2L]))
  expect_identical(report$producer_backend[3L], "fsl")
  expect_identical(report$producer_level[3L], 2L)
  expect_identical(report$producer_artifact[3L], "subject_session_contrasts")
})

test_that("level-specific backend resolution honors level_backends over legacy glm_software", {
  gpa <- create_mock_gpa(include_l1_models = TRUE)
  gpa$glm_software <- c("fsl", "spm", "afni")
  gpa$level_backends <- list(
    l1 = "fsl",
    l2 = "fsl",
    l3 = "afni"
  )
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  expect_identical(fmri.pipeline:::get_level_backend_names(gpa, level = 1L), "fsl")
  expect_identical(fmri.pipeline:::get_level_backend_names(gpa, level = 2L), "fsl")
  expect_identical(fmri.pipeline:::get_level_backend_names(gpa, level = 3L), "afni")
  expect_identical(
    fmri.pipeline:::get_level_backend_names(gpa),
    c("fsl", "afni")
  )

  l1_backends <- fmri.pipeline:::get_glm_backends(gpa, level = 1L)
  l3_backends <- fmri.pipeline:::get_glm_backends(gpa, level = 3L)
  expect_identical(names(l1_backends), "fsl")
  expect_identical(names(l3_backends), "afni")
})

test_that("model-specific execution overrides take precedence over level defaults", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- "fsl"
  gpa$level_backends <- list(l1 = "fsl", l2 = "fsl", l3 = "fsl")
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  execution_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1")
  expect_identical(execution_map$l3_model1, "afni")

  override_map <- fmri.pipeline:::get_effective_model_backends(
    gpa,
    level = 3L,
    model_names = "l3_model1",
    backend_overrides = list(l3 = c(l3_model1 = "spm"))
  )
  expect_identical(override_map$l3_model1, "spm")
})

test_that("3dlmer models default to AFNI execution and FSL producer resolution", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$glm_software <- c("fsl", "afni")
  gpa$level_backends <- list(l1 = "fsl", l2 = "fsl", l3 = "fsl")
  gpa$multi_run <- TRUE
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  execution_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1")
  producer_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1", type = "producer")

  expect_identical(execution_map$l3_model1, "afni")
  expect_identical(producer_map$l3_model1, "fsl")
})

test_that("validate_l3_backend_resolution rejects invalid 3dlmer execution backend", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$l3_models$models$l3_model1$execution_backend <- "fsl"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  execution_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1")
  producer_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1", type = "producer")

  expect_error(
    fmri.pipeline:::validate_l3_backend_resolution(
      gpa = gpa,
      l3_model_names = "l3_model1",
      execution_backend_map = execution_map,
      producer_backend_map = producer_map
    ),
    "must execute only with backend 'afni'"
  )
})

test_that("validate_l3_backend_resolution accepts ordinary AFNI GLMs", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$l3_models$models$l3_model1$l3_input_mode <- "per_session"
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  execution_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1")
  producer_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1", type = "producer")

  expect_no_error(
    fmri.pipeline:::validate_l3_backend_resolution(
      gpa = gpa,
      l3_model_names = "l3_model1",
      execution_backend_map = execution_map,
      producer_backend_map = producer_map
    )
  )
})

test_that("validate_l3_backend_resolution rejects unsupported producer backends for 3dlmer", {
  gpa <- create_mock_gpa(include_l1_models = TRUE, include_l2_models = TRUE, include_l3_models = TRUE)
  gpa$l3_models$models$l3_model1$l3_input_mode <- "3dlmer"
  gpa$l3_models$models$l3_model1$execution_backend <- "afni"
  gpa$l3_models$models$l3_model1$producer_backend <- "spm"
  gpa <- fmri.pipeline:::initialize_glm_backends(gpa)

  execution_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1")
  producer_map <- fmri.pipeline:::get_effective_model_backends(gpa, level = 3L, model_names = "l3_model1", type = "producer")

  expect_error(
    fmri.pipeline:::validate_l3_backend_resolution(
      gpa = gpa,
      l3_model_names = "l3_model1",
      execution_backend_map = execution_map,
      producer_backend_map = producer_map
    ),
    "supports producer_backend='fsl' or 'afni'"
  )
})

test_that("normalize_backend_override_config handles explicit execution/producer shape", {
  # Explicit shape: list(execution = list(l3 = ...), producer = list(l3 = ...))
  cfg <- fmri.pipeline:::normalize_backend_override_config(list(
    execution = list(l3 = c(l3_model1 = "afni")),
    producer = list(l3 = c(l3_model1 = "fsl"))
  ))

  expect_identical(cfg$execution$l3$l3_model1, "afni")
  expect_identical(cfg$producer$l3$l3_model1, "fsl")
  # Shorthand levels should remain empty
expected_empty <- list()
  expect_identical(cfg$execution$l1, expected_empty)
  expect_identical(cfg$execution$l2, expected_empty)
  expect_identical(cfg$producer$l1, expected_empty)
  expect_identical(cfg$producer$l2, expected_empty)
})

test_that("normalize_backend_override_config handles mixed shorthand and explicit shapes", {
  # Shorthand at l1, explicit execution+producer at l3
  cfg <- fmri.pipeline:::normalize_backend_override_config(list(
    l1 = c(model1 = "spm"),
    execution = list(l3 = c(l3_model1 = "afni")),
    producer = list(l3 = c(l3_model1 = "fsl"))
  ))

  # Shorthand goes to execution
  expect_identical(cfg$execution$l1$model1, "spm")
  # Explicit goes to their respective slots
  expect_identical(cfg$execution$l3$l3_model1, "afni")
  expect_identical(cfg$producer$l3$l3_model1, "fsl")
})

test_that("normalize_backend_override_config handles per-level execution/producer sub-lists", {
  # Per-level shape: list(l3 = list(execution = ..., producer = ...))
  cfg <- fmri.pipeline:::normalize_backend_override_config(list(
    l3 = list(
      execution = c(l3_model1 = "afni"),
      producer = c(l3_model1 = "fsl")
    )
  ))

  expect_identical(cfg$execution$l3$l3_model1, "afni")
  expect_identical(cfg$producer$l3$l3_model1, "fsl")
})
