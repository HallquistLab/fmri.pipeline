test_that("custom quantitative AFNI contrasts preserve unit slopes and explicit offsets", {
  dt <- data.frame(age = 21:28, score = c(2, 5, 1, 4, 3, 7, 6, 8), y = c(2, 3, 1, 5, 4, 6, 8, 7))
  mobj <- structure(list(
    lmfit = stats::lm(y ~ age + score, data = dt),
    contrast_spec = list(diagonal = FALSE, custom = list(
      age_slope = c(age = 1, score = 0), score_slope = c(age = 0, score = 1)
    ))
  ), class = "hi_model_spec")

  # Infer quantitative variables from the fit as well as from explicit qVars.
  for (qvars in list(NULL, c("age", "score"))) {
    glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, dt, qVars = qvars)
    expect_equal(glts$label_raw, c("age_slope", "score_slope"))
    expect_equal(glts$code, c("age :", "score :"))
    expect_true(all(glts$valid))
  }

  # A one-column matrix must retain its coefficient name when its row drops dimensions.
  mobj$contrast_spec$custom <- matrix(1, dimnames = list("age_slope", "age"))
  glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, dt)
  expect_equal(glts$code, "age :")
  expect_equal(glts$label_raw, "age_slope")

  # Raw GLTs deliberately use values after ':' as covariate evaluation offsets.
  mobj$contrast_spec <- list(diagonal = FALSE)
  raw_codes <- list(below_center = "age : -1", above_center = "age : 2")
  glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, dt, raw_glt_codes = raw_codes)
  expect_equal(glts$code, unname(unlist(raw_codes)))
  expect_true(all(glts$source == "user"))
  expect_true(all(glts$valid))
})

test_that("custom quantitative AFNI contrasts reject unsupported weights and combinations", {
  dt <- data.frame(age = 21:28, score = c(2, 5, 1, 4, 3, 7, 6, 8), y = c(2, 3, 1, 5, 4, 6, 8, 7))
  mobj <- structure(list(
    lmfit = stats::lm(y ~ age + score, data = dt),
    contrast_spec = list(diagonal = FALSE)
  ), class = "hi_model_spec")
  unsupported <- list(
    negative_age = c(age = -1),
    doubled_age = c(age = 2),
    half_age = c(age = 0.5),
    negative_half_age = c(age = -0.5),
    slope_sum = c(age = 1, score = 1),
    slope_difference = c(age = 1, score = -1)
  )

  for (name in names(unsupported)) {
    custom <- stats::setNames(list(unsupported[[name]]), name)
    custom_matrix <- matrix(unsupported[[name]], nrow = 1L,
                            dimnames = list(name, names(unsupported[[name]])))
    for (input in list(custom, custom_matrix)) {
      # Both current and legacy custom-contrast storage must use the same guard.
      for (legacy in c(FALSE, TRUE)) {
        mobj$contrast_spec$custom <- if (legacy) NULL else input
        mobj$contrast_list <- if (legacy) list(custom = input) else NULL
        expect_error(
          fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, dt, qVars = c("age", "score")),
          paste0("custom contrast '", name, "'.*single quantitative coefficient.*weight \\+1.*covariate offsets")
        )
      }
    }
  }
})
