test_that("L2 refits preserve transformations and reference levels within each scope", {
  data <- expand.grid(run_number = 1:6, session = 1:2, id = c("s1", "s2"))
  data$x <- data$run_number * 10 + data$session * 100 + as.integer(data$id) * 1000
  data$condition <- rep(c("a", "b"), nrow(data) / 2)
  original <- data
  transforms <- list(
    mean = function(x) x - mean(x),
    zscore = function(x) as.vector(scale(x)),
    min = function(x) x - min(x),
    max = function(x) x - max(x)
  )

  for (transform in names(transforms)) {
    mobj <- structure(list(
      level = 2L, covariate_transform = c(x = transform),
      reference_level = c(condition = "b"), contrast_spec = list(diagonal = TRUE)
    ), class = c("hi_model_spec", "list"))
    mobj <- fmri.pipeline:::mobj_fit_lm(
      mobj, model_formula = "~ x + condition", data = data,
      id_cols = c("id", "session", "run_number")
    )

    # A missing run changes the transformation parameters used by the refit.
    retained <- data[data$run_number != 6L, ]
    for (scope in c("id_session", "id")) {
      split_on <- if (scope == "id") "id" else c("id", "session")
      result <- fmri.pipeline:::respecify_l2_models_by_subject(mobj, retained, split_on)
      rows <- result$by_subject
      for (i in seq_len(nrow(rows))) {
        subset <- retained[retained$id == rows$id[i], ]
        if (scope == "id_session") subset <- subset[subset$session == rows$session[i], ]
        subset$x <- transforms[[transform]](subset$x)
        subset$condition <- relevel(factor(subset$condition), ref = "b")
        expected <- model.matrix(~ x + condition, subset)
        actual <- rows$model_matrix[[i]]

        expect_equal(unname(actual), unname(expected), ignore_attr = TRUE,
          info = paste(transform, scope, i))
        expect_identical(colnames(actual), colnames(expected))
        expect_identical(colnames(rows$contrasts[[i]]), colnames(expected))
        expect_identical(rows$grouping_scope[i], scope)
      }
    }
  }
  expect_identical(data, original)
})

test_that("L2 centering uses only complete model rows", {
  data <- data.frame(id = "s1", session = 1L, run_number = 1:6,
    x = c(10, 20, 30, 40, 50, 1000), z = c(2, 0, 1, 3, 2, NA_real_))
  mobj <- structure(list(level = 2L, covariate_transform = c(x = "mean"),
    contrast_spec = list(diagonal = TRUE)), class = c("hi_model_spec", "list"))
  mobj <- fmri.pipeline:::mobj_fit_lm(mobj, model_formula = "~ x + z", data = data,
    id_cols = c("id", "session", "run_number"))
  result <- fmri.pipeline:::respecify_l2_models_by_subject(mobj, data)
  expect_equal(unname(result$by_subject$model_matrix[[1]][, "x"]), c(-20, -10, 0, 10, 20))
})
