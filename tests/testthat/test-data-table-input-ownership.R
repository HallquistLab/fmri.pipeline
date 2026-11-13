test_that("fill_atlas_with_stats does not convert its input by reference", {
  atlas_file <- tempfile(fileext = ".nii.gz")
  file.create(atlas_file)
  stat_df <- data.frame(atlas_value = 1:2, t = c(0.5, 1.5))
  stat_before <- stat_df

  expect_message(
    expect_error(
      fill_atlas_with_stats(
        atlas_file,
        stat_df,
        stat_cols = "missing_stat"
      ),
      "missing_stat"
    ),
    "Coercing stat_dt"
  )

  expect_identical(stat_df, stat_before)
  expect_false(data.table::is.data.table(stat_df))
})
