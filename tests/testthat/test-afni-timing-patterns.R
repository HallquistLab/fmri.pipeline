test_that("AFNI dmBLOCK groups only identical event schedules across runs", {
  dir <- withr::local_tempdir()
  event <- function(onset, value = rep(1, length(onset))) {
    cbind(onset = onset, duration = rep(1, length(onset)), value = value)
  }
  design <- array(list(event(c(2, 5)), event(3),
                       event(c(2, 5), c(8, 9)), event(3, 7),
                       event(2), event(numeric()),
                       event(c(2, 6)), event(3)),
                  dim = c(2, 4), dimnames = list(NULL, c("task", "mod", "short", "shift")))
  fmri.pipeline:::write_afni_timing_files(design, dir)
  expect_setequal(list.files(dir), paste0(c("task_mod", "short", "shift"), "_dmBLOCK.txt"))
  expect_identical(readLines(file.path(dir, "task_mod_dmBLOCK.txt")), c("2*8:1 5*9:1", "3*7:1"))
  expect_identical(readLines(file.path(dir, "short_dmBLOCK.txt")), c("2:1", "*"))
})
