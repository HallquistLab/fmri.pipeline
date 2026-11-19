#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L || args[1L] != "--request" || args[3L] != "--worker") {
  stop("Usage: fwe_task.R --request task.rds --worker worker.R [worker arguments]")
}
Sys.unsetenv("R_TESTS")
tryCatch(fmri.pipeline:::execute_fwe_task(args[2L], args[4L], args[-seq_len(4L)]),
  fwe_worker_error = function(e) {
    message(conditionMessage(e))
    quit(save = "no", status = e$status)
  })
