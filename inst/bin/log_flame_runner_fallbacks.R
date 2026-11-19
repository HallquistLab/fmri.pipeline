#!/usr/bin/env Rscript

print_help <- function() {
  cat(paste(
    "Log FLAME12-to-FLAME1 slice fallback records produced by flame_runner.",
    "Options:",
    "  --feat_dir <path>: FEAT output directory containing flame_runner_fallbacks.tsv.",
    "  --level <n>: FEAT level, used in the logger name.",
    "  --log_txt <path|NULL>: Text lgr log file.",
    "  --log_json <path|NULL>: JSON lgr log file.",
    "  --threshold <level>: lgr threshold.",
    "  --help: Print this help menu.",
    sep = "\n"
  ))
}

args <- commandArgs(trailingOnly = TRUE)

feat_dir <- NULL
level <- NA_integer_
log_txt <- NULL
log_json <- NULL
threshold <- "info"

argpos <- 1L
while (argpos <= length(args)) {
  if (args[[argpos]] == "--feat_dir") {
    feat_dir <- args[[argpos + 1L]]
    argpos <- argpos + 2L
  } else if (args[[argpos]] == "--level") {
    level <- as.integer(args[[argpos + 1L]])
    argpos <- argpos + 2L
  } else if (args[[argpos]] == "--log_txt") {
    log_txt <- args[[argpos + 1L]]
    argpos <- argpos + 2L
  } else if (args[[argpos]] == "--log_json") {
    log_json <- args[[argpos + 1L]]
    argpos <- argpos + 2L
  } else if (args[[argpos]] == "--threshold") {
    threshold <- args[[argpos + 1L]]
    argpos <- argpos + 2L
  } else if (args[[argpos]] == "--help") {
    print_help()
    quit(save = "no", status = 0L)
  } else {
    stop("Unknown argument: ", args[[argpos]])
  }
}

nullish <- function(x) {
  is.null(x) || !nzchar(x) || identical(x, "NULL") || identical(x, "NA")
}

if (nullish(feat_dir)) {
  stop("--feat_dir is required")
}

fit_dir <- if (endsWith(feat_dir, ".gfeat")) file.path(feat_dir, "cope1.feat") else feat_dir
fallback_file <- file.path(fit_dir, "flame_runner_fallbacks.tsv")
if (!file.exists(fallback_file)) {
  quit(save = "no", status = 0L)
}

fallbacks <- tryCatch(
  utils::read.delim(fallback_file, stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) {
    warning("Could not read flame_runner fallback log: ", fallback_file, ": ", conditionMessage(e))
    NULL
  }
)

if (is.null(fallbacks) || nrow(fallbacks) == 0L) {
  quit(save = "no", status = 0L)
}

if (!requireNamespace("lgr", quietly = TRUE)) {
  warning("Package 'lgr' is required to log flame_runner fallback records.")
  quit(save = "no", status = 0L)
}

logger_name <- if (!is.na(level)) {
  paste0("glm_pipeline/l", level, "_estimation")
} else {
  "glm_pipeline/fsl_estimation"
}

lg <- lgr::get_logger(logger_name)
lg$set_threshold(threshold)

if (!nullish(log_txt)) {
  dir.create(dirname(log_txt), recursive = TRUE, showWarnings = FALSE)
  if (!"flame_runner_fallback_txt" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderFile$new(log_txt), name = "flame_runner_fallback_txt")
  }
}

if (!nullish(log_json)) {
  dir.create(dirname(log_json), recursive = TRUE, showWarnings = FALSE)
  if (!"flame_runner_fallback_json" %in% names(lg$appenders)) {
    lg$add_appender(lgr::AppenderJson$new(log_json), name = "flame_runner_fallback_json")
  }
}

for (ii in seq_len(nrow(fallbacks))) {
  rec <- fallbacks[ii, , drop = FALSE]
  outcome <- if (is.null(rec$fallback_status)) {
    "rerunning it with FLAME1 (legacy record, outcome unknown)."
  } else if (rec$fallback_status == 0L) {
    "recovered by rerunning it with FLAME1."
  } else {
    sprintf("FLAME1 retry failed with status %s.", rec$fallback_status)
  }
  lg$warn(
    paste(
      "FSL FLAME12 failed for slice output '%s' in '%s' with status %s;",
      "%s",
      "Original failed log directory: '%s'. Fallback command: %s"
    ),
    rec$ld_dir,
    feat_dir,
    rec$flame12_status,
    outcome,
    rec$failed_ld_dir,
    rec$fallback_command
  )
}
