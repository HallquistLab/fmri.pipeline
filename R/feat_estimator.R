# Higher-level FEAT run modes, using the names accepted by flameo.
feat_requested_estimator <- function(syntax = character()) {
  modes <- c("0" = "ols", "1" = "flame12", "2" = "flame1", "3" = "fe", "4" = "randomise")
  mode <- feat_fsf_value(syntax, "fmri(mixed_yn)") %||% NA_character_
  robust <- feat_fsf_value(syntax, "fmri(robust_yn)") %||% NA_character_
  data.frame(
    feat_requested_estimator = unname(modes[mode]),
    feat_actual_estimator = NA_character_,
    feat_requested_outlier_deweighting = unname(c("0" = FALSE, "1" = TRUE)[robust]),
    feat_actual_outlier_deweighting = NA,
    feat_flame1_fallback_count = NA_integer_
  )
}

# Read only the retained attempt. Its .flame contains the actual FLAME options
# (including --io), even when the outlier retry used a different FSF.
feat_estimator_result <- function(setup) {
  result <- feat_requested_estimator(setup$syntax)
  if (is.null(setup)) return(result)
  fit <- file.path(setup$output_dir, "cope1.feat")
  plan <- file.path(fit, ".flame")
  if (!file.exists(plan)) return(result)
  commands <- trimws(readLines(plan, warn = FALSE))
  commands <- commands[nzchar(commands)]
  if (!length(commands)) return(result)
  if (all(grepl("(^|[/ ])randomise[[:space:]]", commands))) {
    result$feat_flame1_fallback_count <- 0L
    result$feat_actual_estimator <- "randomise"
    result$feat_actual_outlier_deweighting <- FALSE
    return(result)
  }

  mode_pattern <- ".*[[:space:]]--runmode=(ols|fe|flame1|flame12)([[:space:]].*|$)"
  if (!all(grepl(mode_pattern, commands))) return(result)
  modes <- sub(mode_pattern, "\\1", commands)
  robust <- unique(grepl("[[:space:]]--io([[:space:]]|$)", commands))
  if (length(robust) == 1L) result$feat_actual_outlier_deweighting <- robust
  if (setup$level == 3L && any(modes == "flame12")) {
    # Stock FEAT at L2 does not run the pipeline's slice retry wrapper.
    audit <- file.path(fit, "flame_runner_fallbacks.tsv")
    if (!file.exists(audit)) return(result)
    records <- tryCatch(utils::read.delim(audit, quote = "", stringsAsFactors = FALSE),
      error = function(e) NULL)
    if (!"ld_dir" %in% names(records)) return(result)
    if ("fallback_status" %in% names(records)) {
      if (anyNA(records$fallback_status) || any(records$fallback_status != 0L)) return(result)
    } else if (!all(c("flame12_status", "fallback_command") %in% names(records))) {
      return(result)
    }
    # The older audit recorded attempts. A validated fit implies its worker
    # returned success for every slice, including the recorded retries.
    outputs <- sub(".*[[:space:]]--ld=([^[:space:]]+).*", "\\1", commands)
    indices <- match(records$ld_dir, outputs)
    if (anyNA(indices) || anyDuplicated(outputs) || anyDuplicated(indices) ||
        any(modes[indices] != "flame12")) return(result)
    modes[indices] <- "flame1"
    result$feat_flame1_fallback_count <- length(indices)
  } else {
    result$feat_flame1_fallback_count <- 0L
  }
  result$feat_actual_estimator <- paste(intersect(c("ols", "fe", "flame12", "flame1"), modes), collapse = "+")
  result
}
