#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: feat_integrity.R <prepare|validate> <fsf> [folder action]")
if (args[1L] == "prepare") {
  fmri.pipeline:::feat_prepare_run(args[2L], args[3L])
} else if (args[1L] == "validate") {
  fmri.pipeline:::feat_validate_outputs(args[2L])
} else stop("Unknown FEAT integrity action")
