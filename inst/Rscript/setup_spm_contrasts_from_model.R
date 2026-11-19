#!/usr/bin/env Rscript

# This script generates SPM contrasts from an fmri.pipeline L1 model contrast matrix.

args <- commandArgs(trailingOnly = FALSE)

scriptpath <- dirname(sub("--file=", "", grep("--file=", args, fixed = TRUE, value = TRUE), fixed = TRUE))
argpos <- grep("--args", args, fixed = TRUE)
if (length(argpos) > 0L) {
  args <- args[(argpos + 1):length(args)]
} else {
  args <- c()
}

matfile <- NULL
contrast_rds <- NULL
spm_path <- NULL
average_across_runs <- TRUE

argpos <- 1
while (argpos <= length(args)) {
  if (args[argpos] == "-mat_file") {
    matfile <- args[argpos + 1]
    stopifnot(file.exists(matfile))
    argpos <- argpos + 2
  } else if (args[argpos] == "-contrast_rds") {
    contrast_rds <- args[argpos + 1]
    stopifnot(file.exists(contrast_rds))
    argpos <- argpos + 2
  } else if (args[argpos] == "-spm_path") {
    spm_path <- args[argpos + 1]
    stopifnot(dir.exists(spm_path))
    argpos <- argpos + 2
  } else if (args[argpos] == "-average_across_runs") {
    average_across_runs <- as.logical(args[argpos + 1])
    argpos <- argpos + 2
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

if (is.null(matfile) || is.null(contrast_rds)) {
  stop("Both -mat_file and -contrast_rds are required.")
}

library(R.matlab)

spec <- readRDS(contrast_rds)
cmat <- spec$contrast_matrix
if (is.null(cmat) || !inherits(cmat, "matrix") || nrow(cmat) == 0L) {
  stop("contrast_matrix must be a non-empty matrix in contrast_rds.")
}

if (!is.numeric(cmat) || any(!is.finite(cmat))) stop("SPM contrast weights must be finite numbers.")

regressors <- colnames(cmat)
if (is.null(regressors) || length(regressors) != ncol(cmat)) {
  stop("contrast_matrix must have column names matching regressor names.")
}

projection_interaction_terms <- spec$projection_interaction_terms
if (is.null(projection_interaction_terms)) projection_interaction_terms <- character(0)
projection_interaction_terms <- unique(projection_interaction_terms[!is.na(projection_interaction_terms) & nzchar(projection_interaction_terms)])

normalize_projection_modes <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(character(0))
  lower <- tolower(x)
  if (length(lower) == 1L && lower %in% c("none", "off", "false", "no", "0")) return(character(0))
  if (any(lower %in% c("all", "true", "yes", "1"))) {
    return(c("pooled", "session_specific", "session_differences"))
  }
  out <- vapply(lower, function(val) {
    if (val %in% c("pooled", "condition_means", "condition-means", "main_effect")) return("pooled")
    if (val %in% c("session_specific", "session-specific", "cell_means", "cell-means", "by_session")) return("session_specific")
    if (val %in% c("session_difference", "session-difference", "session_differences", "session-differences", "pairwise_diff", "pairwise_differences")) return("session_differences")
    val
  }, character(1))
  bad <- setdiff(out, c("pooled", "session_specific", "session_differences"))
  if (length(bad) > 0L) stop("Unknown projection_interaction_contrast_modes: ", paste(bad, collapse = ", "))
  unique(out)
}

projection_interaction_contrast_modes <- normalize_projection_modes(spec$projection_interaction_contrast_modes)
projection_interaction_run_labels <- spec$projection_interaction_run_labels
if (is.null(projection_interaction_run_labels)) projection_interaction_run_labels <- numeric(0)

projection_main_effect_terms <- spec$projection_main_effect_terms
if (is.null(projection_main_effect_terms)) projection_main_effect_terms <- character(0)
projection_main_effect_terms <- unique(projection_main_effect_terms[!is.na(projection_main_effect_terms) & nzchar(projection_main_effect_terms)])

projection_main_effect_weights <- spec$projection_main_effect_weights
if (!is.null(projection_main_effect_weights)) {
  if (is.matrix(projection_main_effect_weights)) {
    projection_main_effect_weights <- as.data.frame(
      projection_main_effect_weights,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  if (!is.data.frame(projection_main_effect_weights)) {
    stop("projection_main_effect_weights must be a data.frame or matrix when provided.")
  }
  projection_main_effect_weights[] <- lapply(projection_main_effect_weights, as.numeric)
}

spm_info <- readMat(matfile, fixNames = FALSE)
mnames <- unlist(spm_info$mnames)
cpos <- as.vector(spm_info$cpos)

output_dir <- normalizePath(dirname(matfile))

spm_preamble <- c(
  ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
  "spm('defaults', 'fmri');",
  "spm_jobman('initcfg');",
  ""
)

# Only consider columns of interest
col_names <- mnames[cpos]
col_clean <- sub("Sn\\(\\d+\\)\\s+", "", col_names, perl = TRUE)
col_has_bf <- grepl("\\*bf\\(\\d+\\)$", col_clean)
col_base <- sub("\\*bf\\(\\d+\\)$", "", col_clean, perl = TRUE)
session_tokens <- regmatches(col_names, regexec("^Sn\\((\\d+)\\)\\s+", col_names, perl = TRUE))
session_index <- vapply(session_tokens, function(x) {
  if (length(x) >= 2L) as.integer(x[2]) else 1L
}, integer(1))
session_index_full <- rep(NA_integer_, length(mnames))
session_index_full[cpos] <- session_index

# Parse condition and pmod names from SPM labels.
# SPM can encode pmods as either:
# 1) "cond x pmod-pmod^1" (with spaces), or
# 2) "condxpmod-pmod^1" (no spaces).
parse_cond_pmod <- function(x) {
  out <- list(cond = x, pmod = NA_character_)

  if (grepl("\\sx\\s", x, perl = TRUE)) {
    out$cond <- trimws(sub("\\sx\\s.*", "", x, perl = TRUE))
    out$pmod <- trimws(sub(".*\\sx\\s", "", x, perl = TRUE))
  } else {
    parsed <- regexec("^(.*)x(.*-pmod(?:\\^\\d+)?)$", x, perl = TRUE)
    parts <- regmatches(x, parsed)[[1]]
    if (length(parts) == 3L) {
      out$cond <- trimws(parts[2])
      out$pmod <- trimws(parts[3])
    }
  }

  out$pmod <- sub("\\^\\d+$", "", out$pmod, perl = TRUE)
  out
}

parsed_labels <- lapply(col_base, parse_cond_pmod)
cond_name <- vapply(parsed_labels, function(x) x$cond, character(1))
pmod_name <- vapply(parsed_labels, function(x) x$pmod, character(1))

match_columns <- function(reg) {
  idx <- which(
    col_base == reg |
      cond_name == reg |
      pmod_name == reg |
      pmod_name == paste0(reg, "-pmod") |
      col_base == paste0(reg, "-pmod")
  )
  if (length(idx) == 0L) return(integer(0))

  # Prefer canonical basis function if present
  idx_bf1 <- idx[grepl("\\*bf\\(1\\)$", col_clean[idx])]
  if (length(idx_bf1) > 0L) idx <- idx_bf1
  idx
}

full_cmat <- matrix(0, nrow = nrow(cmat), ncol = length(mnames))
rownames(full_cmat) <- rownames(cmat)
colnames(full_cmat) <- mnames
full_cmat_raw <- full_cmat

missing_regs <- character(0)
for (ii in seq_len(nrow(cmat))) {
  for (jj in seq_len(ncol(cmat))) {
    w <- cmat[ii, jj]
    if (is.na(w) || w == 0) next
    reg <- regressors[jj]
    idx <- match_columns(reg)
    if (length(idx) == 0L) {
      missing_regs <- union(missing_regs, reg)
      next
    }
    full_cmat_raw[ii, cpos[idx]] <- full_cmat_raw[ii, cpos[idx]] + w
    if (isTRUE(average_across_runs)) {
      w <- w / length(idx)
    }
    full_cmat[ii, cpos[idx]] <- full_cmat[ii, cpos[idx]] + w
  }
}

if (length(missing_regs) > 0L) {
  stop("Cannot estimate the requested contrasts: no SPM columns matched these regressors: ", paste(missing_regs, collapse = ", "))
}

extra_rows <- list()
extra_names <- character(0)

clean_label <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x, perl = TRUE)
}

add_extra_row <- function(name, vec) {
  if (!any(abs(vec) > 1e-12)) return(NULL)
  name <- clean_label(name)
  if (name %in% c(rownames(full_cmat), extra_names)) {
    suffix <- 2L
    while (paste0(name, "_", suffix) %in% c(rownames(full_cmat), extra_names)) {
      suffix <- suffix + 1L
    }
    name <- paste0(name, "_", suffix)
  }
  extra_rows[[length(extra_rows) + 1L]] <<- vec
  extra_names <<- c(extra_names, name)
}

if (length(projection_interaction_terms) > 0L &&
    any(c("session_specific", "session_differences") %in% projection_interaction_contrast_modes)) {
  all_sessions <- sort(unique(session_index))
  if (length(projection_interaction_run_labels) == length(all_sessions)) {
    label_map <- as.character(projection_interaction_run_labels)
  } else {
    label_map <- as.character(all_sessions)
  }
  names(label_map) <- as.character(all_sessions)

  for (term in projection_interaction_terms) {
    idx_all <- match_columns(term)
    if (length(idx_all) == 0L) next
    term_sessions <- sort(unique(session_index[idx_all]))
    if (length(term_sessions) <= 1L) next

    if ("session_specific" %in% projection_interaction_contrast_modes) {
      for (ss in term_sessions) {
        idx_ss <- idx_all[session_index[idx_all] == ss]
        if (length(idx_ss) == 0L) next
        vec <- rep(0, length(mnames))
        vec[cpos[idx_ss]] <- 1 / length(idx_ss)
        run_label <- label_map[[as.character(ss)]]
        add_extra_row(paste0("proj_int_", term, "_run", run_label), vec)
      }
    }

    if ("session_differences" %in% projection_interaction_contrast_modes && length(term_sessions) >= 2L) {
      for (jj in 2:length(term_sessions)) {
        for (ii in 1:(jj - 1L)) {
          ss_lo <- term_sessions[ii]
          ss_hi <- term_sessions[jj]
          idx_lo <- idx_all[session_index[idx_all] == ss_lo]
          idx_hi <- idx_all[session_index[idx_all] == ss_hi]
          if (length(idx_lo) == 0L || length(idx_hi) == 0L) next
          vec <- rep(0, length(mnames))
          vec[cpos[idx_hi]] <- 1 / length(idx_hi)
          vec[cpos[idx_lo]] <- -1 / length(idx_lo)
          lo_label <- label_map[[as.character(ss_lo)]]
          hi_label <- label_map[[as.character(ss_hi)]]
          add_extra_row(paste0("proj_int_", term, "_run", hi_label, "_minus_run", lo_label), vec)
        }
      }
    }
  }
}

if (length(projection_main_effect_terms) > 0L &&
    !is.null(projection_main_effect_weights) &&
    is.data.frame(projection_main_effect_weights)) {
  if (nrow(projection_main_effect_weights) > 0L) {
    all_sessions <- sort(unique(session_index))
    max_session <- max(all_sessions)

    for (term in projection_main_effect_terms) {
      if (!term %in% names(projection_main_effect_weights)) next

      term_vals <- as.numeric(projection_main_effect_weights[[term]])
      if (length(term_vals) < max_session) next
      term_vals <- term_vals[all_sessions]

      for (ii in seq_len(nrow(full_cmat_raw))) {
        base_row <- full_cmat_raw[ii, ]
        base_name <- rownames(cmat)[ii]
        if (is.na(base_name) || !nzchar(base_name)) base_name <- paste0("contrast_", ii)
        active_sessions <- all_sessions[vapply(
          all_sessions,
          function(ss) any(abs(base_row[session_index_full == ss]) > 1e-12, na.rm = TRUE),
          logical(1)
        )]
        if (length(active_sessions) <= 1L) next

        centered_vals <- term_vals
        centered_vals[active_sessions] <- centered_vals[active_sessions] - mean(centered_vals[active_sessions])
        if (all(abs(centered_vals[active_sessions]) <= 1e-12)) next

        vec <- base_row
        for (ss in active_sessions) {
          idx_ss <- which(session_index_full == ss)
          vec[idx_ss] <- vec[idx_ss] * centered_vals[which(all_sessions == ss)]
        }

        add_extra_row(
          paste0("proj_", term, "_x_", base_name),
          vec
        )
      }
    }
  }
}

if (length(extra_rows) > 0L) {
  extra_mat <- do.call(rbind, extra_rows)
  colnames(extra_mat) <- colnames(full_cmat)
  rownames(extra_mat) <- extra_names
  full_cmat <- rbind(full_cmat, extra_mat)
}

# Drop contrasts with all zero weights
nonzero_rows <- apply(full_cmat, 1, function(x) any(abs(x) > 1e-12))
if (!any(nonzero_rows)) {
  stop("All contrasts are empty after mapping to SPM columns.")
}
if (any(!nonzero_rows)) stop("Zero-weight SPM contrasts: ", paste(rownames(full_cmat)[!nonzero_rows], collapse = ", "))

if (!is.null(spm_info$design_matrix)) {
  X <- as.matrix(spm_info$design_matrix)
  if (ncol(X) != ncol(full_cmat)) stop("SPM design and contrast columns disagree.")
  # Match SPM's singular-value rank tolerance, including small-scale EVs.
  basis <- svd(X, nu = 0L)
  tol <- max(dim(X)) * max(basis$d) * .Machine$double.eps
  Q <- basis$v[, basis$d > tol, drop = FALSE]
  projected <- full_cmat %*% Q %*% t(Q)
  invalid <- apply(abs(full_cmat - projected), 1L, max) > 1e-7 * pmax(1, apply(abs(full_cmat), 1L, max))
  if (any(invalid)) stop("Non-estimable SPM contrasts: ", paste(rownames(full_cmat)[invalid], collapse = ", "))
}
if (anyDuplicated(rownames(full_cmat)) || any(grepl("[\t\r\n]", rownames(full_cmat)))) {
  stop("SPM contrast names must be unique and contain no tabs or newlines.")
}
writeLines(rownames(full_cmat), file.path(output_dir, "spm_expected_contrasts.txt"))

cn <- 1
m_string <- c(
  spm_preamble,
  "% SETUP BATCH JOB STRUCTURE",
  "contrast = struct;",
  "% spmmat",
  paste0("contrast.matlabbatch{1}.spm.stats.con.spmmat = { ['", output_dir, "' filesep 'SPM.mat']};"),
  ""
)

for (ii in seq_len(nrow(full_cmat))) {
  cname <- rownames(full_cmat)[ii]
  cname <- gsub("'", "''", cname, fixed = TRUE)
  m_string <- c(
    m_string,
    paste0("% consess ", cn),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.tcon.name = '", cname, "';"),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.tcon.convec = [ ", paste(full_cmat[ii, ], collapse = ", "), " ];"),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.tcon.sessrep = 'none';")
  )
  cn <- cn + 1
}

m_string <- c(
  m_string,
  "",
  "% delete?",
  "contrast.matlabbatch{1}.spm.stats.con.delete = 1;",
  "% RUN BATCH JOB",
  "spm_jobman('run',contrast.matlabbatch);"
)

cat(m_string, file = file.path(output_dir, "estimate_glm_contrasts.m"), sep = "\n")
