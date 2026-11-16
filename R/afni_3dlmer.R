format_3dlmer_datatable_keys <- function(df, id_col = "id", session_col = "session", max_items = 8L) {
  if (nrow(df) == 0L) return("<none>")

  keys <- paste0(
    "id=", as.character(df[[id_col]]),
    ", session=", as.character(df[[session_col]])
  )
  keys <- unique(keys)
  if (length(keys) > max_items) {
    keys <- c(keys[seq_len(max_items)], sprintf("... (%d more)", length(keys) - max_items))
  }

  paste(keys, collapse = "; ")
}

format_3dlmer_datatable_paths <- function(paths, max_items = 8L) {
  paths <- unique(as.character(paths))
  if (length(paths) == 0L) return("<none>")
  if (length(paths) > max_items) {
    paths <- c(paths[seq_len(max_items)], sprintf("... (%d more)", length(paths) - max_items))
  }

  paste(paths, collapse = "; ")
}

is_missing_3dlmer_field <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x)))
}

find_3dlmer_rows_with_missing_fields <- function(df, fields) {
  if (length(fields) == 0L) return(logical(nrow(df)))

  missing_matrix <- vapply(
    fields,
    function(field) is_missing_3dlmer_field(df[[field]]),
    logical(nrow(df))
  )
  if (is.null(dim(missing_matrix))) return(missing_matrix)

  rowSums(missing_matrix) > 0L
}

find_3dlmer_duplicate_keys <- function(df, id_col = "id", session_col = "session") {
  keys <- paste(as.character(df[[id_col]]), as.character(df[[session_col]]), sep = "\r")
  duplicated(keys) | duplicated(keys, fromLast = TRUE)
}

#' Internal function to build the data table for AFNI 3dLMEr
#'
#' @param subject_data a data.frame containing subject/session-level covariates with unique id/session rows
#' @param input_files a data.frame with columns 'id', 'session', and 'InputFile'
#' @param model_variables a character vector of variables that must be included in the table
#'
#' @return a data.frame formatted for 3dLMEr -dataTable
#' @keywords internal
build_3dlmer_datatable <- function(subject_data, input_files, model_variables) {
  checkmate::assert_data_frame(subject_data)
  checkmate::assert_data_frame(input_files)
  checkmate::assert_character(model_variables)
  checkmate::assert_subset(c("id", "session", model_variables), names(subject_data))
  checkmate::assert_subset(c("id", "session", "InputFile"), names(input_files))

  subject_data <- subject_data %>%
    dplyr::select(dplyr::all_of(c("id", "session", model_variables)))

  input_files <- input_files %>%
    dplyr::select("id", "session", "InputFile")

  if (nrow(input_files) == 0L) {
    stop("AFNI 3dLMEr dataTable has no input rows.", call. = FALSE)
  }
  if (nrow(subject_data) == 0L) {
    stop("AFNI 3dLMEr subject_data has no rows.", call. = FALSE)
  }

  missing_input_rows <- find_3dlmer_rows_with_missing_fields(input_files, c("id", "session", "InputFile"))
  if (any(missing_input_rows)) {
    stop(
      sprintf(
        "AFNI 3dLMEr input_files has missing id, session, or InputFile values for rows: %s.",
        paste(which(missing_input_rows), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_subject_keys <- find_3dlmer_rows_with_missing_fields(subject_data, c("id", "session"))
  if (any(missing_subject_keys)) {
    stop(
      sprintf(
        "AFNI 3dLMEr subject_data has missing id or session values for rows: %s.",
        paste(which(missing_subject_keys), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing_covariates <- find_3dlmer_rows_with_missing_fields(subject_data, model_variables)
  if (any(missing_covariates)) {
    stop(
      sprintf(
        "AFNI 3dLMEr subject_data has missing model variable values for id/session keys: %s.",
        format_3dlmer_datatable_keys(subject_data[missing_covariates, , drop = FALSE])
      ),
      call. = FALSE
    )
  }

  duplicate_input_keys <- find_3dlmer_duplicate_keys(input_files)
  if (any(duplicate_input_keys)) {
    stop(
      sprintf(
        "Duplicate AFNI 3dLMEr input rows for id/session keys: %s. Each contrast dataTable must have exactly one InputFile per id/session.",
        format_3dlmer_datatable_keys(input_files[duplicate_input_keys, , drop = FALSE])
      ),
      call. = FALSE
    )
  }

  duplicate_subject_keys <- find_3dlmer_duplicate_keys(subject_data)
  if (any(duplicate_subject_keys)) {
    stop(
      sprintf(
        "Duplicate AFNI 3dLMEr subject_data rows for id/session keys: %s. Covariates must be one-to-one with AFNI input rows.",
        format_3dlmer_datatable_keys(subject_data[duplicate_subject_keys, , drop = FALSE])
      ),
      call. = FALSE
    )
  }

  input_paths <- as.character(input_files$InputFile)
  missing_files <- input_paths[!file.exists(input_paths)]
  if (length(missing_files) > 0L) {
    stop(
      sprintf(
        "AFNI 3dLMEr InputFile entries do not exist: %s.",
        format_3dlmer_datatable_paths(missing_files)
      ),
      call. = FALSE
    )
  }

  missing_subject_rows <- input_files %>%
    dplyr::anti_join(subject_data %>% dplyr::select("id", "session"), by = c("id", "session"))
  if (nrow(missing_subject_rows) > 0L) {
    stop(
      sprintf(
        "AFNI 3dLMEr dataTable cannot be built because subject_data is missing id/session rows for input keys: %s.",
        format_3dlmer_datatable_keys(missing_subject_rows)
      ),
      call. = FALSE
    )
  }

  # Merge input files with subject data to get covariates
  # subject_data should contain 'id' and 'session' for longitudinal
  dt <- input_files %>%
    dplyr::left_join(subject_data, by = c("id", "session"))

  if (nrow(dt) != nrow(input_files)) {
    stop(
      sprintf(
        "AFNI 3dLMEr dataTable row count changed during merge: %d input row(s), %d merged row(s).",
        nrow(input_files), nrow(dt)
      ),
      call. = FALSE
    )
  }

  # Keep only necessary columns: Subj, model_variables, InputFile
  # 3dLMEr requires the column name 'Subj' for the subject ID
  dt <- dt %>%
    dplyr::rename(Subj = "id") %>%
    dplyr::select(dplyr::all_of(c("Subj", "session", model_variables, "InputFile")))

  duplicate_final_keys <- find_3dlmer_duplicate_keys(dt, id_col = "Subj")
  if (any(duplicate_final_keys)) {
    stop(
      sprintf(
        "AFNI 3dLMEr dataTable contains duplicate Subj/session rows after merge: %s.",
        format_3dlmer_datatable_keys(dt[duplicate_final_keys, , drop = FALSE], id_col = "Subj")
      ),
      call. = FALSE
    )
  }

  return(dt)
}

build_3dlmer_intersection_mask <- function(mask_files, output_file, lg = NULL) {
  checkmate::assert_character(mask_files, min.len = 1L)
  checkmate::assert_file_exists(mask_files)
  checkmate::assert_string(output_file)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  mask_files <- unique(mask_files)
  if (is.null(lg)) lg <- lgr::get_logger("fmri.pipeline/afni_3dlmer/mask")
  if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)

  mask_imgs <- lapply(mask_files, RNifti::readNifti)
  mask_dims <- vapply(mask_imgs, function(img) paste(dim(img), collapse = "x"), character(1))
  if (length(unique(mask_dims)) != 1L) {
    stop(
      sprintf(
        "Cannot build AFNI 3dLMEr intersection mask from mismatched FSL masks: %s",
        paste(unique(mask_dims), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  mask_4d <- do.call(abind::abind, c(mask_imgs, list(along = 4)))
  mask_sum <- rowSums(mask_4d, dims = 3)
  intersect_img <- 1L * (mask_sum == dim(mask_4d)[4L])

  header <- RNifti::niftiHeader(mask_files[1L])
  intersect_nifti <- RNifti::asNifti(intersect_img, header)
  RNifti::writeNifti(intersect_nifti, file = output_file)

  lg$info(
    "Computed AFNI 3dLMEr intersection mask from %d FSL L2 masks: %s",
    length(mask_files), output_file
  )

  output_file
}

resolve_3dlmer_mask <- function(target_dir, harvested_inputs, explicit_mask = NULL, requires_l2 = FALSE, lg = NULL) {
  checkmate::assert_string(target_dir)
  checkmate::assert_data_frame(harvested_inputs)
  checkmate::assert_string(explicit_mask, null.ok = TRUE)
  checkmate::assert_flag(requires_l2)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)

  if (!is.null(explicit_mask)) return(explicit_mask)
  if (!"feat_dir" %in% names(harvested_inputs)) return(NULL)

  mask_files <- unique(vapply(as.character(harvested_inputs$feat_dir), function(fd) {
    if (!isTRUE(requires_l2)) {
      reg_std_mask <- file.path(fd, "reg_standard", "mask.nii.gz")
      if (file.exists(reg_std_mask)) return(reg_std_mask)
    }
    file.path(fd, "mask.nii.gz")
  }, character(1)))
  if (length(mask_files) == 0L || anyNA(mask_files) || any(!file.exists(mask_files))) {
    if (!is.null(lg)) {
      lg$warn("Could not derive AFNI 3dLMEr intersection mask because one or more FSL masks were missing.")
    }
    return(NULL)
  }

  build_3dlmer_intersection_mask(
    mask_files = mask_files,
    output_file = file.path(target_dir, "intersection_mask.nii.gz"),
    lg = lg
  )
}

#' Normalize an AFNI 3dLMEr model formula
#'
#' @param model_formula fixed-effects formula string.
#'
#' @return A character string with whitespace removed.
#' @keywords internal
# Collapse whitespace because AFNI counts spaces inside -model formulas.
sanitize_3dlmer_model_formula <- function(model_formula) {
  checkmate::assert_string(model_formula)
  gsub("[[:space:]]+", "", model_formula)
}

validate_3dlmer_formula_datatable <- function(model_formula, datatable, context = NULL) {
  checkmate::assert_string(model_formula)
  checkmate::assert_data_frame(datatable)
  checkmate::assert_string(context, null.ok = TRUE)

  formula_text <- trimws(model_formula)
  if (!grepl("~", formula_text, fixed = TRUE)) {
    formula_text <- paste("~", formula_text)
  }

  parsed_formula <- tryCatch(
    stats::as.formula(formula_text, env = baseenv()),
    error = function(e) {
      stop(
        sprintf(
          "Could not parse AFNI 3dLMEr lmer_formula%s: %s",
          if (!is.null(context)) paste0(" for ", context) else "",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  formula_vars <- unique(all.vars(parsed_formula))
  missing_vars <- setdiff(formula_vars, names(datatable))
  if (length(missing_vars) > 0L) {
    msg <- sprintf(
      "AFNI 3dLMEr lmer_formula%s references variable(s) not present in dataTable: %s. Available dataTable columns: %s.",
      if (!is.null(context)) paste0(" for ", context) else "",
      paste(missing_vars, collapse = ", "),
      paste(names(datatable), collapse = ", ")
    )
    if ("id" %in% missing_vars && "Subj" %in% names(datatable)) {
      msg <- paste(
        msg,
        "Use 'Subj' instead of 'id' in 3dLMEr random-effects terms because the AFNI dataTable renames subject IDs to Subj."
      )
    }
    stop(msg, call. = FALSE)
  }

  invisible(formula_vars)
}

sanitize_3dlmer_glt_names <- function(glt_names) {
  checkmate::assert_character(glt_names, any.missing = TRUE)

  glt_names[is.na(glt_names)] <- ""
  glt_names <- trimws(glt_names)
  glt_names <- gsub("[^A-Za-z0-9_]+", "_", glt_names)
  glt_names <- gsub("_+", "_", glt_names)
  glt_names <- gsub("^_+|_+$", "", glt_names)

  empty_names <- !nzchar(glt_names)
  if (any(empty_names)) {
    glt_names[empty_names] <- paste0("glt", which(empty_names))
  }

  starts_badly <- !grepl("^[A-Za-z]", glt_names)
  glt_names[starts_badly] <- paste0("glt_", glt_names[starts_badly])
  make.unique(glt_names, sep = "_")
}

empty_3dlmer_glt_table <- function() {
  data.frame(
    label_raw = character(0),
    label_afni = character(0),
    code = character(0),
    source = character(0),
    contrast_type = character(0),
    valid = logical(0),
    message = character(0),
    stringsAsFactors = FALSE
  )
}

as_3dlmer_glt_table <- function(glt_codes, source = "user", contrast_type = "raw") {
  if (is.null(glt_codes) || length(glt_codes) == 0L) return(empty_3dlmer_glt_table())

  if (is.data.frame(glt_codes)) {
    required <- c("label_raw", "code")
    missing_cols <- setdiff(required, names(glt_codes))
    if (length(missing_cols) > 0L) {
      stop("3dLMEr GLT table is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    out <- as.data.frame(glt_codes, stringsAsFactors = FALSE)
    if (!"label_afni" %in% names(out)) out$label_afni <- sanitize_3dlmer_glt_names(out$label_raw)
    if (!"source" %in% names(out)) out$source <- source
    if (!"contrast_type" %in% names(out)) out$contrast_type <- contrast_type
    if (!"valid" %in% names(out)) out$valid <- TRUE
    if (!"message" %in% names(out)) out$message <- ""
    return(out[, names(empty_3dlmer_glt_table()), drop = FALSE])
  }

  if (is.character(glt_codes)) {
    glt_codes <- as.list(glt_codes)
  }
  checkmate::assert_list(glt_codes)
  labels <- names(glt_codes)
  if (is.null(labels)) labels <- paste0("glt", seq_along(glt_codes))
  labels[is.na(labels) | !nzchar(labels)] <- paste0("glt", which(is.na(labels) | !nzchar(labels)))

  data.frame(
    label_raw = labels,
    label_afni = sanitize_3dlmer_glt_names(labels),
    code = as.character(unlist(glt_codes, use.names = FALSE)),
    source = source,
    contrast_type = contrast_type,
    valid = TRUE,
    message = "",
    stringsAsFactors = FALSE
  )
}

parse_3dlmer_glt_blocks <- function(code, datatable_names) {
  checkmate::assert_string(code)
  checkmate::assert_character(datatable_names, min.len = 1L)

  var_pattern <- paste(gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", datatable_names), collapse = "|")
  matches <- gregexpr(paste0("(?<![A-Za-z0-9_.])(", var_pattern, ")\\s*:"), code, perl = TRUE)[[1L]]
  if (identical(matches[1L], -1L)) {
    return(list(blocks = list(), error = "No 'variable :' block found."))
  }

  pre <- substr(code, 1L, matches[1L] - 1L)
  if (nzchar(trimws(pre))) {
    return(list(blocks = list(), error = sprintf("Unexpected text before first GLT variable block: %s", trimws(pre))))
  }

  lengths <- attr(matches, "match.length")
  headers <- regmatches(code, list(matches))[[1L]]
  vars <- trimws(sub(":\\s*$", "", headers))
  blocks <- vector("list", length(matches))
  for (ii in seq_along(matches)) {
    start <- matches[ii] + lengths[ii]
    end <- if (ii < length(matches)) matches[ii + 1L] - 1L else nchar(code)
    blocks[[ii]] <- list(var = vars[ii], expr = trimws(substr(code, start, end)))
  }

  list(blocks = blocks, error = "")
}

validate_3dlmer_glt_table <- function(glt_table, datatable, qVars = NULL, context = NULL) {
  checkmate::assert_data_frame(datatable)
  checkmate::assert_character(qVars, null.ok = TRUE)
  checkmate::assert_string(context, null.ok = TRUE)

  glt_table <- as_3dlmer_glt_table(glt_table)
  if (nrow(glt_table) == 0L) return(glt_table)

  qVars <- qVars %||% character(0)
  datatable_names <- names(datatable)
  numeric_pattern <- "[+-]?(?:(?:[0-9]+(?:\\.[0-9]*)?)|(?:\\.[0-9]+))(?:[eE][+-]?[0-9]+)?"
  errors <- character(nrow(glt_table))

  for (ii in seq_len(nrow(glt_table))) {
    code <- glt_table$code[ii]
    row_errors <- character(0)

    if (is.na(code) || !nzchar(trimws(code))) {
      row_errors <- c(row_errors, "GLT code is empty.")
    } else {
      parsed <- parse_3dlmer_glt_blocks(code, datatable_names)
      if (nzchar(parsed$error)) {
        row_errors <- c(row_errors, parsed$error)
      } else {
        cat_weight_list <- list()
        for (block in parsed$blocks) {
          var <- block$var
          expr <- block$expr
          if (!var %in% datatable_names) {
            row_errors <- c(row_errors, sprintf("Variable '%s' is not in dataTable.", var))
            next
          }

          is_qvar <- var %in% qVars || (is.numeric(datatable[[var]]) && !is.factor(datatable[[var]]))
          tokens <- strsplit(expr, "\\s+")[[1L]]
          tokens <- tokens[nzchar(tokens)]

          if (is_qvar) {
            # In AFNI 3dLMEr, an empty expression after ':' (e.g., 'Age :' or 'Group : 1*A -1*B Age :')
            # is the required syntax to test the linear slope/effect of a quantitative variable, whereas
            # 'Age : <val>' conditions categorical estimates at an explicit covariate offset.
            if (length(tokens) > 0L) {
              bad <- !grepl(paste0("^", numeric_pattern, "$"), tokens, perl = TRUE)
              if (any(bad) || length(tokens) != 1L) {
                row_errors <- c(row_errors, sprintf(
                  "Quantitative variable '%s' must use '%s :' (slope) or one numeric value, e.g. '%s : 0'.",
                  var, var, var
                ))
              }
            }
          } else {
            if (length(tokens) == 0L) {
              row_errors <- c(row_errors, sprintf("Variable '%s' has no GLT weight expression.", var))
              next
            }
            bad_syntax <- !grepl(paste0("^", numeric_pattern, "\\*.+$"), tokens, perl = TRUE)
            if (any(bad_syntax)) {
              row_errors <- c(row_errors, sprintf(
                "Categorical variable '%s' has invalid level weight token(s): %s.",
                var, paste(tokens[bad_syntax], collapse = ", ")
              ))
              next
            }
            levels_seen <- sub(paste0("^", numeric_pattern, "\\*"), "", tokens, perl = TRUE)
            known_levels <- unique(as.character(datatable[[var]]))
            missing_levels <- setdiff(levels_seen, known_levels)
            if (length(missing_levels) > 0L) {
              row_errors <- c(row_errors, sprintf(
                "Categorical variable '%s' references level(s) not present in dataTable: %s.",
                var, paste(missing_levels, collapse = ", ")
              ))
            }
            cat_weight_list[[var]] <- as.numeric(sub("\\*.*$", "", tokens, perl = TRUE))
          }
        }

        if (identical(glt_table$contrast_type[ii], "pairwise_diff") && length(cat_weight_list) > 0L) {
          is_zero_sum <- vapply(cat_weight_list, function(w) {
            length(w) >= 2L && isTRUE(all.equal(sum(w), 0, tolerance = 1e-5))
          }, logical(1))
          is_cond_slice <- vapply(cat_weight_list, function(w) {
            all(w >= -1e-8) && isTRUE(all.equal(sum(w), 1, tolerance = 1e-5))
          }, logical(1))

          if (!any(is_zero_sum)) {
            bad_var <- names(cat_weight_list)[1L]
            row_errors <- c(row_errors, sprintf(
              "Pairwise AFNI 3dLMEr GLT '%s' must use explicit cell-level weights that sum to zero for '%s' (for example, '%s : 1*A -1*B'), not treatment-coded coefficient weights.",
              glt_table$label_raw[ii], bad_var, bad_var
            ))
          } else if (any(!(is_zero_sum | is_cond_slice))) {
            bad_vars <- names(cat_weight_list)[!(is_zero_sum | is_cond_slice)]
            for (bv in bad_vars) {
              row_errors <- c(row_errors, sprintf(
                "Pairwise AFNI 3dLMEr GLT '%s' must use explicit cell-level weights that sum to zero for '%s' (for example, '%s : 1*A -1*B'), not treatment-coded coefficient weights.",
                glt_table$label_raw[ii], bv, bv
              ))
            }
          }
        }
      }
    }

    errors[ii] <- paste(row_errors, collapse = " ")
  }

  glt_table$valid <- !nzchar(errors)
  glt_table$message <- errors
  if (any(!glt_table$valid)) {
    bad <- glt_table[!glt_table$valid, , drop = FALSE]
    stop(
      sprintf(
        "Invalid AFNI 3dLMEr GLT code%s: %s",
        if (!is.null(context)) paste0(" for ", context) else "",
        paste(sprintf("%s: %s", bad$label_raw, bad$message), collapse = " | ")
      ),
      call. = FALSE
    )
  }

  glt_table
}

#' Construct an AFNI 3dLMEr command string
#'
#' @param prefix Output prefix for 3dLMEr.
#' @param model_formula Fixed-effects formula string.
#' @param qVars Character vector of quantitative variables.
#' @param qVarCenters Optional numeric vector of centering values corresponding to \code{qVars}.
#' @param glt_codes List, character vector, or data frame of named GLT code strings.
#' @param data_table_file Path to the data-table text file.
#' @param mask Optional path to the brain-mask file.
#' @param njobs Number of parallel jobs for 3dLMEr.
#' @param ss_type Sum-of-squares type (1--3; defaults to 3).
#'
#' @return A character string containing the 3dLMEr command.
#' @keywords internal
build_3dlmer_command <- function(prefix, model_formula, qVars = NULL, qVarCenters = NULL, glt_codes = NULL,
                                 data_table_file, mask = NULL, njobs = 1, ss_type = 3) {
  checkmate::assert_string(prefix)
  checkmate::assert_string(model_formula)
  checkmate::assert_character(qVars, null.ok = TRUE)
  checkmate::assert_numeric(qVarCenters, null.ok = TRUE, any.missing = FALSE)
  if (!is.null(qVarCenters) && length(qVarCenters) > 0L && length(qVarCenters) != length(qVars)) {
    stop("qVarCenters must have the same length as qVars.", call. = FALSE)
  }
  if (!is.null(glt_codes) &&
      !checkmate::test_list(glt_codes) &&
      !checkmate::test_data_frame(glt_codes) &&
      !checkmate::test_character(glt_codes)) {
    stop("glt_codes must be a list, character vector, data.frame, or NULL.", call. = FALSE)
  }
  checkmate::assert_string(data_table_file)
  checkmate::assert_string(mask, null.ok = TRUE)
  checkmate::assert_integerish(njobs, lower = 1)
  checkmate::assert_integerish(ss_type, lower = 1, upper = 3)

  model_formula <- sanitize_3dlmer_model_formula(model_formula)

  cmd <- glue::glue(
    "3dLMEr -prefix {shQuote(prefix, type = 'sh')} -model {shQuote(model_formula, type = 'sh')}"
  )

  if (!is.null(qVars) && length(qVars) > 0) {
    cmd <- paste(cmd, glue::glue("-qVars {shQuote(paste(qVars, collapse=\",\"), type = 'sh')}"))
    if (!is.null(qVarCenters) && length(qVarCenters) == length(qVars)) {
      centers_str <- paste(format(as.numeric(qVarCenters), trim = TRUE, scientific = FALSE), collapse = ",")
      cmd <- paste(cmd, glue::glue("-qVarCenters {shQuote(centers_str, type = 'sh')}"))
    }
  }

  if (!is.null(mask)) {
    cmd <- paste(cmd, glue::glue("-mask {shQuote(mask, type = 'sh')}"))
  }

  if (njobs > 1) {
    cmd <- paste(cmd, glue::glue("-jobs {njobs}"))
  }

  cmd <- paste(cmd, glue::glue("-SS_type {ss_type}"))

  if (!is.null(glt_codes) && length(glt_codes) > 0) {
    glt_table <- as_3dlmer_glt_table(glt_codes)
    for (i in seq_len(nrow(glt_table))) {
      glt_name <- glt_table$label_afni[i]
      glt_code <- glt_table$code[i]
      # glt_code should already be in AFNI format: "Factor : 1*Level1 -1*Level2"
      cmd <- paste(cmd, glue::glue("-gltCode {glt_name} {shQuote(glt_code, type = 'sh')}"))
    }
  }

  cmd <- paste(cmd, glue::glue("-dataTable {shQuote(paste0('@', data_table_file), type = 'sh')}"))

  return(cmd)
}

#' Helper to write 3dLMEr script and data table
#'
#' @param output_dir directory where files should be written
#' @param dt the data table data.frame
#' @param cmd the 3dLMEr command string
#' @param script_name name of the shell script to write
#'
#' @return list with paths to script and data table
#' @keywords internal
write_3dlmer_files <- function(output_dir, dt, cmd, script_name = "run_3dlmer.sh") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  dt_out <- as.data.frame(dt, stringsAsFactors = FALSE)
  for (cn in setdiff(names(dt_out), "InputFile")) {
    if (is.factor(dt_out[[cn]]) || is.character(dt_out[[cn]])) {
      vals <- as.character(dt_out[[cn]])
      if (any(grepl("\\s", vals, perl = TRUE), na.rm = TRUE)) {
        warning(
          sprintf("Column '%s' in 3dLMEr dataTable contains whitespace; replacing whitespace with '_' for AFNI table parsing.", cn),
          call. = FALSE
        )
        vals <- gsub("\\s+", "_", trimws(vals), perl = TRUE)
      }
      dt_out[[cn]] <- vals
    }
  }
  if ("InputFile" %in% names(dt_out)) {
    in_paths <- as.character(dt_out$InputFile)
    has_ws <- grepl("\\s", in_paths, perl = TRUE)
    if (any(has_ws, na.rm = TRUE)) {
      in_paths[has_ws] <- paste0("\"", in_paths[has_ws], "\"")
      dt_out$InputFile <- in_paths
    }
  }

  dt_file <- file.path(output_dir, "dataTable.txt")
  # write.table with proper formatting for AFNI (tab separated, header)
  write.table(dt_out, file = dt_file, row.names = FALSE, quote = FALSE, sep = "\t")

  script_file <- file.path(output_dir, script_name)
  # Basic script structure
  script_content <- c(
    "#!/bin/bash",
    "set -uo pipefail",
    "",
    "script_dir=\"$(cd \"$(dirname \"$0\")\" && pwd)\"",
    "cd \"$script_dir\"",
    "",
    "complete_file=\".afni_complete\"",
    "fail_file=\".afni_fail\"",
    "start_file=\".afni_start\"",
    "log_file=\"3dLMEr.log\"",
    "rm -f \"$complete_file\" \"$fail_file\"",
    "start_time=\"$(date -Is)\"",
    "printf '%s\\n' \"$start_time\" > \"$start_file\"",
    "",
    "{",
    "  printf '%s\\n' \"[$start_time] Starting AFNI 3dLMEr\"",
    paste0("  printf '%s\\n' ", shQuote(paste("Command:", cmd), type = "sh")),
    paste0("  ", cmd),
    "} >> \"$log_file\" 2>&1",
    "exit_code=$?",
    "end_time=\"$(date -Is)\"",
    "",
    "if [ $exit_code -eq 0 ]; then",
    "  printf '%s\\n%s\\n' \"$start_time\" \"$end_time\" > \"$complete_file\"",
    "  rm -f \"$fail_file\"",
    "else",
    "  printf '%s\\n%s\\n' \"$start_time\" \"$end_time\" > \"$fail_file\"",
    "fi",
    "",
    "exit $exit_code"
  )
  writeLines(script_content, script_file)
  Sys.chmod(script_file, "0755")

  return(list(script = script_file, data_table = dt_file))
}

#' Status checker for 3dLMEr
#'
#' @param output_file expected output file from 3dLMEr
#' @param lg optional logger
#' @param prefix optional prefix for returned column names
#' @return data.frame indicating output existence and completion status
#' @keywords internal
get_3dlmer_status <- function(output_file, lg = NULL, prefix = NULL) {
  checkmate::assert_string(output_file)
  checkmate::assert_class(lg, "Logger", null.ok = TRUE)
  checkmate::assert_string(prefix, null.ok = TRUE)
  if (is.null(lg)) lg <- lgr::get_logger()

  resolved_output <- output_file
  candidates <- unique(c(
    output_file,
    paste0(output_file, ".HEAD"),
    paste0(output_file, ".BRIK"),
    paste0(output_file, "+tlrc.HEAD"),
    paste0(output_file, "+orig.HEAD")
  ))

  if (grepl("\\+tlrc$|\\+orig$", output_file)) {
    candidates <- unique(c(
      output_file,
      paste0(output_file, ".HEAD"),
      paste0(output_file, ".BRIK")
    ))
  } else if (grepl("\\.(nii|nii\\.gz|HEAD|BRIK)$", output_file, perl = TRUE)) {
    sans_ext <- file_sans_ext(output_file)
    if (!is.na(sans_ext)) {
      candidates <- unique(c(
        output_file,
        paste0(sans_ext, "+tlrc.HEAD"),
        paste0(sans_ext, "+orig.HEAD")
      ))
    }
  }

  existing <- candidates[file.exists(candidates)]
  output_exists <- length(existing) > 0L
  if (output_exists) {
    resolved_output <- existing[1L]
  } else {
    lg$debug("No AFNI 3dLMEr output found for '%s'. Checked: %s", output_file, paste(candidates, collapse = ", "))
  }

  output_dir <- dirname(output_file)
  start_file <- file.path(output_dir, ".afni_start")
  complete_file <- file.path(output_dir, ".afni_complete")
  fail_file <- file.path(output_dir, ".afni_fail")
  log_file <- file.path(output_dir, "3dLMEr.log")
  start_file_exists <- file.exists(start_file)
  complete_file_exists <- file.exists(complete_file)
  fail_file_exists <- file.exists(fail_file)

  execution_start <- as.POSIXct(NA)
  execution_end <- as.POSIXct(NA)
  execution_min <- NA_real_
  timing_file <- NULL

  if (complete_file_exists) {
    timing_file <- complete_file
    if (fail_file_exists) {
      lg$warn("Both .afni_complete and .afni_fail files exist in %s", output_dir)
      lg$warn(
        "Treating this AFNI 3dLMEr run as %s based on output files.",
        ifelse(isTRUE(output_exists), "complete", "failed")
      )
    }
    if (!isTRUE(output_exists)) {
      lg$warn("Found .afni_complete in %s, but expected AFNI 3dLMEr output is missing.", output_dir)
    }
    afni_complete <- isTRUE(output_exists)
    afni_failed <- !isTRUE(output_exists)
  } else if (fail_file_exists) {
    lg$debug("Detected AFNI 3dLMEr failure marker in: %s", output_dir)
    timing_file <- fail_file
    afni_complete <- FALSE
    afni_failed <- TRUE
  } else {
    afni_complete <- output_exists
    afni_failed <- if (output_exists) FALSE else NA
    if (start_file_exists) timing_file <- start_file
  }

  if (!is.null(timing_file)) {
    timing <- readLines(timing_file, warn = FALSE)
    if (length(timing) > 0L) {
      timing <- anytime::anytime(timing)
      execution_start <- timing[1L]
      if (length(timing) >= 2L) {
        execution_end <- timing[2L]
        execution_min <- as.numeric(difftime(timing[2L], timing[1L], units = "mins"))
      } else if (!identical(timing_file, start_file)) {
        lg$warn("Did not find two timing entries in %s.", timing_file)
        lg$warn("File contents: %s", paste(timing, collapse = "; "))
      }
    }
  }

  out <- data.frame(
    output_file = output_file,
    afni_output_file_resolved = resolved_output,
    afni_output_file_exists = output_exists,
    afni_start_file = start_file,
    afni_start_file_exists = start_file_exists,
    afni_complete_file = complete_file,
    afni_complete_file_exists = complete_file_exists,
    afni_fail_file = fail_file,
    afni_fail_file_exists = fail_file_exists,
    afni_log_file = log_file,
    afni_log_file_exists = file.exists(log_file),
    afni_execution_start = execution_start,
    afni_execution_end = execution_end,
    afni_execution_min = execution_min,
    afni_complete = afni_complete,
    afni_failed = afni_failed,
    stringsAsFactors = FALSE
  )
  if (!is.null(prefix)) names(out) <- paste0(prefix, names(out))
  out
}

#' Translator from emmeans-style contrast spec to 3dLMEr -gltCode
#'
#' @param mobj high-level model specification object (hi_model_spec)
#' @param data the data.frame used for the model (to check factor levels)
#'
#' @return a list of 3dLMEr gltCode strings
#' @keywords internal
emmeans_to_3dlmer_glt <- function(mobj, data, qVars = NULL, raw_glt_codes = NULL, context = NULL) {
  checkmate::assert_class(mobj, "hi_model_spec")
  checkmate::assert_data_frame(data)
  checkmate::assert_character(qVars, null.ok = TRUE)
  checkmate::assert_string(context, null.ok = TRUE)
  lg <- lgr::get_logger("fmri.pipeline/afni_3dlmer/emmeans_to_3dlmer_glt")

  generated <- empty_3dlmer_glt_table()
  if (is.null(mobj$lmfit)) {
    lg$warn("No lmfit in mobj. Cannot generate emmeans-to-3dlmer GLTs.")
  }

  spec <- mobj$contrast_spec
  if (!is.null(spec) && !is.null(mobj$lmfit)) {
    add_generated <- function(label, code, contrast_type) {
      generated <<- rbind(
        generated,
        data.frame(
          label_raw = label,
          label_afni = label,
          code = code,
          source = "generated",
          contrast_type = contrast_type,
          valid = TRUE,
          message = "",
          stringsAsFactors = FALSE
        )
      )
    }

    eff_weights <- spec$weights %||% "cells"
    vtypes <- extract_lm_var_types(mobj$lmfit)
    cat_cols <- intersect(if (length(spec$cat_vars) > 0L) spec$cat_vars else vtypes$cat_vars, names(data))
    cont_cols <- unique(intersect(c(qVars, vtypes$cont_vars), names(data)))

    get_levels_ordered <- function(col_nm) {
      if (!col_nm %in% names(data)) {
        stop(sprintf("Cannot generate AFNI 3dLMEr GLTs for '%s' because it is not in the dataTable.", col_nm), call. = FALSE)
      }
      x <- data[[col_nm]]
      obs <- unique(as.character(x[!is.na(x)]))
      if (is.factor(x)) {
        intersect(levels(x), obs)
      } else if (!is.null(mobj$factor_levels[[col_nm]])) {
        intersect(as.character(mobj$factor_levels[[col_nm]]), obs)
      } else {
        obs
      }
    }

    # Check if marginalizing across other categorical factors under non-equal weights spans any 0-count cell
    slice_touches_empty_cell <- function(fixed_vals_named, label_desc) {
      if (identical(eff_weights, "equal") || length(cat_cols) <= 1L) return(FALSE)
      marginalized <- setdiff(cat_cols, names(fixed_vals_named))
      if (length(marginalized) == 0L) return(FALSE)

      full_levs <- stats::setNames(lapply(cat_cols, get_levels_ordered), cat_cols)
      full_grid <- expand.grid(full_levs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      if (length(fixed_vals_named) > 0L) {
        keep <- rep(TRUE, nrow(full_grid))
        for (fn in names(fixed_vals_named)) {
          if (fn %in% names(full_grid)) {
            keep <- keep & (full_grid[[fn]] == as.character(fixed_vals_named[[fn]]))
          }
        }
        full_grid <- full_grid[keep, , drop = FALSE]
      }
      if (nrow(full_grid) == 0L) return(TRUE)

      data_keys <- do.call(paste, c(lapply(data[, cat_cols, drop = FALSE], as.character), sep = "\r"))
      grid_keys <- do.call(paste, c(lapply(full_grid[, cat_cols, drop = FALSE], as.character), sep = "\r"))
      has_zero <- any(!grid_keys %in% data_keys)
      if (isTRUE(has_zero)) {
        warning(
          sprintf(
            "Dropping marginal 3dLMEr GLT '%s' under weights='%s' because one or more marginalized cells (%s) have zero observations.",
            label_desc, eff_weights, paste(marginalized, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      has_zero
    }

    build_cell_code <- function(var_vec, val_vec) {
      paste(paste0(var_vec, " : 1*", val_vec), collapse = " ")
    }

    # 1. Diagonal GLTs (quantitative slopes 'qVar :' and, if no other categorical GLTs requested, reference-vs-level diffs)
    if (isTRUE(spec$diagonal)) {
      for (qv in cont_cols) {
        add_generated(label = paste0("EV_", qv), code = paste0(qv, " :"), contrast_type = "diagonal")
      }
      has_other_cat_glts <- length(spec$cond_means) > 0L || isTRUE(spec$cell_means) ||
        length(spec$pairwise_diffs) > 0L || length(spec$interaction_diffs) > 0L
      if (!isTRUE(has_other_cat_glts) && length(cat_cols) > 0L) {
        for (cv in cat_cols) {
          levs <- get_levels_ordered(cv)
          if (length(levs) >= 2L) {
            ref_lev <- levs[1L]
            for (lv in levs[-1L]) {
              lbl <- paste0("EV_", cv, lv)
              if (slice_touches_empty_cell(stats::setNames(lv, cv), lbl) ||
                  slice_touches_empty_cell(stats::setNames(ref_lev, cv), lbl)) {
                next
              }
              add_generated(
                label = lbl,
                code = paste0(cv, " : 1*", lv, " -1*", ref_lev),
                contrast_type = "diagonal"
              )
            }
          }
        }
      }
    }

    # 2. Condition means (supports 'A', 'A:B', and 'A | B')
    if (length(spec$cond_means) > 0L) {
      for (vv in spec$cond_means) {
        ps <- parse_emm_spec_string(vv)
        all_v <- unique(c(ps$by_vars, ps$focal_vars))
        missing_v <- setdiff(all_v, names(data))
        if (length(missing_v) > 0L) {
          stop(sprintf("Cannot generate AFNI 3dLMEr GLTs for '%s' because variable(s) are not in the dataTable: %s.", vv, paste(missing_v, collapse = ", ")), call. = FALSE)
        }
        lev_grid <- expand.grid(stats::setNames(lapply(all_v, get_levels_ordered), all_v), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        # Keep only observed combinations when multiple factors are conditioned
        if (length(all_v) > 1L) {
          obs_keys <- unique(do.call(paste, c(lapply(data[, all_v, drop = FALSE], as.character), sep = "\r")))
          grid_keys <- do.call(paste, c(lapply(lev_grid[, all_v, drop = FALSE], as.character), sep = "\r"))
          lev_grid <- lev_grid[grid_keys %in% obs_keys, , drop = FALSE]
        }
        for (rr in seq_len(nrow(lev_grid))) {
          row_vals <- as.character(unlist(lev_grid[rr, all_v, drop = TRUE]))
          names(row_vals) <- all_v
          if (length(all_v) == 1L) {
            lbl <- paste(all_v, row_vals[1L], sep = ".")
          } else if (!is.null(ps$by_vars) && length(ps$by_vars) > 0L) {
            by_tag <- paste(ps$by_vars, row_vals[ps$by_vars], sep = ".", collapse = "_")
            foc_tag <- paste(paste(ps$focal_vars, collapse = "_"), paste(row_vals[ps$focal_vars], collapse = "."), sep = ".")
            lbl <- paste0(by_tag, "_", foc_tag)
          } else {
            lbl <- paste(paste(ps$focal_vars, collapse = "_"), paste(row_vals[ps$focal_vars], collapse = "."), sep = ".")
          }
          if (slice_touches_empty_cell(row_vals, lbl)) next
          add_generated(label = lbl, code = build_cell_code(all_v, row_vals), contrast_type = "cond_mean")
        }
      }
    }

    # 3. Cell means across all categorical factors
    if (isTRUE(spec$cell_means) && length(cat_cols) > 0L) {
      lev_grid <- expand.grid(stats::setNames(lapply(cat_cols, get_levels_ordered), cat_cols), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      obs_keys <- unique(do.call(paste, c(lapply(data[, cat_cols, drop = FALSE], as.character), sep = "\r")))
      grid_keys <- do.call(paste, c(lapply(lev_grid[, cat_cols, drop = FALSE], as.character), sep = "\r"))
      lev_grid <- lev_grid[grid_keys %in% obs_keys, , drop = FALSE]
      for (rr in seq_len(nrow(lev_grid))) {
        row_vals <- as.character(unlist(lev_grid[rr, cat_cols, drop = TRUE]))
        lbl <- paste(cat_cols, row_vals, sep = ".", collapse = "_")
        add_generated(label = lbl, code = build_cell_code(cat_cols, row_vals), contrast_type = "cell_mean")
      }
    }

    # 4. Overall average response across categorical factors
    if (isTRUE(spec$overall_response) && length(cat_cols) > 0L) {
      if (!slice_touches_empty_cell(character(0), "overall")) {
        parts <- vapply(cat_cols, function(cv) {
          levs <- get_levels_ordered(cv)
          if (identical(eff_weights, "equal")) {
            w <- rep(1 / length(levs), length(levs))
          } else {
            tab <- table(factor(as.character(data[[cv]]), levels = levs))
            w <- as.numeric(tab) / sum(tab)
          }
          paste0(cv, " : ", paste0(round(w, 6), "*", levs, collapse = " "))
        }, character(1))
        add_generated(label = "overall", code = paste(parts, collapse = " "), contrast_type = "overall")
      }
    }

    # 5. Pairwise differences (supports 'A' and conditioned simple effects 'A | B'; rejects all-cell 'A:B' with guidance)
    if (length(spec$pairwise_diffs) > 0L) {
      for (vv in spec$pairwise_diffs) {
        ps <- parse_emm_spec_string(vv)
        if (length(ps$focal_vars) > 1L) {
          stop(
            sprintf(
              paste0(
                "Automatic AFNI 3dLMEr GLT generation supports only one-factor pairwise differences (or conditioned simple effects like 'A | B'). ",
                "Contrast '%s'%s is an all-pairs interaction containing non-separable cross-diagonal comparisons; ",
                "use interaction_diffs = '%s' for factorial difference-of-differences, pairwise_diffs = 'A | B' for simple effects, or provide it explicitly with lmer_glt_codes."
              ),
              vv,
              if (!is.null(context)) paste0(" for ", context) else "",
              vv
            ),
            call. = FALSE
          )
        }
        fv <- ps$focal_vars[1L]
        all_v <- unique(c(ps$by_vars, fv))
        missing_v <- setdiff(all_v, names(data))
        if (length(missing_v) > 0L) {
          stop(sprintf("Cannot generate AFNI 3dLMEr GLTs for '%s' because variable(s) are not in the dataTable: %s.", vv, paste(missing_v, collapse = ", ")), call. = FALSE)
        }
        levs <- get_levels_ordered(fv)
        if (length(levs) < 2L) next
        pair_list <- utils::combn(levs, 2L, simplify = FALSE)

        if (is.null(ps$by_vars) || length(ps$by_vars) == 0L) {
          for (pp in pair_list) {
            lbl <- paste0(fv, ".pw.", pp[1L], ".vs.", pp[2L])
            if (slice_touches_empty_cell(stats::setNames(pp[1L], fv), lbl) ||
                slice_touches_empty_cell(stats::setNames(pp[2L], fv), lbl)) {
              next
            }
            add_generated(
              label = lbl,
              code = paste0(fv, " : 1*", pp[1L], " -1*", pp[2L]),
              contrast_type = "pairwise_diff"
            )
          }
        } else {
          by_grid <- expand.grid(stats::setNames(lapply(ps$by_vars, get_levels_ordered), ps$by_vars), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
          for (br in seq_len(nrow(by_grid))) {
            by_vals <- as.character(unlist(by_grid[br, ps$by_vars, drop = TRUE]))
            names(by_vals) <- ps$by_vars
            by_code <- build_cell_code(ps$by_vars, by_vals)
            by_tag <- paste(ps$by_vars, by_vals, sep = ".", collapse = "_")
            for (pp in pair_list) {
              # Check both cells exist in data for this `by` slice
              sub_mask1 <- rep(TRUE, nrow(data))
              sub_mask2 <- rep(TRUE, nrow(data))
              for (bn in ps$by_vars) {
                sub_mask1 <- sub_mask1 & (as.character(data[[bn]]) == by_vals[[bn]])
                sub_mask2 <- sub_mask2 & (as.character(data[[bn]]) == by_vals[[bn]])
              }
              sub_mask1 <- sub_mask1 & (as.character(data[[fv]]) == pp[1L])
              sub_mask2 <- sub_mask2 & (as.character(data[[fv]]) == pp[2L])
              if (!any(sub_mask1, na.rm = TRUE) || !any(sub_mask2, na.rm = TRUE)) next

              lbl <- paste0(by_tag, "_", fv, ".pw.", pp[1L], ".vs.", pp[2L])
              if (slice_touches_empty_cell(c(by_vals, stats::setNames(pp[1L], fv)), lbl) ||
                  slice_touches_empty_cell(c(by_vals, stats::setNames(pp[2L], fv)), lbl)) {
                next
              }
              add_generated(
                label = lbl,
                code = paste0(by_code, " ", fv, " : 1*", pp[1L], " -1*", pp[2L]),
                contrast_type = "pairwise_diff"
              )
            }
          }
        }
      }
    }

    # 6. Factorial interaction (difference-of-differences) contrasts: 'Group : 1*A -1*B Session : 1*pre -1*post'
    if (length(spec$interaction_diffs) > 0L) {
      for (vv in spec$interaction_diffs) {
        ps <- parse_emm_spec_string(vv)
        if (length(ps$focal_vars) < 2L) next
        all_v <- unique(c(ps$by_vars, ps$focal_vars))
        missing_v <- setdiff(all_v, names(data))
        if (length(missing_v) > 0L) {
          stop(sprintf("Cannot generate AFNI 3dLMEr interaction GLTs for '%s' because variable(s) are not in dataTable: %s.", vv, paste(missing_v, collapse = ", ")), call. = FALSE)
        }
        fac_pairs <- lapply(ps$focal_vars, function(fv) {
          levs <- get_levels_ordered(fv)
          if (length(levs) < 2L) return(list())
          utils::combn(levs, 2L, simplify = FALSE)
        })
        if (any(lengths(fac_pairs) == 0L)) next

        idx_grid <- expand.grid(lapply(fac_pairs, seq_along), KEEP.OUT.ATTRS = FALSE)
        by_grid <- if (!is.null(ps$by_vars) && length(ps$by_vars) > 0L) {
          expand.grid(stats::setNames(lapply(ps$by_vars, get_levels_ordered), ps$by_vars), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        } else {
          data.frame(.dummy = 1L)
        }

        for (br in seq_len(nrow(by_grid))) {
          by_prefix_lbl <- ""
          by_prefix_code <- ""
          by_vals <- character(0)
          if (!is.null(ps$by_vars) && length(ps$by_vars) > 0L) {
            by_vals <- as.character(unlist(by_grid[br, ps$by_vars, drop = TRUE]))
            names(by_vals) <- ps$by_vars
            by_prefix_lbl <- paste0(paste(ps$by_vars, by_vals, sep = ".", collapse = "_"), "_")
            by_prefix_code <- paste0(build_cell_code(ps$by_vars, by_vals), " ")
          }

          for (ir in seq_len(nrow(idx_grid))) {
            chosen_pairs <- lapply(seq_along(ps$focal_vars), function(k) fac_pairs[[k]][[idx_grid[ir, k]]])
            names(chosen_pairs) <- ps$focal_vars

            # Verify all involved cells exist in data
            cell_grid <- expand.grid(c(as.list(by_vals), chosen_pairs), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
            obs_keys <- unique(do.call(paste, c(lapply(data[, colnames(cell_grid), drop = FALSE], as.character), sep = "\r")))
            c_keys <- do.call(paste, c(lapply(cell_grid, as.character), sep = "\r"))
            if (!all(c_keys %in% obs_keys)) next

            code_blocks <- vapply(seq_along(ps$focal_vars), function(k) {
              fv <- ps$focal_vars[k]
              pp <- chosen_pairs[[k]]
              paste0(fv, " : 1*", pp[1L], " -1*", pp[2L])
            }, character(1))
            lbl_blocks <- vapply(chosen_pairs, function(pp) paste0(pp[1L], "_M_", pp[2L]), character(1))
            lbl <- paste0(by_prefix_lbl, paste(ps$focal_vars, collapse = "_x_"), ".", paste(lbl_blocks, collapse = "_x_"))

            any_empty_slice <- any(vapply(seq_len(nrow(cell_grid)), function(cr) {
              rv <- as.character(unlist(cell_grid[cr, , drop = TRUE]))
              names(rv) <- colnames(cell_grid)
              slice_touches_empty_cell(rv, lbl)
            }, logical(1)))
            if (isTRUE(any_empty_slice)) next

            add_generated(
              label = lbl,
              code = paste0(by_prefix_code, paste(code_blocks, collapse = " ")),
              contrast_type = "interaction_diff"
            )
          }
        }
      }
    }

    # 7. Simple slopes ('Group : 1*A Age :') and slope differences ('Group : 1*A -1*B Age :')
    if (length(spec$simple_slopes) > 0L) {
      for (vv in seq_along(spec$simple_slopes)) {
        trend_var <- names(spec$simple_slopes)[vv]
        if (!trend_var %in% names(data)) {
          stop(sprintf("Cannot generate AFNI 3dLMEr simple_slopes for '%s' because it is not in dataTable.", trend_var), call. = FALSE)
        }
        comb_list <- spec$simple_slopes[[vv]]
        if (is.character(comb_list)) comb_list <- as.list(comb_list)
        for (comb in comb_list) {
          comb_vars <- trimws(unlist(strsplit(as.character(comb), "[:*+]", perl = TRUE)))
          comb_vars <- comb_vars[nzchar(comb_vars)]
          lev_grid <- expand.grid(stats::setNames(lapply(comb_vars, get_levels_ordered), comb_vars), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
          obs_keys <- unique(do.call(paste, c(lapply(data[, comb_vars, drop = FALSE], as.character), sep = "\r")))
          g_keys <- do.call(paste, c(lapply(lev_grid[, comb_vars, drop = FALSE], as.character), sep = "\r"))
          lev_grid <- lev_grid[g_keys %in% obs_keys, , drop = FALSE]
          for (rr in seq_len(nrow(lev_grid))) {
            rv <- as.character(unlist(lev_grid[rr, comb_vars, drop = TRUE]))
            names(rv) <- comb_vars
            lbl <- paste(trend_var, paste(rv, collapse = "."), sep = ".")
            if (slice_touches_empty_cell(rv, lbl)) next
            add_generated(
              label = lbl,
              code = paste0(build_cell_code(comb_vars, rv), " ", trend_var, " :"),
              contrast_type = "simple_slope"
            )
          }
        }
      }
    }

    if (length(spec$slope_diffs) > 0L) {
      for (vv in seq_along(spec$slope_diffs)) {
        trend_var <- names(spec$slope_diffs)[vv]
        if (!trend_var %in% names(data)) {
          stop(sprintf("Cannot generate AFNI 3dLMEr slope_diffs for '%s' because it is not in dataTable.", trend_var), call. = FALSE)
        }
        comb_list <- spec$slope_diffs[[vv]]
        if (is.character(comb_list)) comb_list <- as.list(comb_list)
        for (comb in comb_list) {
          ps <- parse_emm_spec_string(paste(as.character(comb), collapse = ":"))
          fac_pairs <- lapply(ps$focal_vars, function(fv) utils::combn(get_levels_ordered(fv), 2L, simplify = FALSE))
          if (any(lengths(fac_pairs) == 0L)) next
          idx_grid <- expand.grid(lapply(fac_pairs, seq_along), KEEP.OUT.ATTRS = FALSE)
          for (ir in seq_len(nrow(idx_grid))) {
            chosen <- lapply(seq_along(ps$focal_vars), function(k) fac_pairs[[k]][[idx_grid[ir, k]]])
            code_blocks <- vapply(seq_along(ps$focal_vars), function(k) {
              paste0(ps$focal_vars[k], " : 1*", chosen[[k]][1L], " -1*", chosen[[k]][2L])
            }, character(1))
            lbl_blocks <- vapply(chosen, function(pp) paste0(pp[1L], "_M_", pp[2L]), character(1))
            lbl <- paste(trend_var, paste(ps$focal_vars, collapse = "_x_"), paste(lbl_blocks, collapse = "_x_"), sep = ".")
            add_generated(
              label = lbl,
              code = paste0(paste(code_blocks, collapse = " "), " ", trend_var, " :"),
              contrast_type = "slope_diff"
            )
          }
        }
      }
    }

    # 8. Custom contrasts (preserved in contrast_spec$custom or contrast_list$custom)
    cmat <- spec$custom %||% mobj$contrast_list$custom
    if (!is.null(cmat)) {
      cmat <- parse_custom_contrast_spec(cmat, reg_names = names(stats::coef(mobj$lmfit)))
    }
    if (!is.null(cmat) && is.matrix(cmat) && nrow(cmat) > 0L) {
      for (i in seq_len(nrow(cmat))) {
        con_name <- rownames(cmat)[i]
        weights <- cmat[i, ]
        non_zero <- which(abs(weights) > 1e-8)
        vars <- names(weights)[non_zero]
        if (!all(vars %in% names(data))) {
          stop(
            sprintf(
              "Cannot convert custom contrast '%s' to AFNI 3dLMEr GLT because regressor(s) are not direct dataTable variables: %s. Provide lmer_glt_codes manually.",
              con_name, paste(setdiff(vars, names(data)), collapse = ", ")
            ),
            call. = FALSE
          )
        }
        code_parts <- vapply(seq_along(vars), function(vi) {
          vn <- vars[vi]
          wv <- as.numeric(weights[non_zero[vi]])
          if ((vn %in% cont_cols) && abs(wv - 1) <= 1e-8) {
            paste0(vn, " :")
          } else {
            paste0(vn, " : ", wv)
          }
        }, character(1))
        add_generated(con_name, paste(code_parts, collapse = " "), "custom")
      }
    }

    # 9. Apply spec$delete to generated GLTs
    if (!is.null(spec$delete) && is.character(spec$delete) && length(spec$delete) > 0L && nrow(generated) > 0L) {
      del_sanitized <- sanitize_3dlmer_glt_names(spec$delete)
      gen_sanitized <- sanitize_3dlmer_glt_names(generated$label_raw)
      keep_gen <- !(generated$label_raw %in% spec$delete | gen_sanitized %in% del_sanitized)
      generated <- generated[keep_gen, , drop = FALSE]
    }
  }

  raw_glt_codes <- raw_glt_codes %||% mobj$lmer_glt_codes
  raw <- as_3dlmer_glt_table(raw_glt_codes, source = "user", contrast_type = "raw")
  out <- rbind(generated, raw)
  if (nrow(out) == 0L) {
    return(out)
  }
  out$label_afni <- sanitize_3dlmer_glt_names(out$label_raw)

  validate_3dlmer_glt_table(out, datatable = data, qVars = qVars, context = context)
}
