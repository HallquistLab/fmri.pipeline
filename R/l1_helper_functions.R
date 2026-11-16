# helper function to return a logical vector of trials included in subset
get_trial_set_from_signal <- function(ss, trial_data=NULL) {
  checkmate::assert_data_frame(trial_data)
  # if user has a cached trial expression, need to evaluate the trial set before getting value df
  if (isTRUE(ss$trial_subset)) {
    # calculate trial set from cached expression
    trial_set <- with(trial_data, eval(parse(text = ss$trial_subset_expression)))
  } else {
    # when FALSE: keep all trials. This is a local variable that is calculated every time this function is called (e.g., modification)
    trial_set <- rep(TRUE, nrow(trial_data))
  }
  return(trial_set)
}

# helper function to compute statistics for how much data is excluded by trial subset
get_trial_subset_stats <- function(signal, trial_data, trial_set=NULL) {
  if (is.null(trial_set)) {
    trial_set <- get_trial_set_from_signal(signal, trial_data)
  }

  overall <- sum(trial_set == TRUE) / length(trial_set)
  tmp <- trial_data %>% bind_cols(trial_set = trial_set)
  by_id <- tmp %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(pct_true = sum(trial_set == TRUE) / n(), .groups = "drop") %>%
    dplyr::summarize(
      mean = mean(pct_true, na.rm = T), sd = stats::sd(pct_true, na.rm = T),
      min = min(pct_true, na.rm = T), max = max(pct_true, na.rm = T)
    ) %>%
    mutate(overall = overall) %>%
      unlist()

  signal$trial_subset_statistics <- by_id
  return(signal)
}

# helper function to figure out expected regressor columns for a given signal based on whether
# the signal has within-subject factors, is a beta series signal, and/or includes a temporal derivative
get_regressors_from_signal <- function(sig) {
  if (!is.null(sig$wi_factors) && is.null(sig$wi_model)) {
    sig <- fit_wi_model(sig)
  }

  # in terms of design, always add derivative columns en bloc after the corresponding non-derivative columns
  # handle expansion of PPI signals slightly differently -- expand w/i factors within psych and PPI parts if needed
  if (checkmate::test_string(sig$ts_multiplier)) {
    if (!is.null(sig$wi_factors)) {
      cols <- paste(sig$name, names(stats::coef(sig$wi_model)), sep = ".")
      cols <- c(cols, sig$ts_multiplier, sub(paste0("^", sig$name), paste(sig$name, "ppi", sep = "."), cols))
    } else {
      cols <- c(sig$name, sig$ts_multiplier, paste(sig$name, "ppi", sep = "."))
    }
  } else if (!is.null(sig$wi_factors)) {
    # use the lmfit object in the signal to determine the columns that will be included
    # for within-subject factor modulation, always include the signal as the prefix on the column names
    cols <- paste(sig$name, names(stats::coef(sig$wi_model)), sep = ".")
  } else if (isTRUE(sig$beta_series)) {
    if (is.data.frame(sig$value)) {
      trials <- sort(unique(sig$value$trial)) # vector of trials for parametric signal
    } else {
      trials <- sort(unique(event_list[[sig$event]]$trial))
    }

    cols <- paste(sig$name, sprintf("%03d", trials), sep = "_t") # signal_t001 etc.
  } else {
    cols <- sig$name # nothing special
  }

  if (isTRUE(sig$add_deriv)) {
    cols <- c(cols, paste0("d_", cols)) # add temporal derivative for each column
  }
  return(cols)
}

build_wi_trial_matrix <- function(df_trials, wi_ev_formula, wi_factors, wi_factor_levels, factor_coding = NULL) {
  ev_rhs <- stats::delete.response(stats::terms(stats::as.formula(wi_ev_formula)))
  cat_wi <- intersect(wi_factors, names(wi_factor_levels %||% list()))
  cont_wi <- setdiff(wi_factors, cat_wi)

  if (length(cat_wi) > 0L) {
    grid_df <- expand.grid(wi_factor_levels[cat_wi], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    for (fv in cat_wi) {
      grid_df[[fv]] <- factor(grid_df[[fv]], levels = wi_factor_levels[[fv]])
      if (!is.null(factor_coding) && length(factor_coding) > 0L) {
        code_this <- if (!is.null(names(factor_coding)) && fv %in% names(factor_coding)) {
          factor_coding[[fv]]
        } else if (is.null(names(factor_coding)) && length(factor_coding) == 1L) {
          factor_coding[[1L]]
        } else {
          NULL
        }
        grid_df[[fv]] <- apply_factor_coding_to_col(grid_df[[fv]], code_this)
      }
    }
  } else {
    grid_df <- data.frame(.row_placeholder = 1L)
  }
  for (cv in cont_wi) {
    grid_df[[cv]] <- 0
  }
  grid_df$.row_placeholder <- NULL
  M_grid <- stats::model.matrix(ev_rhs, data = grid_df)

  if (nrow(df_trials) == 0L) {
    out <- matrix(numeric(0), nrow = 0L, ncol = ncol(M_grid), dimnames = list(NULL, colnames(M_grid)))
    return(out)
  }

  fac_only <- df_trials[, wi_factors, drop = FALSE]
  for (fv in cat_wi) {
    fac_only[[fv]] <- factor(as.character(fac_only[[fv]]), levels = wi_factor_levels[[fv]])
  }
  for (cv in cont_wi) {
    fac_only[[cv]] <- as.numeric(fac_only[[cv]])
  }
  combined <- rbind(fac_only, grid_df[, wi_factors, drop = FALSE])
  for (fv in cat_wi) {
    combined[[fv]] <- factor(as.character(combined[[fv]]), levels = wi_factor_levels[[fv]])
    if (!is.null(factor_coding) && length(factor_coding) > 0L) {
      code_this <- if (!is.null(names(factor_coding)) && fv %in% names(factor_coding)) {
        factor_coding[[fv]]
      } else if (is.null(names(factor_coding)) && length(factor_coding) == 1L) {
        factor_coding[[1L]]
      } else {
        NULL
      }
      combined[[fv]] <- apply_factor_coding_to_col(combined[[fv]], code_this)
    }
  }

  M_all <- stats::model.matrix(ev_rhs, data = combined)
  M_all[seq_len(nrow(df_trials)), , drop = FALSE]
}

# cf. get_regressors_from_signal -- these have to align!
expand_signal <- function(sig) {
  if (checkmate::test_string(sig$ts_multiplier) && is.null(sig$ts_expanded)) {
    sig$ts_expanded <- TRUE
    sig_psych <- sig
    sig_psych$ts_multiplier <- NULL
    sig_psych <- expand_signal(sig_psych)

    sig_physio <- sig
    sig_physio$name <- sig_physio$ts_multiplier
    sig_physio$physio_only <- TRUE
    sig_physio <- list(sig_physio)
    names(sig_physio) <- sig_physio[[1]]$name

    sig_ppi <- sig
    sig_ppi$name <- paste(sig_ppi$name, "ppi", sep = ".")
    sig_ppi <- expand_signal(sig_ppi)

    s_list <- c(sig_psych, sig_physio, sig_ppi)
  } else if (is.null(sig$wi_factors) && isFALSE(sig$beta_series)) {
    sig$is_parametric <- identical(sig$value_type, "parametric")
    s_list <- list(sig)
    names(s_list) <- sig$name
  } else if (!is.null(sig$wi_factors) && isTRUE(sig$beta_series)) {
    stop("In expand_signal, cannot sort out what to do when beta series is enabled and there is a within-subject factor!")
  } else if (!is.null(sig$wi_factors)) {
    if (is.null(sig$wi_model) || is.null(sig$wi_ev_formula)) {
      sig <- fit_wi_model(sig)
    }

    checkmate::assert_data_frame(sig$value)
    miss_fac <- if (nrow(sig$value) > 0L) {
      !stats::complete.cases(sig$value[, sig$wi_factors, drop = FALSE])
    } else {
      logical(0)
    }
    if (any(miss_fac)) {
      warning(
        sprintf(
          "Signal '%s' has %d trial(s) with missing (NA) within-subject factor values (%s); omitting those trials before expanding EVs.",
          sig$name, sum(miss_fac), paste(sig$wi_factors, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    val_clean <- sig$value[!miss_fac, , drop = FALSE]

    model_mat <- build_wi_trial_matrix(
      df_trials = val_clean,
      wi_ev_formula = sig$wi_ev_formula,
      wi_factors = sig$wi_factors,
      wi_factor_levels = sig$wi_factor_levels,
      factor_coding = sig$factor_coding
    )

    coefs <- colnames(model_mat)
    exp_names <- paste(sig$name, coefs, sep = ".")
    is_formula_mode <- identical(sig$wi_expansion_resolved, "formula") &&
      !(all(abs(model_mat) < 1e-8 | abs(model_mat - 1) < 1e-8) && all(rowSums(abs(model_mat)) <= 1 + 1e-8))
    is_param <- identical(sig$value_type, "parametric")

    # Record per-run trial-level model matrices for estimability diagnostics in formula mode
    run_trial_mats <- NULL
    if (identical(sig$wi_expansion_resolved, "formula") && "run_number" %in% names(val_clean)) {
      run_nums <- sort(unique(val_clean$run_number))
      run_trial_mats <- stats::setNames(lapply(run_nums, function(rn) {
        r_mask <- val_clean$run_number == rn
        M_r <- model_mat[r_mask, , drop = FALSE]
        colnames(M_r) <- exp_names
        if ("trial" %in% names(val_clean)) {
          attr(M_r, "trial_ids") <- val_clean$trial[r_mask]
        }
        M_r
      }), paste0("run_number", run_nums))
    }

    base_val_df <- val_clean %>% dplyr::select(-dplyr::all_of(sig$wi_factors))

    s_list <- lapply(seq_along(exp_names), function(ee) {
      sig_copy <- sig
      sig_copy$name <- exp_names[ee]
      sig_copy$is_parametric <- is_param
      sig_copy$parent_signal_name <- sig$name
      sig_copy$wi_expansion_resolved <- sig$wi_expansion_resolved
      if (ee == 1L && !is.null(run_trial_mats)) {
        sig_copy$wi_run_trial_matrices <- run_trial_mats
      }

      w <- if (nrow(model_mat) > 0L) as.numeric(model_mat[, ee]) else numeric(0)
      nz <- which(abs(w) > 1e-8)
      sub_df <- base_val_df[nz, , drop = FALSE]
      if (nrow(sub_df) > 0L) {
        if (isTRUE(is_formula_mode)) {
          # Store raw modulator and contrast weight separately so build_design_matrix can center
          # parametric values per run BEFORE multiplying by signed/overlapping formula weights.
          sub_df$wi_weight <- w[nz]
          if (!isTRUE(is_param)) {
            sub_df$value <- sub_df$value * w[nz]
          }
          sig_copy$formula_weighted <- TRUE
        } else {
          sub_df$value <- sub_df$value * w[nz]
          sig_copy$formula_weighted <- FALSE
        }
      }
      sig_copy$value <- sub_df

      # remove residual within-subject modulation settings
      sig_copy$wi_formula <- sig_copy$wi_ev_formula <- sig_copy$wi_model <- sig_copy$wi_factors <- NULL
      return(sig_copy)
    })
  } else if (isTRUE(sig$beta_series)) {
    checkmate::assert_data_frame(sig$value)
    if (!"trial" %in% names(sig$value)) {
      stop("Beta-series signals require a trial column in their value data.")
    }
    if (anyNA(sig$value$trial)) {
      stop("Beta-series signals require non-missing trial identifiers.")
    }
    trials <- sort(unique(sig$value$trial))
    exp_names <- paste(sig$name, sprintf("%03d", trials), sep = "_t")
    s_list <- lapply(seq_along(exp_names), function(ee) {
      sig_copy <- sig
      sig_copy$name <- exp_names[ee]
      sig_copy$value <- sig$value %>% filter(trial == !!trials[ee])
      sig_copy$beta_series <- FALSE
      sig_copy$beta_series_trial <- trials[ee]
      sig_copy$is_parametric <- FALSE
      return(sig_copy)
    })
  }

  # create derivative signals en bloc (marked with return_deriv = TRUE so place_dmat_on_time_grid
  # replaces them with orthogonalized first-difference temporal derivatives rather than duplicates)
  if (isTRUE(sig$add_deriv)) {
    s_list <- lapply(s_list, function(ss) {
      ss$add_deriv <- FALSE
      return(ss)
    })
    s_list_deriv <- lapply(s_list, function(ss) {
      ss$name <- paste0("d_", ss$name)
      ss$add_deriv <- FALSE
      ss$return_deriv <- TRUE
      return(ss)
    })
    s_list <- c(s_list, s_list_deriv)
  }

  return(s_list)
}

fit_wi_model <- function(sobj) {
  if (is.null(sobj$wi_factors) || length(sobj$wi_factors) == 0L) {
    return(sobj)
  }

  if (is.null(sobj$wi_formula) || !nzchar(trimws(as.character(sobj$wi_formula)[1L]))) {
    sobj$wi_formula <- paste("~ 0 +", paste(sobj$wi_factors, collapse = ":"))
  } else {
    wi_form_obj <- stats::as.formula(sobj$wi_formula)
    wi_exp_mode <- tolower(sobj$wi_expansion %||% "auto")
    if (!identical(wi_exp_mode, "formula") || identical(attr(stats::terms(wi_form_obj), "intercept"), 0L)) {
      wi_form_obj <- stats::update.formula(wi_form_obj, ~ . - 1)
    }
    sobj$wi_formula <- as.character.formula(wi_form_obj)
  }

  checkmate::assert_data_frame(sobj$value)
  checkmate::assert_subset(sobj$wi_factors, names(sobj$value))

  # Establish or enforce canonical factor levels across subjects for categorical wi_factors
  # while preserving numeric continuous modulators as numeric predictors.
  if (is.null(sobj$wi_factor_levels)) {
    lev_map <- list()
    for (fv in sobj$wi_factors) {
      col_v <- sobj$value[[fv]]
      force_cat <- (!is.null(sobj$reference_level) && fv %in% names(sobj$reference_level)) ||
        (!is.null(sobj$factor_coding) && fv %in% names(sobj$factor_coding))
      if (is.numeric(col_v) && !is.factor(col_v) && !isTRUE(force_cat)) {
        next
      }
      raw_f <- factor(col_v)
      if (nlevels(raw_f) < 2L) {
        stop(
          sprintf(
            "Within-subject factor '%s' on signal '%s' must have at least 2 non-missing levels in trial_data (found: %s).",
            fv, sobj$name %||% "<unnamed>", paste(levels(raw_f), collapse = ", ")
          ),
          call. = FALSE
        )
      }
      if (!is.null(sobj$reference_level) && fv %in% names(sobj$reference_level)) {
        ref_val <- as.character(sobj$reference_level[[fv]])
        if (ref_val %in% levels(raw_f)) {
          raw_f <- stats::relevel(raw_f, ref = ref_val)
        } else {
          warning(
            sprintf("Requested reference_level '%s' for wi_factor '%s' not found; keeping '%s'.", ref_val, fv, levels(raw_f)[1L]),
            call. = FALSE
          )
        }
      }
      lev_map[[fv]] <- levels(raw_f)
    }
    sobj$wi_factor_levels <- lev_map
  }

  cat_wi <- intersect(sobj$wi_factors, names(sobj$wi_factor_levels))
  cont_wi <- setdiff(sobj$wi_factors, cat_wi)

  wi_exp <- tolower(trimws(as.character(sobj$wi_expansion %||% "auto")[1L]))
  checkmate::assert_subset(wi_exp, c("auto", "cell_means", "formula"))
  sobj$wi_expansion <- wi_exp

  if (length(cont_wi) > 0L && identical(wi_exp, "cell_means")) {
    stop(
      sprintf(
        "Signal '%s' requested wi_expansion='cell_means', which applies only to categorical wi_factors, but continuous predictor(s) are present: %s.",
        sobj$name %||% "<unnamed>", paste(cont_wi, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Determine whether wi_formula is saturated over the Cartesian categorical factor grid
  if (length(cat_wi) > 0L) {
    grid_df <- expand.grid(sobj$wi_factor_levels[cat_wi], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    for (fv in cat_wi) {
      grid_df[[fv]] <- factor(grid_df[[fv]], levels = sobj$wi_factor_levels[[fv]])
      if (!is.null(sobj$factor_coding) && length(sobj$factor_coding) > 0L) {
        code_this <- if (!is.null(names(sobj$factor_coding)) && fv %in% names(sobj$factor_coding)) {
          sobj$factor_coding[[fv]]
        } else if (is.null(names(sobj$factor_coding)) && length(sobj$factor_coding) == 1L) {
          sobj$factor_coding[[1L]]
        } else {
          NULL
        }
        grid_df[[fv]] <- apply_factor_coding_to_col(grid_df[[fv]], code_this)
      }
    }
  } else {
    grid_df <- data.frame(.row_placeholder = 1L)
  }
  for (cv in cont_wi) grid_df[[cv]] <- 0
  grid_df$.row_placeholder <- NULL

  is_saturated <- FALSE
  if (length(cont_wi) == 0L && length(cat_wi) > 0L) {
    M_grid <- stats::model.matrix(stats::as.formula(sobj$wi_formula), data = grid_df)
    is_saturated <- (ncol(M_grid) == nrow(grid_df)) && (qr(M_grid)$rank == nrow(grid_df))
  }

  if (identical(wi_exp, "cell_means") || (identical(wi_exp, "auto") && isTRUE(is_saturated) && length(cat_wi) > 1L && length(cont_wi) == 0L)) {
    sobj$wi_expansion_resolved <- "cell_means"
    sobj$wi_ev_formula <- paste("~ 0 +", paste(cat_wi, collapse = ":"))
  } else {
    sobj$wi_expansion_resolved <- "formula"
    sobj$wi_ev_formula <- sobj$wi_formula
  }

  complete_mask <- if (nrow(sobj$value) > 0L) {
    stats::complete.cases(sobj$value[, sobj$wi_factors, drop = FALSE])
  } else {
    logical(0)
  }
  wi_df <- sobj$value[complete_mask, sobj$wi_factors, drop = FALSE]
  if (nrow(wi_df) == 0L) {
    wi_df <- grid_df[, sobj$wi_factors, drop = FALSE]
  }
  for (fv in cat_wi) {
    wi_df[[fv]] <- factor(as.character(wi_df[[fv]]), levels = sobj$wi_factor_levels[[fv]])
    if (!is.null(sobj$factor_coding) && length(sobj$factor_coding) > 0L) {
      code_this <- if (!is.null(names(sobj$factor_coding)) && fv %in% names(sobj$factor_coding)) {
        sobj$factor_coding[[fv]]
      } else if (is.null(names(sobj$factor_coding)) && length(sobj$factor_coding) == 1L) {
        sobj$factor_coding[[1L]]
      } else {
        NULL
      }
      wi_df[[fv]] <- apply_factor_coding_to_col(wi_df[[fv]], code_this)
    }
  }
  for (cv in cont_wi) {
    wi_df[[cv]] <- as.numeric(wi_df[[cv]])
  }
  wi_df$dummy <- seq_len(nrow(wi_df))

  ffit <- stats::update.formula(stats::as.formula(sobj$wi_ev_formula), dummy ~ .)
  fit_res <- fit_lm_with_canonical_levels(
    model_formula = ffit,
    model_data = wi_df,
    canonical_levels = sobj$wi_factor_levels,
    factor_coding = sobj$factor_coding
  )
  sobj$wi_model <- fit_res$lmfit
  return(sobj)
}
  
#' S3 method to convert formula to character
#' @param x a formula to convert to character
#' @param ... not used
#' @details copied from formula.tools to avoid NAMESPACE import
#' @export
as.character.formula <- function(x, ...) { 
  form <- paste( deparse(x), collapse=" " )
  form <- gsub( "\\s+", " ", form, perl=FALSE ) # remove multiple spaces
  return(form)
}

#obsolete
# # helper function to expand beta series
# expand_signals <- function(orig) {
#   # use levels of wi_factors here to build out multiple signals, one per level... maybe just copy-paste these inside object, subsetting value?
#   signals_expanded <- lapply(orig, function(ss) {
#     if (is.null(ss$wi_factors)) {
#       return(ss)
#     } else {
#       # extract within-subject regression model matrix (has dummy coding)
#       model_df <- as.data.frame(model.matrix(ss$wi_model))

#       # for each column of the model matrix, pull the signal value data.frame where
#       # the dummy code is 1 (i.e., rows to which the dummy pertains). This divides up the value
#       # data.frame into smaller data.frames, one per dummy combination
#       dummy_value <- lapply(model_df, function(x) {
#         which_1 <- which(x == 1)
#         ss$value %>%
#           dplyr::slice(!!which_1) %>%
#           dplyr::select(-dplyr::all_of(ss$wi_factors))
#       })

#       # create a list of signals, one per dummy level, by copying the master signal
#       # settings, then amending the value field to include only relevant trials/rows
#       signal_list <- lapply(seq_along(dummy_value), function(x) {
#         sig_copy <- ss
#         # update signal name to reflect level of within-subject factor(s)
#         sig_copy$name <- paste(sig_copy$name, make.names(names(dummy_value)[x]), sep = "_")

#         # copy across relevant data.frame for this level
#         sig_copy$value <- dummy_value[[x]]

#         # remove residual within-subject modulation settings
#         sig_copy$wi_formula <- sig_copy$wi_model <- sig_copy$wi_factors <- NULL
#         sig_copy
#       })
#     }
#   })

#   # now handle beta series
# }
