extract_lm_var_types <- function(lmfit) {
  if (is.null(lmfit)) {
    return(list(cat_vars = character(0), cont_vars = character(0)))
  }
  tt <- stats::terms(lmfit)
  rhs_vars <- all.vars(stats::delete.response(tt))
  dclass <- attr(tt, "dataClasses")
  if (is.null(dclass) || length(rhs_vars) == 0L) {
    return(list(cat_vars = character(0), cont_vars = character(0)))
  }
  present <- intersect(rhs_vars, names(dclass))
  cat_vars <- present[dclass[present] %in% c("factor", "character", "ordered")]
  cont_vars <- present[dclass[present] %in% c("integer", "numeric")]
  list(cat_vars = cat_vars, cont_vars = cont_vars)
}

coerce_yaml_char_vec <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.list(x)) {
    if (length(x) == 0L) return(character(0))
    x <- unlist(x, recursive = TRUE, use.names = FALSE)
  }
  as.character(x)
}

parse_custom_contrast_spec <- function(custom_spec, reg_names = NULL, lg = NULL) {
  if (is.null(custom_spec) || length(custom_spec) == 0L) return(NULL)
  if (is.null(lg)) lg <- lgr::get_logger()

  if (is.matrix(custom_spec)) {
    if (nrow(custom_spec) == 0L) return(NULL)
    storage.mode(custom_spec) <- "double"
    if (is.null(colnames(custom_spec)) && !is.null(reg_names) && ncol(custom_spec) == length(reg_names)) {
      colnames(custom_spec) <- reg_names
    }
    return(custom_spec)
  }

  if (is.data.frame(custom_spec)) {
    return(parse_custom_contrast_spec(as.matrix(custom_spec), reg_names = reg_names, lg = lg))
  }

  if (!is.list(custom_spec)) {
    stop("custom contrasts must be a matrix or a named list of contrast vectors.", call. = FALSE)
  }

  cnames <- names(custom_spec)
  if (is.null(cnames)) {
    cnames <- paste0("custom_", seq_along(custom_spec))
  } else {
    empty_nm <- !nzchar(cnames) | is.na(cnames)
    if (any(empty_nm)) cnames[empty_nm] <- paste0("custom_", which(empty_nm))
  }

  # Collect union of named columns if reg_names is not yet known
  if (is.null(reg_names)) {
    named_cols <- unique(unlist(lapply(custom_spec, function(entry) {
      if (is.list(entry)) names(entry) else names(entry)
    })))
    if (length(named_cols) > 0L) reg_names <- named_cols
  }

  rows <- list()
  kept_names <- character(0)
  for (ii in seq_along(custom_spec)) {
    entry <- custom_spec[[ii]]
    if (is.list(entry) && !is.null(names(entry))) {
      entry <- unlist(entry)
    }
    if (!is.numeric(entry)) {
      entry_num <- suppressWarnings(as.numeric(entry))
      if (!is.null(names(entry))) names(entry_num) <- names(entry)
      entry <- entry_num
    }

    if (!is.null(names(entry)) && any(nzchar(names(entry)))) {
      if (is.null(reg_names)) {
        cvec <- entry
      } else {
        unknown_cols <- setdiff(names(entry), reg_names)
        if (length(unknown_cols) > 0L) {
          lg$warn(
            "Custom contrast '%s' references unknown/dropped regressor(s) (%s); dropping this custom contrast.",
            cnames[ii], paste(unknown_cols, collapse = ", ")
          )
          next
        }
        cvec <- stats::setNames(rep(0, length(reg_names)), reg_names)
        cvec[names(entry)] <- as.numeric(entry)
      }
    } else {
      if (!is.null(reg_names) && length(entry) != length(reg_names)) {
        lg$warn(
          "Custom contrast '%s' has length %d, expected %d (%s); dropping this custom contrast.",
          cnames[ii], length(entry), length(reg_names), paste(reg_names, collapse = ", ")
        )
        next
      }
      cvec <- as.numeric(entry)
      if (!is.null(reg_names)) names(cvec) <- reg_names
    }

    rows[[length(rows) + 1L]] <- cvec
    kept_names <- c(kept_names, cnames[ii])
  }

  if (length(rows) == 0L) return(NULL)
  out <- do.call(rbind, rows)
  rownames(out) <- kept_names
  if (!is.null(reg_names) && ncol(out) == length(reg_names)) {
    colnames(out) <- reg_names
  }
  out
}

#' Internal function to add, edit, rename, delete contrasts
#'
#' @param mobj \code{l1_model_spec} or \code{hi_model_spec} object
#' @param signals a list of signals that provide information about regressor specification
#' @param spec_list a list of settings for this model based on an external YAML/JSON specification file
#' @param preselected_pairwise Optional character vector of factor names preselected for pairwise differences
#'   in interactive higher-level model setup.
#' @keywords internal
#' @author Michael Hallquist
#' @importFrom magrittr %>%
#' @importFrom emmeans emmeans emtrends
#' @importFrom lgr get_logger
#' @importFrom checkmate assert_multi_class assert_integerish
#' @importFrom stats coef model.matrix setNames terms
#' @importFrom utils menu modifyList
specify_contrasts <- function(mobj = NULL, signals = NULL, spec_list = NULL, preselected_pairwise = NULL) {
  checkmate::assert_multi_class(mobj, c("l1_model_spec", "l1_wi_spec", "hi_model_spec")) # verify that we have an object of known structure
  checkmate::assert_list(signals, null.ok = TRUE)

  checkmate::assert_integerish(mobj$level, lower = 1L, upper = 3L, len = 1L)
  lg <- lgr::get_logger(paste0("glm_pipeline/l", mobj$level, "_setup"))

  # ideally, we should expect a fitted model object to be passed in to support condition means and differences
  # default values for different contrast selections
  cond_means <- character(0) # factors for which we want model-predicted condition means
  pairwise_diffs <- if (!is.null(preselected_pairwise)) as.character(preselected_pairwise) else character(0) # factors for which we want pairwise differences
  interaction_diffs <- character(0) # factorial interaction (difference-of-differences) contrasts
  cell_means <- FALSE # whether to include cell means across all factors in design
  overall_response <- FALSE # whether to include predicted average response across all predictors
  simple_slopes <- list() # continuous variables over which to test slopes
  slope_diffs <- list() # pairwise/interaction differences in slopes across factor levels
  custom_contrasts <- mobj$contrast_spec$custom %||% mobj$contrast_list$custom
  weights <- "cells"
  factor_coding <- mobj$factor_coding %||% NULL
  cat_vars <- character(0)
  cont_vars <- character(0)
  delete <- character(0) # contrasts to be deleted based on user specification or redundancy

  # subfunction to populate contrast_spec and create $contrast_list from spec
  # uses lexical scope of specify_contrasts
  populate_mobj_contrasts <- function(mobj) {
    mobj$contrast_spec$diagonal <- include_diagonal
    mobj$contrast_spec$cond_means <- cond_means
    mobj$contrast_spec$pairwise_diffs <- pairwise_diffs
    mobj$contrast_spec$interaction_diffs <- interaction_diffs
    mobj$contrast_spec$cell_means <- cell_means
    mobj$contrast_spec$overall_response <- overall_response
    mobj$contrast_spec$simple_slopes <- simple_slopes
    mobj$contrast_spec$slope_diffs <- slope_diffs
    mobj$contrast_spec$custom <- custom_contrasts
    mobj$contrast_spec$weights <- weights
    mobj$contrast_spec$factor_coding <- factor_coding
    mobj$contrast_spec$cat_vars <- cat_vars
    mobj$contrast_spec$regressors <- mobj$regressors
    mobj$contrast_spec$delete <- delete
    if (!is.null(custom_contrasts)) {
      mobj$contrast_list$custom <- custom_contrasts
    }
    mobj <- get_contrasts_from_spec(mobj)
  }

  # list used internally for caching generated contrasts
  if (is.null(mobj$contrast_list)) {
    mobj$contrast_list <- list()
  }

  normalize_contrast_keys <- function(x, lg, allow_noncontrast = character(0)) {
    if (is.null(x) || !checkmate::test_list(x)) return(x)
    nm <- names(x)
    if (is.null(nm)) return(x)

    valid_keys <- c(
      "diagonal", "cond_means", "pairwise_diffs", "interaction_diffs", "cell_means",
      "overall_response", "simple_slopes", "slope_diffs", "custom",
      "weights", "factor_coding", "delete", "wi_contrasts"
    )
    alias_map <- c(
      include_diagonal = "diagonal",
      pairwise_diff = "pairwise_diffs",
      interaction_diff = "interaction_diffs",
      interaction_contrasts = "interaction_diffs",
      cond_mean = "cond_means",
      cell_mean = "cell_means",
      overall_responses = "overall_response",
      simple_slope = "simple_slopes",
      slope_diff = "slope_diffs",
      custom_contrasts = "custom",
      wi_contrast = "wi_contrasts"
    )

    for (alias in names(alias_map)) {
      if (alias %in% nm) {
        target <- alias_map[[alias]]
        if (target %in% nm) {
          lg$warn("Contrast spec includes both '%s' and '%s'; using '%s'.", alias, target, target)
          x[[alias]] <- NULL
        } else {
          lg$warn("Contrast spec uses legacy/alias '%s'; treating as '%s'.", alias, target)
          x[[target]] <- x[[alias]]
          x[[alias]] <- NULL
        }
      }
    }

    nm <- names(x)
    unknown <- setdiff(nm, c(valid_keys, allow_noncontrast))
    if (length(unknown) > 0L) {
      stop(
        "Unknown contrast field(s) in spec: ", paste(unknown, collapse = ", "),
        ". Allowed: ", paste(valid_keys, collapse = ", "),
        call. = FALSE
      )
    }

    return(x)
  }

  # Allow YAML-style nesting under `contrasts:` and normalize legacy aliases.
  if (!is.null(spec_list) && checkmate::test_list(spec_list)) {
    if (!is.null(spec_list$contrasts) && checkmate::test_list(spec_list$contrasts)) {
      spec_list$contrasts <- normalize_contrast_keys(spec_list$contrasts, lg)
      spec_list <- utils::modifyList(spec_list, spec_list$contrasts)
    }

    spec_list <- normalize_contrast_keys(
      spec_list,
      lg,
      allow_noncontrast = c(
        "signals", "name", "level", "model_formula", "num2fac",
        "covariate_transform", "reference_level", "factor_coding", "contrasts", "fsl_outlier_deweighting",
        "l2_scope", "l3_input_mode", "execution_backend", "producer_backend", "inference",
        "random_effects", "lmer_formula", "lmer_mask", "lmer_njobs", "lmer_glt_codes",
        "wi_expansion", "wi_factor_levels"
      )
    )

    spec_list$contrasts <- NULL
  }

  # Auto-add diagonal for intercept-only L2/L3 models to avoid interactive prompts.
  if (is.null(mobj$contrast_spec) && mobj$level %in% c(2L, 3L) && !is.null(mobj$lmfit)) {
    mm <- tryCatch(model.matrix(mobj$lmfit), error = function(e) NULL)
    if (!is.null(mm) && ncol(mm) == 1L && colnames(mm)[1] == "(Intercept)") {
      contrast_keys <- c(
        "diagonal", "cond_means", "pairwise_diffs", "interaction_diffs", "cell_means",
        "overall_response", "simple_slopes", "slope_diffs", "custom",
        "weights", "factor_coding", "delete", "wi_contrasts"
      )
      present_keys <- character(0)
      if (!is.null(spec_list) && checkmate::test_list(spec_list)) {
        present_keys <- intersect(names(spec_list), contrast_keys)
      }
      ignored_keys <- setdiff(present_keys, "diagonal")
      if (length(ignored_keys) > 0L || isTRUE(isFALSE(spec_list$diagonal))) {
        lg$warn(
          "Intercept-only model detected; ignoring contrast fields (%s) and using diagonal contrast only.",
          if (length(ignored_keys) > 0L) paste(ignored_keys, collapse = ", ") else "diagonal = FALSE"
        )
      }
      lg$info("Intercept-only model detected at level %d; adding diagonal contrast automatically.", mobj$level)
      mobj$regressors <- names(coef(mobj$lmfit))
      mobj$contrast_spec <- list(
        diagonal = TRUE,
        cond_means = character(0),
        pairwise_diffs = character(0),
        interaction_diffs = character(0),
        cell_means = FALSE,
        overall_response = FALSE,
        simple_slopes = list(),
        slope_diffs = list(),
        custom = NULL,
        weights = "cells",
        factor_coding = factor_coding,
        cat_vars = character(0),
        regressors = mobj$regressors,
        delete = character(0)
      )
      reg <- mobj$regressors
      cmat <- diag(length(reg))
      rownames(cmat) <- paste0("EV_", reg)
      colnames(cmat) <- reg
      rownames(cmat) <- gsub("(Intercept)", "Intercept", rownames(cmat), fixed = TRUE)
      mobj$contrast_list <- list(diagonal = cmat)
      mobj$contrasts <- cmat
      return(mobj)
    }
  }

  prompt_contrasts <- FALSE
  has_only_lmer_glts <- FALSE
  if (!is.null(spec_list)) {
    # If spec_list is malformed/empty, fall back to interactive prompts.
    if (!checkmate::test_list(spec_list)) {
      lg$warn("contrast spec_list is not a list; falling back to interactive contrast prompts.")
      spec_list <- NULL
    } else {
      has_any_contrast_spec <- any(c(
        "diagonal", "cond_means", "pairwise_diffs", "interaction_diffs", "cell_means",
        "overall_response", "simple_slopes", "slope_diffs", "custom",
        "weights", "factor_coding", "delete", "wi_contrasts"
      ) %in% names(spec_list))
      has_only_lmer_glts <- !isTRUE(has_any_contrast_spec) &&
        !is.null(spec_list$lmer_glt_codes) && length(spec_list$lmer_glt_codes) > 0L
      if (!isTRUE(has_any_contrast_spec) && !isTRUE(has_only_lmer_glts)) {
        if (rlang::is_interactive()) {
          lg$warn("contrast spec_list is empty or missing contrast fields; falling back to interactive contrast prompts.")
          spec_list <- NULL
        } else {
          lg$warn("contrast spec_list has no explicit contrast fields in non-interactive mode; using default contrast settings.")
        }
      }
    }
  }

  if (!is.null(spec_list)) {
    if (is.null(mobj$lmfit)) {
      emmeans_fields <- c(
        "cond_means", "pairwise_diffs", "interaction_diffs",
        "cell_means", "overall_response", "simple_slopes", "slope_diffs"
      )
      has_emmeans_fields <- any(vapply(emmeans_fields, function(fld) {
        val <- spec_list[[fld]]
        if (is.null(val)) return(FALSE)
        if (is.logical(val)) return(isTRUE(val))
        length(val) > 0L
      }, logical(1)))
      if (isTRUE(has_emmeans_fields)) {
        present_emm <- emmeans_fields[vapply(emmeans_fields, function(fld) {
          val <- spec_list[[fld]]
          if (is.null(val)) return(FALSE)
          if (is.logical(val)) return(isTRUE(val))
          length(val) > 0L
        }, logical(1))]
        lg$warn(
          "Contrast spec includes emmeans-based fields (%s), but no model is available. Dropping them. For L1, use wi_contrasts on within-subject signals.",
          paste(present_emm, collapse = ", ")
        )
      }
      spec_list$cond_means <- NULL
      spec_list$pairwise_diffs <- NULL
      spec_list$interaction_diffs <- NULL
      spec_list$cell_means <- NULL
      spec_list$overall_response <- NULL
      spec_list$simple_slopes <- NULL
      spec_list$slope_diffs <- NULL
    }

    if (!is.null(spec_list$diagonal)) {
      checkmate::assert_logical(spec_list$diagonal, len = 1L)
      include_diagonal <- spec_list$diagonal
    } else if (inherits(mobj, "l1_wi_spec") || isTRUE(has_only_lmer_glts)) {
      include_diagonal <- FALSE
    } else {
      lg$debug("In contrast specification from spec file, including diagonal contrasts by default")
      include_diagonal <- TRUE
    }

    if (!is.null(spec_list$cond_means)) {
      cond_means <- coerce_yaml_char_vec(spec_list$cond_means)
      checkmate::assert_character(cond_means)
    }

    if (!is.null(spec_list$pairwise_diffs)) {
      pairwise_diffs <- coerce_yaml_char_vec(spec_list$pairwise_diffs)
      checkmate::assert_character(pairwise_diffs)
    }

    if (!is.null(spec_list$interaction_diffs)) {
      interaction_diffs <- coerce_yaml_char_vec(spec_list$interaction_diffs)
      checkmate::assert_character(interaction_diffs)
    }

    if (!is.null(spec_list$cell_means)) {
      checkmate::assert_logical(spec_list$cell_means, len = 1L)
      cell_means <- spec_list$cell_means
    }

    if (!is.null(spec_list$overall_response)) {
      checkmate::assert_logical(spec_list$overall_response, len = 1L)
      overall_response <- spec_list$overall_response
    }

    if (!is.null(spec_list$simple_slopes)) {
      simple_slopes <- spec_list$simple_slopes
    }

    if (!is.null(spec_list$slope_diffs)) {
      slope_diffs <- spec_list$slope_diffs
    }

    if (!is.null(spec_list$factor_coding)) {
      factor_coding <- unlist(spec_list$factor_coding)
      mobj$factor_coding <- factor_coding
    }

    if (!is.null(spec_list$custom)) {
      reg_for_custom <- if (!is.null(mobj$lmfit)) {
        prefix_tmp <- ifelse(inherits(mobj, "l1_wi_spec"), paste0(mobj$signal_name, "."), "")
        paste0(prefix_tmp, names(coef(mobj$lmfit)))
      } else {
        mobj$regressors
      }
      custom_contrasts <- parse_custom_contrast_spec(spec_list$custom, reg_names = reg_for_custom, lg = lg)
    }

    if (!is.null(spec_list$weights)) {
      checkmate::assert_string(spec_list$weights)
      checkmate::assert_subset(spec_list$weights, c("equal", "cells", "proportional", "flat"))
      weights <- spec_list$weights
    } else {
      lg$debug("In contrast specficiation from spec file, setting emmeans weights to 'cells' by default")
      weights <- "cells"
    }

    if (!is.null(spec_list$delete)) {
      delete <- coerce_yaml_char_vec(spec_list$delete)
      checkmate::assert_character(delete)
    }

    # populate categorical vars for within-subject contrast and higher-level models
    if (inherits(mobj, c("l1_wi_spec", "hi_model_spec"))) {
      vtypes <- extract_lm_var_types(mobj$lmfit)
      cat_vars <- vtypes$cat_vars
      cont_vars <- vtypes$cont_vars
    }

    if (length(cat_vars) == 1L && isTRUE(spec_list$cell_means)) {
      lg$info(
        "cell_means requested with a single factor (%s); dropping because it duplicates cond_means.",
        cat_vars
      )
      spec_list$cell_means <- FALSE
      cell_means <- FALSE
    }

    # now handle population of within-subject contrasts (wi_models)

    if (!is.null(spec_list$wi_contrasts)) {
      checkmate::assert_list(spec_list$wi_contrasts)
      # Validate that wi_contrasts map to known signals with wi_factors.
      if (is.null(signals)) {
        stop("wi_contrasts were provided, but no signals list was passed to specify_contrasts.")
      }
      missing_signals <- setdiff(names(spec_list$wi_contrasts), names(signals))
      if (length(missing_signals) > 0L) {
        stop("wi_contrasts provided for unknown signals: ", paste(missing_signals, collapse = ", "))
      }
      no_wi_factors <- names(spec_list$wi_contrasts)[vapply(
        names(spec_list$wi_contrasts),
        function(sig) is.null(signals[[sig]]$wi_factors),
        logical(1)
      )]
      if (length(no_wi_factors) > 0L) {
        stop("wi_contrasts provided for signals without wi_factors: ", paste(no_wi_factors, collapse = ", "))
      }

      for (ff in seq_along(spec_list$wi_contrasts)) {
        sig <- names(spec_list$wi_contrasts)[ff]

        wi_obj <- list(
          level = 1L, lmfit = signals[[sig]]$wi_model, signal_name = sig,
          wi_factors = signals[[sig]]$wi_factors, wi_formula = signals[[sig]]$wi_formula,
          wi_expansion = signals[[sig]]$wi_expansion,
          wi_expansion_resolved = signals[[sig]]$wi_expansion_resolved
        )
        class(wi_obj) <- c("list", "l1_wi_spec")

        mobj$wi_models[[sig]] <- specify_contrasts(wi_obj, spec_list = spec_list$wi_contrasts[[ff]])
      }
    }

    mobj <- populate_mobj_contrasts(mobj)

  } else if (is.null(mobj$contrast_spec)) {
    mobj$contrast_spec <- list()
    prompt_contrasts <- TRUE
  } else {
    cat("\nCurrent model contrasts:\n\n")
    print(mobj$contrasts)
    cat("\n--------\n\n")

    res <- menu(c("No", "Yes"), title = "Do you want to modify model contrasts?")
    if (res == 2L) {
      prompt_contrasts <- TRUE
      mobj$contrast_list <- mobj$contrast_spec <- list() # clear prior contrasts that were calculated from the specification
    }
  }

  # if the user does not want to modify existing contrasts, return object unchanged
  if (isFALSE(prompt_contrasts)) {
    lg$debug("Returning model object without changing contrasts")
    return(mobj)
  }

  # handle diagonal contrasts and L1 within-factor modulation of signals
  if (inherits(mobj, "l1_wi_spec")) {
    include_diagonal <- FALSE # for within-subject factor, no need to re-ask about diagonal contrasts for a given factor
    cat("\n\n---------\n")
    cat(glue("The signal {mobj$signal_name} is modulated by the following within-subject factors: {c_string(mobj$wi_factors)}.\n", .trim=FALSE))
    cat(glue("  Please specify contrasts for this signal based on this model: {deparse(mobj$wi_formula, width.cutoff=500)}.\n\n", .trim=FALSE))
  } else if (mobj$level == 1L && inherits(mobj, "l1_model_spec")) {
    has_wi_contrasts <- FALSE # default overridden below if relevant
    # in case of level 1 object, look for any within-subject factors and prompt for relevant contrasts of each
    if (is.null(signals)) {
      lg$warn("An L1 model has been passed to specify_contrasts without a signals list. Cannot check for within-subject factors.")
    } else {
      wi_factors <- sapply(mobj$signals, function(x) !is.null(signals[[x]]$wi_factors))
      if (any(wi_factors)) {
        # setup models for each wi factor
        for (ss in names(wi_factors)[wi_factors]) {
          wi_obj <- list(
            level = 1L, lmfit = signals[[ss]]$wi_model, signal_name = ss,
            wi_factors = signals[[ss]]$wi_factors, wi_formula = signals[[ss]]$wi_formula
          )
          class(wi_obj) <- c("list", "l1_wi_spec")
          mobj$wi_models[[ss]] <- specify_contrasts(wi_obj)
        }
        has_wi_contrasts <- sapply(mobj$wi_models, function(m) {
          if (is.null(m$contrasts) || nrow(m$contrasts) == 0L) FALSE else TRUE
        })
      }
    }

    # add L1 diagonal contrasts automatically if there are no wi_factors or the user didn't ask for wi_contrasts
    if (!any(has_wi_contrasts)) {
      lg$info("Automatically adding diagonal contrasts")
      include_diagonal <- TRUE
    } else {
      # if we have wi contrasts, then prompt whether we also want overall diagonal contrasts
      include_diagonal <- menu(c("Yes", "No"), title = "Do you want to include diagonal contrasts for each regressor?")
      include_diagonal <- ifelse(include_diagonal == 1L, TRUE, FALSE)
    }
  } else {
    # at l2 and l3, ask about diagonal contrasts    
    include_diagonal <- menu(c("Yes", "No"), title = "Do you want to include diagonal contrasts for each regressor?")
    include_diagonal <- ifelse(include_diagonal == 1L, TRUE, FALSE)
  }

  format_prompt <- function(vars, signal_name=NULL) {
    if (is.null(signal_name)) {
      return(vars)
    } else {
      return(paste(vars, "modulation of", signal_name))
    }
  }

  # if we are working from a model specification in lm(), ask about factor means and pairwise comparisons
  if (is.null(mobj$lmfit)) {
    lg$debug("No $lmfit object present in mobj within specify_contrasts.")
  } else {
    # these apply to L2 and L3 model specification
    term_labels <- attr(terms(mobj$lmfit), "term.labels") # names of regressors
    int_terms <- attr(terms(mobj$lmfit), "order") # interactions in model
    vtypes <- extract_lm_var_types(mobj$lmfit)
    cat_vars <- vtypes$cat_vars
    cont_vars <- vtypes$cont_vars

    cat(
      "First, please decide on your weighting scheme for computing contrasts.",
      "For details, see the weights section of ?emmeans.",
      "In general, we recommend 'cell' weights in unbalanced designs to",
      "  approximate ordinary marginal means.\n",
      sep = "\n"
    )
    ii <- menu(c(
      "Equal (all cells weighted the same)",
      "Cell (frequencies of cells being averaged)",
      "Proportional (frequencies of factor combinations)",
      "Flat (frequencies of cells being averaged, ignoring empty cells)"
    ),
    title = paste0("How should contrasts be computed by emmeans?")
    )
    weights <- switch(EXPR = ii,
      `1` = "equal",
      `2` = "cells",
      `3` = "proportional",
      `4` = "flat"
    )
    if (is.null(weights)) weights <- "cells"

    # handle model-predicted means for each level of each factor (marginalizing over any other factors)
    for (vv in seq_along(cat_vars)) {
      title_str <- paste0(
        "Do you want to include model-predicted means for ",
        ifelse(inherits(mobj, "l1_wi_spec"), paste(cat_vars[vv], "modulation of", mobj$signal_name), cat_vars[vv]), "?"
      )

      ii <- menu(c("Yes", "No"), title = title_str)
      ii <- ifelse(ii == 1L, TRUE, FALSE)
      if (isTRUE(ii)) {
        cond_means <- c(cond_means, cat_vars[vv])
      }
    }

    # for categorical variables and their interactions, ask about pairwise differences, computed by emmeans
    emm_vars <- cat_vars
    cat_int_names <- character(0)
    if (any(int_terms > 1)) {
      cat_int <- sapply(term_labels[which(int_terms > 1)], function(tt) {
        term_split <- strsplit(tt, ":")[[1]] # split terms on colon
        all(term_split %in% cat_vars)
      })

      if (any(cat_int)) {
        cat_int_names <- names(cat_int)[cat_int]
        emm_vars <- c(emm_vars, cat_int_names)
      }
    }
    
    # handle pairwise differences for all factors (and interactions, if specified)
    for (vv in seq_along(emm_vars)) {
      if (emm_vars[vv] %in% pairwise_diffs) next
      title_str <- paste0(
        "Do you want to include pairwise differences for ",
        ifelse(inherits(mobj, "l1_wi_spec"), paste(emm_vars[vv], "modulation of", mobj$signal_name), emm_vars[vv]), "?"
      )
      ii <- menu(c("Yes", "No"), title = title_str)
      ii <- ifelse(ii == 1L, TRUE, FALSE)
      if (isTRUE(ii)) {
        pairwise_diffs <- unique(c(pairwise_diffs, emm_vars[vv]))
      }
    }

    # handle difference-of-differences interaction contrasts for categorical interactions
    if (length(cat_int_names) > 0L) {
      for (int_nm in cat_int_names) {
        title_str <- paste0(
          "Do you want to include factorial interaction (difference-of-differences) contrasts for ",
          ifelse(inherits(mobj, "l1_wi_spec"), paste(int_nm, "modulation of", mobj$signal_name), int_nm), "?"
        )
        ii <- menu(c("Yes", "No"), title = title_str)
        if (identical(ii, 1L)) {
          interaction_diffs <- unique(c(interaction_diffs, int_nm))
        }
      }
    }

    # handle model-predicted means for each *cell* of a factorial design
    if (length(cat_vars) > 1L) {
      title_str <- glue(
        ifelse(inherits(mobj, "l1_wi_spec"), "For {mobj$signal_name}, do ", "Do "),
        "you want to include model-predicted cell means across the factors: {c_string(cat_vars)}?"
      )
      ii <- menu(c("Yes", "No"), title = title_str)
      cell_means <- ifelse(ii == 1L, TRUE, FALSE)
    }

    title_str <- glue(
      "Do you want to include the overall average response",
      ifelse(inherits(mobj, "l1_wi_spec"), " across the levels of {c_string(cat_vars)}?", "?")
    )
    ii <- menu(c("Yes", "No"), title = title_str)
    overall_response <- ifelse(ii == 1L, TRUE, FALSE)

    # Handle model-predicted slopes across levels of the factors.
    # This should only fire if we have categorical x continuous interactions and
    # we wish to calculate simple slopes.
    if (length(cont_vars) > 0L && any(int_terms > 1)) {
      which_int <- int_terms > 1 # which terms are interactions (order > 1)
      term_split <- strsplit(term_labels[which_int], ":") # split terms on colon
      for (vv in cont_vars) {
        has_vv <- which(sapply(term_split, function(x) any(x == vv))) # does this term have the target var
        for (tt in seq_along(has_vv)) {
          this_term <- term_split[[has_vv[tt]]]
          fac_part <- this_term[!this_term == vv]
          if (length(fac_part) == 0L || !all(fac_part %in% cat_vars)) next
          ii <- menu(c("Yes", "No"), title = paste0(
            "Do you want to include model-predicted simple slopes of ", vv, " across levels of ",
            paste(fac_part, collapse = ", "), "?"
          ))
          if (identical(ii, 1L)) { # build list where top-level elements are continuous variables
            simple_slopes[[vv]] <- c(simple_slopes[[vv]], list(fac_part))
          }
          ii_diff <- menu(c("Yes", "No"), title = paste0(
            "Do you want to include pairwise slope differences of ", vv, " across levels of ",
            paste(fac_part, collapse = ", "), "?"
          ))
          if (identical(ii_diff, 1L)) {
            slope_diffs[[vv]] <- c(slope_diffs[[vv]], list(fac_part))
          }
        }
      }
    }
  }

  # populate contrast specification object
  mobj <- populate_mobj_contrasts(mobj)

  cmat_base <- matrix(numeric(0), nrow = 0, ncol = length(mobj$regressors), dimnames = list(NULL, mobj$regressors))
  cmat <- mobj$contrasts

  # N.B. For now, hand off full model matrix and contrasts to fMRI fitting function (e.g., feat) and
  # let it handle the rank degeneracy. We may need to come back here and drop bad contrasts and use
  # the reduced model matrix in the output. The challenge is that the contrasts calculated by emmeans
  # may be wrong for any column involving aliasing
  if (!is.null(mobj$aliased_terms)) {
    cat(
      "Aliased terms in model.",
      "We will drop aliased columns from the contrast matrix to match the model matrix.",
      "This may invalidate some contrast vectors. Check these manually, check these carefully!",
      "You will probably need to fix the contrasts by hand!",
      sep = "\n"
    )

    has_alias <- which(apply(cmat[, intersect(mobj$aliased_terms, colnames(cmat)), drop=F], 1, function(x) {
      any(abs(x - 0) > 1e-5)
    }))

    if (length(has_alias) > 0L) {
      cat("These contrasts used at least one aliased column. Check them!\n")
      cat(strwrap(paste(names(has_alias), collapse = ", "), 70, exdent = 5), sep = "\n")
    }

    # Not used at present
    cmat_reduce <- cmat[, !colnames(cmat) %in% mobj$aliased_terms, drop=F]
    # write.csv(cmat_reduce, file = "cmat_test_reduce.csv")
  }

  # for a within-subject factor model, no need to request contrast editor since this should be handled at overall model level
  if (inherits(mobj, "l1_wi_spec")) {
    return(mobj)
  }

  if (!rlang::is_interactive()) {
    if (!is.null(mobj$aliased_terms)) mobj$contrasts_noalias <- cmat_reduce
    return(mobj)
  }

  # syntax of contrast
  cat("\n\n-----\nContrast editor\n\n")
  cat("When entering contrasts, you have two options:\n\n")
  cat("  1) enter a vector of numeric values, one for each regressor\n     Example: 1 0 -1 0 0\n\n")
  cat("  2) enter name=value pairs for relevant coefficients\n     Example: pe_1h = 1 pe_2h = -1\n\n")

  add_more <- 1
  while (add_more != 4) {
    # refresh model object with each iteration since some things may change
    mobj <- get_contrasts_from_spec(mobj)

    add_more <- menu(c("Add contrast", "Show contrasts", "Delete contrast", "Done with contrast setup"),
      title = "Contrast setup menu"
    )
    if (identical(add_more, 0L)) break

    if (add_more == 1L) {
      cname <- readline("Enter contrast name: ")
      cat("Columns of contrast matrix (one per regressor)\n")
      cat(strwrap(paste(colnames(cmat_base), collapse = ", "), 70, exdent = 5), sep = "\n")
      centry <- readline("Enter contrast syntax: ")

      # Support column = value named syntax or numeric vector syntax
      if (grepl("=", centry)) {
        # split into pairs
        centry <- gsub("([^\\s]+)\\s*=\\s*([^\\s]+)", "\\1=\\2", centry, perl = TRUE)
        csplit <- do.call(rbind, strsplit(strsplit(centry, "\\s+")[[1L]], "=")) # first element is name, second is value
        cvec <- rep(0, ncol(cmat_base)) %>% setNames(colnames(cmat_base))

        for (ii in seq_len(nrow(csplit))) {
          if (!csplit[ii, 1] %in% colnames(cmat_base)) {
            lg$warn("Cannot find column: %s in contrast matrix. Ignoring contrast input.", csplit[ii, 1])
            cvec <- NULL # do not add
            break
          } else {
            cvec[csplit[ii, 1]] <- as.numeric(csplit[ii, 2])
          }
        }
      } else {
        cvec <- sapply(strsplit(centry, "(,\\s*|\\s+)")[[1L]], as.numeric)
        if (length(cvec) != ncol(cmat_base)) {
          lg$warn("Number of elements in entered contrast does not match number of columns in design. Not adding contrast.")
          cvec <- NULL # do not add
        } else {
          names(cvec) <- colnames(cmat_base)
        }
      }

      if (!is.null(cvec)) {
        mobj$contrast_list$custom <- rbind(mobj$contrast_list$custom, cvec)
        rownames(mobj$contrast_list$custom)[nrow(mobj$contrast_list$custom)] <- cname
        mobj$contrast_spec$custom <- mobj$contrast_list$custom
      }
    } else if (add_more == 2L) {
      summarize_contrasts(mobj$contrasts)
    } else if (add_more == 3L) {
      whichdel <- menu(rownames(mobj$contrasts), title = "Which contrast should be deleted?")
      if (whichdel != 0) {
        mobj$contrast_spec$delete <- c(mobj$contrast_spec$delete, rownames(mobj$contrasts)[whichdel])
      }
    } else if (add_more == 4L) {
      cat("Exiting contrast setup\n")
    }
  }

  if (!is.null(mobj$aliased_terms)) mobj$contrasts_noalias <- cmat_reduce
  return(mobj)
}
