.mixed_by_normalize_pairs_args <- function(pairs_args, field_name) {
  if (is.null(pairs_args) || identical(pairs_args, FALSE)) {
    return(NULL)
  }
  if (isTRUE(pairs_args)) {
    return(list())
  }
  if (is.character(pairs_args) && length(pairs_args) == 1L) {
    return(list(adjust = pairs_args))
  }
  if (!is.list(pairs_args)) {
    stop(field_name, " must be TRUE/FALSE, a single adjustment method, or a named list of pairs() arguments.", call. = FALSE)
  }
  pairs_args
}

.mixed_by_merge_pair_args <- function(base_args, override_args) {
  if (length(override_args) == 0L) {
    return(base_args)
  }
  override_names <- names(override_args)
  if (!is.null(override_names)) {
    named_overrides <- override_names[nzchar(override_names)]
    base_args[named_overrides] <- NULL
  }
  c(base_args, override_args)
}

.mixed_by_split_posthoc_spec <- function(spec, pairs_field) {
  spec[c("outcome", "model_name", "pairs", pairs_field, paste0(pairs_field, "_args"))] <- NULL
  spec
}

.mixed_by_tidy_emm_result <- function(result, component = NULL) {
  if (inherits(result, "emm_list")) {
    if (!is.null(component) && component %in% names(result)) {
      return(broom.mixed::tidy(result[[component]]))
    }

    components <- setdiff(names(result), "contrasts")
    if (length(components) == 1L) {
      return(broom.mixed::tidy(result[[components]]))
    }

    return(data.table::rbindlist(lapply(components, function(component) {
      td <- broom.mixed::tidy(result[[component]])
      td$emm_component <- component
      td
    }), fill = TRUE))
  }
  broom.mixed::tidy(result)
}

.mixed_by_get_pairs_args <- function(spec, pairs_field, pairs_enabled, pairs_args) {
  has_spec_pairs <- pairs_field %in% names(spec) || "pairs" %in% names(spec)

  if (pairs_field %in% names(spec)) {
    spec_pairs <- spec[[pairs_field]]
  } else if ("pairs" %in% names(spec)) {
    spec_pairs <- spec[["pairs"]]
  } else {
    spec_pairs <- NULL
  }

  if (has_spec_pairs) {
    out <- .mixed_by_normalize_pairs_args(spec_pairs, pairs_field)
  } else if (isTRUE(pairs_enabled)) {
    out <- pairs_args
  } else {
    out <- NULL
  }

  spec_args_field <- paste0(pairs_field, "_args")
  if (!is.null(out) && spec_args_field %in% names(spec)) {
    out <- .mixed_by_merge_pair_args(out, spec[[spec_args_field]])
  }

  out
}

.mixed_by_tidy_emm_pairs <- function(result, pairs_args, main_component, pairs_field) {
  if (inherits(result, "emm_list")) {
    if (main_component %in% names(result)) {
      result <- result[[main_component]]
    } else if ("contrasts" %in% names(result)) {
      warning(
        "Using precomputed contrasts from emm_list because the main ", main_component,
        " component is unavailable; ", pairs_field, "_args cannot be reapplied."
      )
      return(broom.mixed::tidy(result[["contrasts"]]))
    } else {
      stop("Cannot apply ", pairs_field, " to this emmeans result list.", call. = FALSE)
    }
  }

  pairs_result <- do.call(utils::getS3method("pairs", "emmGrid"), c(list(x = result), pairs_args))
  broom.mixed::tidy(pairs_result)
}

.mixed_by_run_posthoc_spec <- function(posthoc_fun, posthoc_model, spec, main_component,
                                       pairs_field, pairs_enabled, pairs_args) {
  call_args <- .mixed_by_split_posthoc_spec(spec, pairs_field)
  result <- do.call(posthoc_fun, c(list(object = posthoc_model), call_args))
  these_pairs_args <- .mixed_by_get_pairs_args(spec, pairs_field, pairs_enabled, pairs_args)

  list(
    main = .mixed_by_tidy_emm_result(result, component = main_component),
    pairs = if (is.null(these_pairs_args)) NULL else .mixed_by_tidy_emm_pairs(result, these_pairs_args, main_component, pairs_field)
  )
}

.mixed_by_extract_posthoc_list <- function(nested_list, data_element, spec, value_col, split_on, extract_fun) {
  if (is.null(spec)) {
    return(NULL)
  }

  posthoc_data <- extract_fun(nested_list, data_element, split_on)
  if (nrow(posthoc_data) == 0L) {
    return(NULL)
  }

  out <- lapply(seq_along(spec), function(aa) {
    posthoc_sub <- subset(posthoc_data, outcome == spec[[aa]]$outcome & model_name == spec[[aa]]$model_name)
    other_keys <- setdiff(names(posthoc_sub), value_col)
    posthoc_name <- names(spec)[aa]
    posthoc_sub[, get(value_col)[[1]][[posthoc_name]], by = other_keys]
  })
  names(out) <- names(spec)
  out
}

#' Mixed by runs a set of mixed-effects models for each combination of a set of factors. Its primary use is to
#'   run the same model on different splits of the data.
#'
#' @param data A data.frame or data.table object containing stacked data for each combination of the \code{split_on}
#'   variables. The function will run separate mixed-effect models for each combination. Alternatively, a vector of
#'   filenames can be passed, which will be read in sequentially and fit (.rds, .csv, .dat, and .txt
#'   supported at present).
#' @param outcomes A character vector of outcome variables to be analyzed
#' @param rhs_model_formulae A named list of lme4-format formula specifying the exact model to be run for each data split.
#' @param model_formulae Alternative to the outcome + rhs_model_formulae approach. This is a list of lme4-format
#'   formulae that includes the outcome on the left-hand side. This is useful if the outcomes change from one
#'   model to the next, but you don't want the Cartesian product (combinations) of outcomes and rhs_model_formulae
#' @param split_on A character vector of columns in \code{data} used to split the analyses into separate models.
#' @param external_df An optional data.frame/data.table containing external data that should be joined with \code{data}
#'   prior to model fitting. Useful if \code{data} contains external time series data and \code{external_df} is a
#'   dataset of behavioral variables that do not vary by neural sensor/region. Optionally, this can be a single
#'   filename to a .rds file if you want to pass in the filename, not the data itself.
#' @param external_merge_by A character vector specifying which columns of \code{external_df} and \code{data} should
#'   be used for creating a combined dataset.
#' @param padjust_by A character vector or list consisting of one or more variables over which an adjusted p-value
#'   should be calculated. This defaults to adjusting by each term (fixed effect) in the model, which will adjust for
#'   all tests for a given term across variables in \code{split_on}.
#' @param padjust_method The adjustment method (see \code{?p.adjust}) for adjusting p-values. Multiple values
#'   can be passed as a character vector, in which case multiple corrections will be added as distinct columns.
#' @param outcome_transform A vectorized function that will be applied to the outcome variable prior to running
#'   the model. For example, \code{outcome_transform=function(x) { as.vector(scale(x)) } }.
#' @param scale_predictors An optional vector of predictor names that should be z-scored (unit normalized) prior
#'   to entering into the lmer model. Note that this scaling will be applied to the predictor regardless of which
#'   rhs model is being run as long as that predictor is in the model.
#' @param ncores The number of compute cores to be used in the computation. Defaults to 1.
#' @param cl An optional external cl object (created by a variant of makeCluster) used for computation. Can
#'   save the overhead of starting and stopping many workers in a loop context.
#' @param refit_on_nonconvergence The number of times a model should be refit if it does not converge. Final estimates
#'   from one iteration are used as starting values for the next.
#' @param tidy_args A list of arguments passed to tidy.merMod for creating the coefficient data.frame. By default,
#'   the function only returns the fixed effects and computes confidence intervals using the Wald method.
#' @param lmer_control An lmerControl object specifying any optimization settings to be passed to lmer()
#' @param engine Estimation engine to use. \code{"lme4"} uses \code{lmerTest::lmer()} for
#'   compatibility with the historical output, \code{"rstanarm"} uses \code{rstanarm::stan_lmer()},
#'   and \code{"jlmer"} uses the package's \code{jlmer()} function.
#' @param engine_args A named list of additional arguments passed to the selected engine.
#' @param calculate A character vector specifying what calculations should be returned by the function. 
#'   The options are: "parameter_estimates_reml", "parameter_estimates_ml", "fit_statistics", "fitted", and "residuals".
#' @param emmeans_spec A named list of emmeans calls to be run for each model to obtain model predictions.
#'   Any arguments that are valid for emmeans can be passed through the list structure.
#' @param emmeans_pairs A boolean indicating whether to also return pairwise comparisons (via
#'   \code{pairs()} on the \code{emmeans} result) for each requested element of \code{emmeans_spec}.
#' @param emmeans_pairs_args A named list of arguments passed to \code{pairs()} (e.g.,
#'   \code{list(adjust = "tukey")}).
#' @param emtrends_spec A named list of emtrends calls to be run for each model to obtain model-predicted slopes.
#'   Any arguments that are valid for emtrends can be passed through the list structure.
#' @param emtrends_pairs A boolean indicating whether to also return pairwise comparisons (via
#'   \code{pairs()} on the \code{emtrends} result) for each requested element of \code{emtrends_spec}.
#' @param emtrends_pairs_args A named list of arguments passed to \code{pairs()} (e.g.,
#'   \code{list(adjust = "tukey")}).
#' @param return_models A boolean indicating whether to return fitted model objects, which can be used for
#'   post hoc contrasts, visualization and statistics. Note that model objects can get very large, so be
#'   careful with this option since it could generate a massive data object.
#' @param weights Optional character scalar naming a numeric column in \code{data}
#'   (or \code{external_df}) containing positive observation/trial precision
#'   weights (e.g., \code{"precision_weight"} from \code{decon_epoch_regression()})
#'   to be normalized within each split and passed to the mixed-model engine.
#'   
#' @details In general, restricted maximum likelihood (REML) should be used for making inferences about parameter
#'   estimates, whereas ML should be used for model comparisons based on log-likelihood (e.g., AIC). 
#'   If "parameter_estimates_reml" are requested in \code{calculate}, then models will
#'   be fitted with REML. If "parameter_estimates_ml" and/or "fit_statistics" are requested, models will be
#'   fitted with ML. Note that if both REML and ML are requested, each model is fit twice since the estimators
#'   each have advantages and disadvantages noted above.
#'
#'   For \code{engine = "rstanarm"}, coefficient tables include \code{pd}, the
#'   posterior probability of direction from \code{bayestestR::p_direction()},
#'   and \code{p.value}, its two-sided p-like conversion using
#'   \code{as_p = TRUE}. These are posterior-derived indices rather than
#'   frequentist df-based p-values.
#'   
#' Example of emmeans_spec usage:
#'   mixed_by(data, emmeans_spec=list(
#'     em1=list(specs = ~ memory | noise_level, adjust = "sidak", weights = "cells"),
#'     em2=list(specs = ~ memory * noise_level, weights = "equal"),
#'     em3=list(specs = ~ memory)
#'   ))
#'
#' @return A list containing coefficient tables, fit statistics, optional post-hoc outputs
#'   (\code{emmeans_list}, \code{emmeans_pairs_list}, \code{emtrends_list},
#'   \code{emtrends_pairs_list}), residuals, fitted values, and optionally fitted models.
#'
#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
#' @importFrom checkmate assert_data_frame assert_character assert_subset assert_formula
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom iterators iter
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom broom.mixed tidy
#' @importFrom data.table fread setDT setkeyv setattr copy
#' @importFrom stats fitted na.exclude p.adjust update update.formula
#' @export
mixed_by <- function(data, outcomes = NULL, rhs_model_formulae = NULL, model_formulae = NULL, split_on = NULL,
                     external_df = NULL, external_merge_by = NULL,
                     padjust_by = "term", padjust_method = "BY", outcome_transform = NULL, scale_predictors = NULL,
                     ncores = 1L, cl = NULL, refit_on_nonconvergence = 3,
                     tidy_args = list(effects = "fixed", conf.int = TRUE),
                     lmer_control = lmerControl(optimizer = "nloptwrap"),
                     engine = c("lme4", "rstanarm", "jlmer"), engine_args = list(),
                     calculate=c("parameter_estimates_reml", "parameter_estimates_ml", "fit_statistics"),
                     return_models = FALSE, emmeans_spec = NULL, emmeans_pairs = FALSE,
                     emmeans_pairs_args = list(), emtrends_spec=NULL, emtrends_pairs = FALSE,
                     emtrends_pairs_args = list(), weights = NULL) {

  engine <- match.arg(engine)
  checkmate::assert_list(engine_args, names = "unique")
  if (is.null(weights) && !is.null(engine_args$weights) &&
      checkmate::test_string(engine_args$weights)) {
    weights <- engine_args$weights
    engine_args$weights <- NULL
  }
  checkmate::assert_string(weights, null.ok = TRUE)
  if (engine == "rstanarm" && !requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Package \"rstanarm\" must be installed to use engine = \"rstanarm\".", call. = FALSE)
  }
  if (engine == "rstanarm" && !requireNamespace("bayestestR", quietly = TRUE)) {
    stop("Package \"bayestestR\" must be installed to use engine = \"rstanarm\".", call. = FALSE)
  }
  if (engine == "jlmer" && !is.null(engine_args$jlmer_fun)) {
    stop(
      "engine_args$jlmer_fun is no longer supported; engine = \"jlmer\" always uses fmri.pipeline::jlmer().",
      call. = FALSE
    )
  }

  ## VALIDATE INPUTS
  # support data.frame input for single dataset execution or a vector of files that are imported and fit sequentially
  if (checkmate::test_data_frame(data)) {
    single_df <- TRUE
  } else {
    if (!checkmate::test_character(data)) {
      stop("If data is not a data.frame, it should be a set of file names.")
    }
    checkmate::assert_file_exists(data)
    single_df <- FALSE
  }

  checkmate::assert_character(split_on, null.ok = TRUE, unique = TRUE)
  if (!is.null(padjust_by)) {
    if (!is.list(padjust_by)) {
      padjust_by <- list(padjust_by)
    } # convert to list for consistency
    sapply(padjust_by, checkmate::assert_character, null.ok = TRUE)
  }
  checkmate::assert_string(padjust_method)
  checkmate::assert_integerish(ncores, lower = 1L)
  checkmate::assert_class(cl, "cluster", null.ok = TRUE)
  checkmate::assert_subset(calculate, c("parameter_estimates_reml", "parameter_estimates_ml", "fit_statistics", "residuals", "fitted"))
  checkmate::assert_logical(return_models, len = 1L)
  checkmate::assert_character(scale_predictors, null.ok = TRUE)
  checkmate::assert_list(emmeans_spec, null.ok = TRUE)
  checkmate::assert_logical(emmeans_pairs, len = 1L)
  checkmate::assert_list(emmeans_pairs_args, null.ok = FALSE)
  checkmate::assert_list(emtrends_spec, null.ok = TRUE)
  checkmate::assert_logical(emtrends_pairs, len = 1L)
  checkmate::assert_list(emtrends_pairs_args, null.ok = FALSE)

  # turn off refitting if user specifies 'FALSE'
  if (is.logical(refit_on_nonconvergence) && isFALSE(refit_on_nonconvergence)) {
    refit_on_nonconvergence <- 0L
  }
  checkmate::assert_integerish(refit_on_nonconvergence, null.ok = FALSE)

  validate_form_list <- function(ll) {
    if (inherits(ll, "formula")) {
      ll <- list(ll)  # wrap single formula as list
    } else {
      checkmate::assert_list(ll)
    }
    
    lapply(ll, checkmate::assert_formula)
    if (is.null(names(ll))) {
      nm <- paste0("model", seq_along(ll))
      message("Using default model names of ", paste(nm, collapse=", "))
      names(ll) <- nm
    }
    return(ll)
  }
  
  if (!is.null(model_formulae)) {
    message("Using model_formulae to setup model_set")
    model_formulae <- validate_form_list(model_formulae)
    
    # don't expand combinations, just form outcomes, rhs, and formulae
    model_set <- tibble::tibble(outcome = sapply(model_formulae, formula.tools::lhs.vars), model_name = names(model_formulae))
    model_set$rhs <- lapply(model_formulae, function(x) { update.formula(x, "NULL ~ .") })
    model_set$form <- model_formulae

  } else if (!is.null(rhs_model_formulae)) {
    checkmate::assert_character(outcomes, null.ok = FALSE) # must provide valid outcomes
    rhs_model_formulae <- validate_form_list(rhs_model_formulae)
    
    model_set <- expand.grid(outcome = outcomes, rhs = rhs_model_formulae, stringsAsFactors = FALSE)
    model_set$form <- lapply(seq_len(nrow(model_set)), function(ii) {
      update.formula(model_set$rhs[[ii]], paste(model_set$outcome[[ii]], "~ ."))
    })
    model_set <- tibble::as_tibble(model_set)
    model_set$model_name <- names(model_set$rhs)
  }
  
  # handle external_df
  if (!is.null(external_df)) {
    if (checkmate::test_string(external_df)) {
      checkmate::assert_file_exists(external_df)
      external_df <- readRDS(external_df) # only supports .rds at the moment
    } else {
      checkmate::assert_data_frame(external_df)
      external_df <- if (is.data.table(external_df)) {
        data.table::copy(external_df)
      } else {
        data.table::as.data.table(external_df)
      }
      checkmate::assert_subset(external_merge_by, names(external_df))
      data.table::setkeyv(external_df, external_merge_by) # key local external data by merge columns
    }
  }

  # worker subfunction to fit a given model to a data split
  model_worker <- function(data, model_formula, lmer_control, outcome_transform = NULL, scale_predictors = NULL,
                           REML = TRUE, engine = "lme4", engine_args = list(), weights = NULL) {
    if (!is.null(outcome_transform)) { # apply transformation to outcome
      lhs <- all.vars(model_formula)[1]
      data[[lhs]] <- outcome_transform(data[[lhs]])
    }

    #apply z-scoring to specified predictors
    if (!is.null(scale_predictors)) {
      # technically this will match the outcome, too -- hopefully the user uses outcome_transform in that case
      model_terms <- all.vars(model_formula)
      if (any(model_terms %in% scale_predictors)) {
        for (aa in intersect(scale_predictors, model_terms)) {
          data[[aa]] <- as.vector(scale(data[[aa]]))
        }
      }
    }

    if (!is.null(weights)) {
      w_vec <- as.numeric(data[[weights]])
      bad_w <- !is.finite(w_vec) | (w_vec <= 0)
      if (any(bad_w)) {
        w_vec[bad_w] <- NA_real_
      }
      ok_w <- is.finite(w_vec) & (w_vec > 0)
      if (any(ok_w)) {
        w_vec[ok_w] <- w_vec[ok_w] / mean(w_vec[ok_w])
      }
      if (engine == "jlmer") {
        lhs <- all.vars(model_formula)[1]
        if (any(bad_w)) data[[lhs]][bad_w] <- NA_real_
        vv_jl <- all.vars(model_formula)
        comp_jl <- stats::complete.cases(as.data.frame(data)[, vv_jl, drop = FALSE])
        engine_args$wts <- w_vec[comp_jl]
      } else {
        engine_args$weights <- w_vec
      }
    }
    
    if (engine == "lme4") {
      fit_args <- c(
        list(formula = model_formula, data = data, control = lmer_control, REML = REML, na.action = na.exclude),
        engine_args
      )
      md <- do.call(lmerTest::lmer, fit_args)
    } else if (engine == "rstanarm") {
      fit_args <- c(list(formula = model_formula, data = data, na.action = na.exclude), engine_args)
      md <- do.call(rstanarm::stan_lmer, fit_args)
    } else if (engine == "jlmer") {
      fit_args <- c(list(formula = model_formula, data = data, REML = REML), engine_args)
      md <- do.call(jlmer, fit_args)
    } else {
      stop("Unsupported mixed model engine: ", engine, call. = FALSE)
    }

    if (engine == "lme4" && refit_on_nonconvergence > 0L) {
      rfc <- 0
      while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages)) && rfc < refit_on_nonconvergence) {
        # print(md@optinfo$conv$lme4$conv)
        ss <- lme4::getME(md, "theta")
        lmod <- lmer_control
        lmod$optimizer <- "bobyqa" # produces convergence more reliably
        md <- update(md, start = ss, control = lmod)
        rfc <- rfc + 1 # increment refit counter
      }
    }
    return(md)
  }

  # convert to data.table and nest
  if (isTRUE(single_df)) {
    # Keep all subsequent by-reference operations local to this function.
    data <- if (is.data.table(data)) {
      data.table::copy(data)
    } else {
      data.table::as.data.table(data)
    }
    df_set <- c("internal")
  } else {
    df_set <- data
    # even though we aren't splitting on this, adding it will propagate the filename to the output structure
    if (!".filename" %in% split_on) {
      split_on <- c(".filename", split_on)
    }
  }

  # process emmeans specification setup
  if (!is.null(emmeans_spec)) {
    #has the outcome and model name to lookup whether to run a given emmeans specification
    if (is.null(names(emmeans_spec))) {
      names(emmeans_spec) <- paste("emm", seq_along(emmeans_spec), sep="_")
    } else {
      empty_names <- which(names(emmeans_spec) == "")
      names(emmeans_spec)[empty_names] <- paste("emm", empty_names, sep="_")
    }
    
    emm_metadata <- data.table::rbindlist(lapply(emmeans_spec, function(ee) { data.frame(ee[c("outcome", "model_name")])}))
    emm_metadata$emm_label <- names(emmeans_spec)
    
    bad_outcomes <- !emm_metadata$outcome %in% model_set$outcome
    if (any(bad_outcomes)) {
      cat("Cannot have an outcome in emmeans that is not present in the models.\nIgnoring this specification!\n")
      print(emm_metadata[bad_outcomes, ])
    }
    
    bad_models <- !emm_metadata$model_name %in% model_set$model_name
    if (any(bad_models)) {
      cat("Cannot have model names in emmeans that are not present in the model formula names.\nIgnoring this specification!\n")
      print(emm_metadata[bad_models, ])
    }
    
    # logical and of outcomes and models to keep
    keep_models <- !bad_outcomes & !bad_models
    
    if (sum(keep_models) == 0L) {
      emmeans_spec <- NULL # reset to NULL because no emts specified overlap with the models and outcomes
    } else {
      # subset metadata and trends list to only relevant models
      emm_metadata <- emm_metadata[keep_models,,drop=FALSE]
      emm_metadata$emm_number <- 1:nrow(emm_metadata)
      emmeans_spec <- emmeans_spec[keep_models]
    }
    
  }
  
  # process emtrends specification setup
  if (!is.null(emtrends_spec)) {
    #has the outcome and model name to lookup whether to run a given emtrends specification
    if (is.null(names(emtrends_spec))) {
      names(emtrends_spec) <- paste("emt", seq_along(emtrends_spec), sep="_")
    } else {
      empty_names <- which(names(emtrends_spec) == "")
      names(emtrends_spec)[empty_names] <- paste("emt", empty_names, sep="_")
    }

    emt_metadata <- data.table::rbindlist(lapply(emtrends_spec, function(ee) { data.frame(ee[c("outcome", "model_name")])}))
    emt_metadata$emt_label <- names(emtrends_spec)
    
    bad_outcomes <- !emt_metadata$outcome %in% model_set$outcome
    if (any(bad_outcomes)) {
      cat("Cannot have an outcome in emtrends that is not present in the models.\nIgnoring this specification!\n")
      print(emt_metadata[bad_outcomes, ])
    }
    
    bad_models <- !emt_metadata$model_name %in% model_set$model_name
    if (any(bad_models)) {
      cat("Cannot have model names in emtrends that are not present in the model formula names.\nIgnoring this specification!\n")
      print(emt_metadata[bad_models, ])
    }
    
    # logical and of outcomes and models to keep
    keep_models <- !bad_outcomes & !bad_models
    
    if (sum(keep_models) == 0L) {
      emtrends_spec <- NULL # reset to NULL because no emts specified overlap with the models and outcomes
    } else {
      # subset metadata and trends list to only relevant models
      emt_metadata <- emt_metadata[keep_models,,drop=FALSE]
      emt_metadata$emt_number <- 1:nrow(emt_metadata)
      emtrends_spec <- emtrends_spec[keep_models]
    }
    
  }

  # setup parallel compute
  if (ncores > 1L && is.null(cl)) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(try(stopCluster(cl)))
  } else if (!is.null(cl)) {
    message("Using external cl object for parallelism. mixed_by will not stop the cluster upon completion.")
    message("Use stopCluster() yourself when done with all computation.")
    registerDoParallel(cl, cores = ncores)
  } else {
    registerDoSEQ()
  }

  mresults <- vector(mode = "list", length(df_set)) # preallocate list
  posthoc_helper_exports <- c(
    ".mixed_by_normalize_pairs_args",
    ".mixed_by_merge_pair_args",
    ".mixed_by_split_posthoc_spec",
    ".mixed_by_tidy_emm_result",
    ".mixed_by_get_pairs_args",
    ".mixed_by_tidy_emm_pairs",
    ".mixed_by_run_posthoc_spec"
  )

  # foreach workers need local bindings for unexported package helpers.
  .mixed_by_normalize_pairs_args <- .mixed_by_normalize_pairs_args
  .mixed_by_merge_pair_args <- .mixed_by_merge_pair_args
  .mixed_by_split_posthoc_spec <- .mixed_by_split_posthoc_spec
  .mixed_by_tidy_emm_result <- .mixed_by_tidy_emm_result
  .mixed_by_get_pairs_args <- .mixed_by_get_pairs_args
  .mixed_by_tidy_emm_pairs <- .mixed_by_tidy_emm_pairs
  .mixed_by_run_posthoc_spec <- .mixed_by_run_posthoc_spec
  environment(.mixed_by_get_pairs_args) <- environment()
  environment(.mixed_by_run_posthoc_spec) <- environment()

  # loop over each dataset to be fit
  for (i in seq_along(df_set)) {
    df_i <- df_set[i]

    # read each dataset if operating in multiple data scenario
    if (isFALSE(single_df)) {
      message("Reading file: ", df_i)
      if (grepl(".rds$", df_i, ignore.case = TRUE)) {
        data <- readRDS(df_i)
        if (!is.data.table(data)) {
          data.table::setDT(data)
        }
      } else if (grepl("(.csv|.csv.gz|.csv.bz2|.dat|.txt|.txt.gz|.txt.bz2)", df_i, ignore.case = TRUE, perl = TRUE)) {
        data <- data.table::fread(df_i, data.table = TRUE)
      } else {
        stop("Unable to sort out this data input: ", df_i)
      }
      data[, .filename := basename(df_i)] # add filename
    }

    # validate structure of data against models to be fit
    if (is.null(split_on)) {
      split_on <- "split" # dummy split to make code function consistently
      data[, split := factor(1)]
      has_split <- FALSE
    } else {
      has_split <- TRUE
    }

    # handle external_df, if requested
    if (!is.null(external_df)) {
      checkmate::assert_subset(external_merge_by, names(data))
      data <- merge(data, external_df, by = external_merge_by)
    }

    # verify that outcomes and split variables are present in data
    checkmate::assert_subset(model_set$outcome, names(data))
    checkmate::assert_subset(split_on, names(data))
    if (!is.null(weights)) {
      checkmate::assert_subset(weights, names(data))
      checkmate::assert_numeric(data[[weights]])
    }

    # nest data.tables for each combination of split factors
    data.table::setkeyv(data, split_on)
    data <- data[, .(dt = list(.SD)), by = split_on]

    # for each split and each outcome + rhs, examine whether there are 0 non-NA cases
    # for (rr in seq_len(nrow(data))) {
    #   for (mm in seq_len(nrow(model_set))) {
    #     ff <- model_set$form[[mm]]
    #     vv <- all.vars(ff)
    #     miss_data <- data[rr, dt[[1]]] %>%
    #       dplyr::select(dplyr::all_of(vv)) %>% # just keep model-relevant variables
    #       mutate(any_miss = rowSums(is.na(dplyr::select(., dplyr::any_of(vv)))) > 0)
    #     n_present <- miss_data %>%
    #       dplyr::filter(any_miss == FALSE) %>%
    #       nrow()
    # 
    #     if (n_present == 0) {
    #       ss <- paste(names(data[rr, ..split_on]), data[rr, ..split_on], sep="=", collapse=", ")
    #       ss <- paste(ss, "model=", model_set$model_name[mm], "outcome=", model_set$outcome[mm])
    #       message("No non-missing observations for split combination: ", ss)
    #     }
    #   }
    # }

    # loop over outcomes and rhs formulae within each chunk to maximize compute time by chunk (reduce worker overhead)
    message("Starting processing of data splits")
    foreach_packages <- c("lme4", "lmerTest", "data.table", "dplyr", "broom.mixed", "emmeans")
    if (engine == "rstanarm") {
      foreach_packages <- c(foreach_packages, "rstanarm", "bayestestR")
    }
    mresults[[i]] <- foreach(
      dt_split = iter(data, by = "row"), .packages = foreach_packages,
      .export = posthoc_helper_exports, .noexport = "data", .inorder = FALSE
    ) %dopar% {
      if (nrow(dt_split$dt[[1L]]) == 0L) {
        warning("No rows found in split. Skipping")
        return(NULL)
      }

      split_results <- lapply(seq_len(nrow(model_set)), function(mm) {
        ff <- model_set$form[[mm]]
        ret <- data.table::copy(dt_split)
        ret[, outcome := model_set$outcome[[mm]]]
        ret[, model_name := names(model_set$rhs)[mm]]
        ret[, rhs := as.character(model_set$rhs[mm])]
        ret[, engine := engine]
        
        # defaults
        ret[, coef_df_reml := list()]
        ret[, coef_df_ml := list()]
        ret[, fit_df := list()]
        ret[, residuals := list()]
        ret[, fitted := list()]
        ret[, emm := list()] 
        ret[, emp := list()]
        ret[, emt := list()]
        ret[, etp := list()]
        if (isTRUE(return_models)) ret[, model := list()]
        
        # check missingness
        vv <- all.vars(ff)
        
        miss_data <- as.data.frame(ret$dt[[1]])[, vv, drop = FALSE]
        valid_rows <- stats::complete.cases(miss_data)
        if (!is.null(weights)) {
          w_chk <- as.numeric(ret$dt[[1]][[weights]])
          valid_rows <- valid_rows & is.finite(w_chk) & (w_chk > 0)
        }
        n_present <- sum(valid_rows)
        
        if (n_present == 0) {
          ret[, dt := NULL]
          return(ret)
        }
        
        ##

        thism <- NULL
        thism_ml <- NULL
        estimator <- if (engine == "rstanarm") "bayes" else NA_character_

        glance_model <- function(model, estimator) {
          gd <- broom.mixed::glance(model)
          gd$estimator <- estimator
          gd
        }

        add_bayes_p_values <- function(coef_df, model) {
          pd_df <- as.data.frame(bayestestR::p_direction(model, effects = "fixed", as_p = FALSE))
          p_df <- as.data.frame(bayestestR::p_direction(model, effects = "fixed", as_p = TRUE))

          stopifnot(all(c("term") %in% names(coef_df)))
          stopifnot(all(c("Parameter", "pd") %in% names(pd_df)))
          stopifnot(all(c("Parameter", "p") %in% names(p_df)))

          coef_df$pd <- pd_df$pd[match(coef_df$term, pd_df$Parameter)]
          coef_df$p.value <- p_df$p[match(coef_df$term, p_df$Parameter)]
          coef_df
        }

        if (engine == "rstanarm") {
          needs_model <- any(c(
            "parameter_estimates_reml", "parameter_estimates_ml", "fit_statistics",
            "residuals", "fitted"
          ) %in% calculate) || !is.null(emmeans_spec) || !is.null(emtrends_spec) || isTRUE(return_models)

          if (isTRUE(needs_model)) {
            thism <- model_worker(ret$dt[[1]], ff, lmer_control, outcome_transform, scale_predictors,
                                  REML = NA, engine = engine, engine_args = engine_args, weights = weights)
          }

          if ("parameter_estimates_reml" %in% calculate) {
            coef_reml <- do.call(broom.mixed::tidy, c(list(x = thism), tidy_args))
            coef_reml <- add_bayes_p_values(coef_reml, thism)
            coef_reml$estimator <- estimator
            ret[, coef_df_reml := list(coef_reml)]
          }

          if ("parameter_estimates_ml" %in% calculate) {
            coef_ml <- do.call(broom.mixed::tidy, c(list(x = thism), tidy_args))
            coef_ml <- add_bayes_p_values(coef_ml, thism)
            coef_ml$estimator <- estimator
            ret[, coef_df_ml := list(coef_ml)]
          }

          if ("fit_statistics" %in% calculate) {
            ret[, fit_df := list(glance_model(thism, estimator))]
          }
        } else {
          if ("parameter_estimates_reml" %in% calculate) {
            thism <- model_worker(ret$dt[[1]], ff, lmer_control, outcome_transform, scale_predictors,
                                  REML = TRUE, engine = engine, engine_args = engine_args, weights = weights)
            coef_reml <- do.call(broom.mixed::tidy, c(list(x = thism), tidy_args))
            coef_reml$estimator <- "REML"
            ret[, coef_df_reml := list(coef_reml)]
          }

          if (any(c("parameter_estimates_ml", "fit_statistics") %in% calculate)) {
            thism_ml <- model_worker(ret$dt[[1]], ff, lmer_control, outcome_transform, scale_predictors,
                                     REML = FALSE, engine = engine, engine_args = engine_args, weights = weights) #refit with ML for AIC/BIC
            if ("fit_statistics" %in% calculate) {
              ret[, fit_df := list(glance_model(thism_ml, "ML"))]
            }
            if ("parameter_estimates_ml" %in% calculate) {
              coef_ml <- do.call(broom.mixed::tidy, c(list(x = thism_ml), tidy_args))
              coef_ml$estimator <- "ML"
              ret[, coef_df_ml := list(coef_ml)]
            }
          }

          needs_posthoc <- any(c("residuals", "fitted") %in% calculate) ||
            !is.null(emmeans_spec) || !is.null(emtrends_spec) || isTRUE(return_models)
          if (isTRUE(needs_posthoc) && is.null(thism) && is.null(thism_ml)) {
            thism <- model_worker(ret$dt[[1]], ff, lmer_control, outcome_transform, scale_predictors,
                                  REML = TRUE, engine = engine, engine_args = engine_args, weights = weights)
          }
        }

        posthoc_model <- if (!is.null(thism)) thism else thism_ml
        
        if ("residuals" %in% calculate) {
          ret[, residuals := list(resid = residuals(posthoc_model))]
        }
        
        if ("fitted" %in% calculate) {
          ret[, fitted := list(fitted = fitted(posthoc_model))]
        }
        
        if (isTRUE(return_models)) {
          model_payload <- if (engine == "rstanarm") {
            list(bayes = posthoc_model)
          } else {
            list(reml = thism, ml = thism_ml)
          }
          ret[, model := list(list(model_payload))]
        }
        
        # process emmeans
        if (!is.null(emmeans_spec)) {
          emm_torun <- emm_metadata[
            emm_metadata$outcome == model_set$outcome[[mm]] & emm_metadata$model_name == names(model_set$rhs)[mm],
            ,
            drop = FALSE
          ]
          
          if (nrow(emm_torun) > 0L) {
            this_emmspec <- emmeans_spec[emm_torun$emm_number] #subset list to relevant elements
            emm_runs <- lapply(seq_along(this_emmspec), function(emm_i) {
              .mixed_by_run_posthoc_spec(
                emmeans::emmeans, posthoc_model, this_emmspec[[emm_i]], "emmeans",
                "emmeans_pairs", emmeans_pairs, emmeans_pairs_args
              )
            })
            emms <- lapply(seq_along(emm_runs), function(emm_i) {
              cbind(
                emm_runs[[emm_i]]$main,
                emm_torun[emm_i, c("emm_number", "emm_label"), drop = FALSE]
              ) # don't add model and outcome, which will double at the unnest
            })
            names(emms) <- names(this_emmspec)
            
            #N.B. Need to double list wrap for data.table to keep list class for singleton list
            ret[, emm := list(list(emms))] 

            emps <- lapply(seq_along(emm_runs), function(emm_i) {
              if (is.null(emm_runs[[emm_i]]$pairs)) {
                return(NULL)
              }
              cbind(
                emm_runs[[emm_i]]$pairs,
                emm_torun[emm_i, c("emm_number", "emm_label"), drop = FALSE]
              )
            })
            if (any(!vapply(emps, is.null, logical(1)))) {
              names(emps) <- names(this_emmspec)
              ret[, emp := list(list(emps))]
            }
          }
        }

        # process emtrends
        if (!is.null(emtrends_spec)) {
          emt_torun <- emt_metadata[
            emt_metadata$outcome == model_set$outcome[[mm]] & emt_metadata$model_name == names(model_set$rhs)[mm],
            ,
            drop = FALSE
          ]

          if (nrow(emt_torun) > 0L) {
            this_emtspec <- emtrends_spec[emt_torun$emt_number] #subset list to relevant elements
            emt_runs <- lapply(seq_along(this_emtspec), function(emt_i) {
              .mixed_by_run_posthoc_spec(
                emmeans::emtrends, posthoc_model, this_emtspec[[emt_i]], "emtrends",
                "emtrends_pairs", emtrends_pairs, emtrends_pairs_args
              )
            })
            emts <- lapply(seq_along(emt_runs), function(emt_i) {
              cbind(
                emt_runs[[emt_i]]$main,
                emt_torun[emt_i, c("emt_number", "emt_label"), drop = FALSE]
              ) # don't add model and outcome, which will double at the unnest
            })
            names(emts) <- names(this_emtspec)

            #N.B. Need to double list wrap for data.table to keep list class for singleton list
            ret[, emt := list(list(emts))] 

            etps <- lapply(seq_along(emt_runs), function(emt_i) {
              if (is.null(emt_runs[[emt_i]]$pairs)) {
                return(NULL)
              }
              cbind(
                emt_runs[[emt_i]]$pairs,
                emt_torun[emt_i, c("emt_number", "emt_label"), drop = FALSE]
              )
            })
            if (any(!vapply(etps, is.null, logical(1)))) {
              names(etps) <- names(this_emtspec)
              ret[, etp := list(list(etps))]
            }
          }
        }

        ret[, dt := NULL] # drop original data.table for this split from data

        return(ret)
      })

      split_results <- data.table::rbindlist(split_results, fill=TRUE)

      #need to unnest

      coef_df_reml <- NULL
      result_keys <- c("outcome", "model_name", "rhs", "engine")
      if ("parameter_estimates_reml" %in% calculate) {
        coef_df_reml <- split_results[, coef_df_reml[[1]], by = result_keys] # unnest coefficients
        coef_df_reml <- cbind(dt_split[, ..split_on], coef_df_reml) # add back metadata for this split
      }

      coef_df_ml <- NULL
      if ("parameter_estimates_ml" %in% calculate) {
        coef_df_ml <- split_results[, coef_df_ml[[1]], by = result_keys] # unnest coefficients
        coef_df_ml <- cbind(dt_split[, ..split_on], coef_df_ml) # add back metadata for this split
      }

      fit_df <- NULL
      if ("fit_statistics" %in% calculate) {
        fit_df <- split_results[, fit_df[[1]], by = result_keys] # unnest coefficients
        fit_df <- cbind(dt_split[, ..split_on], fit_df) # add back metadata for this split
      }
      
      emm_data <- NULL
      if (!is.null(emmeans_spec) && "emm" %in% names(split_results)) {
        emm_data <- subset(split_results, sapply(emm, function(x) !is.null(x)), 
                           select=c(split_on, "outcome", "model_name", "rhs", "emm"))
      }

      emp_data <- NULL
      if (!is.null(emmeans_spec) && "emp" %in% names(split_results)) {
        emp_data <- subset(split_results, sapply(emp, function(x) !is.null(x)),
                           select=c(split_on, "outcome", "model_name", "rhs", "emp"))
      }

      emt_data <- NULL
      if (!is.null(emtrends_spec) && "emt" %in% names(split_results)) {
        emt_data <- subset(split_results, sapply(emt, function(x) !is.null(x)), 
                           select=c(split_on, "outcome", "model_name", "rhs", "emt"))
      }

      etp_data <- NULL
      if (!is.null(emtrends_spec) && "etp" %in% names(split_results)) {
        etp_data <- subset(split_results, sapply(etp, function(x) !is.null(x)),
                           select=c(split_on, "outcome", "model_name", "rhs", "etp"))
      }
      
      residual_data <- NULL
      if ("residuals" %in% calculate) {
        residual_data <- split_results[, residuals[[1]], by = .(outcome, model_name, rhs)] # unnest coefficients
        residual_data <- cbind(dt_split[, ..split_on], residual_data) # add back metadata for this split
      }
      
      fitted_data <- NULL
      if ("fitted" %in% calculate) {
        fitted_data <- split_results[, fitted[[1]], by = .(outcome, model_name, rhs)] # unnest coefficients
        fitted_data <- cbind(dt_split[, ..split_on], fitted_data) # add back metadata for this split
      }

      model_data <- NULL
      if (isTRUE(return_models) && "model" %in% names(split_results)) {
        model_data <- subset(split_results, select=c(split_on, "outcome", "model_name", "rhs", "engine", "model"))
      }

      list(coef_df_reml=coef_df_reml, coef_df_ml=coef_df_ml, fit_df=fit_df,
           emm_data=emm_data, emp_data=emp_data, emt_data=emt_data, etp_data=etp_data,
           residuals = residual_data, fitted = fitted_data,
           models = model_data) # return list
    }
    
    rm(data) # cleanup big datasets and force memory release before next iteration
    gc()
  }

  # need to put this above, but trying to avoid tmp objects
  # data[, filename:=df_i] #tag for later
  
  #helper subfunction to pull out and rbind element from doubly nested structure above
  extract_df <- function(nested_list, element, split_on) {
    result_df <- data.table::rbindlist(lapply(nested_list, function(df_set) {
      data.table::rbindlist(lapply(df_set, "[[", element))
    })) # combine results from each df (in the multiple df case)
    
    if (nrow(result_df) > 0L) { # only set keys if data were extracted (setorderv fails on 0-row df)
      data.table::setattr(result_df, "split_on", split_on) #tag split variables for secondary analysis
      data.table::setorderv(result_df, split_on) # since we allow out-of-order foreach, reorder coefs here.
    }
    return(result_df)
  }
  
  coef_results_reml <- NULL
  coef_results_ml <- NULL
  fit_results <- NULL
  fitted_data <- NULL
  residual_data <- NULL
  model_data <- NULL
  if ("parameter_estimates_reml" %in% calculate) { coef_results_reml <- extract_df(mresults, "coef_df_reml", split_on) }
  if ("parameter_estimates_ml" %in% calculate) { coef_results_ml <- extract_df(mresults, "coef_df_ml", split_on) }
  if ("fit_statistics" %in% calculate) { fit_results <- extract_df(mresults, "fit_df", split_on) }
  if ("fitted" %in% calculate) { fitted_data <- extract_df(mresults, "fitted", split_on) }
  if ("residuals" %in% calculate) { residual_data <- extract_df(mresults, "residuals", split_on) }
  if (isTRUE(return_models)) { model_data <- extract_df(mresults, "models", split_on) }
  
  emmeans_list <- .mixed_by_extract_posthoc_list(mresults, "emm_data", emmeans_spec, "emm", split_on, extract_df)
  emmeans_pairs_list <- .mixed_by_extract_posthoc_list(mresults, "emp_data", emmeans_spec, "emp", split_on, extract_df)
  emtrends_list <- .mixed_by_extract_posthoc_list(mresults, "emt_data", emtrends_spec, "emt", split_on, extract_df)
  emtrends_pairs_list <- .mixed_by_extract_posthoc_list(mresults, "etp_data", emtrends_spec, "etp", split_on, extract_df)
  
  #helper subfunction to adjust p-values
  adjust_dt <- function(dt, padjust_by) {
    if (!"p.value" %in% names(dt)) {
      warning("Skipping p-value adjustment because coefficient table has no p.value column.")
      return(dt)
    }
    for (ff in padjust_by) {
      checkmate::assert_subset(ff, names(dt))
      cname <- paste0("padj_", padjust_method, "_", paste(ff, collapse = "_"))
      dt <- dt[, (cname) := p.adjust(p.value, method = padjust_method), by = ff]
      if (isTRUE(all.equal(dt[["p.value"]], dt[[cname]]))) {
        warning("p-value adjustment: ", paste(ff, collapse=", "), 
                " yields the same result as the uncorrected p-value. Setting ",
                cname, " to NA.")
        dt[[cname]] <- NA_real_
      }
    }
    return(dt)
  }
  
  # compute adjusted p values
  if (!is.null(padjust_by)) {
    checkmate::assert_subset(padjust_method, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"))
    if (!is.null(coef_results_reml)) coef_results_reml <- adjust_dt(coef_results_reml, padjust_by)
    if (!is.null(coef_results_ml)) coef_results_ml <- adjust_dt(coef_results_ml, padjust_by)
  }

  # drop off dummy split if irrelevant
  if (isFALSE(has_split)) {
    if (!is.null(coef_results_reml)) { coef_results_reml[, split := NULL] }
    if (!is.null(coef_results_ml)) { coef_results_ml[, split := NULL] }
    if (!is.null(fit_results)) { fit_results[, split := NULL] }
    if (!is.null(model_data)) { model_data[, split := NULL] }
  }

  return(list(coef_df_reml=coef_results_reml, coef_df_ml=coef_results_ml, fit_df=fit_results, 
              emmeans_list=emmeans_list, emmeans_pairs_list=emmeans_pairs_list,
              emtrends_list=emtrends_list, emtrends_pairs_list=emtrends_pairs_list,
              residuals = residual_data,
              fitted = fitted_data, models = model_data))
}
