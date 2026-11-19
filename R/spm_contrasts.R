#' This function reads an SPM.mat file and generate contrasts based on the design matrix specification
#'
#' @param output_dir location for SPM outputs and scripts for estimating contrasts. Must contain SPM.mat already
#' @param condition_contrasts see generate_spm_mat
#' @param unit_contrasts see generate_spm_mat
#' @param effects_of_interest_F see generate_spm_mat
#' @param spm_path see generate_spm_mat
#' @param execute whether to run contrast setup. This depends on SPM.mat having been created already. Default: FALSE
#' @param matlab_cmd see generate_spm_mat
#' @param matlab_args see generate_spm_mat
#' 
#' @importFrom R.matlab readMat
#' @author Michael Hallquist
#' @export
generate_spm_contrasts <- function(output_dir, condition_contrasts=TRUE, unit_contrasts=TRUE, effects_of_interest_F=TRUE,
                                   spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12", execute=FALSE,
                                   matlab_cmd="matlab", matlab_args="-batch") {

  if (missing(output_dir)) { stop("No output_dir provided. This must be the folder containing the SPM.mat file") }
  if (execute && !file.exists(file.path(output_dir, "SPM.mat"))) { stop("No SPM.mat file found in: ", output_dir, ". This must be setup prior to estimating contrasts.") }

  if (condition_contrasts==FALSE && unit_contrasts==FALSE && effects_of_interest_F==FALSE) {
    message("No contrasts were requested by call to generate_spm_contrasts. As a result, nothing will occur here.")
    return(NULL)
  }

  spm_preamble <- c(
    ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
    "spm('defaults', 'fmri');",
    "spm_jobman('initcfg');",
    ""    
  )
  
  #get the names of all regressors
  mtmp <- file.path(output_dir, "extract_design_columns.m")
  mattmp <- file.path(output_dir, "design_columns.mat")
  mnames <- c(
    spm_preamble,
    paste0("load([ '", output_dir, "' filesep 'SPM.mat']);"),
    paste0("mnames=SPM.xX.name(:);"),
    paste0("cpos=SPM.xX.iC; %contrasts of design"),
    paste0("bpos=SPM.xX.iB; %block/run regressors"),
    paste0("npos=SPM.xX.iG; %nuisance regressors"),        
    "design_matrix = SPM.xX.X;",
    "if isfield(SPM.xX, 'xKXs'), design_matrix = SPM.xX.xKXs.X; end;",
    paste0("save('", mattmp, "', 'mnames', 'cpos', 'bpos', 'npos', 'design_matrix', '-v6');")
  )
  cat(mnames, file=mtmp, sep="\n")

  matlab_expression <- paste0("try; run('", gsub("'", "''", mtmp, fixed = TRUE),
    "'); catch ME; disp(getReport(ME,'extended')); exit(1); end; exit(0);")
  extract_cmd <- paste(matlab_cmd, matlab_args, shQuote(matlab_expression))
  if (execute) {
    if (system(extract_cmd) != 0L) stop("SPM design extraction failed.")
  }

  setup_script <- system.file("Rscript", "setup_spm_contrasts.R", package="fmri.pipeline")
  stopifnot(file.exists(setup_script))

  rscript_cmd <- paste0(shQuote(rscript_executable()),
    " --no-save --no-restore ", shQuote(setup_script),
    " -mat_file ", shQuote(mattmp),
    " -condition_contrasts ", as.character(condition_contrasts),
    " -unit_contrasts ", as.character(unit_contrasts),
    " -effects_of_interest_F ", as.character(effects_of_interest_F),
    if (!is.null(spm_path)) paste0(" -spm_path ", shQuote(spm_path)) else "")

  if (execute) {
    if (system(rscript_cmd) != 0L) stop("SPM contrast mapping failed.")
  }

  return(list(extract_cmd=extract_cmd, setup_cmd=rscript_cmd))
}

#' Generate SPM contrasts from an fmri.pipeline L1 model contrast matrix
#'
#' @param output_dir location for SPM outputs and scripts for estimating contrasts. Must contain SPM.mat already
#' @param mobj l1_model_spec object containing $contrasts and regressor names
#' @param spm_path see generate_spm_mat
#' @param execute whether to run contrast setup. This depends on SPM.mat having been created already. Default: FALSE
#' @param matlab_cmd see generate_spm_mat
#' @param matlab_args see generate_spm_mat
#' @param average_across_runs whether to average regressor weights across run-specific columns in SPM.xX
#' @param projection_interaction_terms optional character vector of projected L2-by-L1 interaction term names
#' @param projection_interaction_contrast_modes optional character vector. Allowed values are
#'   'pooled', 'session_specific', and 'session_differences'
#' @param projection_interaction_run_labels optional run labels aligned to SPM sessions for naming
#' @param projection_main_effect_terms optional character vector of projected L2 term names to
#'   estimate in contrast-space against each L1 contrast row
#' @param projection_main_effect_weights optional data.frame or matrix containing projected
#'   L2 term values aligned to SPM session/run order
#' 
#' @export
generate_spm_contrasts_from_model <- function(output_dir, mobj,
                                              spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12",
                                              execute=FALSE, matlab_cmd="matlab", matlab_args="-batch",
                                              average_across_runs=TRUE,
                                              projection_interaction_terms=character(0),
                                              projection_interaction_contrast_modes=character(0),
                                              projection_interaction_run_labels=NULL,
                                              projection_main_effect_terms=character(0),
                                              projection_main_effect_weights=NULL) {
  if (missing(output_dir)) { stop("No output_dir provided. This must be the folder containing the SPM.mat file") }
  if (execute && !file.exists(file.path(output_dir, "SPM.mat"))) { stop("No SPM.mat file found in: ", output_dir, ". This must be setup prior to estimating contrasts.") }

  checkmate::assert_multi_class(mobj, c("l1_model_spec", "l1_wi_spec", "hi_model_spec"))
  if (is.null(mobj$contrasts) || !inherits(mobj$contrasts, "matrix") || nrow(mobj$contrasts) == 0L) {
    stop("Model contrasts are missing or empty for SPM contrast generation.")
  }

  cmat <- mobj$contrasts
  if (is.null(colnames(cmat))) {
    stop("Model contrast matrix is missing regressor column names.")
  }

  if (!is.null(projection_main_effect_weights) &&
      !is.data.frame(projection_main_effect_weights) &&
      !is.matrix(projection_main_effect_weights)) {
    stop("projection_main_effect_weights must be NULL, a data.frame, or a matrix.")
  }

  if (is.matrix(projection_main_effect_weights)) {
    projection_main_effect_weights <- as.data.frame(
      projection_main_effect_weights,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (is.data.frame(projection_main_effect_weights)) {
    projection_main_effect_weights[] <- lapply(projection_main_effect_weights, as.numeric)
  }

  contrast_spec <- list(
    contrast_matrix = cmat,
    regressors = colnames(cmat),
    average_across_runs = average_across_runs,
    projection_interaction_terms = unique(as.character(projection_interaction_terms)),
    projection_interaction_contrast_modes = unique(as.character(projection_interaction_contrast_modes)),
    projection_interaction_run_labels = projection_interaction_run_labels,
    projection_main_effect_terms = unique(as.character(projection_main_effect_terms)),
    projection_main_effect_weights = projection_main_effect_weights
  )
  spec_path <- file.path(output_dir, "spm_contrast_spec.rds")
  saveRDS(contrast_spec, file = spec_path)

  spm_preamble <- c(
    ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
    "spm('defaults', 'fmri');",
    "spm_jobman('initcfg');",
    ""
  )

  mtmp <- file.path(output_dir, "extract_design_columns.m")
  mattmp <- file.path(output_dir, "design_columns.mat")
  mnames <- c(
    spm_preamble,
    paste0("load([ '", output_dir, "' filesep 'SPM.mat']);"),
    paste0("mnames=SPM.xX.name(:);"),
    paste0("cpos=SPM.xX.iC; %contrasts of design"),
    paste0("bpos=SPM.xX.iB; %block/run regressors"),
    paste0("npos=SPM.xX.iG; %nuisance regressors"),
    "design_matrix = SPM.xX.X;",
    "if isfield(SPM.xX, 'xKXs'), design_matrix = SPM.xX.xKXs.X; end;",
    paste0("save('", mattmp, "', 'mnames', 'cpos', 'bpos', 'npos', 'design_matrix', '-v6');")
  )
  cat(mnames, file = mtmp, sep = "\n")

  matlab_expression <- paste0("try; run('", gsub("'", "''", mtmp, fixed = TRUE),
    "'); catch ME; disp(getReport(ME,'extended')); exit(1); end; exit(0);")
  extract_cmd <- paste(matlab_cmd, matlab_args, shQuote(matlab_expression))

  setup_script <- system.file("Rscript", "setup_spm_contrasts_from_model.R", package="fmri.pipeline")
  stopifnot(file.exists(setup_script))

  rscript_cmd <- paste0(shQuote(rscript_executable()),
    " --no-save --no-restore ", shQuote(setup_script),
    " -mat_file ", shQuote(mattmp),
    " -contrast_rds ", shQuote(spec_path),
    " -average_across_runs ", as.character(average_across_runs),
    if (!is.null(spm_path)) paste0(" -spm_path ", shQuote(spm_path)) else "")

  if (execute) {
    if (system(extract_cmd) != 0L) stop("SPM design extraction failed.")
    if (system(rscript_cmd) != 0L) stop("SPM contrast mapping failed.")
  }

  return(list(extract_cmd=extract_cmd, setup_cmd=rscript_cmd))
}
