#' @details Start with \code{\link{setup_glm_pipeline}} to select the GLM engine
#'   and provide the analysis tables. Use \code{\link{run_glm_pipeline}} to submit
#'   the selected models. For a worked comparison of how each suite combines
#'   runs and supplies group contrasts, see \code{vignette("glm-engines", package = "fmri.pipeline")}.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @useDynLib fmri.pipeline, .registration = TRUE
## usethis namespace: end
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ":=", ".", "..split_on", ".filename", ".I", ".SD",
    "afni_complete", "boundary", "child_level",
    "dd", "dt_split", "emm", "emt", "emp", "etp", "event_list",
    "evt_onset", "evt_time", "expect_complete",
    "expect_func_dir", "expect_func_file", "expect_mprage_bet",
    "expect_mprage_dir", "expect_mprage_file", "expect_mprage_warpcoef",
    "feat_complete", "feat_dir_exists",
    "fname", "Freq", "grouping_scope",
    "i", "id", "ii", "img", "img_chunk", "img_exists", "img_stats", "is_complete",
    "isi", "iti_post", "l1_model", "l2_model", "l3_model", "label",
    "lg", "mask_value",
    "model", "model_matrix",
    "model_name", "motion_params_file", "name", "outcome", "p.value",
    "pct_true", "rhs",
    "run", "run_nifti", "run_nifti_present", "sequence_id",
    "session", "si", "STAT", "term", "time_started", "to_run", "trial",
    "value", "volume", "vox_ts", "x", "y", "z"
  ))
}
