# Package index

## All functions

- [`R_batch_job`](https://hallquistlab.github.io/fmri.pipeline/reference/R_batch_job.md)
  : Description of R_batch_job R6 class

- [`R_batch_sequence`](https://hallquistlab.github.io/fmri.pipeline/reference/R_batch_sequence.md)
  : Description of R_batch_sequence R6 class

- [`add_tracked_job_parent()`](https://hallquistlab.github.io/fmri.pipeline/reference/add_tracked_job_parent.md)
  : Add parent/child id relationship to tracking database

- [`afni_3dclusterize`](https://hallquistlab.github.io/fmri.pipeline/reference/afni_3dclusterize.md)
  : wrapper class for 3dClusterize

- [`afni_3dclustsim`](https://hallquistlab.github.io/fmri.pipeline/reference/afni_3dclustsim.md)
  : R6 class for 3dClustSim automation

- [`afni_3dfwhmx`](https://hallquistlab.github.io/fmri.pipeline/reference/afni_3dfwhmx.md)
  : R6 class for running 3dFWHMx on a single input file based on user
  specification

- [`afni_3dfwhmx_list`](https://hallquistlab.github.io/fmri.pipeline/reference/afni_3dfwhmx_list.md)
  : R6 class for running 3dFWHMx on a group of input files using a
  scheduler/cluster

- [`afni_whereami`](https://hallquistlab.github.io/fmri.pipeline/reference/afni_whereami.md)
  : wrapper class for AFNI whereami

- [`as.character(`*`<formula>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/as.character.formula.md)
  : S3 method to convert formula to character

- [`build_design_matrix()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_design_matrix.md)
  : Creates an fmri design matrix, including timing files for for AFNI
  or FSL.

- [`build_fwe_correction()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_fwe_correction.md)
  : Function to walk user through setting up FWE corrections for model
  outputs

- [`build_l1_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_l1_models.md)
  : Interactive function to build an l1 model specification for
  setup_glm_pipeline

- [`build_l2_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_l2_models.md)
  : Interactive function to build an l2 model specification for
  setup_glm_pipeline

- [`build_l3_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_l3_models.md)
  : Interactive function to build an l3 model specification for
  setup_glm_pipeline

- [`cleanup_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/cleanup_glm_pipeline.md)
  : helper function to refresh l3 model status and save gpa object from
  batch pipeline back to its cache

- [`cluster_job_submit()`](https://hallquistlab.github.io/fmri.pipeline/reference/cluster_job_submit.md)
  : This function submits a single script to a high-performance cluster
  using a scheduler (Slurm or TORQUE). It accepts a vector of arguments
  to be passed to the scheduler and a vector of environment variables
  that should be passed to the compute node at job execution.

- [`cluster_submit_shell_jobs()`](https://hallquistlab.github.io/fmri.pipeline/reference/cluster_submit_shell_jobs.md)
  : helper function to submit a set of shell jobs that are independent
  of one another

- [`collect_fwe_results()`](https://hallquistlab.github.io/fmri.pipeline/reference/collect_fwe_results.md)
  : Collect FWE results and persist the artifact manifest

- [`collinearity_diagnostics`](https://hallquistlab.github.io/fmri.pipeline/reference/collinearity_diagnostics.md)
  : Functions for extracting and summarizing collinearity diagnostics
  across the GLM pipeline

- [`combine_feat_l3_to_afni()`](https://hallquistlab.github.io/fmri.pipeline/reference/combine_feat_l3_to_afni.md)
  : Function to combine Feat L3 analyses into AFNI BRIK+HEAD files
  according to user-specified combinations

- [`combine_spm_l3_to_afni()`](https://hallquistlab.github.io/fmri.pipeline/reference/combine_spm_l3_to_afni.md)
  : Combine SPM L3 outputs into AFNI BRIK+HEAD files for visualization

- [`compress_mts_pca()`](https://hallquistlab.github.io/fmri.pipeline/reference/compress_mts_pca.md)
  : Core function for computing multivariate time series compression
  scores by principal components

- [`concat_design_runs()`](https://hallquistlab.github.io/fmri.pipeline/reference/concat_design_runs.md)
  : Concatenate design matrices for each run to form a single design
  with unique baselines per run (ala AFNI)

- [`convolve_cpp`](https://hallquistlab.github.io/fmri.pipeline/reference/convolve_cpp.md)
  : Internal function to convolve two vectors using nested for loop

- [`convolve_double_gamma`](https://hallquistlab.github.io/fmri.pipeline/reference/convolve_double_gamma.md)
  : This function convolves a stimulus vector with the double-gamma hrf

- [`create_fwe_spec()`](https://hallquistlab.github.io/fmri.pipeline/reference/create_fwe_spec.md)
  : helper function to create an FWE object that specifies the type of
  FWE correction and the model outputs to which it is applied

- [`deconvolve_nlreg`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_nlreg.md)
  : C++ port of Bush and Cisler 2013, Magnetic Resonance Imaging Adapted
  from the original provided by Keith Bush as well as C++ code from
  Jiang Bian

- [`deconvolve_nlreg_resample()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_nlreg_resample.md)
  : This function deconvolves the BOLD signal using Bush 2011 method,
  augmented by the resampling approach of Bush 2015.

- [`deconvolve_reglin()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_reglin.md)
  : Regularized linear deconvolution of BOLD time series

- [`detrendts()`](https://hallquistlab.github.io/fmri.pipeline/reference/detrendts.md)
  : Detrend a time series up to quadratic trend. Used by fir1Bandpass
  prior to filtering

- [`diagnose_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/diagnose_pipeline.md)
  : Function for diagnosing errors in a run of the pipeline

- [`do_convolve`](https://hallquistlab.github.io/fmri.pipeline/reference/do_convolve.md)
  : Internal port of fsl do_convolve

- [`event_lock_ts()`](https://hallquistlab.github.io/fmri.pipeline/reference/event_lock_ts.md)
  : function to get interpolated event locked data

- [`evt_align_decon_files()`](https://hallquistlab.github.io/fmri.pipeline/reference/evt_align_decon_files.md)
  : Align a set of deconvolved time series files to an event of interest
  function.

- [`export_collinearity_to_csv()`](https://hallquistlab.github.io/fmri.pipeline/reference/export_collinearity_to_csv.md)
  : Export collinearity diagnostics to CSV files

- [`export_glm_config()`](https://hallquistlab.github.io/fmri.pipeline/reference/export_glm_config.md)
  : Output GLM configuration to a '.yaml' or 'json' file

- [`extract_glm_betas_in_mask()`](https://hallquistlab.github.io/fmri.pipeline/reference/extract_glm_betas_in_mask.md)
  : Extract FSL GLM statistics within spatially aligned masks

- [`extract_l1_collinearity()`](https://hallquistlab.github.io/fmri.pipeline/reference/extract_l1_collinearity.md)
  : Extract collinearity diagnostics from all L1 models

- [`fields_from_spec()`](https://hallquistlab.github.io/fmri.pipeline/reference/fields_from_spec.md)
  : internal function to populate l1_model_set onsets from spec
  specification

- [`fill_atlas_with_stats()`](https://hallquistlab.github.io/fmri.pipeline/reference/fill_atlas_with_stats.md)
  : Populates parcel-wise statistics from an atlas or set of clusters
  back into NIfTIs in the original space

- [`finalize_pipeline_configuration()`](https://hallquistlab.github.io/fmri.pipeline/reference/finalize_pipeline_configuration.md)
  : This is a small helper function to validate the glm_model_arguments
  list structure. It adds a few details such as the output directory to
  make it less burdensome for to setup a pipeline N.B. gpa is a
  shorthand abbreviation for glm_model_arguments, to save typing

- [`fmri.stimulus()`](https://hallquistlab.github.io/fmri.pipeline/reference/fmri.stimulus.md)
  : Convolve a regressor with a hemodynamic response function for fMRI
  analysis.

- [`fmri_ts`](https://hallquistlab.github.io/fmri.pipeline/reference/fmri_ts.md)
  : R6 class representing a multivariate time series object for fMRI
  analysis

- [`fsl_generate_fsf_contrast_syntax()`](https://hallquistlab.github.io/fmri.pipeline/reference/fsl_generate_fsf_contrast_syntax.md)
  : This function generates syntax for FSL Feat .fsf files for the
  contrasts tab of an fMRI analysis. It accepts a numeric contrast
  matrix whose rownames correspond to the name of the contrast.

- [`fsl_generate_fsf_ev_syntax()`](https://hallquistlab.github.io/fmri.pipeline/reference/fsl_generate_fsf_ev_syntax.md)
  : This function generates syntax for FSL Feat .fsf files for the EVs
  tab of a higher-level fMRI analysis. It accepts a numeric design
  matrix whose colum names correspond to individual EVs in the model

- [`fsl_generate_fsf_lvl1_ev_syntax()`](https://hallquistlab.github.io/fmri.pipeline/reference/fsl_generate_fsf_lvl1_ev_syntax.md)
  : This function generates syntax for FSL Feat .fsf files for the EVs
  tab of a first-level fMRI analysis. It accepts a numeric design matrix
  whose colum names correspond to individual EVs in the model

- [`fsl_l1_model()`](https://hallquistlab.github.io/fmri.pipeline/reference/fsl_l1_model.md)
  : This function specifies a first-level GLM model in FSL based on a
  build_design_matrix bdm object.

- [`fwe_result_artifacts()`](https://hallquistlab.github.io/fmri.pipeline/reference/fwe_result_artifacts.md)
  : Select artifacts from a collected FWE result

- [`fwe_spec()`](https://hallquistlab.github.io/fmri.pipeline/reference/fwe_spec.md)
  : Create a familywise-error correction specification

- [`generate_feature`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_feature.md)
  : This function creates K shifts of a neural events vector according
  to the kernel length, K.

- [`generate_feature_armadillo`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_feature_armadillo.md)
  : This function creates K shifts of a neural events vector according
  to the kernel length, K.

- [`generate_mask_diagnostics()`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_mask_diagnostics.md)
  : Generate a set of diagnostic images and plots for brain masks
  relative to one another and a template

- [`generate_run_data_from_bids()`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_run_data_from_bids.md)
  : Function to generate a run_data object from a BIDS-compliant folder

- [`generate_spm_contrasts()`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_spm_contrasts.md)
  : This function reads an SPM.mat file and generate contrasts based on
  the design matrix specification

- [`generate_spm_contrasts_from_model()`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_spm_contrasts_from_model.md)
  : Generate SPM contrasts from an fmri.pipeline L1 model contrast
  matrix

- [`generate_spm_mat()`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_spm_mat.md)
  : Uses an object from build_design_matrix to generate a corresponding
  SPM GLM design matrix.

- [`generate_subject_data_from_bids()`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_subject_data_from_bids.md)
  : Function to generate a subject_data object from a BIDS-compliant
  folder

- [`generate_trial_data_from_bids()`](https://hallquistlab.github.io/fmri.pipeline/reference/generate_trial_data_from_bids.md)
  : Function to generate a trial_data object from a BIDS-compliant
  folder

- [`get_medusa_compression_score()`](https://hallquistlab.github.io/fmri.pipeline/reference/get_medusa_compression_score.md)
  : front-end function for taking a list of windowed time series by mask
  value, interpolating them onto a time grid, and (optionally) averaging
  across voxels/units within a value to derive the mean interpolated
  time series

- [`get_medusa_interpolated_ts()`](https://hallquistlab.github.io/fmri.pipeline/reference/get_medusa_interpolated_ts.md)
  : front-end function for taking a list of windowed time series by mask
  value, interpolating them onto a time grid, and (optionally) averaging
  across voxels/units within a value to derive the mean interpolated
  time series

- [`get_output_directory()`](https://hallquistlab.github.io/fmri.pipeline/reference/get_output_directory.md)
  : small helper function to return the location of an l1 directory
  based on id, session, and run number

- [`get_spm_status()`](https://hallquistlab.github.io/fmri.pipeline/reference/get_spm_status.md)
  : Helper to check whether expected SPM outputs exist

- [`get_tracked_job_status()`](https://hallquistlab.github.io/fmri.pipeline/reference/get_tracked_job_status.md)
  : Query job status in tracking SQLite database

- [`insert_df_sqlite()`](https://hallquistlab.github.io/fmri.pipeline/reference/insert_df_sqlite.md)
  : helper function to insert a keyed data.frame into the sqlite storage
  database

- [`jlmer()`](https://hallquistlab.github.io/fmri.pipeline/reference/jlmer.md)
  : Fit a linear mixed model with Julia MixedModels

- [`lookup_feat_outputs()`](https://hallquistlab.github.io/fmri.pipeline/reference/lookup_feat_outputs.md)
  [`as.data.frame(`*`<feat_output_lookup>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/lookup_feat_outputs.md)
  [`print(`*`<feat_output_lookup>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/lookup_feat_outputs.md)
  : Lookup FSL FEAT output images by model and contrast

- [`meta_mixed_by()`](https://hallquistlab.github.io/fmri.pipeline/reference/meta_mixed_by.md)
  : Function to run Bayesian random-effects meta-regression on
  coefficients from mixed_by

- [`mixed_by()`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md)
  : Mixed by runs a set of mixed-effects models for each combination of
  a set of factors. Its primary use is to run the same model on
  different splits of the data.

- [`plan_fwe_correction()`](https://hallquistlab.github.io/fmri.pipeline/reference/plan_fwe_correction.md)
  : Compile an FWE specification into an executable plan

- [`plot_collinearity()`](https://hallquistlab.github.io/fmri.pipeline/reference/plot_collinearity.md)
  : Create a visualization of collinearity diagnostics

- [`populate_defaults()`](https://hallquistlab.github.io/fmri.pipeline/reference/populate_defaults.md)
  : helper function to copy any missing fields in target

- [`print(`*`<l1_collinearity_summary>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/print.l1_collinearity_summary.md)
  : Print method for l1_collinearity_summary objects

- [`ptfce_spec`](https://hallquistlab.github.io/fmri.pipeline/reference/ptfce_spec.md)
  : R6 class for pTFCE specification for one or more z-statistic images

- [`read_confounds_motion_parameters()`](https://hallquistlab.github.io/fmri.pipeline/reference/read_confounds_motion_parameters.md)
  : helper function to read confounds and motion parameters and combine
  them

- [`read_feat_dir()`](https://hallquistlab.github.io/fmri.pipeline/reference/read_feat_dir.md)
  : helper function to look up core stats outputs from a .feat folder

- [`read_feat_inputs()`](https://hallquistlab.github.io/fmri.pipeline/reference/read_feat_inputs.md)
  : This is a small helper function that returns the inputs provided in
  the feat_files field for a set of .fsf files.

- [`read_fwe_result()`](https://hallquistlab.github.io/fmri.pipeline/reference/read_fwe_result.md)
  [`read_fwe_artifact_manifest()`](https://hallquistlab.github.io/fmri.pipeline/reference/read_fwe_result.md)
  : Read a collected FWE result or artifact manifest

- [`read_gfeat_dir()`](https://hallquistlab.github.io/fmri.pipeline/reference/read_gfeat_dir.md)
  : helper function to look up core stats outputs from a .gfeat folder

- [`refresh_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/refresh_fwe_plan.md)
  : Refresh the state of a compiled FWE plan

- [`refresh_glm_status()`](https://hallquistlab.github.io/fmri.pipeline/reference/refresh_glm_status.md)
  : Refresh GLM backend status for a given analysis level

- [`resolve_fwe_targets()`](https://hallquistlab.github.io/fmri.pipeline/reference/resolve_fwe_targets.md)
  : Resolve an FWE specification to concrete group-analysis targets

- [`rle_dt`](https://hallquistlab.github.io/fmri.pipeline/reference/rle_dt.md)
  : R6 class for a keyed data.table object that uses run-length encoding
  to reduce storage size

- [`run_3dlmer_sepjobs()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_3dlmer_sepjobs.md)
  : function to submit 3dLMEr jobs on a cluster

- [`run_afni_command()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_afni_command.md)
  : Wrapper for running an AFNI command safely within R

- [`run_decon_alignment()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_decon_alignment.md)
  :

  Align a set of deconvolved time series files from
  `voxelwise_deconvolution` to a set of events in a task-related design

- [`run_feat_sepjobs()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_feat_sepjobs.md)
  : function to submit a set of jobs on a cluster to estimate many Feat
  level 1 models

- [`run_fsl_command()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_fsl_command.md)
  : Wrapper for running an FSL command safely within R

- [`run_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_fwe_plan.md)
  : Execute a compiled FWE plan

- [`run_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_glm_pipeline.md)
  : primary function for running a GLM analysis pipeline

- [`run_spm_sepjobs()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_spm_sepjobs.md)
  : Submit SPM jobs to the cluster for unattended MATLAB/Octave
  execution

- [`setup_compute_environment()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_compute_environment.md)
  : Helper function for setting up and modifying the compute environment
  for scripts generated by the pipeline

- [`setup_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md)
  : Main worker function for setting up an analysis pipeline

- [`setup_l1_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_l1_models.md)
  : Function for setting up level 1 model specifications and
  corresponding files

- [`setup_l2_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_l2_models.md)
  : This function generates the inputs for an FSL level 2 analysis,
  where multiple runs for a subject are combined using fixed effects
  estimation.

- [`setup_l3_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_l3_models.md)
  : This function generates the inputs for level 3 analyses, where
  multi-subject data are analyzed in group analyses.

- [`spm_extract_anatomical_rois()`](https://hallquistlab.github.io/fmri.pipeline/reference/spm_extract_anatomical_rois.md)
  : This is a wrapper around the spm_extract_anatomical_rois.m script in
  the inst directory

- [`spm_l1_model()`](https://hallquistlab.github.io/fmri.pipeline/reference/spm_l1_model.md)
  : Setup a first-level SPM model using a build_design_matrix object

- [`spm_l3_model()`](https://hallquistlab.github.io/fmri.pipeline/reference/spm_l3_model.md)
  : Setup a third-level SPM model using lower-level con images

- [`subset_nifti_volumes`](https://hallquistlab.github.io/fmri.pipeline/reference/subset_nifti_volumes.md)
  : Subset Timepoints from a 4D NIfTI Image

- [`summarize_correlations_by_pair()`](https://hallquistlab.github.io/fmri.pipeline/reference/summarize_correlations_by_pair.md)
  : Summarize correlations by regressor pair across all subjects

- [`summarize_vif_by_regressor()`](https://hallquistlab.github.io/fmri.pipeline/reference/summarize_vif_by_regressor.md)
  : Summarize VIF by regressor across all subjects

- [`summary(`*`<bg_block_data>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/summary.bg_block_data.md)
  : Summarize the contents of the subject data

- [`summary(`*`<bg_subject_data>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/summary.bg_subject_data.md)
  : Summarize the contents of the subject data

- [`summary(`*`<glm_pipeline_arguments>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/summary.glm_pipeline_arguments.md)
  : Helper function to better summarize GPA object

- [`summary(`*`<l1_model_set>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/summary.l1_model_set.md)
  : Helper function to better summarize GPA object's L1 models

- [`summary(`*`<l1_model_set_events>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/summary.l1_model_set_events.md)
  : Helper function to better summarize GPA object's l1_models' events

- [`summary(`*`<l1_model_set_models>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/summary.l1_model_set_models.md)
  : Helper function to better summarize GPA object's l1_models' models

- [`summary(`*`<l1_model_set_signals>`*`)`](https://hallquistlab.github.io/fmri.pipeline/reference/summary.l1_model_set_signals.md)
  : Helper function to better summarize GPA object's l1_models' signals

- [`update_tracked_job_status()`](https://hallquistlab.github.io/fmri.pipeline/reference/update_tracked_job_status.md)
  : Update job status in tracking SQLite database

- [`validate_fwe_spec()`](https://hallquistlab.github.io/fmri.pipeline/reference/validate_fwe_spec.md)
  : Validate an FWE correction specification

- [`view_log()`](https://hallquistlab.github.io/fmri.pipeline/reference/view_log.md)
  : Function for quick viewing of log files in the command line

- [`visualize_design_matrix()`](https://hallquistlab.github.io/fmri.pipeline/reference/visualize_design_matrix.md)
  : Visualize design matrix, including event onset times and run
  boundaries

- [`voxelwise_deconvolution()`](https://hallquistlab.github.io/fmri.pipeline/reference/voxelwise_deconvolution.md)
  : Function to perform voxelwise deconvolution on an fMRI dataset using
  the fMRI model arguments object

- [`wait_for_job()`](https://hallquistlab.github.io/fmri.pipeline/reference/wait_for_job.md)
  : This function pauses execution of an R script while a scheduled qsub
  job is not yet complete.

- [`write_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/write_fwe_plan.md)
  [`read_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/write_fwe_plan.md)
  : Write or read a compiled FWE plan

- [`write_fwe_spec()`](https://hallquistlab.github.io/fmri.pipeline/reference/write_fwe_spec.md)
  [`read_fwe_spec()`](https://hallquistlab.github.io/fmri.pipeline/reference/write_fwe_spec.md)
  : Write or read an FWE correction specification as YAML
