# Familywise error correction for voxelwise GLMs in fmri.pipeline

The [whole-brain
walkthrough](https://hallquistlab.github.io/fmri.pipeline/articles/whole-brain.md)
produces a group z-statistic map for the participant-level
faces-minus-houses contrast. This guide corrects that map for the many
spatial tests inside its model mask. You can register pTFCE when
specifying the GLM, run it after estimation, or use the standalone AFNI
tools to make cluster-level inferences.

| Workflow | Inputs | What to inspect |
|----|----|----|
| pTFCE, automatically after L3 or from a saved analysis | Completed FSL L3 z-statistic, model mask, and FEAT smoothness file | Enhanced statistic, FWE threshold table, and a signed thresholded image for each requested alpha |
| AFNI spatial ACF and `3dClustSim` | First-level residuals and corresponding masks, group statistic, and group search mask | Cluster extent threshold, surviving cluster table, and integer cluster labels |
| AFNI residual sign randomization and `3dClustSim` | Residual maps from an appropriate group test and its search mask | Cluster thresholds calibrated from randomized group residuals |

Specify the search region, alpha, and sidedness before looking at the
result. A whole-brain correction covers the voxels in the model mask.
Each correction controls the spatial family for one map. If you make a
single claim across several contrasts, models, or visits, account for
those additional tests in the analysis plan. pTFCE uses spatial support
to enhance the statistic without a single cluster-forming threshold.
AFNI extent correction uses a specified voxelwise threshold and makes an
inference about the extent of the resulting clusters. Report the method
that produced the result. [Spisák et
al. (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6834440/), [Cox et
al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5399747/)

## Correct a completed group analysis

The declarative interface uses the completed `gpa` from the [batch
cache](https://hallquistlab.github.io/fmri.pipeline/articles/whole-brain.html#run-and-reload-the-analysis).
The package currently supports planned and automatic pTFCE for FSL L3.
Install `fmri.pipeline` and the HallquistLab pTFCE fork (version 0.4.0
or later) on the host that runs correction. The tested revision is
pinned in the package’s `DESCRIPTION` under `Remotes`. Verify that the R
interpreter used by your job can load both packages.

### Specify the target and search region

Select the map by session, model, and contrast at each level. Here the
two higher-level `~ 1` models carry the lower-level comparison forward.
The statistical question is the group mean face-minus-house effect.

``` r

library(fmri.pipeline)

spec <- fwe_spec(
  name = "facehouse_ptfce",
  method = "ptfce",
  targets = list(
    session = 1L,
    l1_model = "facehouse",
    l1_cope_name = "faces_minus_houses",
    l2_model = "run_average", l2_cope_name = "EV_Intercept",
    l3_model = "group_mean", l3_cope_name = "EV_Intercept"
  ),
  correction_mask = "model",
  fwe_alpha = c(0.05, 0.01),
  two_sided = TRUE,
  write_threshold_images = TRUE
)
spec
```

With `two_sided = TRUE`, pTFCE enhances both tails on a shared threshold
grid and uses half of the requested alpha per tail. Positive surviving
values favor faces in this example and negative values favor houses. Use
`two_sided = FALSE` for a pre-specified positive-tail test. A planned
test in the opposite direction can use an explicitly reversed contrast.

For a pre-specified small-volume correction, set `correction_mask` to
one nonempty binary NIfTI that is inside and aligned with each target’s
model mask. pTFCE reuses the model smoothness estimate and recalculates
the search volume and number of resels for the smaller region. Report
this as a small-volume result. Atlas masks used to summarize beta
estimates can have multiple labels, whereas the correction mask is
binary.

Keep the correction specification with the analysis script. A saved YAML
specification can be compiled again after refitting.

``` r

write_fwe_spec(spec, "config/facehouse_ptfce.yaml")
spec <- read_fwe_spec("config/facehouse_ptfce.yaml")
```

### Inspect the plan and run it

Planning resolves the named contrasts to files. The default output
directory is `<analysis>/fwe/<spec name>/`. Planning itself creates no
files and runs no jobs.

``` r

fwe_plan <- plan_fwe_correction(gpa, spec, refresh_status = TRUE)
fwe_plan

fwe_plan$targets[, c("target_id", "session", "l1_cope_name", "l3_cope_name")]
fwe_plan$tasks[, c(
  "task_id", "target_id", "zstat_file",
  "correction_mask_file", "smoothness_file", "status"
)]
fwe_plan$artifacts[, c("task_id", "artifact_role", "fwe_alpha", "artifact_file")]
stopifnot(nrow(fwe_plan$tasks) == 1L)
```

Expect one task for the face/house example. For a zero-match error, look
up the maps with `lookup_feat_outputs(gpa, level = 3L)$l3` and compare
the semantic labels and session to `spec$targets`. If you expect the
results but some files are still missing, use `require_inputs = FALSE`
on the plan to inspect the `status` and `missing_inputs` columns. A
blocked task cannot run until its inputs are available.

``` r

preview <- run_fwe_plan(fwe_plan, dry_run = TRUE)
preview$execution[, c("task_id", "execution_status", "command")]

write_fwe_plan(fwe_plan, "config/facehouse_ptfce_plan.rds")

fwe_run <- run_fwe_plan(fwe_plan, scheduler = "local")
```

The dry run lets you check exact files and commands before starting.
Local execution is synchronous. For batch processing, call
[`run_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_fwe_plan.md)
with `scheduler = "slurm"` or `"torque"` and set resources for the
correction jobs, such as `memgb_per_command` and `time_per_command`.
Collect after the jobs finish. In a new R session,
[`read_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/write_fwe_plan.md)
reloads the saved plan.

### Collect and interpret the results

``` r

fwe_result <- collect_fwe_results(fwe_run, require_complete = TRUE)

thresholds <- fwe_result_artifacts(fwe_result, artifact_role = "threshold_table")
stopifnot(nrow(thresholds) == 1L)
readr::read_csv(thresholds$artifact_file, show_col_types = FALSE)

surviving <- fwe_result_artifacts(
  fwe_result,
  artifact_role = "thresholded_statistic",
  fwe_alpha = 0.05
)
surviving[, c("target_id", "l1_cope_name", "l3_cope_name",
              "artifact_file", "can_use_as_mask")]
```

The threshold table records the enhanced statistic cutoff (`z_ptfce`),
requested FWE alpha (`p_value`), search volume in voxels, and number of
surviving voxels. A row with `significant_voxels = 0` is a valid result.
The thresholded image preserves the sign of surviving values and is zero
elsewhere. Use the [beta extraction
guide](https://hallquistlab.github.io/fmri.pipeline/articles/beta-extraction.html#use-surviving-voxels-after-ptfce)
to turn its support into an ROI mask and extract unthresholded COPEs.
The enhanced values in that image describe inference, while COPEs carry
the effect size.

Collecting writes `manifest/artifacts.csv` and `manifest/result.rds`
under the correction output directory. The manifest preserves target
labels, source files, and threshold settings. To review a previously
saved result, call
[`collect_fwe_results()`](https://hallquistlab.github.io/fmri.pipeline/reference/collect_fwe_results.md)
on its directory or result RDS.

If collection fails, run
`collect_fwe_results(fwe_run, require_complete = FALSE)` and inspect the
manifest’s `artifact_status`, `task_status`, and `log_file` columns. If
[`run_fwe_plan()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_fwe_plan.md)
stopped before returning a run object, collect from the plan’s output
directory when a result snapshot exists there. You can also use
`stop_on_error = FALSE` to retain task statuses when rerunning. Check
that all inputs are readable and that the worker can load pTFCE.
Verified completion requires matching input content and settings, a
successful worker receipt, and validated output files. Files from older
runs without a receipt are unverified. A rerun skips verified tasks by
default and recomputes stale or incomplete outputs. Use `force = TRUE`
to rerun verified tasks as well.

## Run pTFCE automatically after L3

To put correction in the original pipeline submission, add the following
`inference` block to the existing `group_mean` model in
`facehouse.yaml`. The [interactive model
builder](https://hallquistlab.github.io/fmri.pipeline/articles/interactive-model-builder.html#level-3-builder)
also offers group inference for compatible models.

``` yaml
    inference:
      method: ptfce
      targets:
        l1_cope_name: faces_minus_houses
      correction_mask: model
      fwe_alpha: [0.05, 0.01]
      two_sided: true
      compute:
        scheduler: inherit
```

Rebuild L3 from the edited specification. Alternatively, attach the spec
created above using `gpa$l3_models$models$group_mean$inference <- spec`.
[`export_glm_config()`](https://hallquistlab.github.io/fmri.pipeline/reference/export_glm_config.md)
preserves declarative inference settings.

``` r

gpa <- build_l3_models(gpa, from_spec_file = "facehouse.yaml")

run_glm_pipeline(
  gpa, l1_model_names = "facehouse",
  l2_model_names = "run_average", l3_model_names = "group_mean"
)
```

The pipeline collects the completed fit metadata and then runs a
separate `run_fwe` stage. Check the group fit and the correction status
in the [diagnosis
guide](https://hallquistlab.github.io/fmri.pipeline/articles/diagnosis.md).
Set `run_fwe = FALSE` in
[`run_glm_pipeline()`](https://hallquistlab.github.io/fmri.pipeline/reference/run_glm_pipeline.md)
to skip correction for one submission.

If only correction needs to run again, reload the completed batch’s
`gpa` and use the correction-only entry point. It requires the model’s
`inference` settings or a declarative specification configured with
[`build_fwe_correction()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_fwe_correction.md).

``` r

gpa <- run_fwe_correction(gpa, l3_model_names = "group_mean")

# For scheduler execution, wait until the submitted correction jobs finish.
gpa <- refresh_fwe_correction(gpa, require_complete = TRUE)
gpa$fwe_status
names(gpa$fwe_runs)
```

For a standalone FEAT directory without a pipeline configuration, supply
the completed `.gfeat` directly. Choose `zstat_numbers` after checking
the contrast names in FEAT.

``` r

ptfce <- ptfce_spec$new(
  gfeat_dir = group_map$analysis_dir,
  zstat_numbers = group_map$l3_cope_number,
  fwe_p = c(0.05, 0.01),
  two_sided = TRUE
)
ptfce$run()
```

The installed `ptfce_zstat.R` script accepts a FEAT group directory or
an individual z-statistic and its mask and smoothness inputs. Run it
with `--help` to see both forms.

## Apply cluster thresholds with AFNI

These standalone helpers run AFNI commands and do not require an
`fwe_spec`. For the ACF route, estimate the spatial autocorrelation from
first-level GLM residuals. Then simulate fields with that spatial
structure to derive cluster extent thresholds for the group search
volume. The assumptions and performance depend on the design and spatial
noise properties. [Cox et
al. (2017)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5399747/)

### Estimate ACF from the contributing runs

Start with the group map and mask selected in the [whole-brain
walkthrough](https://hallquistlab.github.io/fmri.pipeline/articles/whole-brain.html#find-the-group-map-and-its-inputs).
The example below assumes all participants and both of their runs
contributed. If any runs were excluded or a group fit used a subset,
first restrict `run_maps` to the actual contributors. Pair each run’s
residual time series with the mask from that same fit.

``` r

run_maps <- lookup_feat_outputs(
  gpa, level = 1L, what = "cope", refresh_status = TRUE
)$l1 |>
  dplyr::filter(
    l1_model == "facehouse", l1_cope_name == "faces_minus_houses",
    session == 1L, analysis_status == "complete", image_exists
  )

run_dirs <- unique(run_maps$image_feat_dir)
residuals <- file.path(run_dirs, "stats", "res4d.nii.gz")
run_masks <- file.path(run_dirs, "mask.nii.gz")
stopifnot(length(run_dirs) > 0L, all(file.exists(c(residuals, run_masks))))

acf_dir <- file.path(gpa$output_directory, "fwe", "facehouse_acf")
dir.create(acf_dir, recursive = TRUE, showWarnings = FALSE)

run_acf <- vapply(seq_along(residuals), function(i) {
  fwhmx <- afni_3dfwhmx$new(
    input_file = residuals[i], mask_file = run_masks[i]
  )
  fwhmx$run()
  stopifnot(fwhmx$is_complete())
  setNames(fwhmx$get_acf_params()[1:3], c("a", "b", "c"))
}, FUN.VALUE = numeric(3))

acf_summary <- data.frame(residuals, t(run_acf))
acf_summary
stopifnot(
  all(is.finite(run_acf)),
  all(run_acf["a", ] >= 0 & run_acf["a", ] <= 1),
  all(run_acf[c("b", "c"), ] > 0)
)

sim <- afni_3dclustsim$new(
  fwhmx_input_files = residuals,
  fwhmx_mask_files = run_masks,
  clustsim_mask = group_mask,
  out_dir = acf_dir,
  prefix = "facehouse_acf",
  pthr = 0.001, athr = 0.05,
  scheduler = "local", ncpus = 1L, seed = 42L
)
sim$submit()
sim$is_complete()
```

Inspect the per-run ACF estimates before pooling them. `3dFWHMx` saves
its outputs alongside each residual file. `afni_3dclustsim` reads those
estimates and takes the arithmetic mean of the three parameters across
runs. If the checks fail, review the affected residuals, masks, and
`3dFWHMx` output before pooling. If the fit did not retain residuals,
obtain them from the matching GLM before using this route. Wait until
`sim$is_complete()` returns `TRUE` before applying the threshold.

### Apply and record the matching threshold

This example uses voxelwise `p = 0.001`, clusterwise alpha `0.05`, and
NN1 (face adjacency). AFNI’s bi-sided setting finds positive and
negative clusters separately using a two-tail threshold. Use the same
search volume, sidedness, connectivity, and voxelwise p value as in the
selected simulation row. [AFNI 3dClustSim
documentation](https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dClustSim.html)

``` r

stopifnot(sim$is_complete())
sim$get_clustsim_df() |>
  dplyr::filter(nn == 1L, sided == "bi", pthr == 0.001, athr == 0.05)

cobj <- sim$apply_clustsim(
  statistic_nifti = group_map$image_file,
  NN = 1L, sided = "bi", pthr = 0.001, athr = 0.05,
  voxelwise_stat = list(stat_type = "z"),
  output_cluster_mask = TRUE,
  output_thresholded_image = TRUE,
  add_whereami = FALSE
)
cobj$get_clust_df()
cobj$get_outputs()
```

Report the cluster-level threshold, observed extent, and location for
surviving clusters. Cluster-extent inference does not establish a
voxelwise claim for each voxel in a cluster. The integer `cluster_map`
is the input for [extracting one effect per
cluster](https://hallquistlab.github.io/fmri.pipeline/articles/beta-extraction.html#use-a-cluster-map-after-afni-correction).

### Use residual sign randomization for a suitable group test

The `residuals_file` mode of `afni_3dclustsim` generates null statistic
maps with `3dttest++` and passes them to `3dClustSim`. Use it when
independent, sign-exchangeable group residuals represent the intended
one-sample test. This wrapper’s residual sign flips do not model
repeated measures, group-label exchange, or arbitrary covariate designs.
Use a null-generation method justified for your group design. Supply a
group-level 4D residual image with one volume per observation and its
matching 3D search mask.

``` r

# Paths from the one-sample analysis being corrected
residual_dir <- file.path(gpa$output_directory, "fwe", "one_sample_residuals")
dir.create(residual_dir, recursive = TRUE, showWarnings = FALSE)

sim <- afni_3dclustsim$new(
  residuals_file = "/study/one_sample/group_residuals.nii.gz",
  residuals_mask_file = "/study/one_sample/group_mask.nii.gz",
  residuals_njobs = 30L,
  out_dir = residual_dir,
  prefix = "one_sample",
  scheduler = "slurm"
)
sim$submit()
```

After completion, inspect `sim$get_clustsim_df()` and apply the matching
threshold to the statistic from that one-sample analysis. Keep the
simulation settings, mask, and statistic together when recording the
result.
