# Familywise error correction for voxelwise GLMs in fmri.pipeline

This vignette provides an overview of familywise error correction for
voxelwise GLM analyses in fmri.pipeline.

## Automatic correction at the end of L3

FSL L3 models can request pTFCE when they are specified, before any
images exist. The last step of the interactive
[`build_l3_models()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_l3_models.md)
dialog offers group inference. The equivalent YAML configuration is:

``` yaml
l3_models:
  group:
    model_formula: "~ age"
    execution_backend: fsl
    l3_input_mode: subject_rows
    contrasts:
      diagonal: true
    inference:
      method: ptfce
      fwe_alpha: [0.05, 0.01]
      two_sided: true
      correction_mask: model
      compute:
        scheduler: inherit
        memgb_per_command: 8
        time_per_command: "00:30:00"
```

`inference` applies to all contrasts of its model unless `targets`
selects a subset, for example `targets: {l3_cope_name: age}`. Correction
controls FWE separately within each selected map, including both tails
when `two_sided` is true. It does **not** jointly correct across
contrasts or models. Omit `inference` or set it to `false` to leave
correction disabled.

For programmatic configuration, assign the same specification to the
model:

``` r

gpa$l3_models$models$group$inference <- fwe_spec(
  "group_ptfce", targets = list(l3_model = "group"),
  fwe_alpha = c(0.05, 0.01), two_sided = TRUE
)
export_glm_config(gpa, "glm_config.yaml")
run_glm_pipeline(gpa, l3_model_names = "group")
```

The pipeline waits for L3 estimation, collects its output metadata, then
submits a separate `run_fwe` stage. That stage waits for correction
workers and collects their results before completing.
`run_glm_pipeline(..., run_fwe = FALSE)` skips correction for that
invocation. Only implemented, compatible methods are offered; this
release supports automatic pTFCE for FSL L3 outputs. Unsupported
backend/method combinations fail before submission. Automatic correction
uses the pipeline scheduler, or runs synchronously with
`compute: {scheduler: local}`. Use the correction-only entry point to
submit to a different scheduler.

The main YAML export/import preserves model-level inference and any
additional declarative specifications in `gpa$fwe_correction`.
[`build_fwe_correction()`](https://hallquistlab.github.io/fmri.pipeline/reference/build_fwe_correction.md)
can configure these additional specifications before estimation is
complete. Older live R6 correction objects remain usable through their
standalone APIs, but must be replaced with declarative specifications
for automatic execution.

To inspect or retry correction after estimation, reload the saved
pipeline object and use the correction-only entry point:

``` r

gpa <- readRDS(gpa$output_locations$object_cache)
gpa <- run_fwe_correction(gpa, l3_model_names = "group")
# After scheduler jobs finish:
gpa <- refresh_fwe_correction(gpa, require_complete = TRUE)
gpa$fwe_status
gpa$fwe_runs$group_ptfce$result
```

This does not refit GLMs. Local execution completes synchronously;
scheduler execution returns after submission. Correction status is
stored separately from estimation status. A valid result with zero
significant voxels is still complete; `significant_voxels` in the
threshold table records that outcome.

Completion requires a successful worker receipt, matching content hashes
of inputs and inference settings, and intact validated outputs. Old
results without receipts are unverified and will be recomputed. Changing
a mask, statistic image, alpha, or sidedness invalidates completion even
when filenames are unchanged. A forced or failed rerun cannot reuse an
earlier success receipt. Use
[`collect_fwe_results()`](https://hallquistlab.github.io/fmri.pipeline/reference/collect_fwe_results.md)
to refresh a saved result snapshot before reuse.

For a custom ROI, supply a nonempty binary mask aligned with and
contained within the model mask. pTFCE reuses the model’s stationary
smoothness estimate but recalculates search volume and resel count for
the correction region. The worker rejects incompatible geometry and
non-finite values inside the mask. Install `pTFCE` and `fmri.pipeline`
in the R environment used by compute jobs; workers do not install
packages automatically. `min_cluster_voxels` is deprecated because it
never filtered inference outputs; cluster-size filtering belongs to
separate reporting with `ptfce_spec$get_clusters(clust_nvox = ...)`.

## Standalone declarative pTFCE planning and execution

For completed level-3 FSL/FEAT analyses, the preferred workflow is to
describe the intended correction once, compile that description into an
inspectable plan, and execute the plan locally or through the configured
scheduler. The plan resolves model and contrast names through the
pipeline output inventory, rather than requiring you to construct GFEAT
paths or call a correction method for each map yourself.

At present, pTFCE is the FWE method with a plan-and-execution adapter. A
pTFCE plan needs, for every selected level-3 statistic map, the `zstat`
image, the correction mask, and the FEAT `stats/smoothness` file. By
default, the analysis-specific FEAT mask is used for correction. The R
interpreter used on the execution host must also be able to load the
`pTFCE` package and its image I/O dependencies.

### 1. Describe the correction

Create an `fwe_spec` from semantic selectors. For example, the following
specification applies pTFCE to every level-3 z-statistic map from
`l3_group`. Narrow the scope by adding selectors such as `l3_cope_name`,
`l2_model`, or `l1_model`; selector values may be vectors.

``` r

all_l3_spec <- fwe_spec(
  name = "l3_group_ptfce",
  targets = list(l3_model = "l3_group"),
  method = "ptfce",
  fwe_alpha = c(0.05, 0.01),
  two_sided = TRUE,
  write_threshold_images = TRUE,
  correction_mask = "model",
  compute = list(scheduler = "slurm")
)

all_l3_spec
```

`fwe_alpha` controls the FWE thresholds reported and, when
`write_threshold_images = TRUE`, the thresholded pTFCE images written
for each selected map. With `two_sided = TRUE`, pTFCE corrects positive
and negative tails separately and uses half of each requested FWE alpha
for the per-tail threshold. Set `two_sided = FALSE` when a one-sided
positive-tail analysis is the intended inferential test.

Use `correction_mask = "model"` for the FEAT mask associated with each
target, or supply one explicit NIfTI mask when the same small-volume
mask should be used for every selected map. Specifications are ordinary
serializable objects, so they can be kept under version control as YAML.

``` r

write_fwe_spec(all_l3_spec, "config/fwe/l3_group_ptfce.yaml")

# In a later R session:
all_l3_spec <- read_fwe_spec("config/fwe/l3_group_ptfce.yaml")
```

### 2. Compile and inspect the plan

Compile the specification against a `glm_pipeline_arguments` object
after the selected level-3 FEAT analyses have completed. Planning is
side-effect free: it does not create output directories, submit jobs, or
run pTFCE.

``` r

fwe_plan <- plan_fwe_correction(gpa, all_l3_spec)
fwe_plan

# One ready task per selected level-3 z-statistic map.
fwe_plan$tasks[, c(
  "task_id", "l3_cope_name", "zstat_file",
  "correction_mask_file", "smoothness_file", "status"
)]

# Expected pTFCE image, threshold table, and thresholded maps.
fwe_plan$artifacts[, c("task_id", "artifact_role", "fwe_alpha", "artifact_file")]
```

If you are planning before all inputs exist, use a non-strict preview.
Blocked tasks identify the missing input instead of running with an
incomplete setup.

``` r

fwe_plan <- plan_fwe_correction(
  gpa, all_l3_spec, require_inputs = FALSE
)
fwe_plan$tasks[, c("task_id", "l3_cope_name", "status", "missing_inputs")]
```

Once the plan is ready, it may also be saved and executed from another R
session or submission script.

``` r

write_fwe_plan(fwe_plan, "derivatives/fwe/l3_group_ptfce/plan.rds")
```

### 3. Preview, execute, and collect results

First use a dry run to inspect the generated command for each task
without creating files or submitting work.

``` r

dry_run <- run_fwe_plan(fwe_plan, dry_run = TRUE)
dry_run$execution[, c("task_id", "execution_status", "command")]
```

For a small analysis, run pTFCE synchronously on the current host.
Otherwise, omit `scheduler` to use the scheduler compiled into the plan,
or set it explicitly to `"slurm"` or `"torque"` and provide
task-appropriate resources.

``` r

# Local execution
fwe_run <- run_fwe_plan(fwe_plan, scheduler = "local")

# Example scheduler submission (use this instead of the local call above).
# fwe_run <- run_fwe_plan(
#   fwe_plan,
#   scheduler = "slurm",
#   commands_per_cpu = 1L,
#   cpus_per_job = 1L,
#   memgb_per_command = 8,
#   time_per_command = "00:30"
# )
```

Local execution returns after pTFCE has finished. Scheduler execution
returns after submission, so refresh and collect after those jobs have
completed. Real runs automatically persist a result snapshot and
`manifest/artifacts.csv` in the correction output directory; collecting
again updates those canonical files with the current artifact state.

``` r

fwe_plan <- refresh_fwe_plan(fwe_run$plan)
fwe_result <- collect_fwe_results(
  fwe_plan,
  require_complete = TRUE
)

# Select thresholded pTFCE maps by their semantic contrast and FWE alpha.
fwe_masks <- fwe_result_artifacts(
  fwe_result,
  artifact_role = "thresholded_statistic",
  l3_cope_name = "age",
  fwe_alpha = 0.05
)
fwe_masks[, c("target_id", "artifact_file", "can_use_as_mask")]
```

The artifact manifest retains the target, method, threshold, and
source-map metadata alongside each output path. It is therefore the
stable handoff point for downstream analysis: select the desired
`artifact_file` rather than reconstructing a path from a GFEAT
directory. The legacy sections below remain available for manually
configuring and running the earlier correction tools.

## AFNI 3dClustSim

As outlined in Cox et al. (2017, Brain Connectivity), the conventional
Gaussian random field theory approach to cluster correction in fMRI
analysis is anticonservative in many cases because the empirical
autocorrelation function of noise in fMRI has a much longer tail than
the Gaussian distribution would predict. This was highlighted in Woo &
Wager (2014) and elaborated by Eklund (2016, to greater effect).

There are a couple of improvements proposed and tested by Cox to control
the familywise error rate to (nearly) .05. The package supports both of
these methods. Let’s do a quick walkthrough.

### Method 1: AFNI 3dClustSim + 3dFWHMx -ACF

The basic approach to 3dClustSim (using -fwhm or -fwhmxyz) assumes a
Gaussian-distributed noise structure. As noted above, this is not
realistic in fMRI data. The first correction attempts to overcome the
assumption of Gaussianity by fitting a long-tailed distribution to the
noise structure of empirical data. Specifically, one passes in the
first-level GLM residuals (one volume per timepoint) as inputs to
3dFWHMx and includes the -ACF flag. This flag will estimate the
autocorrelation function from the data by fitting parameters of a
mono-exponential + Gaussian mixture that can better accommodate the
long-tailed spatial autocorrelations in fMRI data.

3dFWHMx produces ACF parameters for every volume in each run’s residuals
file. The usual approach (per AFNI) is to take the average of these
parameters per run as the best estimate of the spatial structure of
noise in that run. Then, to estimate cluster correction thresholds in
3dClustSim, one further averages these run-averaged ACF parameters
across all runs to achieve an overall ACF for the entire dataset. In Cox
et al., 2017, the median ACF parameters across datasets are input to
3dClustSim – I suppose the mean could work, too, depending on whether
there are outliers that would sway it. The estimated values of this
averaged 3-parameter ACF are passed to 3dClustSim with the -acf
argument.

Based on the -acf input, 3dClustSim now estimates null datasets that
have the spatial distribution of noise derived from the empirical data,
which better controls for FWE relative to the standard Gaussian
assumption/distribution. Cox and colleagues show that at reasonable
voxelwise thresholds (p \< .005 or p \< .001), the -ACF approach is
close to the nominal .05 alpha level (4 - 8%) when a) the spatial
smoothing kernel is \> 4mm FWHM, and b) the test/contrast comes from an
event-related design with many repeats of relatively brief stimuli. The
false positive rate (FPR) for block designs and small smoothing kernels
(4mm FWHM) are worse, sometimes \> 15%. Also, note that a voxelwise
threshold of .01 or more (e.g., .02) is far too generous for cluster
correction and should be avoided in general (see Woo and Wager for more
info). The above summary is depicted in Cox et al., 2017, Fig 1 G, H, I.

#### Implementing the -ACF correction in fmri.pipeline

##### Example 1: Setting up a 3dFWHMx object yourself

If you would like to run 3dFWHMx yourself using the afni_3dfwhmx wrapper
class, here is an example:

``` r

# example FEAT LVL1 output
path <- "/proj/mnhallqlab/studies/MMClock/MR_Proc/10637_20140304/mni_5mm_aroma/sceptic_vchosen_ventropy_dauc_pemax_vtime_preconvolve/FEAT_LVL1_run1.feat"
res4d_file <- file.path(path, "stats", "res4d.nii.gz")
fwhmx_mask_file <- file.path(path, "mask.nii.gz")

test <- afni_3dfwhmx$new(input_file = res4d_file, mask_file = fwhmx_mask_file, average = "geometric")
test$run() # will run 3dFWHMx in the local compute environment (use force = TRUE to re-run completed file)
test$is_complete() # TRUE if 3dFWHMx has already completed and has generated its output files
test$get_acf_params() # average ACF params (typically combined with other datasets, then averaged, for 3dClustSim)
test$get_fwhm_by_volume() %>% head() # FWHM for each sub-brik (volume)
test$get_acf_by_radius() %>% head() # FWHM for each sub-brik (volume)

# commented out to avoid deleting files of interest
#test$delete_outputs(prompt = TRUE) # delete files generated by this input -- it will ask you to confirm each
```

This shows how you can run 3dFWHMx through R and read the outputs of the
ACF fitting and FWHM calculation.

##### Example 2: Running 3dFWHMx for a set of residual files

The afni_3dfwhmx class is a wrapper around 3dFWHMx. In the typical ACF
approach, we pass the average ACF parameters across all datasets to
3dClustSim -acf. That is, we must provide only one set of ACF parameter
estimates to 3dClustSim -acf so that the program can generate null
datasets with this spatial autocorrelation structure. We usually use
3dClustSim to calculate an FWE correction on a group analysis that may
reflect hundreds of first-level datasets (fMRI runs). The residuals from
each first-level GLM must be passed through 3dFWHMx so that the ACF
parameters can inform the overall (group) ACF estimates.

Said differently, the typical 3dClustSim -acf approach requires ACF
parameter averaging at two stages. First, we take the average ACF
parameters across all volumes from the first-level residuals. This is
controlled by the -geom or -arith flag of 3dFWHMx, which specifies how
volume-wise ACF estimates should be combined to yield a run-wise ACF.
Then, we average the run-wise ACF parameters to derive a sample-wise ACF
average that reflects the average spatial autocorrelation across all
datasets in the analysis.

To aid in the estimation of ACF parameters across many runs, we provide
the afni_3dfwhmx_list class. Conceptually, this is a list of
afni_3dfwhmx objects, each of which pertains to a single run of data.
The goal of the afni_3fwhmx_list class is to run 3dFWHMx for all input
datasets and to provide average values for the ACF across the datasets.
This average can then be passed directly to 3dClustSim -acf.

Here’s a quick example:

``` r

# get ACF parameters for 8 first-level residual files for testing
# this example is simpler than typical applications where we have many subject and many runs
path <- "/proj/mnhallqlab/studies/MMClock/MR_Proc/10637_20140304/mni_5mm_aroma/sceptic_vchosen_ventropy_dauc_pemax_vtime_preconvolve"
res4d_files <- list.files(pattern = "res4d.nii.gz", path = path, full.names = T, recursive = T)
fwhmx_mask_files <- list.files(pattern = "mask.nii.gz", path = path, full.names = T, recursive = T)

# these are the files to be passed through 3dFWHMx
print(res4d_files)

# create the list object requires you to provide input and mask files for each residuals dataset of interest.
test <- afni_3dfwhmx_list$new(input_files = res4d_files, mask_files = fwhmx_mask_files, scheduler = "slurm")

# submit 3dFWHMx jobs to cluster -- the object will divide the list into smaller jobs if needed
# if 3dFWHMx has completed for all inputs already, the $submit method will just return immediately and note
# that there is nothing to submit. If you want to force re-estimation, use $submit(force = TRUE).
test$submit()

# has 3dFWHMx completed for all datasets provided to the list?
test$is_complete()

# average a, b, c parameters across datasets -- good for input to 3dClustSim -acf
test$get_acf_average()

# average effective FWHM estimate using -ACF
test$get_effective_fwhm()

# per-dataset ACF summary if you want to see what got averaged
test$get_acf_df()
```

##### Example 3: Setup 3dClustSim + 3dFWHMx -acf approach

Although you are welcome to use the afni_3dfwhmx and afni_3dfwhmx_list
directly, if your purpose is to run 3dClustSim using the ACF approach,
you can ask the pipeline to run and compile the 3dFWHMx results for all
level 1 residuals in a single step.

The afni_3dclustsim object will handle the 3dFWHMx batch estimation as
an initial compute job, then pass the averaged ACF parameters to
3dClustSim as a dependent compute job (i.e., 3dFWHMx has to complete
before 3dClustSim starts).

This is achieved using the `fwhmx_input_files` and `fwhmx_mask_files`
arguments to `afni_3dclustsim$new()`. Here, let’s pass the same level 1
residuals and mask files as above (in the 3dFWHMx list step). Here, we
ask 3dClustSim to use 8 cores for parallel estimation and we use the MNI
2009c brain mask (at 2.3mm) as the volume over which to correct for FWE
(this should match the empirical data).

``` r

clustsim_acf <- afni_3dclustsim$new(
  fwhmx_input_files = res4d_files, fwhmx_mask_files = fwhmx_mask_files,
  scheduler = "slurm", out_dir = "/proj/mnhallqlab/users/michael/fmri.pipeline", # where you want the 3dClustSim output files to go
  clustsim_mask = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", ncpus = 8
)

# Submit the job to the scheduler. If the 3dFWHMx calls have not been completed yet for all fwhmx_input_files, this will occur first via
# a separate job submission
clustsim_acf$submit()

# check whether it's complete
clustsim_acf$is_complete()
```

Once the above job completes, we can look at the number of voxels needed
at different combinations of voxelwise and clusterwise thresholds.

``` r

clustsim_acf$get_clustsim_df() %>% head()
```

If you were curious about a cluster correction result, you can use the
[`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)
function to obtain this. For example, imagine that I want to look at the
whole-brain FWE cluster threshold for all clusters defined by NN=1
(faces touching) and a voxelwise significance level of p = .001. Here’s
the result.

``` r

clustsim_acf$get_clustsim_df() %>% dplyr::filter(nn==1 & pthr==.001)
```

If you would like to apply the 3dClustSim correction to your data to
obtain an FWE-corrected map and corresponding clusters, use the
`$apply_clustsim()` method. Note that this method requires a test
statistic file as its input and you should only give voxelwise statistic
maps that relate to the 3dClustSim that was run. That is, if you run
3dFWHMx on 50 L1 residual files, you should only apply the resulting
3dClustSim results to group analyses that pertain to those 50 files.

For now, I violate this here in the interest of computational efficiency
(don’t want to run 3dFWHMx on a massive set of files for a demo). Still,
these are z-statistic maps from a larger analysis of which these
subjects were a part.

By default, this method clusterizes using a cluster threshold of p = .05
and a voxelwise threshold of p = .001, bisided, NN = 1 (faces touch).

Also, by default, this method has `add_whereami = TRUE`, which looks up
anatomical labels for each cluster that is found.

``` r

cobj <- clustsim_acf$apply_clustsim(
  statistic_nifti = file.path(
    "/proj/mnhallqlab/studies/MMClock/group_analyses/MMClock_aroma_preconvolve_fse_groupfixed",
    "sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept-Age.gfeat/cope4.feat/stats/zstat1.nii.gz"
  )
)

# will include anatomical labels, too, by default
cobj$get_clust_df()
```

### Method 2: AFNI 3dClustSim with 3dttest++ permutations

The second method proposed by Cox involves generating null datasets from
the empirical distribution of residuals in the group analysis. That is,
if we run our group GLM, we can compute a $`Y - \hat{Y}`$ value for
every voxel and every subject. This results in a voxelwise map per
subject that contains the residual from the group analysis.

Rather than attempt to approximate the spatial autocorrelation function
using a given semiparametric distribution, the permutation-based
approach tries to get test statistics on a group analysis in which the
residuals have been shuffled (permuted and/or sign-flipped) to create a
set of null maps that have a spatial correlation structure that is
directly based on the data. These datasets are then thresholded and
clusterized by 3dClustSim directly, rather tha having 3dClustSim
simulate data with a given spatial noise structure.

Here is the direct text from Cox et al., 2017

    1. Compute the residuals of the model at each voxel at the group level.
    2. Generate a null distribution by randomizing among subjects the signs of the residuals in the test (and permuting subject data sets between groups, if doing a two-sample test without covariates), 
      repeat the t-tests (with covariates, if present), and iterate 10,000 times.
    3. Take the 10,000 3D t-statistic maps from the randomization and use those as input to 3dClustSim (with no additional smoothing or random deviate generation): 
      threshold the maps at a large number of p values (e.g., 0.01, 0.005, and 0.001), clusterize them, then record the maximum cluster size for each p value from each t-map, 
      accumulate statistics across the 10,000 t-maps, and then determine the cluster-size threshold to achieve a given FPR, for each p value.

At present, we adopt the potentially flawed, but likely awfully close
and usually correct, approach of just doing sign-flipping on the
residuals – the approach advocated for one-sample t-tests. That is, the
signs of residuals are randomized across subjects, but we do not permute
group labels, nor do we permute the covariates.

``` r

clustsim_perm <- afni_3dclustsim$new(
  residuals_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-pe/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_pe.gfeat/cope1.feat/stats/res4d.nii.gz",
  residuals_mask_file = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-pe/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_pe.gfeat/cope1.feat/mask.nii.gz",
  residuals_njobs = 32,  scheduler = "slurm", prefix = "res4d", iter = 30000,
  clustsim_mask = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", ncpus = 8
)
```

In this approach, the number of null datasets generated is identical to
the number of iterations of 3dClustSim operations to estimate the
cluster thresholds for FWE correction. In the permutation/sign-flipping
approach, we can split up the total number of null datasets (here:
30000) across a number of separate HPC scheduler jobs. This is
controlled by `residuals_njobs` and it’s reasonable to park 1000-3000
iterations/datasets in each job (i.e., `iter / residuals_njobs ~ 10`).

You can submit the 3dttest++ permutation + 3dClustSim jobs using the
`$submit()` method. Note that if the permutations have already completed
(and are saved), this step will not be repeated. Likewise, if 3dClustSim
is already complete, the `$submit()` method will tell you this and not
re-run the simulation unless you use `$submit(force = TRUE)`.

``` r

clustsim_perm$submit()

clustsim_perm$is_complete()
```

Now that we have the 3dClustSim results, we can adopt the same approach
as with the -ACF object above.

``` r

clustsim_perm$get_clustsim_df() %>% dplyr::filter(nn==1 & pthr==.001)
```

And if we want to clusterize and threshold a statistic file that relates
to this analysis, we can use `$apply_clustsim()` to any statistics from
the GLM whose residuals we used for permutation.

What are the contrasts in this analysis? If you don’t recall from
setting this up manually (or otherwise), we can use the handy
`read_feat_dir` function in fmri.pipeline.

``` r

feat_info <- read_feat_dir("/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-pe/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_pe.gfeat/cope1.feat")

feat_info$cope_df %>% select(cope_number, contrast_name, z)
```

Let’s take a look at FWE clusters from zstat6, the overall contrast.

``` r

cobj <- clustsim_perm$apply_clustsim(
  statistic_nifti = "/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-pe/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_pe.gfeat/cope1.feat/stats/zstat6.nii.gz"
)
```

We can get info about the FWE clusters from this map, based on the
3dClustSim permutation-based thresholds (again, default cluster p \< .05
and voxel p \< .001, bisided).

``` r

cobj$get_clust_df()
```

If we want to look at the file outputs of the clusterize operation, you
can obtain them like this:

``` r

cobj$get_outputs()
```

### AFNI 3dClustSim: how to decide what to do

The nice part about the ACF approach is that the 3dClustSim results can
be used for any group analysis that included the lower-level runs whose
ACFs were estimated. This means you can get a lot of mileage from the
3dClustSim results across essentially all group analyses in a dataset.

The bad part about the ACF approach is that it may be slightly
anticonservative, esp. if you use liberal voxelwise thresholds (e.g. p =
.005) or block design effects.

The nice part about the 3dttest++ permutation approach is that it
tightly controls FWE to the nominal level. The bad parts are: how to
handle more complex designs (not just one-sample t-tests) and that in
principle, you need to re-run 3dClustSim for each zstat – or minimally,
FEAT GLM. This is true because each group analysis has a single
residuals while, which is used for FWE.

In practice, within a set of GLM analyses of a dataset, there is often
little difference if you switch from ACF to permutation, or switch from
one permutation 3dClustSim result to another (e.g., a different GLM in
the same dataset).

## Extracting betas (regression coefficients) from FWE-corrected maps

Okay, great, now I have my FWE clusters for a contrast of interest. How
do I get betas (aka regression coefficients) from my whole-brain
analyses? This will get easier soon in the pipeline, but for now, we can
pass the cluster mask from the clusterize operation to the
`extract_glm_betas_in_mask` function in the package. This will ask you
from which models you wish to extract coefficients and will then extract
coefficients from each cluster in the FWE mask for all models you
specify/request in the prompts. Note that this means that you can
extract betas from one model based on FWE clusters drawn from another
model (that was passed through FWE correction). This is handy for
checking representational specificity, testing dissociation of different
effects, and so on.

To get the outputs generated by clusterize step of `$apply_clustsim()`,
you can use the `$get_outputs()` method of the returned clusterize
object.

``` r

cobj$get_outputs()
```

Here, the `cluster_table` is the text file of cluster locations and
sizes from 3dClusterize, `cluster_map` is the integer-valued mask of all
clusters that survived the corrections requested, and
`cluster_masked_data` is the map of interest (here zstat1) that has been
thresholded at the settings requested (essentially anything in that map
is FWE-significant and all non-significant voxels are removed).

For extracting betas, we’d want to hand off the `cluster_map` as the
input to `extract_glm_betas_in_mask` like so:

``` r

extract_glm_betas_in_mask(gpa, mask_files = cobj$get_outputs()$cluster_map)
```

### Subclustering

Oftentimes, if you have a robust map, a whole-brain FWE correction will
produce some remarkably large clusters (e.g., 10,000 voxels). Although
it is technically correct that this entire mass is statistically
significant, it can be deeply dissatisfying as a neuroscientist to see a
cluster that spans many brains regions that are functionally and/or
anatomically distinct. These ‘mega clusters’ likely reflect the
co-activation of potentially distinct neural circuits/systems. If we
extract betas from a mega cluster, however, we will be left with just
one number (usually the mean) for each subject that represents neural
activity across the wide array of regions that span the cluster.

It is possible that the voxelwise activation pattern is indeed
homogeneous within a mega cluster and that the mean is representative.
But I think this is a strong and often untested assumption. Moreover,
inferences about mega clusters are often awkward, leading to a sentence
like, “Individuals who had higher scores on the working memory test also
showed stronger reward prediction activity across the parietal and
occipital cortices.” Perhaps it is true that the between-subjects
association of prediction errors with working memory spans the expanse
of parietal and occipital cortices, yet it is even more plausible that
the strength of that association may vary by region or subregion in ways
that our mean beta coefficient cannot uncover.

When a mega cluster emerges, an intuitive approach is to increase the
statistical threshold relative to the FWE p \< .05 level in order to see
if the mega cluster breaks up into smaller clusters that modulated even
more strongly (e.g,. at FWE p \< .005). If this happens (breaking into
smaller clusters as one raises the threshold), it can lead to greater
anatomical specificity while still preserving Type I error control
(indeed, those smaller clusters survive a more stringent correction).

Building on this intuitive approach, I have implemented a simple
subclustering algorithm that can be applied to any cluster set
identified by 3dClusterize. After setting up an afni_3dclusterize object
and running it to obtain the cluster data.frame, you run subclustering
as follows:

``` r

cobj$generate_subclusters(break_nvox = 800, min_subclust_nvox = 100, max_subclust_nvox = 4500, max_iter = 100, print_progress = TRUE)
```

#### Subclustering algorithm

The subclustering algorithm begins by searching for any clusters that
exceed the user-specified cluster size (in voxels). This is controlled
by the `break_nvox` argument and defaults to `800`. If any cluster is
larger than this threshold, subclustering will be attempted. For each
large cluster, the algorithm begins with the threshold applied on the
overall map in the initial 3dClusterize setup. It then iteratively
increases the threshold by a step size (default .1) and reclusters the
data *within the larger cluster*. That is, the algorithm attempts to
break up each large cluster into 2+ subclusters. The number and size of
the subclusters are controlled by the following arguments:

- `min_subclust_nvox`: The minimum number of voxels that must be present
  in every subcluster. Default: 25.
- `max_subclust_nvox`: The maximum number of voxels allowed in each
  subcluster. Default: unlimited.
- `min_n_subclust`: The lowest number of subclusters that will be
  accepted. Default: 2 (anything less wouldn’t really be subclustering!)
- `max_n_subclust`: The highest number of subclusters that will be
  accepted. Default: unlimited.

At each step, the algorithm checks whether the subclusters identified at
the current thresholds satisfy the size and number criteria provided by
the user. If the solution meets these criteria, the algorithm then
attempts to refine the solution by walking backwards in its thresholding
in smaller steps. For example, if a z-statistic threshold of 4.2 did not
meet the subclustering criteria, whereas a threshold of 4.3 did (one
step in according to a .1 step size), the algorithm will try values in
between, such as 4.28. The lowest threshold that still satisfies the
subclustering criteria will be retained, thereby maximizing the size of
each subcluster.

#### Examining subclusters

Once you have run subclustering, you can see the subclusters (if any)
using the `$get_clust_df()` method, as before. Note that you can
suppress subclusters in this data.frame with
`$get_clust_df(include_subclusters = FALSE)`. In the cluster data.frame,
you’ll see a column called `subroi_num`, which is the subcluster number
within each corresponding `roi_num`. By default, the subcluster maps and
3dClusterize tables are also saved in the same folder as the
input/threshold image, with a suffix like “\_roi1_subclusters.1D”. You
can view the locations of these files with
`$get_outputs()$subcluster_files`.

``` r

cobj$get_clust_df()
cobj$get_outputs()$subcluster_files
```

### Extracting clusters that overlap with a parcellation or atlas

Although subclustering can be quite useful for finding smaller
functional ROIs that may have dissociable statistics, you also tend to
lose many FWE-significant voxels from the thresholding process. One
solution to this is to extract voxelwise statistics from an a priori
atlas, as noted above. This is achieved in the pipeline using
[`extract_glm_betas_in_mask()`](https://hallquistlab.github.io/fmri.pipeline/reference/extract_glm_betas_in_mask.md).

A related approach is to extract coefficients from the atlas only for
parcels that are in the FWE-significant map. That is, we can check
whether each parcel in the atlas overlaps with the clusterized map. If
it does, then we could extract subject- or run-level statistics from
that parcel. For parcels that do not overlaps FWE-significant clusters,
however, we typically have little reason to interrogate them further.

You can generate an atlas subset that overlaps with your cluster map
using the `subset_atlas_against_clusters()` method. This method accepts
an integer-valued brain parcellation map of the same resolution as the
input/threshold map. For each parcel (unique integer value) in the atlas
map, it checks whether the parcel overlaps substantially with the
clusters in the thresholded statistical map. If there is substantial
overlap, the atlas parcel is retained. If not, the parcel is removed.
The amount of overlap is controlled by the `minimum_overlap` argument,
which defaults to 0.8 (i.e., 80%; values can be between .01 and 1.0).

After this process, the atlas subset file is saved in the same folder as
the 3dclusterize input/threshold image, named according to the
combination of the threshold file and the atlas file. You can specify a
unique name with the `output_atlas` argument to the method. Or you can
suppress creation of an atlas overlap file with `output_atlas = FALSE`.
This will just print the retained and rejected parcels, but not generate
any new files.

By default, for retained atlas parcels (i.e., those where there is
substantial overlap), all voxels within the parcel are retained. This is
consistent with the view of an atlas as a relatively independent basis
for functional region definition. You can, however, ask that only voxels
within the retained parcels that overlap with the thresholded and
clusterized map by retained. This is achieved with
`mask_by_overlap = TRUE`, perhaps at the expense of somewhat
circular/data-dependent parcel definition.

``` r

atlas <- "/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_DorsAttn_2.3mm.nii.gz"
cobj$subset_atlas_against_clusters(atlas, minimum_overlap = 0.7, output_atlas = "default")
cobj$get_outputs()
```

You can now use the atlas subset files of interest as an input mask to
[`extract_glm_betas_in_mask()`](https://hallquistlab.github.io/fmri.pipeline/reference/extract_glm_betas_in_mask.md).
