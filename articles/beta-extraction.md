# Extracting contrast estimates from masks and atlases

After a whole-brain fit, you might want one value per participant to
plot the magnitude and direction of an effect within a region. This
guide uses the completed FSL analysis from [whole-brain group
maps](https://hallquistlab.github.io/fmri.pipeline/articles/whole-brain.md).
It covers an independent atlas, clusters from AFNI, and voxel sets from
pTFCE. For the latter two, first run the corresponding method in the
[FWE correction
guide](https://hallquistlab.github.io/fmri.pipeline/articles/fwe.md).

[`extract_glm_betas_in_mask()`](https://hallquistlab.github.io/fmri.pipeline/reference/extract_glm_betas_in_mask.md)
extracts FSL FEAT images. A **COPE** is the estimated contrast of
regression coefficients. In the example, it is the face coefficient
minus the house coefficient. `zstat` and `tstat` quantify evidence
against the null, while `varcope` is the estimated within-image variance
of the contrast at each voxel. The mean of a z-map is a descriptive
summary. Likewise, the mean of the voxelwise VARCOPEs is generally
different from the sampling variance of a regional mean.

For AFNI and SPM fits, use the [image-based extraction
example](#extract-directly-from-afni-or-spm-images) with the fitted
coefficient or contrast images.

## Choose the maps to summarize

| Question | FEAT level to extract | Result |
|----|----|----|
| How did the estimated contrast differ across runs? | L1 | A value per participant, visit, run, contrast, and ROI |
| What is each participant’s combined effect? | L2 in the two-run example | A value per participant, visit, contrast, and ROI |
| What is the fitted sample effect in the ROI? | L3 | A value per group map and ROI |

In a single-run FSL analysis, use L1 for the participant effects. L3
contains group effects, so there is no subject ID on those rows. The L3
session identifies the fitted session or is `0` for a pooled fit.

Work with the completed batch cache as described in [run and reload the
analysis](https://hallquistlab.github.io/fmri.pipeline/articles/whole-brain.html#run-and-reload-the-analysis).
The extraction function reads the completed FEAT setup tables. Refresh
the levels you plan to use if your snapshot predates estimation.

``` r

library(fmri.pipeline)
library(dplyr)

gpa <- refresh_glm_status(gpa, level = 2L, glm_software = "fsl")
feat <- lookup_feat_outputs(gpa, level = 2L, what = "cope", source = "setup")
subject_maps <- feat$l2 |>
  filter(
    l1_model == "facehouse", l1_cope_name == "faces_minus_houses",
    l2_model == "run_average", l2_cope_name == "EV_Intercept", session == 1L
  )

subject_maps |>
  count(analysis_status, is_passthrough, image_exists)
```

Compare this inventory with the participant roster and the inputs of the
specific group fit. Missing files are reported and omitted during
extraction. For a one-run L2 passthrough, the only available logical L2
image is the L1 COPE. Inspect the participant and session counts before
analyzing the table.

## Prepare the mask

A mask is a 3D NIfTI already aligned to **every** statistic image being
summarized. The extractor checks dimensions, voxel sizes, spatial units,
and voxel-to-world transforms. Verify both the template space and the
image grid. If a mask comes from a different grid in the same space,
resample it onto a representative COPE grid with nearest-neighbor
interpolation for discrete labels. If the spaces differ, apply the
appropriate spatial transform first. Overlay the result on
representative COPEs and inspect coverage.

- A binary mask uses `0` for background and `1` for one ROI.
- An integer atlas or cluster mask uses `0` for background and a
  distinct integer per parcel or cluster. Keep a lookup table from label
  to region name.
- A continuous mask returns voxel rows with the continuous values in
  `mask_value`. It does not automatically calculate weighted means.
- With a 4D input, the function uses only the first volume. Supply
  explicit 3D masks when using a multi-volume atlas.

For atlas extraction, start with an atlas in the appropriate template
space. Use distinct mask basenames if extracting several atlases into
one output directory. The basename is used in the table and output
filename.

### Use a cluster map after AFNI correction

The [AFNI correction
example](https://hallquistlab.github.io/fmri.pipeline/articles/fwe.html#apply-cluster-thresholds-with-afni)
creates `cobj`, an `afni_3dclusterize` object. Each surviving cluster
has a different integer label in `cluster_map`. Use `cluster_table` to
identify the cluster IDs and their locations.

``` r

cluster_files <- cobj$get_outputs()
cobj$get_clust_df()
mask_file <- cluster_files$cluster_map
```

If no clusters survive, there may be no usable label map. Record the
null result and skip cluster-based extraction. For descriptive
exploration of large surviving clusters, `generate_subclusters()` raises
the statistic threshold within the parent cluster and searches for
smaller connected components. Alternatively,
`subset_atlas_against_clusters()` can keep parcels that overlap the
cluster mask. In the example below, each retained parcel must have at
least 70% of its voxels inside the cluster mask. The default
`mask_by_overlap = FALSE` preserves the full retained parcels.

``` r

cobj$generate_subclusters(break_nvox = 800, min_subclust_nvox = 100)
cobj$get_outputs()$subcluster_files

cobj$subset_atlas_against_clusters(
  atlas_file = "/study/atlases/aligned_parcels.nii.gz",
  minimum_overlap = 0.70, mask_by_overlap = FALSE
)
cobj$get_outputs()$atlas_files
```

Report any data-driven refinement along with the initial selection rule.
These summaries inherit the selection of the parent corrected map.

### Use surviving voxels after pTFCE

The [pTFCE result
manifest](https://hallquistlab.github.io/fmri.pipeline/articles/fwe.html#collect-and-interpret-the-results)
gives you a thresholded, signed statistic image. Convert it to integer
labels to summarize by direction. The example makes one ROI for all
surviving positive voxels and another for all surviving negative voxels,
including disconnected components. The output image inherits geometry
from the statistic image. The underlying values for extraction still
come from the unthresholded COPEs.

``` r

ptfce_image <- fwe_result_artifacts(
  fwe_result,
  artifact_role = "thresholded_statistic",
  l3_cope_name = "EV_Intercept", fwe_alpha = 0.05
)
stopifnot(nrow(ptfce_image) == 1L, ptfce_image$can_use_as_mask)

signed_stat <- RNifti::readNifti(ptfce_image$artifact_file)
labels <- signed_stat
labels[] <- 0L
labels[signed_stat > 0] <- 1L
labels[signed_stat < 0] <- 2L
table(as.integer(labels))

if (!any(labels != 0L)) {
  stop("No voxels survived at the selected alpha.")
}

mask_file <- file.path(gpa$output_directory, "faces_minus_houses_ptfce_05_signs.nii.gz")
RNifti::writeNifti(labels, mask_file, datatype = "int16")
```

## Extract and check the table

For an independent atlas, set `mask_file` to your aligned integer atlas
instead. If needed, obtain that file by resampling a pre-specified atlas
to the grid of one of the L2 COPEs as described above.

``` r

mask_file <- "/study/atlases/aligned_parcels.nii.gz"
```

Run the following at an interactive R console. For this L2 call, choose
`facehouse` when prompted for the contributing L1 model. For L3
extraction, the prompts also ask for lower-level models. The
`extract_l1`, `extract_l2`, and `extract_l3` arguments each take one
name, `"all"`, `"prompt"`, or `"none"`. Selecting `"all"` at a higher
level still asks about its inputs.

``` r

roi <- extract_glm_betas_in_mask(
  gpa,
  mask_files = mask_file,
  what = "cope",
  extract_l1 = "none",
  extract_l2 = "run_average",
  extract_l3 = "none",
  aggregate = TRUE,
  scheduler = "local",
  out_dir = file.path(gpa$output_directory, "roi_effects")
)

stopifnot(NROW(roi$l2) > 0L)
effect <- roi$l2 |>
  filter(
    l1_cope_name == "faces_minus_houses", l2_cope_name == "EV_Intercept",
    session == 1L
  ) |>
  select(id, session, mask_name, mask_value, value)

effect |>
  group_by(mask_name, mask_value) |>
  summarize(
    n_participants = n_distinct(id),
    n_finite = sum(is.finite(value)),
    .groups = "drop"
  )
```

The function returns a list with `l1`, `l2`, and `l3` elements. By
default it also writes CSVs grouped by L1 model, mask basename,
statistic, and level. The L2 example writes
`roi_effects/L1m-facehouse/<mask_name>_cope_l2.csv.gz`. Keep the
contrast and session columns when joining output tables that contain
several fits. With `write_data = FALSE`, only the returned object is
produced.

For integer masks, `aggregate = TRUE` returns one summary per label and
image, using `mean` unless you set `aggFUN`. Near-zero statistic values
are treated as missing by the extractor to account for voxels outside a
fit’s coverage. An ROI with no usable voxels can have a non-finite
summary, and coverage can vary across participants. The returned `x`,
`y`, and `z` are the label’s centroid in NIfTI world coordinates for
aggregated rows. With `aggregate = FALSE`, the table includes voxel
values, one-based indices `i`, `j`, `k`, and world coordinates for each
voxel. Start with a small selection because a voxelwise table can grow
quickly.

A plot can show the participant values that went into a regional
summary.

``` r

ggplot2::ggplot(effect, ggplot2::aes(x = factor(mask_value), y = value)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey65") +
  ggplot2::geom_jitter(width = 0.12, height = 0, na.rm = TRUE) +
  ggplot2::labs(x = "ROI label", y = "Faces minus houses COPE")
```

If the ROI was selected using the same group contrast and participants,
use this plot descriptively. The selected effect sizes can be inflated
and an ordinary second test of the same effect within the selected
region does not supply independent confirmation. For confirmatory ROI
inference, define regions independently of the tested data, for example
using an a priori atlas, an independent localizer, or held-out
participants. A different contrast from the same data is independent
only under an appropriate design and null hypothesis. [Kriegeskorte et
al. (2009)](https://www.nature.com/articles/nn.2303)

## Extract directly from AFNI or SPM images

In a native AFNI or SPM analysis, the L1 fit supplies one contrast image
per participant and visit by default. The [engine
guide](https://hallquistlab.github.io/fmri.pipeline/articles/glm-engines.md)
shows where to find each fit. Read `contrasts.tsv` from an AFNI L1
directory or `spm_contrasts.tsv` from an SPM L1 directory. Both identify
the effect image in `coef_file`. Select the intended contrast by name
and use the fit’s coverage mask to restrict the summary. The following
example takes a mean per parcel from one fit.

``` r

fit_dir <- "/study/analysis/afni_l1/sub-01/facehouse"
manifest <- read.delim(file.path(fit_dir, "contrasts.tsv"))
coef_file <- manifest$coef_file[manifest$contrast_name == "faces_minus_houses"]
stopifnot(length(coef_file) == 1L, file.exists(coef_file))

estimate <- RNifti::readNifti(coef_file)
atlas <- RNifti::readNifti("/study/atlases/aligned_parcels.nii.gz")
coverage <- RNifti::readNifti(file.path(fit_dir, "mask.nii.gz"))

same_grid <- function(img) {
  identical(dim(img), dim(atlas)) &&
    isTRUE(all.equal(
      RNifti::xform(img, useQuaternionFirst = FALSE),
      RNifti::xform(atlas, useQuaternionFirst = FALSE)
    ))
}
stopifnot(same_grid(estimate), same_grid(coverage))

keep <- atlas > 0 & coverage != 0 & is.finite(estimate)
tapply(estimate[keep], as.integer(atlas[keep]), mean)
```

For SPM, use the SPM fit directory, its manifest, and its model mask
(typically `mask.nii`). Verify the template and voxel units as in [mask
preparation](#prepare-the-mask). Repeat for each fit and keep the
participant, session, model, and contrast with each set of parcel
values.
