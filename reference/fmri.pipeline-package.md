# fmri.pipeline: Hierarchical fMRI GLMs with FSL, AFNI, and SPM

Specify, run, and track hierarchical general linear models (GLMs) for
task fMRI data with FSL, AFNI, and SPM. Generate suite-specific designs
and contrasts from shared trial, run, and participant tables, and submit
time-series and group analyses locally or through a cluster scheduler.

## Details

Start with
[`setup_glm_pipeline`](https://hallquistlab.github.io/fmri.pipeline/reference/setup_glm_pipeline.md)
to select the GLM engine and provide the analysis tables. Use
[`run_glm_pipeline`](https://hallquistlab.github.io/fmri.pipeline/reference/run_glm_pipeline.md)
to submit the selected models. For a worked comparison of how each suite
combines runs and supplies group contrasts, see
[`vignette("glm-engines", package = "fmri.pipeline")`](https://hallquistlab.github.io/fmri.pipeline/articles/glm-engines.md).

## See also

Useful links:

- <https://hallquistlab.github.io/fmri.pipeline/>

- Report bugs at <https://github.com/HallquistLab/fmri.pipeline/issues>

## Author

**Maintainer**: Michael Hallquist <michael.hallquist@gmail.com>

Authors:

- Michael Hallquist <michael.hallquist@gmail.com>
