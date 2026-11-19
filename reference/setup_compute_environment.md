# Configure the compute environment for pipeline jobs

Configure the compute environment for pipeline jobs

## Usage

``` r
setup_compute_environment(gpa, preselect_action = NULL)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object

- preselect_action:

  optional menu action to apply without prompting.

## Value

a modified gpa object containing settings for the compute environment

## Details

Set shell initialization commands in `gpa$parallel$compute_environment`.
Use `global` for commands shared by all jobs and the `fsl`, `afni`, or
`spm` entry for commands specific to an engine. For SPM, make the
configured MATLAB-compatible executable available here and set the SPM
installation path in `gpa$glm_settings$spm$spm_path`.
