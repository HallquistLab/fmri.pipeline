# Create a familywise-error correction specification

Create a declarative, serializable specification describing which
completed group-analysis contrasts should receive a familywise-error
correction and which correction method should be used. The specification
contains no live execution objects and can therefore be stored safely in
a pipeline object or YAML file.

## Usage

``` r
fwe_spec(
  name,
  targets,
  method = "ptfce",
  level = 3L,
  correction_mask = "model",
  compute = list(scheduler = "inherit"),
  schema_version = fwe_spec_schema_version(),
  ...
)
```

## Arguments

- name:

  unique name for this correction.

- targets:

  non-empty named list selecting level-3 model and contrast fields.
  Valid fields are `session`, `l1_model`, `l1_cope_number`,
  `l1_cope_name`, `l2_model`, `l2_cope_number`, `l2_cope_name`,
  `l3_model`, `l3_cope_number`, and `l3_cope_name`. Values may be
  vectors.

- method:

  correction method name or named method configuration. Currently
  supported specifications are `"ptfce"` and
  `"afni_3dclustsim_permutation"`. Only pTFCE currently has a planning
  and execution adapter; the AFNI specification can be saved for future
  use.

- level:

  GLM level to correct. Only level 3 is currently supported.

- correction_mask:

  `"model"` to use each analysis mask, a NIfTI path, or a list with
  `source` equal to `"model"` or `"file"` (and `path` for the latter).

- compute:

  named list of compute settings. By default the scheduler is inherited
  from the pipeline.

- schema_version:

  specification schema version. Normally left at its default.

- ...:

  named method-specific options overriding method defaults.

## Value

an object of class `fwe_spec`.
