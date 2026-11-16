# Internal function to add, edit, rename, delete contrasts

Internal function to add, edit, rename, delete contrasts

## Usage

``` r
specify_contrasts(
  mobj = NULL,
  signals = NULL,
  spec_list = NULL,
  preselected_pairwise = NULL
)
```

## Arguments

- mobj:

  `l1_model_spec` or `hi_model_spec` object

- signals:

  a list of signals that provide information about regressor
  specification

- spec_list:

  a list of settings for this model based on an external YAML/JSON
  specification file

- preselected_pairwise:

  Optional character vector of factor names preselected for pairwise
  differences in interactive higher-level model setup.

## Author

Michael Hallquist
