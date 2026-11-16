# Construct an AFNI 3dLMEr command string

Construct an AFNI 3dLMEr command string

## Usage

``` r
build_3dlmer_command(
  prefix,
  model_formula,
  qVars = NULL,
  qVarCenters = NULL,
  glt_codes = NULL,
  data_table_file,
  mask = NULL,
  njobs = 1,
  ss_type = 3
)
```

## Arguments

- prefix:

  Output prefix for 3dLMEr.

- model_formula:

  Fixed-effects formula string.

- qVars:

  Character vector of quantitative variables.

- qVarCenters:

  Optional numeric vector of centering values corresponding to `qVars`.

- glt_codes:

  List, character vector, or data frame of named GLT code strings.

- data_table_file:

  Path to the data-table text file.

- mask:

  Optional path to the brain-mask file.

- njobs:

  Number of parallel jobs for 3dLMEr.

- ss_type:

  Sum-of-squares type (1–3; defaults to 3).

## Value

A character string containing the 3dLMEr command.
