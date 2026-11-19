# Function to walk user through setting up FWE corrections for model outputs

Function to walk user through setting up FWE corrections for model
outputs

## Usage

``` r
build_fwe_correction(gpa, lg = NULL)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing configured L3 models

- lg:

  a logger object for logging messages

## Value

returns a modified gpa object that contains FWE specifications in
\$fwe_correction

## Details

currently only looks a level 3
