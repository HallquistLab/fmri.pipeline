# helper function to create an FWE object that specifies the type of FWE correction and the model outputs to which it is applied

helper function to create an FWE object that specifies the type of FWE
correction and the model outputs to which it is applied

## Usage

``` r
create_fwe_spec(gpa, level = 3L, lg = NULL)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing configured L3 models;
  estimation need not be complete

- level:

  an integer indicating what modeling level this FWE correction applies
  to (1, 2, or 3)

- lg:

  a logger object for logging messages

## Details

note that this only works for level = 3 for now
