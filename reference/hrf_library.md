# Build a small, constrained library of causal HRFs

Build a small, constrained library of causal HRFs

## Usage

``` r
hrf_library(
  tr,
  peak_delays = c(5, 6, 7),
  dispersions = c(0.8, 1, 1.2),
  undershoot_ratios = 6,
  duration = 32,
  preset = c("custom", "standard", "extended", "dense")
)
```

## Arguments

- tr:

  Acquisition interval in seconds.

- peak_delays:

  Positive response-delay parameters in seconds.

- dispersions:

  Positive response-dispersion parameters.

- undershoot_ratios:

  Positive response-to-undershoot ratios.

- duration:

  Kernel duration in seconds (at least two acquisition intervals).

- preset:

  Optional preset grid configuration: \`"custom"\` (default, uses
  supplied \`peak_delays\`, \`dispersions\`, \`undershoot_ratios\`),
  \`"standard"\` (delays 5, 6, 7s; dispersions 0.8, 1, 1.2),
  \`"extended"\` (delays 4–8s; dispersions 0.8, 1, 1.25), or \`"dense"\`
  (delays 4–8.5s by 0.5s; dispersions 0.8, 1, 1.2, 1.4).

## Value

A unit-sum kernel matrix, with time rows, candidate columns, and a
parameters attribute. The canonical candidate (6, 1, 6) is first if
present.

## Details

Candidates use the SPM difference-of-gammas parameterization with
undershoot delay 16 seconds and undershoot dispersion 1. These are
candidate shapes, not a universally validated range for every population
or protocol.
