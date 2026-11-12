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
  preset = c("custom", "standard", "extended", "canonical_27", "dense", "manifold",
    "enet", "medusa"),
  normalize = c("auto", "signed_sum", "positive_lobe", "peak"),
  n_enet = 45L
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

  Positive response-to-undershoot ratios (`> 1`; `Inf` gives a pure
  single-gamma HRF with zero post-stimulus undershoot).

- duration:

  Kernel duration in seconds (at least two acquisition intervals).

- preset:

  Optional preset grid configuration: `"custom"` (default, uses supplied
  `peak_delays`, `dispersions`, `undershoot_ratios`), `"standard"`
  (delays 5, 6, 7s; dispersions 0.8, 1, 1.2), `"extended"` (delays 4–8s;
  dispersions 0.8, 1, 1.25), `"canonical_27"` (fixed-undershoot 0.5s
  canonical grid: delays 4–8s by 0.5s; dispersions 0.9, 1.0, 1.1; 27
  candidates), `"dense"` (delays 4–8.5s by 0.5s; dispersions 0.8, 1,
  1.2, 1.4), `"manifold"` or `"medusa"` (pruned 3D physiological
  manifold spanning fast/narrow/minimal-undershoot subcortical through
  slow cortical shapes; ~44 candidates), or `"enet"` (maximin cosine
  epsilon-net covering the expanded 3D parameter cube with `n_enet`
  candidates).

- normalize:

  Normalization convention across candidate columns: `"auto"` (uses
  `"positive_lobe"` whenever undershoot varies across candidates,
  includes `Inf`, differs from the canonical fixed ratio 6, or uses
  `"manifold"`, `"medusa"`, or `"enet"` so varying undershoot ratios do
  not inflate deconvolved neural amplitudes; uses `"signed_sum"` for
  legacy fixed-ratio-6 grids), `"signed_sum"` (`sum(h) == 1`),
  `"positive_lobe"` (`sum(pmax(h, 0)) == 1`), or `"peak"`
  (`max(h) == 1`).

- n_enet:

  Number of candidates retained when `preset = "enet"` (default 45).

## Value

A normalized kernel matrix, with time rows, candidate columns, and
`parameters`, `tr`, and `normalize` attributes. The canonical candidate
`(6, 1, 6)` is first if present.

## Details

Candidates use the SPM difference-of-gammas parameterization with
undershoot delay 16 seconds and undershoot dispersion 1. Setting
`undershoot_ratios = Inf` omits the negative gamma lobe (single-gamma
HRF).
