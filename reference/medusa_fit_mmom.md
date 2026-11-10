# Fit a REPAIR modified method-of-moments (M-MOM / M-MOMS) measurement-error model

Estimates a linear regression model when predictors, the outcome, or
both are likelihood-only predictions measured with heteroscedastic
sampling error (Vig & Hallquist, 2026). Delegates directly to
[`repair::fit_mmom()`](https://hallquistlab.github.io/repair/reference/fit_mmom.html)
when the repair package is installed and uses an identical built-in
fallback otherwise.

## Usage

``` r
medusa_fit_mmom(
  outcome,
  predictors,
  err_x,
  err_y = NULL,
  err_xy = NULL,
  regularization = c("step_up", "minimum", "none", "custom"),
  alpha = NULL,
  backend = c("auto", "r", "rcpp")
)
```

## Arguments

- outcome:

  Numeric outcome vector of length \\N\\.

- predictors:

  Numeric predictor matrix (\\N \times q\\, excluding the intercept).

- err_x:

  Numeric 3D array (\\q \times q \times N\\) of per-observation
  predictor sampling-error covariance matrices
  \\\hat{\boldsymbol{\Psi}}\_{x,j}\\.

- err_y:

  Non-negative numeric vector of length \\N\\ giving per-observation
  outcome sampling-error variances \\\hat{\psi}\_{y,j}\\.

- err_xy:

  Optional numeric matrix (\\N \times q\\) of per-observation
  outcome-predictor sampling-error covariances
  \\\hat{\boldsymbol{\psi}}\_{xy,j}^\top\\.

- regularization:

  Regularization policy: `"step_up"` (default; adaptive M-MOMS step-up
  search), `"minimum"` (fixed minimum \\\alpha = q + 2\\), `"none"`, or
  `"custom"`.

- alpha:

  Optional custom regularization parameter when
  `regularization = "custom"`.

- backend:

  Computation backend passed to
  [`repair::fit_mmom()`](https://hallquistlab.github.io/repair/reference/fit_mmom.html)
  when available (`"auto"`, `"r"`, or `"rcpp"`).

## Value

A list with `converged`, `coefficients` (data.frame with `term`,
`estimate`, `se`, `zstat`, `pval`, `ci_low`, `ci_high`), `vcov`,
`preliminary_coef`, `corrected_design_moment`, `analysis_weights`,
`diagnostics`, and `repair_fit` (the raw repair object when repair is
installed).
