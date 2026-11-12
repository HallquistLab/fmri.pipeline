# MEDuSA: from blurred BOLD signals to trial-level questions

## Start with a scientific question

Suppose a learning task contains a cue, a choice, and feedback on every
trial. You want to know whether a reward prediction error (RPE) predicts
activity specifically around feedback, and whether that relationship
differs across people. A conventional parametric GLM asks whether BOLD
covaries with an RPE-weighted, HRF-convolved event regressor. **MEDuSA
(Multilevel Event-locked Deconvolved Signal Analysis) instead supports
trial-level modeling of event-aligned estimates of latent activity
across a range of time lags.** These are complementary ways of asking
about brain–behavior relationships.

The challenge is that feedback does not arrive on an empty stage: BOLD
from the cue and choice is still evolving, and the next trial may start
before that response has ended. Recovering a feedback-specific effect
therefore requires choices about vascular response shape,
regularization, event overlap, the duration of neural activity, and how
trial-level estimation uncertainty is propagated to group mixed models
and errors-in-variables connectivity models. PRs \#42–#64 improve each
stage of this problem; they do not turn deconvolution into a direct
recording of neuronal firing.

This primer assumes familiarity with regression, covariance, and matrix
multiplication. Read it before the [worked
workflow](https://hallquistlab.github.io/fmri.pipeline/articles/medusa-workflows.md).
The [technical
review](https://github.com/HallquistLab/fmri.pipeline/blob/master/medusa_deconvolution_review_prs_42_53.qmd)
preserves earlier derivations and simulation demonstrations; its
filename is retained for existing links, while this vignette and
`medusa-workflows.Rmd` document the full architecture through PR \#64.

## Three quantities that should stay distinct

| Quantity | What it represents | What we can infer |
|----|----|----|
| Observed BOLD, $`y(t)`$ | The measured, preprocessed fMRI signal | A vascularly filtered and noisy consequence of activity |
| Latent activity estimate, $`\hat{s}(t)`$ | A regularized inverse estimate under an assumed HRF | A model-dependent activity proxy on the acquisition grid |
| Event-level estimate, $`\hat{a}_i(\tau)`$ | A readout or fitted source amplitude for trial $`i`$ at lag $`\tau`$ | A quantity to relate to behavior, conditional on the readout/source model |

A finer interpolation grid supplies smoother alignment to events
occurring between acquisitions. It does **not** create independent
sub-TR measurements. Likewise, a sharp peak in an estimated activity
trajectory does not establish millisecond timing or a neuronal firing
rate in physiological units.

### The forward model: why BOLD is blurred

Convolution means that each brief input contributes a delayed, extended
response. If the input is a sequence of impulses, multiply each trial’s
amplitude by a shifted HRF and add the responses. In matrix form,

``` math
\mathbf y = H\mathbf s + X_0\boldsymbol\gamma + \mathbf e.
```

Here $`\mathbf y`$ contains the observed time points; each column of
$`H`$ is the HRF shifted to one possible activity time; $`\mathbf s`$
contains latent activity at those times; $`X_0\boldsymbol\gamma`$
represents nuisance effects; and $`\mathbf e`$ is residual noise.
Ordinary task GLMs constrain $`\mathbf s`$ through event regressors
before fitting. Continuous deconvolution estimates a much more flexible
$`\mathbf s`$ and then aligns it to events.

[Glover (1999)](https://pubmed.ncbi.nlm.nih.gov/10191170/) provides an
early fMRI precedent for deblurring with deconvolution, including
discussion of where linearity breaks down. The relevant connection is
the inverse problem, not a validation of all current MEDuSA options.

## Why we need regularization

Many activity sequences can generate almost the same smooth BOLD signal.
Trying to invert the blur exactly can amplify a small measurement error
into a large oscillation. Regularization makes the problem tractable by
preferring some solutions over others:

``` math
\hat{\mathbf s} = \arg\min_{\mathbf s}
\left\{\underbrace{\|W(\mathbf y-H\mathbf s)\|^2}_{\text{fit the observations}}
+\underbrace{\lambda\mathbf s^\top P\mathbf s}_{\text{discourage implausible activity}}\right\}.
```

$`W`$ accounts for temporal noise covariance, $`P`$ specifies what is
penalized, and $`\lambda`$ controls the tradeoff. The numerical
implementation also adds a small stabilizing diagonal term. Increasing
$`\lambda`$ usually gives a more stable, more biased estimate. The
penalty is part of the scientific measurement model, not merely a
computational setting.

| Penalty | Intuition | Practical consideration |
|----|----|----|
| `ridge` | Discourage large activity values | Can substantially attenuate amplitudes |
| `diff1` | Discourage abrupt changes between neighboring samples | Does not penalize constant offsets |
| `diff2` | Discourage rapid changes in slope | Does not penalize constants or linear ramps |
| `matern1` | Mix amplitude control with first-difference smoothing | Anchors low frequencies as well as smoothing |
| `sobolev` | Mix amplitude, first-difference, and curvature penalties | Can reduce baseline wander and ringing in suitable designs |

The last two are fixed discrete quadratic penalties in this package.
Their names connect them to smoothness and covariance models in
statistics; they do not imply a fitted Gaussian process with an
estimated correlation length. They also do not replace drift removal or
motion handling.

**Generalized cross-validation (GCV)** selects a tradeoff by balancing
reconstruction error against the effective flexibility of the fitted
smoother. [Golub, Heath, and Wahba
(1979)](https://doi.org/10.1080/00401706.1979.10489751) is the
statistical precedent. Both
[`deconvolve_reglin()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_reglin.md)
and
[`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_crossrun_hrf.md)
support automatic prewhitened GCV tuning when `lambda = NULL` (in
[`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_crossrun_hrf.md),
`hrf_tune_decon_lambda()` pools across all usable runs in an
order-invariant manner and balances observation-space GCV against
Stage-3 gain-restored variance and modulation transfer). In MEDuSA, a
good BOLD reconstruction score alone need not identify the best
trial-level RPE estimate. Selecting settings because they maximize the
eventual RPE correlation would bias the analysis toward that result; use
design-matched simulations, independent data, or training-only
selection.

PRs \#50–#51 exploit the finite temporal reach of the HRF and AR(1)
noise model to accelerate matrix calculations. Their main scientific
yield is feasibility: more realistic validation, sensitivity analyses,
and larger datasets. Faster algebra alone does not increase the
information in the scans. The optimized continuous solver and the
cross-run API have different memory costs; the latter can still
construct/cache large inverse operators.

## Learn the HRF without fitting the target effect

### A vascular delay can look like a cognitive delay

If the assumed HRF is too slow, deconvolution may move estimated
activity earlier to explain the observed BOLD. A mismatch can therefore
change which within-trial event appears to carry a signal. [Handwerker,
Ollinger, and D’Esposito
(2004)](https://pubmed.ncbi.nlm.nih.gov/15050587/) documented HRF
variation across subjects and regions and its consequences for
statistical analysis. This motivates checking HRF sensitivity, not
assigning a universally fast or slow HRF to a named brain region.

[`estimate_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/estimate_crossrun_hrf.md)
chooses from a constrained, causal HRF library generated by
[`hrf_library()`](https://hallquistlab.github.io/fmri.pipeline/reference/hrf_library.md),
which offers `"custom"`, `"standard"` (9 candidates), `"extended"` (15
candidates), `"canonical_27"` (27 candidates: fixed $`r = 6`$ undershoot
across a $`0.5\text{ s}`$ delay grid over $`4.0\text{--}8.0\text{ s}`$
with dispersions $`0.9, 1.0, 1.1`$), `"dense"` (40 candidates),
`"medusa"` / `"manifold"` (~44 candidates: pruned 3D physiological
manifold), and `"enet"` (`n_enet = 45` maximin cosine $`\epsilon`$-net)
preset grids. For target run A, it learns the shape using B and C, then
freezes that shape before deconvolving A. Repeat this for each run. This
is **cross-fitting**: the target run cannot choose its own HRF. In
addition to discrete candidate selection (`selection_mode = "hard"`),
[`estimate_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/estimate_crossrun_hrf.md)
supports **soft evidence-weighted kernel averaging**
(`selection_mode = "soft_average"`), which forms a normalized composite
HRF (`selected_kernels`) over the top locally compatible shapes
(`kernel_cos >= 0.85`) weighted by spread-normalized Boltzmann evidence
and shrunk toward `fallback_kernel` by Empirical-Bayes departure
confidence $`\omega \in [0, 1]`$. Soft kernel averaging eliminates
discrete grid quantization steps on sub-grid vascular delays and
out-of-family nonlinear Windkessel-Balloon hemodynamics without
requiring an overly dense ($`0.25\text{ s}`$) candidate grid. The
separation addresses one source of circularity; it does not eliminate
HRF–neural-duration ambiguity, guarantee transfer across runs, or make
folds independent participants. The general concern about selecting and
testing on the same data is explained by [Kriegeskorte et
al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2841687/).

### Why a pruned 3D physiological manifold with positive-lobe normalization?

Flexible nonlinear HRF parameterizations such as the 7-parameter
superposition of three inverse-logit sigmoids ([Lindquist et al.,
2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC3318970/)) can fit
arbitrary trial-averaged waveforms, but in single-voxel or parcel
event-related fMRI their unconstrained rise times and plateau widths
admit non-zero pre-onset response ($`h(0) > 0`$) and absorb neural
duration or latency differences directly into the vascular filter.
Conversely, canonical cortical double-gamma grids with fixed peak delays
($`\ge 4\text{ s}`$), narrow dispersion bands ($`0.9\text{--}1.1`$), and
a fixed post-stimulus undershoot ratio ($`r = 6`$) miss the faster
($`2.5\text{--}3.5\text{ s}`$), narrower
($`\text{FWHM} \approx 2.3\text{--}3.2\text{ s}`$), and minimal- or
absent-undershoot responses documented in subcortical nuclei, thalamus,
brainstem, primary sensory cortex, and deep gray matter ([Lewis et al.,
2018](https://doi.org/10.7554/eLife.35710); [Chen et al.,
2023](https://doi.org/10.1162/imag_a_00018); [Prince et al.,
2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9708069/)).

PR \#66 resolves this tradeoff across three coordinated layers:

1.  **Pruned 3D physiological manifold (`preset = "medusa"` /
    `"manifold"`) and $`\epsilon`$-net (`preset = "enet"`):**
    [`hrf_library()`](https://hallquistlab.github.io/fmri.pipeline/reference/hrf_library.md)
    spans three physiological axes—peak delay
    ($`2.5\text{--}8.0\text{ s}`$ by $`0.5\text{ s}`$), a TR-aware
    dispersion floor ($`\ge 0.60`$ at fast
    $`\text{TR} \le 1.2\text{ s}`$, bounded below by
    $`\min(0.80, 0.45\,\text{TR})`$ at slower TRs to prevent sub-Nyquist
    narrower-than-TR aliasing), and undershoot ratio
    $`r \in \{3, 6, 12, \infty\}`$ (where `undershoot_ratios = Inf`
    omits the negative gamma lobe to form a pure single-gamma HRF).
    Rather than scoring the full $`240`$-candidate Cartesian product,
    `"medusa"` / `"manifold"` prunes biophysically implausible corners
    (such as ultra-fast peaks paired with broad dispersion or deep
    undershoot) and applies deterministic greedy farthest-point
    selection in $`L_2`$ waveform cosine space
    (`hrf_select_maximin_subset()`, seeded at canonical `(6, 1, 6)`) to
    retain 44 non-redundant shapes; `"enet"` constructs a maximin cosine
    $`\epsilon`$-net of size `n_enet` across the full 3D cube, while
    `"canonical_27"` preserves the legacy 27-candidate fixed-undershoot
    ($`r = 6`$) grid.
2.  **Broad cross-undershoot scale protection
    (`normalize = "positive_lobe"`):** Under classical net-integral
    normalization (`"signed_sum"`, $`\sum_t h_t = 1`$), varying the
    post-stimulus undershoot distorts positive peak height: an HRF with
    deep undershoot ($`r = 3`$, ~33% negative-to-positive lobe area) has
    a smaller net sum and therefore an artificially ~1.45$`\times`$
    taller positive peak than a zero-undershoot ($`r = \infty`$)
    subcortical HRF, compressing deconvolved neural amplitudes and
    Stage-3 slopes in deep-undershoot regions by ~33–50%. When
    `normalize = "auto"` (default),
    [`hrf_library()`](https://hallquistlab.github.io/fmri.pipeline/reference/hrf_library.md)
    automatically switches to `"positive_lobe"`
    ($`\sum_t \max(h_t, 0) = 1`$) whenever `preset` is `"medusa"`,
    `"manifold"`, or `"enet"`, or whenever `undershoot_ratios` varies,
    includes `Inf`, or differs from `6`. In addition,
    `build_decon_kernel_grid()` (`infer_kernel_grid_normalize()`)
    detects heterogeneous negative-to-positive lobe ratios
    (`diff(range(neg / pos)) > 0.05`) even on unattributed custom
    matrices, propagating the resolved `kernel_normalize` across
    [`estimate_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/estimate_crossrun_hrf.md),
    [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_crossrun_hrf.md)
    (including `lambda = NULL` GCV tuning and Wiener fits),
    [`deconvolve_reglin()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_reglin.md),
    [`voxelwise_deconvolution()`](https://hallquistlab.github.io/fmri.pipeline/reference/voxelwise_deconvolution.md),
    [`diagnose_crossrun_recovery()`](https://hallquistlab.github.io/fmri.pipeline/reference/diagnose_crossrun_recovery.md),
    and
    [`deconvolution_impulse_diagnostics()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolution_impulse_diagnostics.md).
3.  **Sub-quadratic C++ scoring and native cascaded AR($`p`$)
    whitening:** To score $`K \approx 44\text{--}45`$ candidates without
    slowing Step 1 or Step 2, `src/crossrun_hrf.cpp` orthonormalizes
    tall-skinny design matrices via `arma::qr_econ` + $`p \times p`$ SVD
    rather than $`O(m^3)`$ full SVD, applies exact cascaded AR($`p`$) /
    Toeplitz ACF / spike-floor whitening in $`O(m p)`$ across contiguous
    segments and censor gaps via `ArWhitener` (eliminating dense
    $`m \times m`$ whitening matrices), uses lazy $`O(K)`$
    best-vs-fallback encompassing for $`K > 18`$ with
    spectral-rank-bounded incompatibility multiplicity
    ($`k_{\text{eff}} \le 18`$), quantizes parcel AR strata to at most 6
    representative profiles during large-library Step-1 scoring, and
    screens coupled trial-free condition designs on a 12-anchor maximin
    subset.

### What do the HRF training amplitudes mean?

With `amplitude_model = "condition"`, training uses one amplitude per
condition in each run. With `"trial"`, training allows event-specific
amplitudes while sharing a candidate HRF shape. In dense, closely
coupled multi-condition designs, the implementation compares designs
freeing one condition at a time and keeps the most informative training
design, limiting the ambiguity from fitting all neighboring events
freely at once. Imagine feedback responses proportional to signed RPE:
large positive and negative responses can average close to zero. The
condition model can then have little mean response from which to learn
an HRF, even though trial variation is meaningful. Free trial amplitudes
can help, but require more parameters and enough independent timing
information. The training function never takes RPE values as input.

This distinction resembles a familiar choice between condition-level
GLMs and beta-series models. [Mumford et
al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3251697/) explains
why individual-trial estimates become unstable when event regressors are
highly correlated and compares ways of addressing that problem. MEDuSA’s
free-event HRF fitting is not their least-squares-separate algorithm.

### Why colored noise and conservative fallback matter

A slowly varying noise sequence can resemble a smooth evoked response.
Treating successive fMRI volumes as independent can overstate evidence
for a customized HRF. AR(1) models neighboring errors as correlated;
whitening transforms the fitting problem so that such persistence
contributes less misleading evidence. [Woolrich et
al. (2001)](https://www.fmrib.ox.ac.uk/datasets/techrep/tr01mw1/tr01mw1/node1.html)
gives the fMRI context for temporal autocorrelation modeling.

The current cross-run default uses `noise_model = "ar1"` and
`evidence = "calibrated"`. It estimates noise covariance (including
post-selection 2-pass alternating REML AR($`p`$) noise refitting under
the target-excluded selected HRF, gap-local PACF estimation across
scrubbed frames, exact subset Toeplitz Cholesky whitening across censor
gaps, and heterogeneity-gated SUSAN bilateral noise pooling), compares
candidates, and retains the prespecified fallback when evidence is weak,
inconsistent, or task-free (0-event runs). PR \#52 refined the
candidate/fallback subspace comparison (using principal-angle F-tests
when task rank exceeds 3), added spectral effective-candidate
multiplicity correction ([Galwey,
2009](https://doi.org/10.1002/gepi.20408)) for extended/dense libraries
($`K > 9`$), and refined the design used for free-trial training, making
the comparison more useful in dense designs. A fallback is a recorded
modeling decision, not a failed voxel and not proof that the canonical
HRF is true. The evidence threshold corrects library comparisons, **not
whole-brain multiplicity**, and is not a posterior HRF probability.

### Borrowing information across voxels

An independently defined parcel can provide a more stable HRF reference
than a noisy single voxel. Hard `pooling = "parcel"` shares a selected
shape. Soft voxel regularization uses a parcel anchor while permitting
voxel differences. Parcel anchor admissibility supports both strict
voxel counting and pooled parcel evidence via the Cauchy combination
test ([Liu & Xie, 2020](https://doi.org/10.1080/01621459.2018.1554485)),
which remains valid under arbitrary spatial correlation within a parcel.
PR \#52 also adds Empirical-Bayes concordance weighting ($`c_v^2`$) so
voxels with disagreeing candidate-score profiles borrow less from
neighboring profiles; optional 3D spatial coordinates (`spatial_coords`,
`spatial_fwhm`) make borrowing locally Gaussian-weighted within each
parcel.

This pools information about **HRF shape**, not neural activity values.
Parcel boundaries need not match vascular territories, and hard pooling
can hurt heterogeneous regions. Start with voxel selection and no
spatial penalty as a transparent comparison; add borrowing only with a
justified atlas and sensitivity checks. [Prince et al. (2022),
GLMsingle](https://pmc.ncbi.nlm.nih.gov/articles/PMC9708069/) is a
useful neighboring approach: candidate HRFs and regularized single-trial
estimates have benefits in other fMRI workflows. Its estimator and
empirical results should not be attributed to MEDuSA.

## Put amplitudes on an interpretable scale

Suppose two voxels have the same measurement noise but one responds
twice as strongly to the task. Dividing each voxel by its **total**
standard deviation makes the stronger voxel’s denominator larger too,
compressing the difference you wanted to study. The same issue arises
when runs differ in task modulation variance. Z-scoring is a legitimate
choice for standardized within-run questions, but changes the amplitude
estimand.

`scaling = "noise_sd"` instead divides by an estimated **marginal
observation-noise SD** ($`\sigma_{\text{marg}}`$)—not the AR innovation
SD $`\sigma_u`$ and not total signal-plus-noise SD. In modern cross-run
fits with adequate residual degrees of freedom,
`hrf_estimate_gls_noise_sd()` estimates $`\sigma_{\text{marg}}`$ from
leverage-standardized, Winsorized GLS residuals under the selected HRF
and AR($`p`$) correlation model; AR($`p`$)-adjusted second differences
($`y_t-2y_{t-1}+y_{t-2}`$ via `estimate_decon_noise_sd()`) serve as a
residualized omitted-event guard and standalone/low-df fallback. Because
`"noise_sd"`, `"demean"`, and `"zscore"` center each run on usable
observations, tonic between-run or between-person baseline levels are
removed by design. This remains a working noise estimate: severe
unmodeled neural signals or misspecified covariance can still affect it.

| Mode | Units of input to the inverse | Suitable interpretation |
|----|----|----|
| `zscore` | Total within-run SD units | Standardized activity relative to total signal variation |
| `noise_sd` | Estimated noise SD units | Activity relative to estimated measurement noise |
| `demean` | Original input units after centering | Preserve a meaningful acquisition/preprocessing scale |
| `psc` | Percent of a positive BOLD baseline | Percent-change input; unsuitable for zero-mean residuals |
| `none` | Unchanged input units | Already centered data or controlled simulations |

`restore_scale = TRUE` in epoch regression reverses the recorded divisor
for amplitude reporting. It does not restore an identified neural
baseline, undo preprocessing, or calibrate spikes per second. With a
fixed linear operator, multiplying a whole run by a constant does not
change a within-run Pearson correlation; pooled correlations and
amplitude comparisons across differently scaled runs can change. Keep
this distinct from the benefits of de-mixing.

## Alignment and de-mixing answer different questions

### Alignment asks: which samples support this epoch?

The established
[`voxelwise_deconvolution()`](https://hallquistlab.github.io/fmri.pipeline/reference/voxelwise_deconvolution.md)
→
[`evt_align_decon_files()`](https://hallquistlab.github.io/fmri.pipeline/reference/evt_align_decon_files.md)
route remains useful. Recent changes extract NIfTI values internally,
check spatial geometry, retain invalid voxels with status labels, record
scaling and file identities, and interpolate each unit before averaging.
An interpolation should not bridge a censoring gap or combine
disconnected trial segments.

Overlap policies make a scientific choice explicit: `allow` retains
requested windows, `truncate` respects supplied collision boundaries,
and `exclude` removes epochs whose requested target windows overlap.
Truncation does not remove the already blurred contribution of a
preceding event. Excluding trials can also change the behavioral
distribution being analyzed. Inspect missingness and `epoch_status`, not
just the final plot.

### De-mixing asks: how did the estimator mix neighboring sources?

Even with the correct HRF, a regularized inverse spreads and shrinks a
brief input. Write the fixed inverse as $`A`$, so that
$`\hat{\mathbf s}=A\mathbf y`$. Then, ignoring nuisance terms for this
illustration,

``` math
\mathbb E[\hat{\mathbf s}\mid\mathbf s]=\underbrace{AH}_{R}\mathbf s.
```

$`R`$ is the **resolution matrix**: feed in an isolated source and its
column shows the expected recovered blur. If a true amplitude of 1
produces a readout of 0.2, the readout has gain 0.2. If it also
contributes 0.1 to the following trial, there is cross-talk. A raw slope
of 0.2 can thus reflect an estimator’s attenuation rather than weak
neural coding.

[`decon_epoch_regression()`](https://hallquistlab.github.io/fmri.pipeline/reference/decon_epoch_regression.md)
uses the actual event schedule, frozen HRF, scaling, and inverse
settings to predict and invert these distortions:

1.  **Condition-mean and joint 2nd-order Tikhonov non-target de-mixing
    (`demix_other_events = TRUE`):** When non-target events (such as
    cues or choices) precede the target event, projecting out
    $`M_{\text{other}}`$ individual non-target trial columns consumes
    degrees of freedom and inflates target-trial variance. By default
    (`demix_ridge >= 0.03`, default `0.05`),
    `build_other_demix_operator()` estimates and subtracts the
    **condition-mean** non-target response across events in
    `events[[r]]`, removing schedule-confounded cue bleed without
    inflating single-trial noise. When non-target events carry strong
    trial-varying parametric confounds that correlate with the target
    predictor, setting `demix_ridge < 0.03` (e.g., `0.02`) activates a
    2nd-order iterated Tikhonov joint $`(T, O)`$ de-mixing operator.
2.  **2nd-order Tikhonov target-trial de-mixing, per-trial gain
    restoration, and matched covariate alignment
    (`demix_target_trials = TRUE`):** Forms the $`N_r \times N_r`$
    stratum resolution matrix $`K_\tau = S_\tau D G_\tau`$ and applies a
    2nd-order iterated Tikhonov inverse
    (`solve_duration_demix_operator()`) with per-trial effective gain
    restoration (`eff_gain_vec`) unified across `"impulse"`, `"window"`,
    and `"duration"` source models, while gating unobservable
    boundary/gap trial-lags (`obs_frac_tau < 0.20`) to `NA`. When
    fitting slopes against time-varying covariates, centered predictors
    are passed through the matched post-demixing resolution operator
    $`M_{\text{post}}`$ (`P0-1` in PR \#64) so residual off-diagonal
    target shrinkage at rapid ISIs does not bias estimated slopes.
3.  **Within-trial multi-lag point-spread de-mixing
    (`demix_within_trial = TRUE`):** When $`\ge 3`$ lags
    $`\tau_1, \dots, \tau_L`$ are analyzed simultaneously, temporal
    smoothing and HRF spread cause activity at earlier lags to leak into
    later lags within the same trial. A 2nd-order Tikhonov multi-lag
    inverse (`within_trial_ridge = 0.06`) sharpens descriptive
    within-trial waveform recovery across the requested lag grid while
    preserving cross-trial operator covariances (leave
    `demix_within_trial = FALSE` when passing outputs to Stage-4
    [`medusa_multilag_infer()`](https://hallquistlab.github.io/fmri.pipeline/reference/medusa_multilag_infer.md),
    which performs group-level spectral de-blurring and correlated-SURE
    P-spline smoothing).
4.  **Operator-propagated uncertainty, Empirical-Bayes precision
    weights, and spectrum-continuous ZCA whitening:** Because trials
    preceded by short jitter intervals or rapid neighboring events
    undergo stronger unmixing amplification, their de-mixed estimates
    have larger error variance and non-zero trial-to-trial error
    covariance:
    ``` math
    \Sigma_{\text{raw},r} = \sigma_{e,r}^2 \, T_{\text{raw},r} \Omega_r T_{\text{raw},r}^\top, \qquad
    \Sigma_{\text{demix},r} = G_r^\top \Sigma_{\text{raw},r} G_r.
    ```
    [`decon_epoch_regression()`](https://hallquistlab.github.io/fmri.pipeline/reference/decon_epoch_regression.md)
    extracts exact per-trial standard errors (`se_raw`, `se_demixed`),
    lag-1 trial error autocorrelation (`error_ac1`), within-unit mean-1
    precision weights (`precision_weight`), and unnormalized
    cohort-comparable inverse-variance weights
    (`precision_weight_cohort`, preserving between-participant
    reliability differences for
    `mixed_by(..., weights = "precision_weight_cohort")`). It supports
    precision-weighted WLS (`weighting = "precision"`), ZCA-whitened GLS
    (`weighting = "gls"`), and exported within-run ZCA-whitened panel
    variables (`whiten_panel = TRUE`). `compute_zca_whitening_matrix()`
    preserves run-mean null modes via exact centered-trace scaling
    ($`\operatorname{tr}(Q_m \Phi Q_m)/(M-1) = 1`$), damps near-singular
    modes smoothly via a spectrum-continuous fourth-root floor
    $`(\lambda_+^2 + \epsilon^2)^{-1/4}`$, and scales off-diagonal
    whitening adaptively by the Bartlett lag-1 autocorrelation of the
    trial predictors.

Correction remains conditional on the source model. Events that are
indistinguishable in the design stay difficult to separate. Small gains
amplify noise, and ridge stabilization deliberately retains bias. An
attractive corrected curve is insufficient evidence of accurate
recovery. Simulate both positive effects and active nulls: for example,
strong cue activity with **zero** feedback RPE effect.

## Duration changes the question being estimated

A decision process lasting 2 seconds can create more integrated activity
than one lasting 1 second, even at exactly the same activity rate. If
value predicts faster decisions, an impulse model may mistake shorter
duration for weaker value coding. A nuisance RT column alone is not
equivalent to changing the shape of each event’s source.

PR \#54 adds `source_model = "duration"`. A target source occupies the
interval from `onset + lag` to `onset + lag + duration`. Readout adjusts
to the interval’s midpoint and width, and the correction accounts for a
constant task-related rate whose integrated contribution changes with
duration. [Grinband et
al. (2008)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2654219/) provides
a direct neuroscience precedent for variable-duration event models. RT
is a hypothesized duration marker, not a measurement of all neuronal
processing.

| Source model | Hypothesis at lag $`\tau`$ | Amplitude meaning |
|----|----|----|
| `impulse` | Brief source at onset + $`\tau`$ | Unit-impulse amplitude in the package’s discrete convention |
| `window` | Fixed-width boxcar centered at onset + $`\tau`$ | Height of a local sustained source |
| `duration`, `rate` | Trial-specific interval beginning at onset + $`\tau`$ | Boxcar height per acquisition bin; longer intervals have greater integrated activity |
| `duration`, `integral` | Trial-specific interval normalized to unit summed mass | Total discrete source mass, comparable across durations |

“Rate” is a model coefficient, not calibrated firing rate. At fixed TR,
discrete mass is proportional to a time integral; comparisons across TRs
require the sampling convention to be reconciled. For duration models,
lag zero means an interval **starting at the recorded event**, not a
point sample at onset. Each lag is a separate shifted-source analysis;
the resulting curve is not an unconstrained reconstruction of a single
within-trial firing waveform.

Use the actual decision-start onset with RT duration. A response-time
onset followed by another RT-long interval would describe a different
process. Positive durations in HRF training event tables and in epoch
regression should express compatible hypotheses. `trial_durations`
overrides target durations for the readout stage only; it does not refit
the HRF.

## What the recent changes buy you

| PRs | Problem addressed | Practical MEDuSA yield |
|----|----|----|
| \#42 | Extraction, geometry, missing signals, file identity | More trustworthy mapping from image voxels to aligned observations |
| \#43 | Implicit units and approximate Wiener tuning | Explicit scale records and corrected finite-observation tuning |
| \#44 | Disconnected interpolation and hidden overlap decisions | Auditable epoch support and impulse diagnostics |
| \#45–#46 | Fixed-HRF mismatch and spurious adaptation to colored noise | Optional cross-run shape learning with noise-aware evidence/fallback |
| \#47–#49 | Unchecked coefficient, amplitude, and waveform recovery | Validation of different estimands; fixes to padding, compression, and borrowing |
| \#50–#51 | Expensive dense algebra | Faster core solvers and whitening; larger sensitivity studies become practical |
| \#52 | Weak HRF evidence in signed/dense designs and unreliable anchors | Refined reference tests, stabilized fits, concordance-weighted borrowing |
| \#53 | Scaling compression, ringing, noisy readouts, cross-event leakage | Noise scaling, composite penalties, kernel readouts, operator-based epoch regression |
| \#54 | Duration correlated with predictors | Variable-duration source models distinguishing rate from integrated activity |
| \#55–#62 | Heteroskedastic & correlated trial errors, non-canonical AR($`p`$) noise, and multi-lag/EIV inference | Soft Boltzmann HRF averaging, post-selection 2-pass REML AR($`p`$) + SUSAN noise pooling, spike-robust GLS marginal `noise_sd`, exact subset Toeplitz gapped covariance, parcel pre-aggregation, [`medusa_multilag_infer()`](https://hallquistlab.github.io/fmri.pipeline/reference/medusa_multilag_infer.md), and REPAIR (`M-MOMS`) connectivity/trait disattenuation |
| \#64–#65 | Rapid-ISI covariate mismatch, parcel task-curvature SE inflation, shared-noise PPI bias, and partial-lag/multi-covariate inference | Post-demixing covariate operator (`M_post`), task-free `noise_innovations` + multi-stratum parcel operators, `precision_weight_cohort`, spectrum-continuous ZCA, causal observability gating (`obs_frac_tau >= 0.20`), order-invariant multi-run GCV, joint multi-covariate/Level-2 `medusa_multilag_infer(..., contrast_source = "unsmoothed")`, and `decon_connectivity_repair(..., shared_error = "auto")` |
| \#66 | Cortical-only HRF bounds, cross-undershoot amplitude distortion, and $`O(K^2 + m^3)`$ scaling on expanded libraries | 3D physiological manifold (`"medusa"` / `"manifold"`, ~44 shapes) & maximin $`\epsilon`$-net (`"enet"`) presets, `"canonical_27"`, single-gamma `Inf` undershoot, broad `"positive_lobe"` normalization across deconvolution/GCV/diagnostics, and native C++ cascaded AR($`p`$) whitening + $`O(K)`$ encompassing |

These features are available through **two interfaces**, not one
automatic replacement pipeline. The file route supports newer scaling,
penalties, and readouts. Cross-run HRF learning and operator-based epoch
regression currently require the matrix route. The next vignette shows
both and how to move from NIfTI files to matrices.

## Recommended best-practice configuration (5-stage MEDuSA checklist)

Across the five-pillar Monte Carlo validation suite (evaluating trial
amplitude recovery, within-trial waveform alignment, inter-region
connectivity, sample-size/trial-count power surfaces, and nonlinear
Windkessel-Balloon + `fmrisim` composite scanner realism), the following
five-stage configuration consistently achieves the lowest group slope
RMSE, well-calibrated 95% confidence intervals, and Type I error control
under schedule and duration confounds:

| Stage | Recommended Setting | Scientific Rationale |
|----|----|----|
| **Stage 0: Cross-Run Scaling** | `scaling = "noise_sd"` in [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/fmri.pipeline/reference/deconvolve_crossrun_hrf.md) + `restore_scale = TRUE` in [`decon_epoch_regression()`](https://hallquistlab.github.io/fmri.pipeline/reference/decon_epoch_regression.md) | Normalizes runs by **marginal observation-noise SD** ($`\sigma_{\text{marg}}`$, estimated via leverage-standardized Winsorized GLS residuals rather than AR innovation SD or total signal `"zscore"`), preventing task modulation variance, low-frequency drift, or motion spikes from compressing run amplitudes. |
| **Stage 1: HRF & Noise Estimation** | `hrf_library(tr, preset = "medusa")` + `estimate_crossrun_hrf(..., noise_model = "ar1", evidence = "calibrated", selection_mode = "soft_average", spatial_groups = parcel_ids)` | Learns regional vascular delays, widths, and undershoot depths out-of-sample across runs over a pruned 44-candidate 3D physiological manifold (`"positive_lobe"` normalization preventing cross-undershoot amplitude distortion; use `preset = "canonical_27"` for the fixed-undershoot 27-candidate cortical grid), with spatially-correlated Stouffer/Brown + Cauchy parcel evidence pooling (`hrf_parcel_pooled_evidence()`), local Boltzmann kernel averaging (`selection_mode = "soft_average"`), and post-selection 2-pass REML AR($`p`$) + SUSAN bilateral noise learning with exact gapped C++ cascaded whitening. |
| **Stage 2: Deconvolution** | `deconvolve_crossrun_hrf(..., lambda = NULL, penalty = "sobolev", backend = "cpp")` | Combines amplitude, first-difference, and second-difference smoothness penalties under exact stationary/gapped AR($`p`$)-GLS prewhitening. By default (`lambda = NULL`), `hrf_tune_decon_lambda()` pools across all usable runs in an order-invariant manner to select $`\lambda`$ via gain-calibrated, task-weighted adaptive prewhitened GCV over `lambda_grid`. |
| **Stage 3: Epoch De-Mixing, Parcel Pre-Aggregation & Duration Modeling** | `decon_epoch_regression(..., demix_other_events = TRUE, demix_target_trials = TRUE, restore_scale = TRUE, aggregate_by_parcel = TRUE, tau2_mode = "marginal", whiten_panel = TRUE, return_panel = TRUE)` | Subtracts preceding event bleed (`demix_ridge = 0.05`, or `< 0.03` for trial-varying non-target confounds), restores per-trial target gain with matched $`M_{\text{post}}`$ covariate alignment, gates unobservable boundary/gap trials (`obs_frac_tau < 0.20`), and aggregates parcels using inverse-variance + HRF-evidence weights, exact multi-stratum operators, and task-free Step-1 `noise_innovations` spatial covariance (retaining `shared_noise_loading`). Leave `demix_within_trial = FALSE` when feeding [`medusa_multilag_infer()`](https://hallquistlab.github.io/fmri.pipeline/reference/medusa_multilag_infer.md), and set `source_model = "duration"` (`duration_scale = "rate"`) for reaction-time epochs. |
| **Stage 4: Group Multilevel, Trajectory & Connectivity Inference** | `mixed_by(..., weights = "precision_weight_cohort")`, [`medusa_multilag_infer()`](https://hallquistlab.github.io/fmri.pipeline/reference/medusa_multilag_infer.md), and `decon_connectivity_repair(..., use_whitened = FALSE, shared_error = "auto")` | Uses `precision_weight_cohort` to preserve both per-trial and between-participant precision in [`mixed_by()`](https://hallquistlab.github.io/fmri.pipeline/reference/mixed_by.md) (paired with matched `*_whitened` columns for additive behavioral models); [`medusa_multilag_infer()`](https://hallquistlab.github.io/fmri.pipeline/reference/medusa_multilag_infer.md) for joint multi-covariate/Level-2 trajectory de-blurring, Marra–Wood simultaneous bands, and unbiased window contrasts (`contrast_source = "unsmoothed"`); and unwhitened `decon_demixed` + `se_demixed` + `shared_error = "auto"` in [`decon_connectivity_repair()`](https://hallquistlab.github.io/fmri.pipeline/reference/decon_connectivity_repair.md) for M-MOMS disattenuation and cross-parcel shared-noise correction. |

## Decide whether the result is useful

Evaluate at least four outcomes separately: agreement with known
simulated amplitudes, recovery of timing/waveforms, stability of
brain–behavior slopes, and false positives under null effects with
realistic interfering activity. Better correlation need not mean
unbiased amplitude; better BOLD reconstruction need not mean a better
neural waveform. Numerical speedups and simulation recovery percentages
answer different questions.

For human data, compare a fixed canonical HRF with the learned
alternative, inspect fallback and noise estimates, examine sensitivity
to source duration and readout width, and test repeatability across runs
or sessions. Neural latency and vascular latency remain partly
confounded, HRFs may change across runs, and AR(1) may not capture the
actual noise. There is no established universal percentage gain for real
MEDuSA datasets from these PRs.

Finally, trials, voxels, and neighboring lags are not independent
participants. Use the trial panel or participant-level slopes in a
prespecified multilevel/group analysis; handle dependence and
multiplicity across the tested lag/region family. The convenience
statistics returned by epoch regression are not a substitute for that
inferential design.

## Suggested reading order

1.  **Physical motivation:** [Glover, 1999, NeuroImage
    9:416–429](https://doi.org/10.1006/nimg.1998.0419). Read for
    hemodynamic deblurring and the limits of linearity.
2.  **HRF variation:** [Handwerker et al., 2004, NeuroImage
    21:1639–1651](https://pubmed.ncbi.nlm.nih.gov/15050587/). Read
    before treating a canonical response as universal.
3.  **Temporal noise:** [Woolrich et al., 2001, NeuroImage
    14:1370–1386](https://doi.org/10.1006/nimg.2001.0931). Connect
    prewhitening to familiar fMRI GLM practice.
4.  **Regularization selection:** [Golub et al., 1979, Technometrics
    21:215–223](https://doi.org/10.1080/00401706.1979.10489751). Read
    after the inverse-problem section for GCV’s statistical foundation.
5.  **Trial overlap:** [Mumford et al., 2012, NeuroImage
    59:2636–2643](https://pmc.ncbi.nlm.nih.gov/articles/PMC3251697/).
    Compare trial estimation approaches without assuming they implement
    MEDuSA.
6.  **Duration:** [Grinband et al., 2008, NeuroImage
    43:509–520](https://pmc.ncbi.nlm.nih.gov/articles/PMC2654219/). Read
    before deciding whether RT defines duration or an additional
    predictor.
7.  **A related modern workflow:** [Prince et al., 2022, eLife
    11:e77599](https://pmc.ncbi.nlm.nih.gov/articles/PMC9708069/).
    Compare HRF libraries and regularized single-trial estimation with
    the continuous-then-event-readout approach here.
8.  **Selection and inference:** [Kriegeskorte et al., 2009, Nature
    Neuroscience
    12:535–540](https://pmc.ncbi.nlm.nih.gov/articles/PMC2841687/). Read
    before choosing settings or regions using the effect you intend to
    test.
