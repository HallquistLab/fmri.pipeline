# Coefficient recovery results

See [the study report](../../../COEFFICIENT_RECOVERY_VALIDATION.md) for methods,
findings, assumptions, verification and reproduction.

The study compares onset, window-sum and calibrated-template readouts with
ordinary and transformed-design regression. A direct forward GLS coefficient
fit is a reference, without independent trial-amplitude estimates. Readout and
inverse choices are frozen before evaluation.

`primary.csv` and `primary_intervals.csv` give every primary-lambda result.
`paired_primary.csv` retains within-dataset comparisons. Complete summaries at
all lambdas are in `summary.csv.gz` and `paired.csv.gz`; open with
`read.csv(gzfile(path))`. All 32,760 independent-dataset/method/lambda rows are
retained in `replicate_metrics.rds` for paired reanalysis.

The 840 evaluation datasets are independent of development and the preceding
study. `settings.rds` stores the frozen protocols for four disjoint scenario
groups. `fit_diagnostics.csv` records dataset seeds and HRF error;
`provenance.txt` records runtime/source information, and `checksums.md5`
fingerprints the generated artifacts.

Figures show coefficient estimates, direct-effect null rejection and lambda
sensitivity. Error bars are Monte Carlo uncertainty across independent datasets.
The t intervals audited for coverage are intervals from each simulated dataset.

Transformed regression changes coefficients, not reconstructed neural or
trial-amplitude values. It does not independently improve their correlations.
