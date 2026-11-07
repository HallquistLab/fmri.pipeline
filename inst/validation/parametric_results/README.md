# Parametric recovery study results

See [`PARAMETRIC_RECOVERY_VALIDATION.md`](../../../PARAMETRIC_RECOVERY_VALIDATION.md)
for the generating model, estimands, methods, interpretation and reproduction.

`primary.csv` gives every method/scenario mean at the fixed primary lambda 0.01.
`primary_intervals.csv` includes Monte Carlo uncertainty across independent
datasets. `paired_primary.csv` subtracts the second diagnostic reference from
the first, as specified in the report.

`summary.csv.gz` and `paired.csv.gz` retain all outcomes at all eight lambdas;
use `read.csv(gzfile(path))`. `replicate_metrics.rds` retains the 48,000 individual
dataset/method/lambda rows for paired reanalysis. `fit_diagnostics.csv` includes
dataset seeds and HRF/covariance error. `settings.rds` contains frozen protocols
for the three disjoint scenario groups; `provenance.txt` records their runtime
and source hashes. Evaluation and development seed blocks are disjoint.

The figures show correlations under genuine trial variation, regularization
curves for diagnostic references, and downstream null calibration. Error bars
in the null figure are Monte Carlo intervals across datasets, rather than
coefficient confidence intervals from a single simulated participant.

The report generator verifies complete method/lambda/dataset coverage, exact
seed provenance, compatible protocols and source hashes, and disjoint scenarios
before assembling results. No evaluation truth chooses a method or lambda.
