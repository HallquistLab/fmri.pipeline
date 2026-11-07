# Noise-aware HRF validation artifacts

See [the implementation and validation report](../../../NOISE_AWARE_HRF_RECOVERY.md).

- 480 independent datasets: 12 scenarios, 40 fresh seed indices (4--43) each.
  Development-pilot indices 1--3 are excluded.
- Five methods and three fixed lambdas; lambda 0.01 is primary. No winning
  setting is selected from the recovery results.
- `summary.csv`: all means, Monte Carlo intervals, and finite dataset counts.
- `paired_vs_matched_canonical.csv`: changes versus canonical GLS with identical
  estimated noise, acquisition mask, scaling, and inverse nuisance policy.
- `recovery_comparison.png`: primary-lambda paired changes and Monte Carlo intervals.
- `provenance.txt`: script hashes, seeds, and inference unit.
- `benchmark.csv`: final native implementation, 50,000 voxels, three runs,
  300 acquisitions, nine HRFs, 56 event amplitudes, medians of three repetitions.
- `benchmark_before_projection_optimization.csv`: same inputs before the safe
  projection-norm optimization; legacy OLS has no REML calculation to optimize.

R 4.5.1 / OpenBLAS 0.3.30, one BLAS/OpenMP thread. Raw multirun/dataset metrics,
frozen settings, and environments remain in `/tmp/fmri-hrf-noise/reported/`.
Benchmarks are in `/tmp/fmri-hrf-noise/benchmark-optimized/` and `benchmark/`.
The optimized installed binary reran one dataset in each scenario, matching all
180 method/lambda metric summaries exactly (maximum absolute difference zero).

Source MD5 checksums for the final numerical implementation:

- `R/crossrun_hrf.R`: `d998cf391a84503306b92ea24471e9aa`
- `R/hrf_noise_evidence.R`: `d65390438449f2dec8c4eac5efce2f06`
- `R/diagnose_crossrun_recovery.R`: `68c8dcb102610ab6f99d5740b601c43b`
- `src/crossrun_hrf.cpp`: `80b2264b6499213fbeab2b5aad83b824`
