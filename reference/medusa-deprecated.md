# Deprecated MEDuSA deconvolution, event-alignment, and multilevel functions

All MEDuSA deconvolution, event-locked interpolation, and split-wise
multilevel modeling functions have moved to the standalone \`medusa\`
package. Calling these legacy entry points from \`fmri.pipeline\` emits
a deprecation warning and forwards to \`medusa\` when installed.

## Usage

``` r
voxelwise_deconvolution(...)

deconvolve_nlreg(...)

deconvolve_nlreg_resample(...)

deconvolve_reglin(...)

event_lock_ts(...)

get_medusa_interpolated_ts(...)

get_medusa_compression_score(...)

compress_mts_pca(...)

evt_align_decon_files(...)

run_decon_alignment(...)

mixed_by(...)

meta_mixed_by(...)

jlmer(...)

fmri_ts
```

## Format

An object providing \`\$new(...)\` that emits a deprecation warning and
delegates to \`medusa::fmri_ts\$new(...)\`.

## Arguments

- ...:

  Arguments forwarded directly to the corresponding function in
  \`medusa\`.

## Value

The return value of the corresponding function in \`medusa\`.
