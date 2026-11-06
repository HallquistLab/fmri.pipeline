# This function deconvolves the BOLD signal using Bush 2011 method, augmented by the resampling approach of Bush 2015.

This function deconvolves the BOLD signal using Bush 2011 method,
augmented by the resampling approach of Bush 2015.

## Usage

``` r
deconvolve_nlreg_resample(
  bold_obs,
  kernel,
  nev_lr = 0.01,
  epsilon = 0.005,
  beta = 40,
  n_resample = 30,
  trim_kernel = TRUE,
  normalize = TRUE
)
```

## Arguments

- bold_obs:

  observed BOLD signal

- kernel:

  assumed kernel of the BOLD signal

- nev_lr:

  learning rate for the assignment of neural events. Default: .01

- epsilon:

  relative error change (termination condition). Default: .005

- beta:

  slope of the sigmoid transfer function. Default: 40

- n_resample:

  number of resampling steps for deconvolution (at least two). Default:
  30

- trim_kernel:

  whether to remove the first K time points from the deconvolved vector,
  corresponding to kernel leftovers from convolution. Default: TRUE

- normalize:

  whether to z-score the input before both the base fit and residual
  resampling. Default: TRUE. FALSE is for inputs already standardized by
  the caller; Bush resampling is not calibrated for arbitrary input
  units.

## Value

list containing the following fields: - NEVest - the base neural event
estimate - NEVmean - the mean neural event estimate - NEVstd - the std
dev. of the neural event estimate - NEVcupp - the mean (upper limit 95 -
NEVclow - the mean (lower limit 95 - BLDmean - the mean of the BOLD
estimate - BLDstd - the std dev. of the BOLD estimate - BLDcupp - the
mean BOLD (upper limit 95 - BLDclow - the mean BOLD (lower limit 95

## Details

Author: Keith Bush Institution: University of Arkansas at Little Rock
Date: Aug. 12, 2013

## Author

Keith Bush
