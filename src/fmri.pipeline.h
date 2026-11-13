#ifndef _fmripipeline_FMRIPIPELINE_h
#define _fmripipeline_FMRIPIPELINE_h

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::vec convolve_cpp(const arma::vec &a, const arma::vec &b);
arma::vec convolve_double_gamma(const arma::vec &stimulus, double a1, double a2, double b1, double b2, double cc);

#endif
