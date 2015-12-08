#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (arma::vec pars, arma::mat data, arma::vec nsamples, int nrand, double **rand, arma::imat data_rand, arma::vec logL);

// function for calculating a subset of the log-likelihood
double loglike_sub (arma::vec pars, arma::mat data, arma::vec nsamples, int nrand, double **rand, arma::imat data_rand, int ***randindexes, int **nrandindexes, int randi, int randj, arma::vec logL);

//function to scale proposal variance
double adapt_scale(int nacc, int niter, double desacc, double propscale, int totiter, double maxscale, double niterdim);

// [[Rcpp::export]]
arma::mat logisticMH (arma::mat data, arma::vec nsamples, arma::vec ini_pars, arma::mat priors, int niter, double scale, int nadapt, int nprintsum, double maxscale, double niterdim, int nrand, List randindexesL, arma::imat data_rand);

#endif // __FUNCTIONS__
