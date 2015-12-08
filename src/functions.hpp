#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, arma::vec logL);

// function for calculating a subset of the log-likelihood
double loglike_sub (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, int ***randindexes, int **nrandindexes, int randi, int randj, arma::vec logL);

//function to scale proposal variance
double adapt_scale(int nacc, int niter, double desacc, double propscale, int totiter, double maxscale, double niterdim);

// [[Rcpp::export]]
arma::mat logisticMH (arma::mat dataR, arma::vec nsamplesR, arma::vec ini_pars, arma::mat priors, int niter, double scale, int nadapt, int nprintsum, double maxscale, double niterdim, int nrand, List randindexesL, arma::imat data_randR);

#endif // __FUNCTIONS__
