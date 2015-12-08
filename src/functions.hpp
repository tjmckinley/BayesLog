#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <Rcpp.h>

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, NumericVector logL);

// function for calculating a subset of the log-likelihood
double loglike_sub (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, int ***randindexes, int **nrandindexes, int randi, int randj, NumericVector logL);

//function to scale proposal variance
double adapt_scale(int nacc, int niter, double desacc, double propscale, int totiter, double maxscale, double niterdim);

// [[Rcpp::export]]
NumericMatrix logisticMH (NumericMatrix dataR, NumericVector nsamplesR, NumericVector ini_pars, NumericMatrix priors, int niter, double scale, int nadapt, int nprintsum, double maxscale, double niterdim, int nrand, List randindexesL, IntegerMatrix data_randR);

#endif // __FUNCTIONS__
