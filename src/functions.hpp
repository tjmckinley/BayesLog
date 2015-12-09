#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// function for generating MVN RVs
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);

// function for calculating the log-likelihood
double loglike (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, double *logL);

// function for calculating a subset of the log-likelihood
double loglike_sub (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, int ***randindexes, int **nrandindexes, int randi, int randj, double *logL);

//function to scale proposal variance
double adapt_scale(int nacc, int niter, double desacc, double propscale, int totiter, double maxscale, double niterdim);

//function to calculate means and covariance matrices for adaptive MCMC
void adapt_update(int i, int ninitial, int niter, int npars, double adaptscale, arma::vec *tempmn, arma::mat *meanmat, arma::mat *meanmat1, NumericMatrix posterior, arma::mat *propcov, int subrow, arma::ivec *elements);

// [[Rcpp::export]]
NumericMatrix logisticMH (NumericMatrix dataR, NumericVector nsamplesR, NumericVector ini_pars, NumericMatrix priors, int niter, double scale, int nadapt, int nprintsum, double maxscale, double niterdim, int nrand, List randindexesL, IntegerMatrix data_randR, IntegerVector nblock, List blockR);

#endif // __FUNCTIONS__
