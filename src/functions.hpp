#ifndef __FUNCTIONS__
#define __FUNCTIONS__

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

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
void adapt_update(int i, int ninitial, int niter, int npars, arma::vec *tempmn, arma::mat *meanmat, arma::mat *meanmat1, NumericMatrix posterior, arma::mat *propcov, int subrow, arma::ivec *elements);

// [[Rcpp::export]]
NumericMatrix logisticMH (NumericMatrix dataR, NumericVector nsamplesR, NumericVector ini_pars, NumericMatrix priors, int niter, int ninitial, double scale, int nadapt, int nprintsum, double maxscale, double niterdim, int nrand, List randindexesL, IntegerMatrix data_randR, IntegerVector nblock, List blockR, int printini, List noncentreintR, NumericVector noncentreintRE);

//function for producing sens, spec, ppv and npv from posterior predictions
// [[Rcpp::export]]
List classification (NumericMatrix pred, IntegerVector obs, NumericVector thresh);

#endif // __FUNCTIONS__
