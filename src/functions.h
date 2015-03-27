#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (arma::mat pars, arma::vec indpars, arma::mat data, arma::vec nsamples, arma::ivec randint, arma::vec rand);

// function for calculating the log-likelihood based on change in a single random intercept term
double loglike_randint (arma::mat pars, arma::vec indpars, arma::mat data, arma::vec nsamples, arma::ivec randint, arma::vec rand, arma::vec rand_prop, arma::ivec cumrandindex, int nrand);

//function to calculate posterior mean and variance recursively
void calcMeanVar(int i, int ninitial, NumericVector *tempmn, NumericVector *tempvar, IntegerVector *tempcounts, NumericMatrix posterior, int postelement, int parelement, int indelement);

// [[Rcpp::export]]
List logisticMH (arma::mat data, arma::vec nsamples, int nrandint, arma::ivec randint, arma::ivec cumrandindex, IntegerVector factindex, IntegerVector cumfactindex, NumericVector ini_pars, NumericVector ini_sigma, double ini_sigmarand, int gen_inits, NumericMatrix priors, int niter, int nitertraining, double scale, int orignpars, int varselect, int ninitial, int random, int nprintsum);
