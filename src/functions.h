#include <Rcpp.h>

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (NumericMatrix pars, IntegerVector indpars, NumericMatrix data, IntegerVector nsamples, IntegerVector randint, NumericVector rand);

// function for calculating the log-likelihood based on change in a single random intercept term
double loglike_randint (NumericMatrix pars, IntegerVector indpars, NumericMatrix data, IntegerVector nsamples, IntegerVector randint, NumericVector rand, NumericVector rand_prop, IntegerVector cumrandindex, int nrand);

//function to calculate posterior mean and variance recursively
void calcMeanVar(int i, int ninitial, NumericVector *tempmn, NumericVector *tempvar, IntegerVector *tempcounts, NumericMatrix posterior, int postelement, int parelement, int indelement);

// [[Rcpp::export]]
List logisticMH (NumericMatrix data, IntegerVector nsamples, int nrandint, IntegerVector randint, IntegerVector cumrandindex, IntegerVector factindex, IntegerVector cumfactindex, NumericVector ini_pars, NumericVector ini_sigma, double ini_sigmarand, int gen_inits, NumericMatrix priors, int niter, int nitertraining, double scale, int orignpars, int varselect, int ninitial, int random, int nprintsum);
