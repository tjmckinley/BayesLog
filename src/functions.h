#include <RcppArmadillo.h>
#include <omp.h>

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (arma::mat pars, arma::vec indpars, arma::mat data, arma::vec nsamples, arma::ivec randint, arma::vec rand, arma::vec logL);

// function for calculating the log-likelihood based on change in a single random intercept term
double loglike_randint (arma::mat pars, arma::vec indpars, arma::mat data, arma::vec nsamples, arma::ivec randint, arma::vec rand, arma::vec rand_prop, arma::ivec cumrandindex, int nrand);

//function to calculate posterior mean and variance recursively
void calcMeanVar(int i, int ninitial, arma::vec *tempmn, arma::vec *tempvar, arma::ivec *tempcounts, arma::mat *posterior, int postelement, int parelement, int indelement);

//function to scale proposal variance
double adapt_scale(int nacc, int niter, double desacc, double propscale, int totiter, double maxscale, double niterdim);

// [[Rcpp::export]]
List logisticMH (arma::mat data, arma::vec nsamples, int nrandint, arma::ivec randint, arma::ivec cumrandindex, IntegerVector factindex, IntegerVector cumfactindex, arma::vec ini_pars, arma::vec ini_sigma, double ini_sigmarand, int gen_inits, arma::mat priors, int niter, int nitertraining, double scale, int orignpars, int varselect, int ninitial, int random, int nprintsum, double maxscale, double niterdim);
