#include <Rcpp.h>

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (NumericMatrix pars, IntegerVector indpars, NumericMatrix data);

//function to calculate posterior mean and variance recursively
void calcMeanVar(int i, int ninitial, NumericVector *tempmn, NumericVector *tempvar, IntegerVector *tempcounts, NumericMatrix posterior, int postelement, int parelement, int indelement);
