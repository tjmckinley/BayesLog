#include "functions.hpp"

// function for generating MVN RVs

// code by Ahmadou Dicko from Rcpp Gallery
// http://gallery.rcpp.org/articles/simulate-multivariate-normal/

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) 
{
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
