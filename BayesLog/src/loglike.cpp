#include <Rcpp.h>
#include "functions.h"

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (NumericMatrix pars, IntegerVector indpars, NumericMatrix data)
{
    //'pars' is a vector of parameters (beta0, beta, sigma2)
    //'indpars' is vector of indicators
    //'data' is matrix of data with response in first column
    
    int i, j;
    double LL = 0, nu;
    
    for(i = 0; i < data.nrow(); i++)
    {
        //initialise linear component
        nu = pars(0, 0);
        for(j = 1; j < data.ncol(); j++)
        {
            //add contribution for each covariate
            nu += pars(0, j) * indpars[j] * data(i, j);
        }
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data(i, 0) == 0 ? (1.0 - nu):nu);
        LL += log(nu);
    }
    return LL;
}
