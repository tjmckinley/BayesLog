#include <Rcpp.h>
#include "functions.h"

using namespace Rcpp;

// function for calculating the log-likelihood
double loglike (NumericMatrix pars, IntegerVector indpars, NumericMatrix data, IntegerVector nsamples, IntegerVector randint, NumericVector rand)
{
    //'pars' is a vector of parameters (beta0, beta, sigma2)
    //'indpars' is vector of indicators
    //'data' is matrix of data with response in first column
    //'nsamples' is vector of counts for each row of 'data'
    //'randint' is random intercept indictator
    //'rand' is vector of random intercept parameters
    
    int i, j;
    double LL = 0, nu;
    
    for(i = 0; i < data.nrow(); i++)
    {
        //initialise linear component
        nu = pars(0, 0);
        //add random intercept term
        nu += rand[randint[i]];
        for(j = 1; j < data.ncol(); j++)
        {
            //add contribution for each covariate
            nu += pars(0, j) * indpars[j] * data(i, j);
        }
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data(i, 0) == 0 ? (1.0 - nu):nu);
        LL += nsamples[i] * log(nu);
    }
    return LL;
}

// function for calculating the log-likelihood based on change in a single random intercept term
double loglike_randint (NumericMatrix pars, IntegerVector indpars, NumericMatrix data, IntegerVector nsamples, IntegerVector randint, NumericVector rand, NumericVector rand_prop, IntegerVector cumrandindex, int nrand)
{
    //'pars' is a vector of parameters (beta0, beta, sigma2)
    //'indpars' is vector of indicators
    //'data' is matrix of data with response in first column
    //'nsamples' is vector of counts for each row of 'data'
    //'randint' is random intercept indictator
    //'rand' is vector of random intercept parameters
    //'rand_prop' is vector of proposed random intercept parameters
    //'cumrandindex' is a vector of cumulative indexes
    //'nrand' denotes the random intercept term over which to calculate the 
    //  change in the log-likelihood
    
    int i, j;
    double LL = 0, nu, nuP;
    
    for(i = cumrandindex[nrand]; i < cumrandindex[nrand + 1]; i++)
    {
        //initialise linear component
        nu = pars(0, 0);
        for(j = 1; j < data.ncol(); j++)
        {
            //add contribution for each covariate
            nu += pars(0, j) * indpars[j] * data(i, j);
        }
        nuP = nu;
        //add random intercept term
        nu += rand[randint[i]];
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data(i, 0) == 0 ? (1.0 - nu):nu);
        LL -= nsamples[i] * log(nu);
        
        //add random intercept term
        nuP += rand_prop[randint[i]];
        //convert to correct scale
        nuP = exp(nuP) / (1.0 + exp(nuP));
        //calculate log-likelihood contribution
        nuP = (data(i, 0) == 0 ? (1.0 - nuP):nuP);
        LL += nsamples[i] * log(nuP);
    }
    return LL;
}
