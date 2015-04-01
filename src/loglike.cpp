#include "functions.h"

// function for calculating the log-likelihood
double loglike (arma::mat pars, arma::vec indpars, arma::mat data, arma::vec nsamples, arma::ivec randint, arma::vec rand, arma::vec logL)
{
    //'pars' is a vector of parameters (beta0, beta, sigma2)
    //'indpars' is vector of indicators
    //'data' is matrix of data with response in first column
    //'nsamples' is vector of counts for each row of 'data'
    //'randint' is random intercept indictator
    //'rand' is vector of random intercept parameters
    //'logL' is a vector to record log-likelihood components
    
    int i, j;
    double LL = 0.0, nu;
    
    int nrow = data.n_rows;
    int ncol = data.n_cols;
    
    #pragma omp parallel for private(i, nu, j) schedule(static)
    for(i = 0; i < nrow; i++)
    {
        //initialise linear component
        nu = pars(0, 0);
        //add random intercept term
        nu += rand[randint[i]];
        for(j = 1; j < ncol; j++)
        {
            //add contribution for each covariate
            nu += pars(0, j) * indpars[j] * data(i, j);
        }
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data(i, 0) == 0 ? (1.0 - nu):nu);
        logL[i] = nsamples[i] * log(nu);
    }
    #pragma omp parallel for reduction (+:LL)
    for(i = 0; i < nrow; i++) LL += logL[i];
    return LL;
}

// function for calculating the log-likelihood based on change in a single random intercept term
double loglike_randint (arma::mat pars, arma::vec indpars, arma::mat data, arma::vec nsamples, arma::ivec randint, arma::vec rand, arma::vec rand_prop, arma::ivec cumrandindex, int nrand)
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
        for(j = 1; j < data.n_cols; j++)
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
