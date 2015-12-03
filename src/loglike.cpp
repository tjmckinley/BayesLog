#include "functions.hpp"

// function for calculating the log-likelihood
double loglike (arma::vec pars, arma::mat data, arma::vec nsamples, arma::vec logL)
{
    //'pars' is a vector of regression parameters
    //'data' is matrix of data with response in first column
    //'nsamples' is vector of counts for each row of 'data'
    //'logL' is a vector to record log-likelihood components
    
    int i, j;
    double LL = 0.0, nu;
    
    int nrow = data.n_rows;
    int ncol = data.n_cols;
    
    #pragma omp parallel for private(i, nu, j) schedule(static)
    for(i = 0; i < nrow; i++)
    {
        //initialise linear component
        nu = 0.0;
        for(j = 0; j < (ncol - 1); j++)
        {
            //add contribution for each covariate
            nu += pars(j) * data(i, j + 1);
        }
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data(i, 0) == 0 ? (1.0 - nu):nu);
        logL(i) = nsamples(i) * log(nu);
    }
    #pragma omp parallel for reduction (+:LL)
    for(i = 0; i < nrow; i++) LL += logL(i);
    return LL;
}
