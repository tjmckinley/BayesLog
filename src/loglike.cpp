#include "functions.hpp"

// function for calculating the log-likelihood
double loglike (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, double *logL)
{
    //'pars' is a vector of regression parameters
    //'nrow' and 'ncol' are sizes of data set
    //'data' is array of data with response in first column
    //'nsamples' is vector of counts for each row of 'data'
    //'nrand' is number of random effect terms
    //'rand' is matrix of random effect parameters
    //'data_rand' is matrix matching individuals to correct random effect terms
    //'logL' is a vector to record log-likelihood components
    
    int i, j;
    double LL = 0.0, nu;
    
    #pragma omp parallel for private(i, nu, j) schedule(static)
    for(i = 0; i < nrow; i++)
    {
        //initialise linear component
        nu = 0.0;
        for(j = 0; j < (ncol - 1); j++)
        {
            //add contribution for each covariate
            nu += pars[j] * data[i][j + 1];
        }
        //add contribution from random effect terms
        if(nrand > 0)
        {
            for(j = 0; j < nrand; j++) nu += rand[j][data_rand[i][j]];
        } 
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data[i][0] == 0 ? (1.0 - nu):nu);
        logL[i] = nsamples[i] * log(nu);
    }
    #pragma omp parallel for reduction (+:LL)
    for(i = 0; i < nrow; i++) LL += logL[i];
    return LL;
}

// function for calculating a subset of the log-likelihood
double loglike_sub (double *pars, int nrow, int ncol, double **data, double *nsamples, int nrand, double **rand, int **data_rand, int ***randindexes, int **nrandindexes, int randi, int randj, double *logL)
{
    //'pars' is a vector of regression parameters
    //'data' is matrix of data with response in first column
    //'nsamples' is vector of counts for each row of 'data'
    //'nrand' is number of random effect terms
    //'rand' is matrix of random effect parameters
    //'data_rand' is matrix matching individuals to correct random effect terms
    //'randindexes' is multidimensional arry of indexes relating
    //  individuals to specific random effects
    //'nrandindexes' contains number of indexes for each effect
    //'randi' and 'randj' denote specific element to pull from indexing arrays
    //'logL' is a vector to record log-likelihood components
    
    int i, j, k;
    double LL = 0.0, nu;
    
    #pragma omp parallel for private(i, nu, j, k) schedule(static)
    for(k = 0; k < nrow; k++)
    {
        i = randindexes[randi][randj][k];
        //initialise linear component
        nu = 0.0;
        for(j = 0; j < (ncol - 1); j++)
        {
            //add contribution for each covariate
            nu += pars[j] * data[i][j + 1];
        }
        //add contribution from random effect terms
        for(j = 0; j < nrand; j++) nu += rand[j][data_rand[i][j]];
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data[i][0] == 0 ? (1.0 - nu):nu);
        logL[k] = nsamples[i] * log(nu);
    }
    #pragma omp parallel for reduction (+:LL)
    for(i = 0; i < nrow; i++) LL += logL[i];
    return LL;
}
