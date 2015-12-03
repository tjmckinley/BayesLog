#include "functions.hpp"

// a Metropolis-Hastings algorithm for fitting a logistic regression model

arma::mat logisticMH (arma::mat data, arma::vec nsamples, arma::vec ini_pars, arma::mat priors, int niter, double scale, int nadapt, int nprintsum, double maxscale, double niterdim)
{
    // 'data' is a matrix of data with the first column equal to the response variable
    // 'nsamples' corresponds to aggregated counts for each row of 'data'
    // 'nrandint' corresponds to number of random intercept terms
    // 'ini_pars' is a vector of initial values for the unknown parameters
    // 'priors' is an (npars x 2) matrix containing the mean and var for Normal priors
    // 'niter' is the number of iterations over which to run the chain
    // 'scale' is mixing proportion for adaptive MCMC
    // 'nadapt' is the number of iterations to run between updates of the adaptive proposal
    // 'nprintsum' controls how often run time information is printed to the screen
    // 'maxscale' is the maximum scaling of the adaptive proposal variance at each update
    // 'niterdim' is the iteration at which the diminishing adaptation component kicks in
    
    //initialise indexing variables
    int i, j, k;
    
    //set number of threads
    int ncores = omp_get_num_procs();
    omp_set_num_threads(ncores);
    
    // calculate number of parameters
    int npars = ini_pars.size();
    
    //set terms for parallelising code
    arma::vec logL(data.n_rows);
    
    //print runtime information to the screen    
    i = 0;
    for(j = 0; j < nsamples.size(); j++) i += nsamples[j];
    Rprintf("\nNumber of samples in data set = %d\n", i);
    Rprintf("Number of unique samples in data set = %d\n", data.n_rows);
    Rprintf("Run time information printed to screen every %d iterations\n", nprintsum);
    
    Rprintf("\nNumber of iterations = %d\n", niter);
    Rprintf("Scale for adaptive proposal = %f\n", scale);
    Rprintf("Number of regression parameters = %d\n", npars);
    Rprintf("Adapt every %d iterations\n", nadapt);
    Rprintf("Max scale for adapting = %f\n", maxscale);
    Rprintf("%d iterations before diminishing adaptation kicks in\n", (int) niterdim);
    
    //print prior information to screen
    Rprintf("\nPriors: mean = %f variance = %f\n", priors(0, 0), priors(0, 1));
    
    // set up output vector of length 'niter' to record chain
    // (append extra column for unnormalised posterior)
    int noutput = npars + 1;
    arma::mat output(niter, noutput);
    output.zeros();
    
    // initialise chain and set up vector to hold proposals
    arma::vec pars(npars);
    arma::vec pars_prop(npars);
    pars.zeros(); pars_prop.zeros();
    
    //declare variables
    double LL_curr, LL_prop, acc_curr, acc_prop, acc;
    
    //set initial values
    for(i = 0; i < npars; i++) pars(i) = ini_pars(i);
    
    //check the initial values produce a finite log-posterior
    LL_curr = loglike(pars, data, nsamples, logL);
    // calculate log-likelihood – log-prior
    acc_curr = LL_curr;
    //add priors for regression parameters
    for(j = 0; j < npars; j++) acc_curr += R::dnorm(pars(j), priors(j, 0), sqrt(priors(j, 1)), 1);
    if(R_finite(acc_curr) == 0) Rcpp::stop("Initial values are not viable");
    
    //print initial values to screen
    Rprintf("\nInitial values:\n");
    for(j = 0; j < npars; j++) Rprintf("pars[%d] = %f\n", j, pars(j));
    Rprintf("\n");
    
    //set up vectors for adaptive proposals
    arma::vec propsd(npars);
    propsd.fill(0.1 * 0.1);
    
    //set up vectors to record acceptance rates
    arma::ivec nacc(npars);
    arma::ivec nacc1(npars);
    nacc.zeros(); nacc1.zeros();
    
    arma::ivec nattempt(npars);
    arma::ivec nattempt1(npars);
    nattempt.zeros(); nattempt1.zeros();
    
    double tempacc, minnacc, maxnacc;
    
    // run chain
    Rprintf("Starting run:\n");
    for(i = 0; i < niter; i++)
    {  
        //check for user interruptions
        if (i % 10 == 0) R_CheckUserInterrupt();
        
        //reset proposals 
        for(k = 0; k < npars; k++) pars_prop(k) = pars(k);
        
        //now propose updates for regression terms
        for(k = 0; k < npars; k++)
        {
            // propose new parameter value
            if(R::runif(0.0, 1.0) < scale) pars_prop(k) = R::rnorm(pars(k), 0.1);
            else pars_prop(k) = R::rnorm(pars(k), propsd(k));
            nattempt(k)++;
            
            // calculate log-likelihood
            LL_prop = loglike(pars_prop, data, nsamples, logL);
            acc_prop = LL_prop;
            acc_curr = LL_curr;
            //adjust for prior
            acc_curr += R::dnorm(pars(k), priors(k, 0), sqrt(priors(k, 1)), 1);
            acc_prop += R::dnorm(pars_prop(k), priors(k, 0), sqrt(priors(k, 1)), 1);
            acc = acc_prop - acc_curr;
            //proposals cancel
            
            //accept/reject proposal
            if(R_finite(acc) != 0)
            {
                if(log(R::runif(0.0, 1.0)) < acc)
                {
                    pars(k) = pars_prop(k);
                    nacc(k)++;
                    LL_curr = LL_prop;
                    acc_curr = acc_prop;
                }
                else pars_prop(k) = pars(k);
            }
            else pars_prop(k) = pars(k);
        }
        
        // save current value of chain into output matrix
        for(j = 0; j < npars; j++) output(i, j) = pars(j);
        output(i, npars) = acc_curr;
        
        // record posterior mean and variances for use in proposals
        if ((i + 1) % nadapt == 0)
        {
            for(j = 0; j < npars; j++) propsd(j) = adapt_scale(nacc(j) - nacc1(j), nattempt(j) - nattempt1(j), 0.44, propsd(j), (double) i + 1, maxscale, niterdim);
            
            //update counts
            for(j = 0; j < npars; j++)
            {
                nacc1(j) = nacc(j);
                nattempt1(j) = nattempt(j);
            }
        }  
        
        if((i + 1) % nprintsum == 0)
        {          
            // print some output to screen for book-keeping
            minnacc = ((double) nacc(0)) / ((double) nattempt(0));
            maxnacc = minnacc;
            for(j = 1; j < npars; j++)
            {
                tempacc = ((double) nacc(j)) / ((double) nattempt(j));
                minnacc = (minnacc < tempacc ? minnacc:tempacc);
                maxnacc = (maxnacc > tempacc ? maxnacc:tempacc);
            }
            Rprintf("i = %d minacc = %f maxacc = %f\n", i + 1, minnacc, maxnacc);
               
            //reset counters 
            for(j = 0; j < npars; j++)
            {
                nacc(j) = 0; nacc1(j) = 0;
                nattempt(j) = 0; nattempt1(j) = 0;
            }
        }
    }
    
    //return output
    return(output);
}

