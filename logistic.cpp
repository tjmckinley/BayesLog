#include <Rcpp.h>
using namespace Rcpp;

//function to scale proposal variance
double adapt_scale(int nacc, int niter, double desacc, double propscale)
{
    double propscale1 = propscale;
    double accrate = (double) nacc;
    accrate = accrate / ((double) niter);

    if(niter > 0)
    {
        if(accrate <= desacc) propscale1 = propscale1 / (2.0 - (accrate / desacc));
        else propscale1 = propscale1 * (2.0 - ((1.0 - accrate) / (1.0 - desacc)));
    }
    else propscale1 = propscale;
    return propscale1;
}

// function for calculating the log-likelihood
double loglike (double *pars, int *indpars, NumericMatrix data)
{
    //'pars' is a vector of parameters (beta0, beta, sigma2)
    //'indpars' is vector of indicators
    //'data' is matrix of data with response in first column
    
    int i, j;
    double LL = 0, nu;
    
    for(i = 0; i < data.nrow(); i++)
    {
        //initialise linear component
        nu = pars[0];
        for(j = 1; j < data.ncol(); j++)
        {
            //add contribution for each covariate
            nu += pars[j] * indpars[j - 1] * data(i, j);
        }
//        if(R_finite(nu) == 0) Rprintf("nu %d = %f\n", j, nu);
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
//        if(R_finite(nu) == 0) Rprintf("nu2 %d = %f\n", j, nu);
        //calculate log-likelihood contribution
        nu = (data(i, 0) == 0 ? (1.0 - nu):nu);
//        if(R_finite(nu) == 0) Rprintf("nu3 %d = %f\n", j, nu);
        LL += log(nu);
//        if(R_finite(LL) == 0)
//        {
//            Rprintf("nu4 %d = %f Y = %f\n", j, nu, data(i, 0));
//            getchar();
//        }
    }
    return LL;
}

// a Metropolis-Hastings algorithm for fitting the logistic variable selection model

// [[Rcpp::export]]
NumericMatrix logisticMH (NumericMatrix data, NumericVector ini_pars, int gen_inits, NumericMatrix priors, int niter, double scale)
{
    // 'data' is a matrix of data with the first column equal to the response variable
    // 'ini_pars' is a vector of initial values for the unknown parameters
    // 'gen_inits' is an indicator controlling whether initial values need to be generated
    // 'priors' is an (npars x 2) matrix containing the mean and var for Normal priors
    // 'niter' is the number of iterations over which to run the chain
    // 'scale' is mixing proportion for adaptive MCMC
    
    //initialise indexing variables
    int i, j, k;
    
    // calculate number of parameters
    int nregpars = ini_pars.size() - 2;
    int npars = nregpars + 2;
    int *indpars = (int *) malloc (nregpars * sizeof(int));
    
    Rprintf("\nNumber of iterations = %d\n", niter);
    Rprintf("Scale for adaptive proposal = %f\n", scale);
    Rprintf("Number of regression parameters = %d\n", nregpars);
    
    // set up output vector of length 'niter' to record chain
    // (append extra column for unnormalised posterior)
    NumericMatrix output(niter, nregpars + npars + 1);
    
    // initialise chain and set up vector to hold proposals
    double *pars = (double *) malloc (npars * sizeof(double));
    double *pars_prop = (double *) malloc (npars * sizeof(double));
    for(i = 0; i < npars; i++)
    {
        pars[i] = 0;
        pars_prop[i] = 0;
    }
    
    //declare variables
    double LL_curr, LL_prop, acc_curr, acc_prop, acc;
    
    //generate or read in initial values
    k = 0;
    int ind = 0;
    while(k < 100 && ind == 0)
    {
        //generate initial values if required
        if(gen_inits == 1)
        {
            for(i = 0; i < (npars - 1); i++) pars[i] = rnorm(1, 0.0, 1.0)[0];
            pars[npars - 1] = runif(1, 0.0, 20.0)[0];
        }
        else
        {
            for(i = 0; i < npars; i++) pars[i] = ini_pars[i];
        }
            
        for(i = 0; i < nregpars; i++) indpars[i] = 1;//(int) (runif(1, 0, 1)[0] > 0.5 ? 1:0);
        
        if(pars[npars - 1] < 0.0)
        {
            free(pars); pars = NULL;
            free(pars_prop); pars_prop = NULL;
            free(indpars); indpars = NULL;
            stop("\nVariance hyperparameter not positive");
        }
        
        //check the initial values produce a finite log-posterior
        LL_curr = loglike(pars, indpars, data);
        // calculate log-likelihood – log-prior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr -= 0.5 * pow((pars[0] - priors(0, 0)) / sqrt(priors(0, 1)), 2.0);
        //add priors for regression parameters
        for(j = 0; j < nregpars; j++)
        {
            if(indpars[j] == 1)
            {
                acc_curr -= log(sqrt(pars[npars - 1]));
                acc_curr -= 0.5 * pow((pars[j + 1] - priors(j + 1, 0)) / sqrt(pars[npars - 1]), 2.0);
            }
        }
        //variance component is uniform and so is sucked into normalising constant
        if(R_finite(acc_curr) != 0) ind = 1;
        else
        {
            if(gen_inits == 0) k = 99;
        }
        k++;
    }
    if(k == 100 && R_finite(acc_curr) == 0)
    {    
        free(pars); pars = NULL;
        free(pars_prop); pars_prop = NULL;
        free(indpars); indpars = NULL;
        stop("\nInitial values produce non-finite log-likelihood");
    }
    
    Rprintf("\nInitial values:\n");
    for(j = 0; j < npars; j++) Rprintf("pars[%d] = %f\n", j, pars[j]);
    Rprintf("\nInitial indicators:\n");
    for(j = 0; j < nregpars; j++) Rprintf("indpars[%d] = %d\n", j, indpars[j]);
    Rprintf("\n");
    
    // set up adaptive proposal distribution
    double adaptscale = pow(2.38, 2.0);
    double *tempsd = (double *) malloc (npars * sizeof(double));
    for(j = 0; j < npars; j++) tempsd[j] = 0.1;
    
    //set up vectors to record acceptance rates
    int *nacc = (int *) malloc(npars * sizeof(int));
    for(i = 0; i < npars; i++) nacc[i] = 0;
    int minnacc, maxnacc;
    
    // run chain
    for(i = 0; i < niter; i++)
    {
        for(j = 0; j < npars; j++) pars_prop[j] = pars[j];
        
        // propose new value for intercept
        if((i + 1) <= 100) pars_prop[0] = rnorm(1, pars[0], 0.1)[0];
        else
        {
            if(runif(1, 0.0, 1.0)[0] < scale) pars_prop[0] = rnorm(1, pars[0], 0.1)[0];
            else pars_prop[0] = rnorm(1, pars[0], tempsd[0])[0];
        }
//        Rprintf("pars[%d] = %f pars_prop[%d] = %f\n", 0, pars[0], 0, pars_prop[0]);
//        Rprintf("priors = %f priors = %f\n", priors(0, 0), priors(0, 1));
        // update log-likelihood
        LL_prop = loglike(pars_prop, indpars, data);
//        Rprintf("LL_curr = %f LL_prop = %f\n", LL_curr, LL_prop);
        // calculate log-likelihood – log-prior
        acc_prop = LL_prop;
//        Rprintf("acc_prop = %f\n", acc_prop);
        acc_prop -= (1.0 / (2.0 * sqrt(priors(0, 1)))) * pow(pars_prop[0] - priors(0, 0), 2.0);
//        Rprintf("acc_prop = %f\n", acc_prop);
        acc_curr = LL_curr;
        acc_curr -= (1.0 / (2.0 * sqrt(priors(0, 1)))) * pow(pars[0] - priors(0, 0), 2.0);
        //accept/reject proposal
        acc = acc_prop - acc_curr;
//        Rprintf("acc_curr = %f acc_prop = %f\n", acc_curr, acc_prop);
//        Rprintf("acc = %f acc0 = %f\n", acc, exp(acc));
        if(R_finite(acc) != 0)
        {
            if(log(runif(1, 0.0, 1.0)[0]) < acc)
            {
                pars[0] = pars_prop[0];
                nacc[0]++;
                LL_curr = LL_prop;
            }
            else pars_prop[0] = pars[0];
        }
        else pars_prop[0] = pars[0];
        
        //now propose moves for remaining regression terms
        for(j = 0; j < nregpars; j++)
        {
            if(indpars[j] == 1)
            {
                // propose new value for parameter
                if((i + 1) <= 100) pars_prop[j + 1] = rnorm(1, pars[j + 1], 0.1)[0];
                else
                {
                    if(runif(1, 0.0, 1.0)[0] < scale) pars_prop[j + 1] = rnorm(1, pars[j + 1], 0.1)[0];
                    else pars_prop[j + 1] = rnorm(1, pars[j + 1], tempsd[j + 1])[0];
                }
                // calculate log-likelihood – log-prior
                LL_prop = loglike(pars_prop, indpars, data);
                acc_prop = LL_prop;
                acc_prop -= log(sqrt(pars_prop[npars - 1]));
                acc_prop -= 0.5 * pow((pars_prop[j + 1] - priors(j + 1, 0)) / sqrt(pars_prop[npars - 1]), 2.0);
                acc_curr = LL_curr;
                acc_curr -= log(sqrt(pars[npars - 1]));
                acc_curr -= 0.5 * pow((pars[j + 1] - priors(j + 1, 0)) / sqrt(pars[npars - 1]), 2.0);
                //accept/reject proposal
                acc = acc_prop - acc_curr;
//                Rprintf("acc%d = %f\n", j + 1, exp(acc));
                if(R_finite(acc) != 0)
                {
                    if(log(runif(1, 0.0, 1.0)[0]) < acc)
                    {
                        pars[j + 1] = pars_prop[j + 1];
                        nacc[j + 1]++;
                        LL_curr = LL_prop;
                    }
                    else pars_prop[j + 1] = pars[j + 1];
                }
                else pars_prop[j + 1] = pars[j + 1];
            }
        }
        
        //now propose update for variance hyperparameter
        j = npars - 1;
        if((i + 1) <= 100) pars_prop[j] = rnorm(1, pars[j], 0.1)[0];
        else
        {
            if(runif(1, 0.0, 1.0)[0] < scale) pars_prop[j] = rnorm(1, pars[j], 0.1)[0];
            else pars_prop[j] = rnorm(1, pars[j], tempsd[j])[0];
        }
//        Rprintf("pars_prop_prop[%d] = %f\n", j, pars_prop[j]);
        //check validity
        if(sqrt(pars_prop[j]) > priors(npars - 1, 0) && sqrt(pars_prop[j]) < priors(npars - 1, 1))
        {
//            Rprintf("pars_prop[%d] = %f\n", j, pars_prop[j]);
            // calculate log-prior
            acc_prop = 0.0;
            acc_curr = 0.0;
            for(k = 0; k < nregpars; k++)
            {
                acc_prop -= log(sqrt(pars[npars - 1]));
                acc_prop -= 0.5 * pow((pars[k + 1] - priors(k + 1, 0)) / sqrt(pars[npars - 1]), 2.0);
                acc_curr -= log(sqrt(pars_prop[npars - 1]));
                acc_curr -= 0.5 * pow((pars[k + 1] - priors(k + 1, 0)) / sqrt(pars_prop[npars - 1]), 2.0);
            }
            //accept/reject proposal
            acc = acc_prop - acc_curr;
//            Rprintf("acc%d = %f\n", j, exp(acc));
            if(R_finite(acc) != 0)
            {
                if(log(runif(1, 0.0, 1.0)[0]) < acc)
                {
                    pars[j] = pars_prop[j];
                    nacc[j]++;
                }
                else pars_prop[j] = pars[j];
            }
            else pars_prop[j] = pars[j];
        }
        else pars_prop[j] = pars[j];
        
        //calculate current unnormalised posterior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr -= 0.5 * pow((pars[0] - priors(0, 0)) / sqrt(priors(0, 1)), 2.0);
        //add priors for regression parameters
        for(j = 0; j < nregpars; j++)
        {
            if(indpars[j] == 1)
            {
                acc_curr -= log(sqrt(pars[npars - 1]));
                acc_curr -= 0.5 * pow((pars[j + 1] - priors(j + 1, 0)) / sqrt(pars[npars - 1]), 2.0);
            }
        }
        //variance component is uniform and so is sucked into normalising constant
        
        if(R_finite(acc_curr) == 0)
        {
            free(pars); pars = NULL;
            free(pars_prop); pars_prop = NULL;
            free(tempsd); tempsd = NULL;
            free(indpars); indpars = NULL;
            free(nacc); nacc = NULL;
            stop("Non-finite posterior produced");
        }
        
        // save current value of chain into output matrix
        for(j = 0; j < npars; j++) output(i, j) = pars[j];
        for(j = 0; j < nregpars; j++) output(i, npars + j) = indpars[j];
        output(i, npars + nregpars) = acc_curr;
        
        // calculations for adaptive proposal
        if((i + 1) % 1000 == 0)
        {
            //update proposal variances
            for(j = 0; j < npars; j++)
            {
                tempsd[j] = adapt_scale(nacc[j], 1000, 0.44, tempsd[j]);
//                Rprintf("tempsd[%d] = %f nacc = %d\n", j, tempsd[j], nacc[j]);
            }
            
            // print some output to screen for book-keeping
            minnacc = nacc[0];
            maxnacc = nacc[0];
            for(j = 1; j < npars; j++)
            {
                minnacc = (minnacc < nacc[j] ? minnacc:nacc[j]);
                maxnacc = (maxnacc > nacc[j] ? maxnacc:nacc[j]);
            }
            Rprintf("i = %d min accrate = %f max accrate = %f\n", i + 1, ((double) minnacc) / ((double) 1000), ((double) maxnacc) / ((double) 1000));
            for(j = 0; j < npars; j++) nacc[j] = 0;
        }
    }
    
    //free memory from the heap
    free(pars); pars = NULL;
    free(pars_prop); pars_prop = NULL;
    free(tempsd); tempsd = NULL;
    free(indpars); indpars = NULL;
    free(nacc); nacc = NULL;
    
    return(output);
}

