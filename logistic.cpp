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

//function to calculate posterior mean and variance
void calc_meanvar(int i, int ninitial, int npars, NumericVector *tempmn, NumericVector *tempvar, IntegerVector *tempcounts, IntegerVector *indpars, NumericMatrix posterior)
{
    int j, k, m;

    //update means and variances for adaptive proposal
    if((i + 1) == ninitial)
    {
        //first: means
        for(j = 0; j < npars; j++)
        {
            (*tempmn)[j] = 0;
            for(k = 0; k <= i; k++)
            {
                if((*indpars)[j - 1] == 1)
                {
                    (*tempcounts)[j]++;
                    (*tempmn)[j] += posterior(k, j);
                }
            }
            if((*tempcounts)[j] > 1) (*tempmn)[j] = (*tempmn)[j] / ((double) (*tempcounts)[j]);
            else (*tempmn)[j] = 0.0;
        }
        //second: variances
        for(j = 0; j < npars; j++)
        {
            (*tempvar)[j] = 0;
            for(k = 0; k <= i; k++)
            {
                if((*indpars)[j - 1] == 1) (*tempvar)[j] += pow(posterior(k, j), 2.0);
            }
            if((*tempcounts)[j] > 1)
            {
                (*tempvar)[j] = (*tempvar)[j] - (*tempcounts)[j] * pow((*tempmn)[j], 2.0);
                (*tempvar)[j] = (*tempvar)[j] / ((double) ((*tempcounts)[j] - 1));
            }
            else (*tempvar)[j] = 1.0;
        }
    }
    else
    {
        //start recursively updating variance 
        for(j = 0; j < npars; j++)
        {
            if((*indpars)[j - 1] == 1)
            {
                (*tempcounts)[j]++;
                if((*tempcounts)[j] > 2) (*tempvar)[j] = ((*tempvar)[j] * ((*tempcounts)[j] - 2)) + (((*tempcounts)[j] - 1) * pow((*tempmn)[j], 2.0));
            }
        }
        //recursively update mean and variance
        for(j = 0; j < npars; j++)
        {
            if((*indpars)[j - 1] == 1)
            {
                if((*tempcounts)[j] > 2)
                {
                    (*tempmn)[j] = ((*tempmn)[j] * ((*tempcounts)[j] - 1) + posterior(i, j)) / ((double) (*tempcounts)[j]);
                    (*tempvar)[j] += pow(posterior(i, j), 2.0) - (*tempcounts)[j] * pow((*tempmn)[j], 2.0);
                    (*tempvar)[j] = (*tempvar)[j] / ((double) (*tempcounts)[j] - 1);
                }
                else
                {
                    if((*tempcounts)[j] == 2)
                    {
                        //following should be fine since initialised at zero
                        (*tempmn)[j] = 0;
                        for(k = 0; k <= i; k++) (*tempmn)[j] += posterior(k, j);
                        (*tempmn)[j] = (*tempmn)[j] / ((double) (*tempcounts)[j]);
                        (*tempvar)[j] = 0;
                        for(k = 0; k <= i; k++) (*tempvar)[j] += pow(posterior(k, j), 2.0);
                        (*tempvar)[j] = (*tempvar)[j] - (*tempcounts)[j] * pow((*tempmn)[j], 2.0);
                        (*tempvar)[j] = (*tempvar)[j] / ((double) ((*tempcounts)[j] - 1));
                    }
                    else
                    {
                        (*tempmn)[j] = 0.0;
                        (*tempvar)[j] = 1.0;
                    }
                }
            }
        }
    }
    return;
}

// function for calculating the log-likelihood
double loglike (NumericVector pars, IntegerVector indpars, NumericMatrix data)
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
        //convert to correct scale
        nu = exp(nu) / (1.0 + exp(nu));
        //calculate log-likelihood contribution
        nu = (data(i, 0) == 0 ? (1.0 - nu):nu);
        LL += log(nu);
    }
    return LL;
}

// a Metropolis-Hastings algorithm for fitting the logistic variable selection model

// [[Rcpp::export]]
NumericMatrix logisticMH (NumericMatrix data, IntegerVector factindex, IntegerVector cumfactindex, NumericVector ini_pars, int gen_inits, NumericMatrix priors, int niter, double scale, int orignpars, int varselect, int ninitial, int random)
{
    // 'data' is a matrix of data with the first column equal to the response variable
    // 'factindex' is a vector containing number of levels for each variable
    // 'cumfactindex' is a vector indexing start point for each variable
    //      (accounting for different numbers of levels)
    // 'ini_pars' is a vector of initial values for the unknown parameters
    // 'gen_inits' is an indicator controlling whether initial values need to be generated
    // 'priors' is an (npars x 2) matrix containing the mean and var for Normal priors
    // 'niter' is the number of iterations over which to run the chain
    // 'scale' is mixing proportion for adaptive MCMC
    // 'orignpars' is number of variables (disregarding dummy variables)
    // 'varselect' is an indicator controlling whether variable selection is to be done
    // 'ninitial' is the number of ietrations to run before starting to calculate the posterior
    //  mean and variance for use in proposal steps
    // 'random' is an indicator corresponding to whether a "fixed" or "random" effect required
    
    //initialise indexing variables
    int i, j, k;
    
    // calculate number of parameters
    int nregpars = ini_pars.size() - 2;
    int npars = nregpars + 2;
    IntegerVector indpars(nregpars);
    
    Rprintf("\nNumber of iterations = %d\n", niter);
    Rprintf("Scale for adaptive proposal = %f\n", scale);
    Rprintf("Number of regression parameters = %d\n", nregpars);
    Rprintf("Variable selection (1/0): %d\n", varselect);
    Rprintf("Random effect used (1/0): %d\n", random);
    Rprintf("Number of initial iterations before change in add/rem proposals: %d\n", ninitial);
    
    // set up output vector of length 'niter' to record chain
    // (append extra column for unnormalised posterior)
    NumericMatrix output(niter, nregpars + npars + 1);
    
    // initialise chain and set up vector to hold proposals
    NumericVector pars(npars);
    NumericVector pars_prop(npars);
    
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
        //set variance component if fixed effect
        if(random == 0) pars[npars - 1] = priors(0, 1);
            
        if(varselect == 1)
        {
            for(i = 0; i < nregpars; i++) indpars[i] = (int) (runif(1, 0, 1)[0] > 0.5 ? 1:0);
        }
        else
        {
            for(i = 0; i < nregpars; i++) indpars[i] = 1;
        }
        
        if(pars[npars - 1] < 0.0) stop("\nVariance hyperparameter not positive");
        
        //check the initial values produce a finite log-posterior
        LL_curr = loglike(pars, indpars, data);
        // calculate log-likelihood – log-prior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
        //add priors for regression parameters
        for(j = 0; j < nregpars; j++)
        {
            if(indpars[j] == 1) acc_curr += R::dnorm(pars[j + 1], priors(j + 1, 0), sqrt(pars[npars - 1]), 1);
        }
        //add variance component
        if(random == 1) acc_curr += R::dunif(pars[npars - 1], priors(npars - 1, 0), priors(npars - 1, 1), 1);
        if(R_finite(acc_curr) != 0) ind = 1;
        else
        {
            if(gen_inits == 0) k = 99;
        }
        k++;
    }
    if(k == 100 && R_finite(acc_curr) == 0) stop("\nInitial values produce non-finite log-likelihood");
    
    Rprintf("\nInitial values:\n");
    for(j = 0; j < npars; j++) Rprintf("pars[%d] = %f\n", j, pars[j]);
    Rprintf("\nInitial indicators:\n");
    for(j = 0; j < nregpars; j++) Rprintf("indpars[%d] = %d\n", j, indpars[j]);
    Rprintf("\n");
    
    //set up vectors for recording posterior mean and variances
    NumericVector tempmn(npars);
    NumericVector tempvar(npars, 1.0);
    IntegerVector tempcounts(npars);
    
    // set up adaptive proposal distribution
    double adaptscale = pow(2.38, 2.0);
    NumericVector tempsd(npars, 0.1);
    
    //set up vectors to record acceptance rates
    IntegerVector nacc(npars);
    IntegerVector nattempt(npars);
    IntegerVector nacc_add(npars);
    IntegerVector nattempt_add(npars);
    IntegerVector nacc_del(npars);
    IntegerVector nattempt_del(npars);
    
    double tempacc;
    double minnacc, maxnacc;
    double minnacc_add, maxnacc_add;
    double minnacc_del, maxnacc_del;
    
    //set up proposal probabilities to control variable selection
    //and mixing
    double psamp = (varselect == 1 ? 0.5:1.0);
    
    // run chain
    for(i = 0; i < niter; i++)
    {
        for(j = 0; j < npars; j++) pars_prop[j] = pars[j];
        
        // propose new value for intercept
        if(runif(1, 0.0, 1.0)[0] < scale) pars_prop[0] = rnorm(1, pars[0], 0.1)[0];
        else pars_prop[0] = rnorm(1, pars[0], tempsd[0])[0];
        nattempt[0]++;
        
        // update log-likelihood
        LL_prop = loglike(pars_prop, indpars, data);
        // calculate log-likelihood – log-prior
        acc_prop = LL_prop;
        acc_prop += R::dnorm(pars_prop[0], priors(0, 0), sqrt(priors(0, 1)), 1);
        acc_curr = LL_curr;
        acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
        //proposals cancel
        
        //accept/reject proposal
        acc = acc_prop - acc_curr;
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
        for(k = 0; k < orignpars; k++)
        {
            if(indpars[cumfactindex[k] - 1] == 1)
            {
                //if variable is included, then propose to remove or move
                if(runif(1, 0.0, 1.0)[0] < psamp)
                {
                    //MOVEMENT
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        // propose new value for parameter
                        if(runif(1, 0.0, 1.0)[0] < scale) pars_prop[j] = rnorm(1, pars[j], 0.1)[0];
                        else pars_prop[j] = rnorm(1, pars[j], tempsd[j])[0];
                        nattempt[j]++;
                    }
                    // calculate log-likelihood – log-prior
                    LL_prop = loglike(pars_prop, indpars, data);
                    acc_prop = LL_prop;
                    acc_curr = LL_curr;
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc_prop += R::dnorm(pars_prop[j], priors(j, 0), sqrt(pars[npars - 1]), 1);
                        acc_curr += R::dnorm(pars[j], priors(j, 0), sqrt(pars[npars - 1]), 1);
                    }
                    acc = acc_prop - acc_curr;
                    //proposals cancel
                    
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(runif(1, 0.0, 1.0)[0]) < acc)
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars[j] = pars_prop[j];
                                nacc[j]++;
                            }
                            LL_curr = LL_prop;
                        }
                        else
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) pars_prop[j] = pars[j];
                        }
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) pars_prop[j] = pars[j];
                    }
                }
                else
                {
                    //REMOVAL
                
                    //SET new parameter values
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        pars_prop[j] = 0.0;
                        indpars[j - 1] = 0;
                        nattempt_del[j]++;
                    }
                    
                    // calculate log-likelihood – log-prior
                    LL_prop = loglike(pars_prop, indpars, data);
                    acc_prop = LL_prop;
                    acc_curr = LL_curr;
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) acc_curr += R::dnorm(pars[j], priors(j, 0), sqrt(pars[npars - 1]), 1);
                    acc = acc_prop - acc_curr;
                    
                    //adjust for proposals
                    acc -= log(1.0 - psamp);
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) acc += R::dnorm(pars[j], tempmn[j], sqrt(tempvar[j]), 1);
                    
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(runif(1, 0.0, 1.0)[0]) < acc)
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars[j] = pars_prop[j];
                                nacc_del[j]++;
                            }
                            LL_curr = LL_prop;
                        }
                        else
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars_prop[j] = pars[j];
                                indpars[j - 1] = 1;
                            }
                        }
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars_prop[j] = pars[j];
                            indpars[j - 1] = 1;
                        }
                    }
                }
            }
            else
            {
                //ADDITION
                
                //simulate new parameter values
                for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                {
                    pars_prop[j] = rnorm(1, tempmn[j], sqrt(tempvar[j]))[0];
                    indpars[j - 1] = 1;
                    nattempt_add[j]++;
                }
                
                // calculate log-likelihood – log-prior
                LL_prop = loglike(pars_prop, indpars, data);
                acc_prop = LL_prop;
                acc_curr = LL_curr;
                for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) acc_prop += R::dnorm(pars_prop[j], priors(j, 0), sqrt(pars[npars - 1]), 1);
                acc = acc_prop - acc_curr;
                
                //adjust for proposals
                acc += log(1.0 - psamp);
                for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) acc -= R::dnorm(pars_prop[j], tempmn[j], sqrt(tempvar[j]), 1);
                
                //accept/reject proposal
                if(R_finite(acc) != 0)
                {
                    if(log(runif(1, 0.0, 1.0)[0]) < acc)
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars[j] = pars_prop[j];
                            nacc_add[j]++;
                        }
                        LL_curr = LL_prop;
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars_prop[j] = pars[j];
                            indpars[j - 1] = 0;
                        }
                    }
                }
                else
                {
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        pars_prop[j] = pars[j];
                        indpars[j - 1] = 0;
                    }
                }
            }
        }
        
        if(random == 1)
        {
            //now propose update for variance hyperparameter
            j = npars - 1;
            if(runif(1, 0.0, 1.0)[0] < scale) pars_prop[j] = rnorm(1, pars[j], 0.1)[0];
            else pars_prop[j] = rnorm(1, pars[j], tempsd[j])[0];
            nattempt[j]++;
            
            //check validity
            if(sqrt(pars_prop[j]) > priors(npars - 1, 0) && sqrt(pars_prop[j]) < priors(npars - 1, 1))
            {
                // calculate log-prior
                acc_prop = 0.0;
                acc_curr = 0.0;
                for(k = 0; k < nregpars; k++)
                {
                    if(indpars[k] == 1)
                    {
                        acc_prop += R::dnorm(pars[k + 1], priors(k + 1, 0), sqrt(pars_prop[npars - 1]), 1);
                        acc_curr += R::dnorm(pars[k + 1], priors(k + 1, 0), sqrt(pars[npars - 1]), 1);
                    }
                }
                //proposals cancel
                
                //accept/reject proposal
                acc = acc_prop - acc_curr;
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
        }
        
        //calculate current unnormalised posterior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
        //add priors for regression parameters
        for(j = 0; j < nregpars; j++)
        {
            if(indpars[j] == 1) acc_curr += R::dnorm(pars[j + 1], priors(j + 1, 0), sqrt(pars[npars - 1]), 1);
        }
        //add variance component
        if(random == 1) acc_curr += R::dunif(pars[npars - 1], priors(npars - 1, 0), priors(npars - 1, 1), 1);
        
        if(R_finite(acc_curr) == 0) stop("Non-finite posterior produced");
        
        // save current value of chain into output matrix
        for(j = 0; j < npars; j++) output(i, j) = pars[j];
        for(j = 0; j < nregpars; j++) output(i, npars + j) = indpars[j];
        output(i, npars + nregpars) = acc_curr;
        
        // calculations for adaptive proposal
        if((i + 1) % 1000 == 0)
        {
            //update proposal variances
            for(j = 0; j < (random == 1 ? npars:(npars - 1)); j++) tempsd[j] = adapt_scale(nacc[j], nattempt[j], 0.44, tempsd[j]);
            
            // print some output to screen for book-keeping
            minnacc = (nattempt[0] > 0 ? (((double) nacc[0]) / ((double) nattempt[0])):0.0);
            maxnacc = minnacc;
            for(j = 1; j < (random == 1 ? npars:(npars - 1)); j++)
            {
                tempacc = ((double) nacc[j]) / ((double) nattempt[j]);
                if(R_finite(tempacc) != 0)
                {
                    minnacc = (minnacc < tempacc ? minnacc:tempacc);
                    maxnacc = (maxnacc > tempacc ? maxnacc:tempacc);
                }
            }
            if(varselect == 1)
            {
                minnacc_add = (nattempt_add[1] > 0 ? (((double) nacc_add[1]) / ((double) nattempt_add[1])):0.0);
                maxnacc_add = minnacc_add;
                for(j = 1; j < nregpars; j++)
                {
                    tempacc = ((double) nacc_add[j + 1]) / ((double) nattempt_add[j + 1]);
                    if(R_finite(tempacc) != 0)
                    {
                        minnacc_add = (minnacc_add < tempacc ? minnacc_add:tempacc);
                        maxnacc_add = (maxnacc_add > tempacc ? maxnacc_add:tempacc);
                    }
                }
                minnacc_del = (nattempt_del[1] > 0 ? (((double) nacc_del[1]) / ((double) nattempt_del[1])):0.0);
                maxnacc_del = minnacc_del;
                for(j = 1; j < nregpars; j++)
                {
                    tempacc = ((double) nacc_del[j + 1]) / ((double) nattempt_del[j + 1]);
                    if(R_finite(tempacc) != 0)
                    {
                        minnacc_del = (minnacc_del < tempacc ? minnacc_del:tempacc);
                        maxnacc_del = (maxnacc_del > tempacc ? maxnacc_del:tempacc);
                    }
                }
                Rprintf("i = %d minmove = %f maxmove = %f minadd = %f maxadd = %f mindel = %f maxdel = %f\n", i + 1, minnacc, maxnacc, minnacc_add, maxnacc_add, minnacc_del, maxnacc_del);
            }
            else
            {
                for(j = 0; j < nregpars; j++)
                {
                    if(nattempt_add[j] > 0 || nattempt_del[j] > 0) stop("Incorrect addition/removal move");
                }
                Rprintf("i = %d minmove = %f maxmove = %f\n", i + 1, minnacc, maxnacc);
            }
                
            for(j = 0; j < npars; j++)
            {
                nacc[j] = 0;
                nattempt[j] = 0;
                nacc_add[j] = 0;
                nattempt_add[j] = 0;
                nacc_del[j] = 0;
                nattempt_del[j] = 0;
            }
        }
        
        // record posterior mean and variances for use in addition/removal steps
        if ((i + 1) >= ninitial) calc_meanvar(i, ninitial, npars, &tempmn, &tempvar, &tempcounts, &indpars, output);
    }
    
    return(output);
}

