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
NumericMatrix logisticMH (NumericMatrix data, IntegerVector factindex, IntegerVector cumfactindex, NumericVector ini_pars, NumericVector ini_sigma2, int gen_inits, NumericMatrix priors, int niter, double scale, int orignpars, int varselect, int ninitial, int random)
{
    // 'data' is a matrix of data with the first column equal to the response variable
    // 'factindex' is a vector containing number of levels for each variable
    // 'cumfactindex' is a vector indexing start point for each variable
    //      (accounting for different numbers of levels)
    // 'ini_pars' is a vector of initial values for the unknown parameters
    // 'ini_sigma2' is a vector of initial values for the unknown random effects components
    // 'gen_inits' is an indicator controlling whether initial values need to be generated
    // 'priors' is an (npars x 2) matrix containing the mean and var for Normal priors
    // 'niter' is the number of iterations over which to run the chain
    // 'scale' is mixing proportion for adaptive MCMC
    // 'orignpars' is number of variables (disregarding dummy variables)
    // 'varselect' is an indicator controlling whether variable selection is to be done
    // 'ninitial' is the number of ietrations to run before starting to calculate the posterior
    //  mean and variance for use in proposal steps
    // 'random' is an indicator corresponding to whether a "fixed" (0), 
    //      global "random" (1) or local "random" (2) effect required
    
    //initialise indexing variables
    int i, j, k;
    
    // calculate number of parameters
    int npars = ini_pars.size();
    int nregpars = npars - 1;
    IntegerVector indpars(nregpars);
    
    Rprintf("\nNumber of iterations = %d\n", niter);
    Rprintf("Scale for adaptive proposal = %f\n", scale);
    Rprintf("Number of regression parameters (excluding intercept) = %d\n", nregpars);
    Rprintf("Variable selection (1/0): %d\n", varselect);
    Rprintf("Number of initial iterations before change in add/rem proposals: %d\n", ninitial);
    
    if(random == 0) Rprintf("\nFIXED variance components used\n");
    if(random == 1) Rprintf("\nGLOBAL RANDOM variance components used\n");
    if(random == 2) Rprintf("\nLOCAL RANDOM variance components used\n");
    
    // set up output vector of length 'niter' to record chain
    // (append extra column for unnormalised posterior)
    int noutput;
    if(random == 0) noutput = npars + nregpars + 1;
    if(random == 1) noutput = npars + nregpars + 2;
    if(random == 2) noutput = npars + 2 * nregpars + 1;
    NumericMatrix output(niter, noutput);
    
    // initialise chain and set up vector to hold proposals
    NumericVector pars(npars);
    NumericVector pars_prop(npars);
    NumericVector sigma2(nregpars);
    NumericVector sigma2_prop(nregpars);
    
    //declare variables
    double LL_curr, LL_prop, acc_curr, acc_prop, acc, sigma2full, sigma2full_prop;
    
    //generate or read in initial values
    k = 0;
    int ind = 0;
    while(k < 100 && ind == 0)
    {
        //generate initial values if required
        if(gen_inits == 1)
        {
            for(i = 0; i < npars; i++) pars[i] = rnorm(1, 0.0, 1.0)[0];
            if(random == 1)
            {
                sigma2full = runif(1, 0.0, 20.0)[0];
                sigma2full = pow(sigma2full, 2.0);
                for(i = 0; i < nregpars; i++) sigma2[i] = sigma2full;
            }
            if(random == 2)
            {
                for(i = 0; i < nregpars; i++)
                {
                    sigma2[i] = runif(1, 0.0, 20.0)[0];
                    sigma2[i] = pow(sigma2[i], 2.0);
                }
            }
        }
        else
        {
            for(i = 0; i < npars; i++) pars[i] = ini_pars[i];
            for(i = 0; i < nregpars; i++) sigma2[i] = ini_sigma2[i];
            if(random == 1)
            {
                sigma2full = sigma2[0];
                //reset just in case error in input
                for(i = 0; i < nregpars; i++) sigma2[i] = sigma2full;
            }   
        }
        //set variance component if fixed effect
        if(random == 0)
        {
            for(i = 0; i < nregpars; i++) sigma2[i] = priors(0, 1);
        }
            
        if(varselect == 1)
        {
            for(i = 0; i < nregpars; i++) indpars[i] = (int) (runif(1, 0, 1)[0] > 0.5 ? 1:0);
        }
        else
        {
            for(i = 0; i < nregpars; i++) indpars[i] = 1;
        }
        
        for(i = 0; i < nregpars; i++) if(sigma2[i] < 0.0) stop("\nVariance hyperparameter %d not positive\n", i);
        
        //check the initial values produce a finite log-posterior
        LL_curr = loglike(pars, indpars, data);
        // calculate log-likelihood – log-prior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
        //add priors for regression parameters
        for(j = 0; j < nregpars; j++)
        {
            if(indpars[j] == 1) acc_curr += R::dnorm(pars[j + 1], priors(j + 1, 0), sqrt(sigma2[j]), 1);
        }
        //add variance component
        if(random > 0)
        {
            if(random == 2)
            {
                for(j = 0; j < nregpars; j++) if(indpars[j] == 1) acc_curr += R::dunif(sqrt(sigma2[j]), priors(npars, 0), priors(npars, 1), 1);
            }
            else acc_curr += R::dunif(sqrt(sigma2full), priors(npars, 0), priors(npars, 1), 1);
        }
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
    if(random > 0)
    {
        if(random == 2)
        {
            Rprintf("\nInitial variance components:\n");
            for(j = 0; j < nregpars; j++) Rprintf("sigma2[%d] = %f\n", j, sigma2[j]);
        }
        else Rprintf("\nInitial variance component: %f\n", sigma2full);
    }
    else Rprintf("\nFixed variance component: %f\n", priors(0, 1));
    Rprintf("\n");
    
    //set up vectors for recording posterior mean and variances
    int nadaptpars;
    if(random == 0) nadaptpars = npars;
    if(random == 1) nadaptpars = npars + 1;
    if(random == 2) nadaptpars = npars + nregpars;
    NumericVector tempmn(nadaptpars);
    NumericVector tempvar(nadaptpars, 1.0);
    IntegerVector tempcounts(nadaptpars);
    
    // set up adaptive proposal distribution
    double adaptscale = pow(2.38, 2.0);
    NumericVector tempsd(npars, 0.1);
    NumericVector tempsigmasd(nregpars, 0.1);
    double tempsigmafullsd = 0.1;
    
    //set up vectors to record acceptance rates
    IntegerVector nacc(npars);
    IntegerVector nattempt(npars);
    IntegerVector naccsigma(nregpars);
    IntegerVector nattemptsigma(nregpars);
    IntegerVector nacc_add(npars);
    IntegerVector nattempt_add(npars);
    IntegerVector nacc_del(npars);
    IntegerVector nattempt_del(npars);
    int naccsigmafull = 0, nattemptsigmafull = 0;
    
    double tempacc;
    double minnacc, maxnacc;
    double minnaccsigma, maxnaccsigma;
    double minnacc_add, maxnacc_add;
    double minnacc_del, maxnacc_del;
    
    //set up proposal probabilities to control variable selection
    //and mixing
    double psamp = (varselect == 1 ? 0.5:1.0);
    
    // run chain
    for(i = 0; i < niter; i++)
    {
        for(j = 0; j < npars; j++) pars_prop[j] = pars[j];
        for(j = 0; j < nregpars; j++) sigma2_prop[j] = sigma2[j];
        
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
                        acc_curr += R::dnorm(pars[j], priors(j, 0), sqrt(sigma2[j - 1]), 1);
                        acc_prop += R::dnorm(pars_prop[j], priors(j, 0), sqrt(sigma2[j - 1]), 1);
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
                        sigma2_prop[j - 1] = 0.0;
                        indpars[j - 1] = 0;
                        nattempt_del[j]++;
                    }
                    
                    // calculate log-likelihood – log-prior
                    LL_prop = loglike(pars_prop, indpars, data);
                    acc_prop = LL_prop;
                    acc_curr = LL_curr;
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc_curr += R::dnorm(pars[j], priors(j, 0), sqrt(sigma2[j - 1]), 1);
                        if(random == 2) acc_curr += R::dunif(sqrt(sigma2[j - 1]), priors(npars, 0), priors(npars, 1), 1);
                    }
                    acc = acc_prop - acc_curr;
                    
                    //adjust for proposals
                    acc -= log(1.0 - psamp);
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc += R::dnorm(pars[j], tempmn[j], sqrt(tempvar[j]), 1);
                        if(random == 2) acc += R::dnorm(sqrt(sigma2[j - 1]), tempmn[npars + j], sqrt(tempvar[npars + j]), 1);
                    }
                    
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(runif(1, 0.0, 1.0)[0]) < acc)
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars[j] = pars_prop[j];
                                sigma2[j - 1] = sigma2_prop[j - 1];
                                nacc_del[j]++;
                            }
                            LL_curr = LL_prop;
                        }
                        else
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars_prop[j] = pars[j];
                                sigma2_prop[j - 1] = sigma2[j - 1];
                                indpars[j - 1] = 1;
                            }
                        }
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars_prop[j] = pars[j];
                            sigma2_prop[j - 1] = sigma2[j - 1];
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
                    if(random == 2)
                    {
                        sigma2_prop[j - 1] = rnorm(1, tempmn[npars + j], sqrt(tempvar[npars + j]))[0];
                        sigma2_prop[j - 1] = pow(sigma2_prop[j - 1], 2.0);
                    }
                    if(random == 1) sigma2_prop[j - 1] = sigma2full;
                    indpars[j - 1] = 1;
                    nattempt_add[j]++;
                }
                
                // calculate log-likelihood – log-prior
                LL_prop = loglike(pars_prop, indpars, data);
                acc_prop = LL_prop;
                acc_curr = LL_curr;
                for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                {
                    acc_prop += R::dnorm(pars_prop[j], priors(j, 0), sqrt(sigma2_prop[j - 1]), 1);
                    if(random == 2) acc_prop += R::dunif(sqrt(sigma2_prop[j - 1]), priors(npars, 0), priors(npars, 1), 1);
                }
                acc = acc_prop - acc_curr;
                
                //adjust for proposals
                acc += log(1.0 - psamp);
                for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                {
                    acc -= R::dnorm(pars_prop[j], tempmn[j], sqrt(tempvar[j]), 1);
                    if(random == 2) acc -= R::dnorm(sqrt(sigma2_prop[j - 1]), tempmn[npars + j], sqrt(tempvar[npars + j]), 1);
                }
                
                //accept/reject proposal
                if(R_finite(acc) != 0)
                {
                    if(log(runif(1, 0.0, 1.0)[0]) < acc)
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars[j] = pars_prop[j];
                            sigma2[j - 1] = sigma2_prop[j - 1];
                            nacc_add[j]++;
                        }
                        LL_curr = LL_prop;
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars_prop[j] = pars[j];
                            sigma2_prop[j - 1] = sigma2[j - 1];
                            indpars[j - 1] = 0;
                        }
                    }
                }
                else
                {
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        pars_prop[j] = pars[j];
                        sigma2_prop[j - 1] = sigma2[j - 1];
                        indpars[j - 1] = 0;
                    }
                }
            }
        }
        
        if(random > 0)
        {
            if(random == 1)
            {
                //now propose update for SD hyperparameter
                if(runif(1, 0.0, 1.0)[0] < scale) sigma2full_prop = rnorm(1, sqrt(sigma2full), 0.1)[0];
                else sigma2full_prop = rnorm(1, sqrt(sigma2full), tempsigmafullsd)[0];
                nattemptsigmafull++;
                
                //check validity of SD hyperparameter
                if(sigma2full_prop > priors(npars, 0) && sigma2full_prop < priors(npars, 1))
                {
                    //convert to variance for consistency
                    sigma2full_prop = pow(sigma2full_prop, 2.0);
                    
                    // calculate log-prior for regression terms
                    acc_prop = 0.0;
                    acc_curr = 0.0;
                    for(k = 0; k < nregpars; k++)
                    {
                        if(indpars[k] == 1)
                        {
                            acc_prop += R::dnorm(pars[k + 1], priors(k + 1, 0), sqrt(sigma2full_prop), 1);
                            acc_curr += R::dnorm(pars[k + 1], priors(k + 1, 0), sqrt(sigma2full), 1);
                        }
                    }
                    //priors for SD cancel
                    //proposals cancel
                    
                    //accept/reject proposal
                    acc = acc_prop - acc_curr;
                    if(R_finite(acc) != 0)
                    {
                        if(log(runif(1, 0.0, 1.0)[0]) < acc)
                        {
                            sigma2full = sigma2full_prop;
                            naccsigmafull++;
                        }
                        else sigma2full_prop = sigma2full;
                    }
                    else sigma2full_prop = sigma2full;
                }
                else sigma2full_prop = sigma2full;
                //reset individual sigma2
                for(j = 0; j < nregpars; j++) sigma2[j] = sigma2full;
            }
            else
            {
                for(j = 0; j < nregpars; j++)
                {
                    if(indpars[j] == 1)
                    {
                        //now propose update for SD hyperparameter
                        if(runif(1, 0.0, 1.0)[0] < scale) sigma2_prop[j] = rnorm(1, sqrt(sigma2[j]), 0.1)[0];
                        else sigma2_prop[j] = rnorm(1, sqrt(sigma2[j]), tempsigmasd[j])[0];
                        nattemptsigma[j]++;
                        
                        //check validity of SD hyperparameter
                        if(sigma2_prop[j] > priors(npars, 0) && sigma2_prop[j] < priors(npars, 1))
                        {
                            //convert to variance for consistency
                            sigma2_prop[j] = pow(sigma2_prop[j], 2.0);
                            
                            // calculate log-prior for regression terms
                            acc_prop += R::dnorm(pars[j + 1], priors(j + 1, 0), sqrt(sigma2_prop[j]), 1);
                            acc_curr += R::dnorm(pars[j + 1], priors(j + 1, 0), sqrt(sigma2[j]), 1);
                            //priors for SD cancel
                            //proposals cancel
                            
                            //accept/reject proposal
                            acc = acc_prop - acc_curr;
                            if(R_finite(acc) != 0)
                            {
                                if(log(runif(1, 0.0, 1.0)[0]) < acc)
                                {
                                    sigma2[j] = sigma2_prop[j];
                                    naccsigma[j]++;
                                }
                                else sigma2_prop[j] = sigma2[j];
                            }
                            else sigma2_prop[j] = sigma2[j];
                        }
                        else sigma2_prop[j] = sigma2[j];
                    }
                }
            }
        }
        
        //calculate current unnormalised posterior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
        //add priors for regression parameters
        for(j = 0; j < nregpars; j++)
        {
            if(indpars[j] == 1) acc_curr += R::dnorm(pars[j + 1], priors(j + 1, 0), sqrt(sigma2[j]), 1);
        }
        //add variance component
        if(random > 0)
        {
            if(random == 2)
            {
                for(j = 0; j < nregpars; j++)
                {
                    if(indpars[j] == 1) acc_curr += R::dunif(sqrt(sigma2[j]), priors(npars, 0), priors(npars, 1), 1);
                }
            }
            else acc_curr += R::dunif(sqrt(sigma2full), priors(npars, 0), priors(npars, 1), 1);
        }
        if(R_finite(acc_curr) == 0) stop("Non-finite posterior produced");
        
        // save current value of chain into output matrix
        for(j = 0; j < npars; j++) output(i, j) = pars[j];
        if(random == 1)
        {
            output(i, npars) = sigma2full;
            for(j = 0; j < nregpars; j++) output(i, npars + 1 + j) = indpars[j];
            output(i, npars + nregpars + 1) = acc_curr;
        }
        if(random == 2)
        {
            for(j = 0; j < nregpars; j++) output(i, npars + j) = sigma2[j];
            for(j = 0; j < nregpars; j++) output(i, npars + nregpars + j) = indpars[j];
            output(i, npars + 2 * nregpars) = acc_curr;
        }
        
        // calculations for adaptive proposal
        if((i + 1) % 1000 == 0)
        {
            //update proposal variances
            for(j = 0; j < npars; j++) tempsd[j] = adapt_scale(nacc[j], nattempt[j], 0.44, tempsd[j]);
            if(random == 1) tempsigmafullsd = adapt_scale(naccsigmafull, nattemptsigmafull, 0.44, tempsigmafullsd);
            if(random == 2)
            {
                for(j = 0; j < nregpars; j++) tempsigmasd[j] = adapt_scale(naccsigma[j], nattemptsigma[j], 0.44, tempsigmasd[j]);
            }
            
            // print some output to screen for book-keeping
            minnacc = (nattempt[0] > 0 ? (((double) nacc[0]) / ((double) nattempt[0])):0.0);
            maxnacc = minnacc;
            for(j = 1; j < npars; j++)
            {
                tempacc = ((double) nacc[j]) / ((double) nattempt[j]);
                if(R_finite(tempacc) != 0)
                {
                    minnacc = (minnacc < tempacc ? minnacc:tempacc);
                    maxnacc = (maxnacc > tempacc ? maxnacc:tempacc);
                }
            }
            if(random == 2)
            {
                // print some output to screen for book-keeping
                minnaccsigma = (nattemptsigma[0] > 0 ? (((double) naccsigma[0]) / ((double) nattemptsigma[0])):0.0);
                maxnaccsigma = minnaccsigma;
                for(j = 1; j < nregpars; j++)
                {
                    tempacc = ((double) naccsigma[j]) / ((double) nattemptsigma[j]);
                    if(R_finite(tempacc) != 0)
                    {
                        minnaccsigma = (minnaccsigma < tempacc ? minnaccsigma:tempacc);
                        maxnaccsigma = (maxnaccsigma > tempacc ? maxnaccsigma:tempacc);
                    }
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
                if(random == 0) Rprintf("i = %d minmove = %f maxmove = %f minadd = %f maxadd = %f mindel = %f maxdel = %f\n", i + 1, minnacc, maxnacc, minnacc_add, maxnacc_add, minnacc_del, maxnacc_del);
                else
                {
                    if(random == 1) Rprintf("i = %d minmove = %f maxmove = %f sigma = %f minadd = %f maxadd = %f mindel = %f maxdel = %f\n", i + 1, minnacc, maxnacc, naccsigmafull / ((double) nattemptsigmafull), minnacc_add, maxnacc_add, minnacc_del, maxnacc_del);
                    else Rprintf("i = %d minmove = %f maxmove = %f minsigma = %f maxsigma = %f minadd = %f maxadd = %f mindel = %f maxdel = %f\n", i + 1, minnacc, maxnacc, minnaccsigma, maxnaccsigma, minnacc_add, maxnacc_add, minnacc_del, maxnacc_del);
                }
            }
            else
            {
                for(j = 0; j < nregpars; j++)
                {
                    if(nattempt_add[j] > 0 || nattempt_del[j] > 0) stop("Incorrect addition/removal move");
                }
                if(random == 0) Rprintf("i = %d minmove = %f maxmove = %f\n", i + 1, minnacc, maxnacc);
                else
                {
                    if(random == 1) Rprintf("i = %d minmove = %f maxmove = %f sigma = %f\n", i + 1, minnacc, maxnacc, naccsigmafull / ((double) nattemptsigmafull));
                    else Rprintf("i = %d minmove = %f maxmove = %f minsigma = %f maxsigma = %f\n", i + 1, minnacc, maxnacc, minnaccsigma, maxnaccsigma);
                }
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
            for(j = 0; j < nregpars; j++)
            {
                naccsigma[j] = 0;
                nattemptsigma[j] = 0;
            }
            naccsigmafull = 0;
            nattemptsigmafull = 0;
        }
        
        // record posterior mean and variances for use in addition/removal steps
        if ((i + 1) >= ninitial) calc_meanvar(i, ninitial, nadaptpars, &tempmn, &tempvar, &tempcounts, &indpars, output);
    }
    
    return(output);
}

