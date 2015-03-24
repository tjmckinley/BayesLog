#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);

//function to draw from a multivariate normal
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma)
{
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

//function to calculate multivariate normal density
double mvdnormArma(arma::vec x, arma::vec mu, arma::mat sigma)
{
    int ncols = sigma.n_cols;
    arma::vec x1 = (x - mu);
    arma::mat sigmainv = arma::inv(sigma);
    arma::rowvec xmu = x1.t() * sigmainv;
    double ans = dot(xmu, x1);
    double detsig = arma::det(sigma);
    detsig = (detsig < 0 ? -detsig:detsig);
    return -(ncols / 2.0) * log2pi - 0.5 * detsig - 0.5 * ans;
}

//function to calculate posterior mean and variance
void calc_meanvar(int i, int ninitial, NumericVector *tempmn, NumericVector *tempvar, IntegerVector *tempcounts, NumericMatrix posterior, int postelement, int parelement, int indelement)
{
    int j, k, m;

    //update means and variances for adaptive proposal
    if((i + 1) == ninitial)
    {
        //calculate mean
        for(k = 0; k <= i; k++)
        {
            if(posterior(k, indelement) == 1)
            {
                (*tempcounts)[parelement]++;
                (*tempmn)[parelement] += posterior(k, postelement);
            }
        }
        if((*tempcounts)[parelement] > 1) (*tempmn)[parelement] = (*tempmn)[parelement] / ((double) (*tempcounts)[parelement]);
        else (*tempmn)[parelement] = 0.0;
        
        //calculate variance
        (*tempvar)[parelement] = 0.0;
        for(k = 0; k <= i; k++) if(posterior(k, indelement) == 1) (*tempvar)[parelement] += pow(posterior(k, postelement), 2.0);
        if((*tempcounts)[parelement] > 1)
        {
            (*tempvar)[parelement] -= (*tempcounts)[parelement] * pow((*tempmn)[parelement], 2.0);
            (*tempvar)[parelement] = (*tempvar)[parelement] / ((double) ((*tempcounts)[parelement] - 1.0));
        }
        else (*tempvar)[parelement] = pow(0.1, 2.0);
    }
    else
    {
        //start recursively updating variance
        if(posterior(i, indelement) == 1)
        {
            (*tempcounts)[parelement]++;
            if((*tempcounts)[parelement] > 2) (*tempvar)[parelement] = ((*tempvar)[parelement] * ((*tempcounts)[parelement] - 2.0)) + (((*tempcounts)[parelement] - 1.0) * pow((*tempmn)[parelement], 2.0));
        }
        //now update mean
        if((*tempcounts)[parelement] > 2) (*tempmn)[parelement] = (((*tempcounts)[parelement] - 1.0) * (*tempmn)[parelement] + posterior(i, postelement)) / ((double) (*tempcounts)[parelement]);
        else
        {
            if((*tempcounts)[parelement] == 2)
            {
                //calculate mean
                for(k = 0; k <= i; k++)
                {
                    if(posterior(k, indelement) == 1) (*tempmn)[parelement] += posterior(k, postelement);
                }
                (*tempmn)[parelement] = (*tempmn)[parelement] / ((double) (*tempcounts)[parelement]);
            }
            else (*tempmn)[parelement] = 0.0;
        }
        //now update variance
        if(posterior(i, indelement) == 1)
        {
            if((*tempcounts)[parelement] > 2)
            {
                (*tempvar)[parelement] += pow(posterior(i, postelement), 2.0) - (*tempcounts)[parelement] * pow((*tempmn)[parelement], 2.0);
                (*tempvar)[parelement] = (*tempvar)[parelement] / ((double) (*tempcounts)[parelement] - 1);
            }
            else
            {
                if((*tempcounts)[parelement] == 2)
                {
                    (*tempvar)[parelement] = 0.0;
                    for(k = 0; k <= i; k++) if(posterior(k, indelement) == 1) (*tempvar)[parelement] += pow(posterior(k, postelement), 2.0);
                    (*tempvar)[parelement] = (*tempvar)[parelement] - (*tempcounts)[parelement] * pow((*tempmn)[parelement], 2.0);
                    (*tempvar)[parelement] = (*tempvar)[parelement] / ((double) ((*tempcounts)[parelement] - 1));
                }
                else (*tempvar)[parelement] = pow(0.1, 2.0);
            }
        }
    }
    return;
}

// function for calculating the log-likelihood
double loglike (arma::mat pars, IntegerVector indpars, NumericMatrix data)
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

// a Metropolis-Hastings algorithm for fitting the logistic variable selection model

// [[Rcpp::export]]
NumericMatrix logisticMH (NumericMatrix data, IntegerVector factindex, IntegerVector cumfactindex, NumericVector ini_pars, NumericVector ini_sigma, int gen_inits, NumericMatrix priors, int niter, int nitertraining, double scale, int orignpars, int varselect, int ninitial, int random)
{
    // 'data' is a matrix of data with the first column equal to the response variable
    // 'factindex' is a vector containing number of levels for each variable
    // 'cumfactindex' is a vector indexing start point for each variable
    //      (accounting for different numbers of levels)
    // 'ini_pars' is a vector of initial values for the unknown parameters
    // 'ini_sigma' is a vector of initial values for the unknown random effects components
    // 'gen_inits' is an indicator controlling whether initial values need to be generated
    // 'priors' is an (npars x 2) matrix containing the mean and var for Normal priors
    // 'niter' is the number of iterations over which to run the chain
    // 'nitertraining' is the number of iterations over which to run the chain for the training run
    //      (to generate initial values)
    // 'scale' is mixing proportion for adaptive MCMC
    // 'orignpars' is number of variables (disregarding dummy variables)
    // 'varselect' is an indicator controlling whether variable selection is to be done
    // 'ninitial' is the number of ietrations to run before starting to calculate the posterior
    //  mean and variance for use in proposal steps
    // 'random' is an indicator corresponding to whether a "fixed" (0), 
    //      global "random" (1) or local "random" (2) effect required
    
    //initialise indexing variables
    int i, j, k, m;
    
    // calculate number of parameters
    int npars = ini_pars.size();
    int nregpars = npars - 1;
    IntegerVector indpars(npars);
    
    Rprintf("\nNumber of iterations = %d\n", niter);
    Rprintf("Scale for adaptive proposal = %f\n", scale);
    Rprintf("Number of regression parameters (excluding intercept) = %d\n", nregpars);
    Rprintf("Variable selection (1/0): %d\n", varselect);
    Rprintf("Number of initial iterations before adaptive proposal starts: %d\n", ninitial);
    if(nitertraining > 0) Rprintf("Length of training run = %d\n", nitertraining);
    
    if(random == 0) Rprintf("\nFIXED SD components used\n");
    if(random == 1) Rprintf("\nGLOBAL RANDOM SD component used\n");
    if(random == 2) Rprintf("\nLOCAL RANDOM SD components used\n");
    
    // set up output vector of length 'niter' to record chain
    // (append extra column for unnormalised posterior)
    int noutput;
    if(random == 0) noutput = 2 * npars + 1;
    if(random == 1) noutput = 2 * npars + 2;
    if(random == 2) noutput = 3 * npars + 1;
    NumericMatrix output(niter, noutput);
    
    // initialise chain and set up vector to hold proposals
    arma::mat pars(2, npars);
    pars.zeros();
    arma::mat pars_prop(2, npars);
    pars_prop.zeros();
    
    //declare variables
    double LL_curr, LL_prop, acc_curr, acc_prop, acc, sigmafull, sigmafull_prop;
    
    //generate or read in initial values
    k = 0;
    int ind = 0;
    while(k < 100 && ind == 0)
    {
        //generate initial values if required
        if(gen_inits == 1)
        {
            for(i = 0; i < npars; i++) pars(0, i) = rnorm(1, 0.0, 1.0)[0];
            if(random == 1)
            {
                sigmafull = runif(1, priors(npars, 0), priors(npars, 1))[0];
                for(i = 0; i < nregpars; i++) pars(1, i + 1) = sigmafull;
            }
            if(random == 2)
            {
                for(i = 0; i < npars; i++) pars(1, i) = runif(1, priors(npars, 0), priors(npars, 1))[0];
            }
        }
        else
        {
            for(i = 0; i < npars; i++) pars(0, i) = ini_pars[i];
            if(random == 2)
            {
                for(i = 0; i < npars; i++) pars(1, i) = ini_sigma[i];
            }
            if(random == 1)
            {
                sigmafull = ini_sigma[0];
                for(i = 0; i < nregpars; i++) pars(1, i + 1) = sigmafull;
            }   
        }
        //set variance component if fixed effect
        if(random == 0)
        {
            for(i = 0; i < nregpars; i++) pars(1, i + 1) = sqrt(priors(0, 1));
        }
        //sample initial values for variable selection if required    
        if(varselect == 1 && nitertraining == 0)
        {
            indpars[0] = 1;
            for(i = 1; i < npars; i++) indpars[i] = (int) (runif(1, 0, 1)[0] > 0.5 ? 1:0);
        }
        else
        {
            for(i = 0; i < npars; i++) indpars[i] = 1;
        }
        
        for(i = 0; i < npars; i++) if(pars(1, i) < 0.0) stop("\nSD hyperparameter %d not positive\n", i);
        
        //check the initial values produce a finite log-posterior
        LL_curr = loglike(pars, indpars, data);
        // calculate log-likelihood – log-prior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
        //add priors for regression parameters
        for(j = 0; j < npars; j++)
        {
            if(indpars[j] == 1) acc_curr += R::dnorm(pars(0, j), priors(j, 0), pars(1, j), 1);
        }
        //add variance component
        if(random == 1) acc_curr += R::dunif(sigmafull, priors(npars, 0), priors(npars, 1), 1);
        if(random == 2)
        {
            for(j = 0; j < npars; j++)
            {
                if(indpars[j] == 1) acc_curr += R::dunif(pars(1, j), priors(npars, 0), priors(npars, 1), 1);
            }
        }
        if(R_finite(acc_curr) != 0) ind = 1;
        else
        {
            if(gen_inits == 0) k = 99;
        }
        k++;
    }
    if(k == 100 && R_finite(acc_curr) == 0) stop("\nInitial values produce non-finite log-likelihood");
    
    //print initial values to screen
    Rprintf("\nInitial values:\n");
    for(j = 0; j < npars; j++) Rprintf("pars[%d] = %f\n", j, pars(0, j));
    if(varselect > 0)
    {
        Rprintf("\nInitial indicators:\n");
        for(j = 0; j < nregpars; j++) Rprintf("indpars[%d] = %d\n", j + 1, indpars[j + 1]);
    }
    if(random == 0) Rprintf("\nFixed SD component: %f\n", sqrt(priors(0, 1)));
    if(random == 1) Rprintf("\nInitial SD component: %f\n", sigmafull);
    if(random == 2)
    {
        Rprintf("\nInitial SD components:\n");
        for(j = 0; j < npars; j++) Rprintf("sigma[%d] = %f\n", j, pars(1, j));
    }
    Rprintf("\n");
    
    // set up adaptive proposal distribution
    double adaptscale_sing = pow(2.38, 2.0);
    
    //use cube class to create 3D matrix for adaptive proposal for local RE
    NumericVector tempmn(npars);
    NumericVector tempsigmamn(npars);
    NumericVector tempvar(npars, 0.1 * 0.1);
    NumericVector tempsigmavar(npars, 0.1 * 0.1);
    IntegerVector tempcounts(npars);
    IntegerVector tempsigmacounts(npars);
    double temp = 0.0;
    
//    //print to screen as check
//    for(i = 0; i < npars; i++)
//    {
//        Rprintf("Covariance %d\n", i);
//        for(j = 0; j < 2; j++)
//        {
//            for(k = 0; k < 2; k++) Rprintf("%f\t", tempcovini(j, k, i));
//            Rprintf("\n");
//        }
//        Rprintf("\n");
//    }
//    
//    //print to screen as check
//    for(i = 0; i < npars; i++)
//    {
//        Rprintf("Covariance (slice) %d\n", i);
//        for(j = 0; j < 2; j++)
//        {
//            for(k = 0; k < 2; k++) Rprintf("%f\t", tempcov.slice(i)(j, k));
//            Rprintf("\n");
//        }
//        Rprintf("\n");
//    }
    
    //set up vectors to record acceptance rates
    IntegerVector nacc(npars);
    IntegerVector nattempt(npars);
    IntegerVector nacc_sigma(npars);
    IntegerVector nattempt_sigma(npars);
    IntegerVector nacc_add(npars);
    IntegerVector nattempt_add(npars);
    IntegerVector nacc_del(npars);
    IntegerVector nattempt_del(npars);
    int naccsigmafull = 0, nattemptsigmafull = 0;
    double tempacc;
    double minnacc, maxnacc;
    double minnacc_add, maxnacc_add;
    double minnacc_del, maxnacc_del;
    double minnacc_sigma, maxnacc_sigma;
    
    //set up proposal probabilities to control variable selection
    //and mixing
    double psamp = (varselect == 1 && nitertraining == 0 ? 0.5:1.0);
    
    // run chain
    if(nitertraining > 0) Rprintf("Starting TRAINING run:\n");
    else Rprintf("Starting run:\n");
    for(i = 0; i < niter; i++)
    {
        if(nitertraining > 0)
        {
            if((i + 1) % nitertraining == 0)
            {
                //reset after training run
                psamp = 0.5;
                i = 0;
                
                for(j = 0; j < npars; j++)
                {
                    tempmn[j] = 0.0;
                    tempsigmamn[j] = 0.0;
                    tempvar[j] = 0.1 * 0.1;
                    tempsigmavar[j] = 0.1 * 0.1;
                    tempcounts[j] = 0;
                    tempsigmacounts[j] = 0;
                    
                    nacc[j] = 0;
                    nattempt[j] = 0;
                    nacc_sigma[j] = 0;
                    nattempt_sigma[j] = 0;
                    nacc_add[j] = 0;
                    nattempt_add[j] = 0;
                    nacc_del[j] = 0;
                    nattempt_del[j] = 0;
                }
                naccsigmafull = 0;
                nattemptsigmafull = 0;
                
                nitertraining = 0;
                
                //print initial values to screen
                Rprintf("\nInitial values after training run:\n");
                for(j = 0; j < npars; j++) Rprintf("pars[%d] = %f\n", j, pars(0, j));
                if(random == 0) Rprintf("\nFixed SD component: %f\n", sqrt(priors(0, 1)));
                if(random == 1) Rprintf("\nInitial SD component after training run: %f\n", sigmafull);
                if(random == 2)
                {
                    Rprintf("\nInitial SD components after training run:\n");
                    for(j = 0; j < npars; j++) Rprintf("sigma[%d] = %f\n", j, pars(1, j));
                }
                Rprintf("\n");
                
                Rprintf("Starting MAIN run:\n");
            }
        }
                
        for(j = 0; j < npars; j++) for(k = 0; k < 2; k++) pars_prop(k, j) = pars(k, j);
        
        // propose new value for intercept
        if(runif(1, 0.0, 1.0)[0] < scale) pars_prop(0, 0) = rnorm(1, pars(0, 0), 0.1)[0];
        else pars_prop(0, 0) = rnorm(1, pars(0, 0), sqrt(tempvar[0]))[0];
        nattempt[0]++;
        
        // calculate log-likelihood – log-prior
        LL_prop = loglike(pars_prop, indpars, data);
        acc_prop = LL_prop;
        acc_curr = LL_curr;
        
        acc_curr += R::dnorm(pars(0, 0), priors(0, 0), pars(1, 0), 1);
        acc_prop += R::dnorm(pars_prop(0, 0), priors(0, 0), pars_prop(1, 0), 1);
        //proposals cancel
        
        //accept/reject proposal
        acc = acc_prop - acc_curr;
        if(R_finite(acc) != 0)
        {
            if(log(runif(1, 0.0, 1.0)[0]) < acc)
            {
                pars(0, 0) = pars_prop(0, 0);
                nacc[0]++;
                LL_curr = LL_prop;
            }
            else pars_prop(0, 0) = pars(0, 0);
        }
        else pars_prop(0, 0) = pars(0, 0);
        
        //now propose moves for remaining regression terms
        for(k = 0; k < orignpars; k++)
        {
            if(indpars[cumfactindex[k]] == 1)
            {
                //if variable is included, then propose to remove or move
                if(runif(1, 0.0, 1.0)[0] < psamp)
                {
//                    Rprintf("move\n");
                    //MOVEMENT
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        // propose new value for intercept
                        if(runif(1, 0.0, 1.0)[0] < scale) pars_prop(0, j) = rnorm(1, pars(0, j), 0.1)[0];
                        else pars_prop(0, j) = rnorm(1, pars(0, j), sqrt(tempvar[j]))[0];
                        nattempt[j]++;
                    }
                    // calculate log-likelihood – log-prior
                    LL_prop = loglike(pars_prop, indpars, data);
                    acc_prop = LL_prop;
                    acc_curr = LL_curr;
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc_curr += R::dnorm(pars(0, j), priors(j, 0), pars(1, j), 1);
                        acc_prop += R::dnorm(pars_prop(0, j), priors(j, 0), pars_prop(1, j), 1);
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
                                pars(0, j) = pars_prop(0, j);
                                nacc[j]++;
                            }
                            LL_curr = LL_prop;
                        }
                        else
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) pars_prop(0, j) = pars(0, j);
                        }
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++) pars_prop(0, j) = pars(0, j);
                    }
                }
                else
                {
//                    Rprintf("rem\n");
                    //REMOVAL
                
                    //SET new parameter values
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        pars_prop(0, j) = 0.0;
                        pars_prop(1, j) = 0.0;
                        indpars[j] = 0;
                        nattempt_del[j]++;
                    }
                    
                    // calculate log-likelihood – log-prior
                    LL_prop = loglike(pars_prop, indpars, data);
                    acc_prop = LL_prop;
                    acc_curr = LL_curr;
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc_curr += R::dnorm(pars(0, j), priors(j, 0), pars(1, j), 1);
                        if(random == 2) acc_curr += R::dunif(pars(1, j), priors(npars, 0), priors(npars, 1), 1);
                    }
                    
                    //adjust for proposals
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc_curr -= R::dnorm(pars(0, j), priors(j, 0), pars(1, j), 1);
                        if(random == 2)
                        {
                            temp = (tempcounts[j] < 2 ? sqrt(tempsigmavar[j]):(sqrt(tempsigmavar[j]) / adaptscale_sing));
                            acc_curr -= R::dnorm(pars(1, j), tempsigmamn[j], temp, 1);
                            acc_curr += log(R::pnorm(priors(npars, 1), tempsigmamn[j], temp, 1, 0) - R::pnorm(priors(npars, 0), tempsigmamn[j], temp, 1, 0));
                        }
                    }
                    acc = acc_prop - acc_curr;
//                    acc -= log(1.0 - psamp);
                    
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(runif(1, 0.0, 1.0)[0]) < acc)
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars.col(j) = pars_prop.col(j);
                                nacc_del[j]++;
                            }
                            LL_curr = LL_prop;
                        }
                        else
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars_prop.col(j) = pars.col(j);
                                indpars[j] = 1;
                            }
                        }
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars_prop.col(j) = pars.col(j);
                            indpars[j] = 1;
                        }
                    }
                }
            }
            else
            {
//                Rprintf("add\n");
                //ADDITION
                if(runif(1, 0.0, 1.0)[0] < (1.0 - psamp))
                {
                    //simulate new parameter values
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        if(random == 0) pars_prop(1, j) = sqrt(priors(0, 1));
                        if(random == 1) pars_prop(1, j) = sigmafull;
                        if(random == 2)
                        {
                            //use rejection sampling to sample from truncated normal
                            pars_prop(1, j) = priors(npars, 1) + 1.0;
                            while(pars_prop(1, j) > priors(npars, 1))
                            {
                                pars_prop(1, j) = runif(1, priors(npars, 0), priors(npars, 1))[0];
                                temp = (tempcounts[j] < 2 ? sqrt(tempsigmavar[j]):(sqrt(tempsigmavar[j]) / adaptscale_sing));
                                //set numerator
                                acc_prop = R::pnorm(priors(npars, 1), tempsigmamn[j], temp, 1, 0);
                                acc_prop -= R::pnorm(priors(npars, 0), tempsigmamn[j], temp, 1, 0);
                                acc_prop = R::dnorm(pars_prop(1, j), tempsigmamn[j], temp, 1) - log(acc_prop);
                                //set denominator
                                acc_curr = priors(npars, 1) - priors(npars, 0);
                                acc_curr = (acc_curr + 1.0) / acc_curr;
                                acc_curr = log(acc_curr) + log(tempsigmamn[j]);
                                //accept/reject proposal
                                acc = acc_prop - acc_curr;
                                if(log(runif(1, 0.0, 1.0)[0]) >= acc) pars_prop(1, j) = priors(npars, 1) + 1.0;           
                            }
                        }
                        pars_prop(0, j) = rnorm(1, priors(j, 0), pars_prop(1, j))[0];
                        indpars[j] = 1;
                        nattempt_add[j]++;
                    }
                    
                    // calculate log-likelihood – log-prior
                    LL_prop = loglike(pars_prop, indpars, data);
                    acc_prop = LL_prop;
                    acc_curr = LL_curr;
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc_prop += R::dnorm(pars_prop(0, j), priors(j, 0), pars_prop(1, j), 1);
                        if(random == 2) acc_prop += R::dunif(pars_prop(1, j), priors(npars, 0), priors(npars, 1), 1);
                    }
                    
                    //adjust for proposals
                    for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                    {
                        acc_prop -= R::dnorm(pars_prop(0, j), priors(j, 0), pars_prop(1, j), 1);
                        if(random == 2)
                        {
                            acc_prop -= R::dnorm(pars_prop(1, j), tempsigmamn[j], temp, 1);
                            acc_prop += log(R::pnorm(priors(npars, 1), tempsigmamn[j], temp, 1, 0) - R::pnorm(priors(npars, 0), tempsigmamn[j], temp, 1, 0));
                        }
                    }
                    acc = acc_prop - acc_curr;
//                    acc += log(1.0 - psamp);
                    
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(runif(1, 0.0, 1.0)[0]) < acc)
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars.col(j) = pars_prop.col(j);
                                nacc_add[j]++;
                            }
                            LL_curr = LL_prop;
                        }
                        else
                        {
                            for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                            {
                                pars_prop.col(j) = pars.col(j);
                                indpars[j] = 0;
                            }
                        }
                    }
                    else
                    {
                        for(j = cumfactindex[k]; j < (cumfactindex[k] + factindex[k]); j++)
                        {
                            pars_prop.col(j) = pars.col(j);
                            indpars[j] = 0;
                        }
                    }
                }
            }
        }
        
        if(random == 1)
        {
            //now propose update for SD hyperparameter
            if(runif(1, 0.0, 1.0)[0] < scale) sigmafull_prop = rnorm(1, sigmafull, 0.1)[0];
            else sigmafull_prop = rnorm(1, sigmafull, sqrt(tempsigmavar[0]))[0];
            nattemptsigmafull++;
            
            if(sigmafull_prop > priors(npars, 0) && sigmafull_prop < priors(npars, 1))
            {
                // calculate log-prior for regression terms
                acc_prop = 0.0;
                acc_curr = 0.0;
                for(k = 0; k < nregpars; k++)
                {
                    if(indpars[k + 1] == 1)
                    {
                        acc_prop += R::dnorm(pars(0, k + 1), priors(k + 1, 0), sigmafull_prop, 1);
                        acc_curr += R::dnorm(pars(0, k + 1), priors(k + 1, 0), sigmafull, 1);
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
                        sigmafull = sigmafull_prop;
                        for(k = 0; k < nregpars; k++) pars(1, k + 1) = sigmafull;
                        naccsigmafull++;
                    }
                    else sigmafull_prop = sigmafull;
                }
                else sigmafull_prop = sigmafull;
            }
            else sigmafull_prop = sigmafull;
        }
        
        if(random == 2)
        {
            for(j = 0; j < npars; j++)
            {
                //now propose update for SD hyperparameter
                if(runif(1, 0.0, 1.0)[0] < scale) pars_prop(1, j) = rnorm(1, pars(1, j), 0.1)[0];
                else pars_prop(1, j) = rnorm(1, pars(1, j), sqrt(tempsigmavar[j]))[0];
                nattempt_sigma[j]++;
                
                if(pars_prop(1, j) > priors(npars, 0) && pars_prop(1, j) < priors(npars, 1))
                {
                    // calculate log-prior for regression terms
                    acc_prop = 0.0;
                    acc_curr = 0.0;
                    acc_prop += R::dnorm(pars(0, j), priors(j, 0), pars_prop(1, j), 1);
                    acc_curr += R::dnorm(pars(0, j), priors(j, 0), pars(1, j), 1);
                    //priors for SD cancel
                    //proposals cancel
                    
                    //accept/reject proposal
                    acc = acc_prop - acc_curr;
                    if(R_finite(acc) != 0)
                    {
                        if(log(runif(1, 0.0, 1.0)[0]) < acc)
                        {
                            pars(1, j) = pars_prop(1, j);
                            nacc_sigma[j]++;
                        }
                        else pars_prop(1, j) = pars(1, j);
                    }
                    else pars_prop(1, j) = pars(1, j);
                }
                else pars_prop(1, j) = pars(1, j);
            }
        }
        
        //calculate current unnormalised posterior
        acc_curr = LL_curr;
        //add prior for intercept
        acc_curr += R::dnorm(pars(0, 0), priors(0, 0), sqrt(priors(0, 1)), 1);
        //add priors for regression parameters
        for(j = 0; j < npars; j++)
        {
            if(indpars[j] == 1) acc_curr += R::dnorm(pars(0, j), priors(j, 0), pars(1, j), 1);
        }
        //add variance component
        if(random == 1) acc_curr += R::dunif(sigmafull, priors(npars, 0), priors(npars, 1), 1);
        if(random == 2)
        {
            for(j = 0; j < npars; j++)
            {
                if(indpars[j] == 1) acc_curr += R::dunif(pars(1, j), priors(npars, 0), priors(npars, 1), 1);
            }
        }
        if(R_finite(acc_curr) == 0) stop("Non-finite posterior produced");
        
        // save current value of chain into output matrix
        for(j = 0; j < npars; j++) output(i, j) = pars(0, j);
        if(random == 0)
        {
            for(j = 0; j < npars; j++) output(i, npars + j) = indpars[j];
            output(i, 2 * npars) = acc_curr;
        }
        if(random == 1)
        {
            output(i, npars) = sigmafull;
            for(j = 0; j < npars; j++) output(i, npars + 1 + j) = indpars[j];
            output(i, 2 * npars + 1) = acc_curr;
        }
        if(random == 2)
        {
            for(j = 0; j < npars; j++) output(i, npars + j) = pars(1, j);
            for(j = 0; j < npars; j++) output(i, 2 * npars + j) = indpars[j];
            output(i, 3 * npars) = acc_curr;
        }
        
        if((i + 1) % 1000 == 0)
        {            
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
                minnacc_sigma = (nattempt_sigma[0] > 0 ? (((double) nacc_sigma[0]) / ((double) nattempt_sigma[0])):0.0);
                maxnacc_sigma = minnacc_sigma;
                for(j = 1; j < npars; j++)
                {
                    tempacc = ((double) nacc_sigma[j]) / ((double) nattempt_sigma[j]);
                    if(R_finite(tempacc) != 0)
                    {
                        minnacc_sigma = (minnacc_sigma < tempacc ? minnacc_sigma:tempacc);
                        maxnacc_sigma = (maxnacc_sigma > tempacc ? maxnacc_sigma:tempacc);
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
                if(random == 1) Rprintf("i = %d minmove = %f maxmove = %f sigma = %f minadd = %f maxadd = %f mindel = %f maxdel = %f\n", i + 1, minnacc, maxnacc, naccsigmafull / ((double) nattemptsigmafull), minnacc_add, maxnacc_add, minnacc_del, maxnacc_del);
                if(random == 2) Rprintf("i = %d minmove = %f maxmove = %f minsigma = %f maxsigma = %d minadd = %f maxadd = %f mindel = %f maxdel = %f\n", i + 1, minnacc, maxnacc, minnacc_sigma, maxnacc_sigma, minnacc_add, maxnacc_add, minnacc_del, maxnacc_del);
            }
            else
            {
                for(j = 0; j < nregpars; j++)
                {
                    if(nattempt_add[j] > 0 || nattempt_del[j] > 0) stop("Incorrect addition/removal move");
                }
                if(random == 0) Rprintf("i = %d minmove = %f maxmove = %f\n", i + 1, minnacc, maxnacc);
                if(random == 1) Rprintf("i = %d minmove = %f maxmove = %f sigma = %f\n", i + 1, minnacc, maxnacc, naccsigmafull / ((double) nattemptsigmafull));
                if(random == 2) Rprintf("i = %d minmove = %f maxmove = %f minsigma = %f maxsigma = %d\n", i + 1, minnacc, maxnacc, minnacc_sigma, maxnacc_sigma);
            }
                
            for(j = 0; j < npars; j++)
            {
                nacc[j] = 0;
                nattempt[j] = 0;
                nacc_sigma[j] = 0;
                nattempt_sigma[j] = 0;
                nacc_add[j] = 0;
                nattempt_add[j] = 0;
                nacc_del[j] = 0;
                nattempt_del[j] = 0;
            }
            naccsigmafull = 0;
            nattemptsigmafull = 0;
        }
        
        // record posterior mean and variances for use in proposals
        if ((i + 1) >= ninitial)
        {
            if(random == 0)
            {
                for(j = 0; j < npars; j++)
                {
                    //rescale variances if necessary
                    if(tempcounts[j] >= 2) tempvar[j] = tempvar[j] / adaptscale_sing;
                    calc_meanvar(i, ninitial, &tempmn, &tempvar, &tempcounts, output, j, j, npars + j);
                    if(tempcounts[j] >= 2) tempvar[j] = tempvar[j] * adaptscale_sing;
                 }
            }
            if(random == 1)
            {
                for(j = 0; j < npars; j++)
                {   
                    //rescale variances if necessary
                    if(tempcounts[j] >= 2) tempvar[j] = tempvar[j] / adaptscale_sing; 
                    calc_meanvar(i, ninitial, &tempmn, &tempvar, &tempcounts, output, j, j, npars + 1 + j);
                    if(tempcounts[j] >= 2) tempvar[j] = tempvar[j] * adaptscale_sing;
                }
                if((i + 1) > ninitial) tempsigmavar[0] = tempsigmavar[0] / adaptscale_sing;
                //next line links to intercept for inclusion criteria in order to update sigma (fudge)
                calc_meanvar(i, ninitial, &tempsigmamn, &tempsigmavar, &tempsigmacounts, output, 0, 0, npars + 1);
                if((i + 1) > ninitial) tempsigmavar[0] = tempsigmavar[0] * adaptscale_sing;
            }
            if(random == 2)
            {
                for(j = 0; j < npars; j++)
                {
                    //rescale variances if necessary
                    if(tempcounts[j] >= 2)
                    {
                        tempvar[j] = tempvar[j] / adaptscale_sing;
                        tempsigmavar[j] = tempsigmavar[j] / adaptscale_sing; 
                    }
                    calc_meanvar(i, ninitial, &tempmn, &tempvar, &tempcounts, output, j, j, 2 * npars + j);
                    calc_meanvar(i, ninitial, &tempsigmamn, &tempsigmavar, &tempsigmacounts, output, j + npars, j, 2 * npars + j);
                    if(tempcounts[j] >= 2)
                    {
                        tempvar[j] = tempvar[j] * adaptscale_sing;
                        tempsigmavar[j] = tempsigmavar[j] * adaptscale_sing; 
                    } 
                }
            }
        }
    }
    
    //print estimates of posterior means and variances as a test
    Rprintf("TEST MEANS\n");
    if(random == 2) 
    {
        for(j = 0; j < npars; j++) Rprintf("tempmn[%d] = %f tempvar[%d] = %f counts = %d\n", j, tempmn[j], j, sqrt(tempvar[j] / adaptscale_sing), tempcounts[j]);
        for(j = 0; j < npars; j++) Rprintf("tempsigmamn[%d] = %f tempsigmavar[%d] = %f counts = %d\n", j, tempsigmamn[j], j, sqrt(tempsigmavar[j] / adaptscale_sing), tempcounts[j]);
    }
    else
    {
        for(j = 0; j < npars; j++) Rprintf("tempmn[%d] = %f tempsigma[%d] = %f counts = %d\n", j, tempmn[j], j, sqrt(tempvar[j] / adaptscale_sing), tempcounts[j]);
        if(random == 1) Rprintf("tempmn[%d] = %f tempsigma[%d] = %f counts = %d\n", npars, tempmn[npars], npars, sqrt(tempvar[npars] / adaptscale_sing), tempcounts[npars]);
    }

//    //print to screen as check
//    for(i = 0; i < npars; i++)
//    {
//        Rprintf("Covariance %d\n", i);
//        for(j = 0; j < 2; j++)
//        {
//            for(k = 0; k < 2; k++) Rprintf("%f\t", tempcovini(j, k, i));
//            Rprintf("\n");
//        }
//        Rprintf("\n");
//    }
    
    //print to screen as check
//    if(random == 2)
//    {
//        for(i = 0; i < npars; i++)
//        {
//            Rprintf("Covariance (slice) %d\n", i);
//            for(j = 0; j < 2; j++)
//            {
//                for(k = 0; k < 2; k++) Rprintf("%f\t", tempcov.slice(i)(j, k) / adaptscale);
//                Rprintf("\n");
//            }
//            Rprintf("\n");
//        }
//    }
    
    return(output);
}

