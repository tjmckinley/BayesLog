#include "functions.hpp"

// a Metropolis-Hastings algorithm for fitting a logistic regression model

NumericMatrix logisticMH (NumericMatrix dataR, NumericVector nsamplesR, NumericVector ini_pars, NumericMatrix priors, int niter, int ninitial, double scale, int nadapt, int nprintsum, double maxscale, double niterdim, int nrand, List randindexesL, IntegerMatrix data_randR, IntegerVector nblock, List blockR, int printini, List noncentreintR, NumericVector noncentreintRE)
{
    // 'dataR' is a matrix of data with the first column equal to the response variable
    // 'nsamplesR' corresponds to aggregated counts for each row of 'data'
    // 'nrandint' corresponds to number of random intercept terms
    // 'ini_pars' is a vector of initial values for the unknown parameters
    // 'priors' is an (npars x 2) matrix containing the mean and var for Normal priors
    // 'niter' is the number of iterations over which to run the chain
    // 'scale' is mixing proportion for adaptive MCMC
    // 'nadapt' is the number of iterations to run between updates of the adaptive proposal
    // 'nprintsum' controls how often run time information is printed to the screen
    // 'maxscale' is the maximum scaling of the adaptive proposal variance at each update
    // 'niterdim' is the iteration at which the diminishing adaptation component kicks in
    // 'nrand' is number of random effect terms
    // 'randindexesL' is a list of matrices of indexes for referencing random effects terms
    // 'data_randR' is matrix of data corresponding to random effect terms
    // 'nblock' is vector of block lengths
    // 'blockR' is list of indexes for block updating
    // 'printini' is a binary indicator. If equal to 1 then prints setup information
    //      else it just prints chain information (suppresses extraneous information
    //      when running multiple chains)
    // 'noncentreintR' is list denoting which parameters to non-centre with respect
    // to the intercept
    // 'noncentreintRE' is vector denoting which RE parameters to non-centre with respect
    // to the intercept
    
    //initialise indexing variables
    int i, j, k, m;
    
    //set number of threads
    int ncores = omp_get_num_procs();
    omp_set_num_threads(ncores);
    
    // calculate number of parameters
    int npars = ini_pars.size();
    
    //set terms for parallelising code
    double *logL = (double *) R_alloc(dataR.nrow(), sizeof(double));
    for(i = 0; i < dataR.nrow(); i++) logL[i] = 0.0;
    
    //print runtime information to the screen    
    i = 0;
    for(j = 0; j < nsamplesR.size(); j++) i += nsamplesR[j];
    if(printini == 1)
    {
        Rprintf("\nNumber of samples in data set = %d\n", i);
        Rprintf("Number of unique samples in data set = %d\n", dataR.nrow());
        
        Rprintf("\nRun time information printed to screen every %d iterations\n", nprintsum);
        Rprintf("Number of iterations = %d\n", niter);
        Rprintf("Scale for adaptive proposal = %.2f\n", scale);
        Rprintf("Number of regression parameters = %d\n", npars);
        if(nrand > 0) Rprintf("Number of random intercepts = %d\n", nrand);
        Rprintf("Adapt every %d iterations\n", nadapt);
        Rprintf("Max scale for adapting = %.2f\n", maxscale);
        Rprintf("%d iterations before diminishing adaptation kicks in\n", (int) niterdim);
        
        //print prior information to screen
        Rprintf("\nPriors: mean = %.2f variance = %.2f\n", priors(0, 0), priors(0, 1));
        if(nrand > 0) Rprintf("Priors RE SD: lower = %.2f upper = %.2f\n", priors(npars - 1, 0), priors(npars - 1, 1));
    }
    
    //convert Rcpp objects to native C objects for fast processing
    int data_nrows = dataR.nrow();
    int data_ncols = dataR.ncol();
    double **data = (double **) R_alloc (data_nrows, sizeof(double *));
    for(i = 0; i < data_nrows; i++) data[i] = (double *) R_alloc (data_ncols, sizeof(double));
    for(i = 0; i < data_nrows; i++) for(j = 0; j < data_ncols; j++) data[i][j] = dataR(i, j);
    
    double *nsamples = (double *) R_alloc (data_nrows, sizeof(double));
    for(i = 0; i < data_nrows; i++) nsamples[i] = nsamplesR(i);
    
    int **data_rand = (int **) R_alloc (data_nrows, sizeof(int *));
    for(i = 0; i < data_nrows; i++) data_rand[i] = (int *) R_alloc (nrand, sizeof(int));
    for(i = 0; i < data_nrows; i++) for(j = 0; j < nrand; j++) data_rand[i][j] = data_randR(i, j);
    
    //extract indexes from list object for random effect terms
    int *nrandlevels = (int *) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int));
    int *cumrandlevels = (int *) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int));
    int **nrandindexes = (int **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int *));
    int ***randindexes = (int ***) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int **));
    if(nrand > 0)
    {
        for(j = 0; j < nrand; j++)
        {
            nrandlevels[j] = as<List>(randindexesL[j]).size();
            nrandindexes[j] = (int *) R_alloc (nrandlevels[j], sizeof(int));
            randindexes[j] = (int **) R_alloc (nrandlevels[j], sizeof(int *));
            for(k = 0; k < nrandlevels[j]; k++)
            {
                nrandindexes[j][k] = as<IntegerVector>(as<List>(randindexesL[j])[k]).size();
                randindexes[j][k] = (int *) R_alloc (nrandindexes[j][k], sizeof(int));
                for(m = 0; m < nrandindexes[j][k]; m++) randindexes[j][k][m] = as<IntegerVector>(as<List>(randindexesL[j])[k])[m];
            }
        }
        cumrandlevels[0] = 0;
        for(j = 1; j < nrand; j++) cumrandlevels[j] = cumrandlevels[j - 1] + nrandlevels[j - 1];
//        for(j = 0; j < nrand; j++)
//        {
//            for(k = 0; k < nrandlevels[j]; k++)
//            {
//                for(m = 0; m < nrandindexes[j][k]; m++) Rprintf("rand(%d, %d, %d) = %d ", j, k, m, randindexes[j][k][m]);
//                Rprintf("\n");
//            }
//        }
    }
    else
    {
        for(j = 0; j < nrand; j++)
        {
            cumrandlevels[j] = 0;
            nrandlevels[j] = 0;
        }
    }
    
    // set up output vector of length 'niter' to record chain
    // (append extra column for unnormalised posterior)
    int noutput = npars + cumrandlevels[nrand - 1] + nrandlevels[nrand - 1] + 1;
    NumericMatrix output(niter, noutput);
    
    // initialise chain and set up vector to hold proposals
    double *pars = (double *) R_alloc (npars, sizeof(double));
    double *pars_prop = (double *) R_alloc (npars, sizeof(double));
    for(i = 0; i < npars; i++)
    {
        pars[i] = 0.0;
        pars_prop[i] = 0.0;
    }
    
    double **rand = (double **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(double *));
    double **rand_prop = (double **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(double *));
    if(nrand > 0)
    {
        for(i = 0; i < nrand; i++)
        {
            rand[i] = (double *) R_alloc (nrandlevels[i], sizeof(double));
            rand_prop[i] = (double *) R_alloc (nrandlevels[i], sizeof(double));
            for(j = 0; j < nrandlevels[i]; j++)
            {
                rand[i][j] = 0.0;
                rand_prop[i][j] = 0.0;
            }
        }
    }
    
    //declare variables
    double LL_curr, LL_prop, acc_curr, acc_prop, acc;
    
    //set initial values
    for(i = 0; i < npars; i++) pars[i] = ini_pars(i);
    //check the initial values produce a finite log-posterior
    LL_curr = loglike(pars, data_nrows, data_ncols, data, nsamples, nrand, rand, data_rand, logL);
    acc_curr = LL_curr;
    //add hierarchical terms
    if(nrand > 0)
    {
        for(j = 0; j < nrand; j++)
        {
            for(k = 0; k < nrandlevels[j]; k++) acc_curr += R::dnorm(rand[j][k], 0.0, pars[j + npars - nrand], 1);
        }
    }
    //add priors for regression parameters
    for(j = 0; j < (npars - nrand); j++) acc_curr += R::dnorm(pars[j], priors(j, 0), sqrt(priors(j, 1)), 1);
    if(nrand > 0) for(j = (npars - nrand); j < npars; j++) acc_curr += R::dunif(pars[j], priors(j, 0), priors(j, 1), 1);
    
//    //print initial values to screen
//    Rprintf("\nInitial values:\n");
//    for(j = 0; j < npars; j++) Rprintf("pars[%d] = %.2f\n", j, pars[j]);
//    Rprintf("\n");
    
    //check viability
    if(R_finite(acc_curr) == 0) Rcpp::stop("Initial values are not viable: LL = %f acc = %f", LL_curr, acc_curr);
    
    //set number of non-centering components
    int nnoncentre = noncentreintR.size();
    
    //set up objects for non-centring
    arma::field<arma::ivec> noncentreint (nnoncentre);
    for(m = 0; m < nnoncentre; m++)
    {
        SEXP tempivR = noncentreintR[m];
        arma::ivec tempiv = as<arma::ivec>(tempivR); 
        noncentreint(m) = tempiv;
    }
    
    double *propsdnon = (double *) R_alloc (nnoncentre, sizeof(double));
    int *naccnon = (int *) R_alloc (nnoncentre, sizeof(int));
    int *naccnon1 = (int *) R_alloc (nnoncentre, sizeof(int));
    int *nattemptnon = (int *) R_alloc (nnoncentre, sizeof(int));
    int *nattemptnon1 = (int *) R_alloc (nnoncentre, sizeof(int));
    
    //initialise objects if required
    for(i = 0; i < nnoncentre; i++)
    {
        propsdnon[i] = 0.1;
        naccnon[i] = 0;
        naccnon1[i] = 0;
        nattemptnon[i] = 0;
        nattemptnon1[i] = 0;
    }
    
    //check whether non-centring really required
    if(nnoncentre == 1)
    {
        j = 0;
        for(i = 0; i < (npars - nrand); i++) j += noncentreint(0)(i);
        if(j == 0) nnoncentre = 0;
    }
    
    //set number of non-centering components
    int nnoncentreRE = noncentreintRE.size();
    
    double *propsdnonRE = (double *) R_alloc (nnoncentreRE, sizeof(double));
    int *naccnonRE = (int *) R_alloc (nnoncentreRE, sizeof(int));
    int *naccnonRE1 = (int *) R_alloc (nnoncentreRE, sizeof(int));
    int *nattemptnonRE = (int *) R_alloc (nnoncentreRE, sizeof(int));
    int *nattemptnonRE1 = (int *) R_alloc (nnoncentreRE, sizeof(int));
    
    //initialise objects if required
    for(i = 0; i < nnoncentreRE; i++)
    {
        propsdnonRE[i] = 0.1;
        naccnonRE[i] = 0;
        naccnonRE1[i] = 0;
        nattemptnonRE[i] = 0;
        nattemptnonRE1[i] = 0;
    }
    
    //check whether non-centring really required
    if(nnoncentreRE == 1)
    {
        j = 0;
        for(i = 0; i < nnoncentreRE; i++) j += noncentreintRE(i);
        if(j == 0) nnoncentreRE = 0;
    }
    
    if(nnoncentreRE > 0) for(i = 0; i < nnoncentreRE; i++) noncentreintRE(i)--;
    
    //extract number of blocks
    int nblocks = nblock.size();
    
    //set up vectors for adaptive proposals
    double *propsd = (double *) R_alloc (nblocks + nrand, sizeof(double));
    for(i = 0; i < (nblocks + nrand); i++) propsd[i] = 0.1;
    for(i = 0; i < nblocks; i++) propsd[i] = (nblock(i) > 1 ? 1.0:0.1);
    
    double **propsd_rand = (double **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(double *));
    if(nrand > 0)
    {
        for(i = 0; i < nrand; i++)
        {
            propsd_rand[i] = (double *) R_alloc (nrandlevels[i], sizeof(double));
            for(j = 0; j < nrandlevels[i]; j++) propsd_rand[i][j] = 0.1;
        }
    }
    
    //set up covariance matrices for block-updating if required
    arma::field<arma::ivec> block (nblocks);
    arma::field<arma::vec> tempmn (nblocks);
    arma::field<arma::vec> pars_a (nblocks);
    arma::field<arma::vec> pars_prop_a (nblocks);
    arma::field<arma::mat> meanmat (nblocks);
    arma::field<arma::mat> meanmat1 (nblocks);
    arma::field<arma::mat> propcov (nblocks);
    arma::field<arma::mat> propcovI (nblocks);
    arma::vec tempv (1);
    arma::mat tempm (1, 1);
    for(m = 0; m < nblocks; m++)
    {
        SEXP tempivR = blockR[m];
        arma::ivec tempiv = as<arma::ivec>(tempivR); 
        block(m) = tempiv;
        
        tempv.set_size(nblock(m));
        tempv.zeros();
        tempm.set_size(nblock(m), nblock(m));
        tempm.zeros();
        tempmn(m) = tempv;
        pars_a(m) = tempv;
        pars_prop_a(m) = tempv;
        meanmat(m) = tempm;
        meanmat1(m) = tempm;
        
        tempm.eye();
        tempm.diag() *= pow(0.1, 2.0) / ((double) nblock(m));
        propcov(m) = tempm;
        propcovI(m) = tempm;
    }
    
    //set up vectors to record acceptance rates
    int *nacc = (int *) R_alloc (nblocks + nrand, sizeof(int));
    int *nacc1 = (int *) R_alloc (nblocks + nrand, sizeof(int));
    int *nattempt = (int *) R_alloc (nblocks + nrand, sizeof(int));
    int *nattempt1 = (int *) R_alloc (nblocks + nrand, sizeof(int));
    for(i = 0; i < (nblocks + nrand); i++)
    {
        nacc[i] = 0;
        nacc1[i] = 0;
        nattempt[i] = 0;
        nattempt1[i] = 0;
    }
    
    int **nacc_rand = (int **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int *));
    int **nattempt_rand = (int **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int *));
    int **nacc1_rand = (int **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int *));
    int **nattempt1_rand = (int **) R_alloc ((nrand == 0 ? 1:nrand), sizeof(int *));
    if(nrand > 0)
    {
        for(i = 0; i < nrand; i++)
        {
            nacc_rand[i] = (int *) R_alloc (nrandlevels[i], sizeof(int));
            nattempt_rand[i] = (int *) R_alloc (nrandlevels[i], sizeof(int));
            nacc1_rand[i] = (int *) R_alloc (nrandlevels[i], sizeof(int));
            nattempt1_rand[i] = (int *) R_alloc (nrandlevels[i], sizeof(int));
            for(j = 0; j < nrandlevels[i]; j++)
            {
                nacc_rand[i][j] = 0;
                nattempt_rand[i][j] = 0;
                nacc1_rand[i][j] = 0;
                nattempt1_rand[i][j] = 0;
            }
        }
    }
    
    double tempacc, minacc, maxacc;
    double minaccnon = 0.0, maxaccnon = 0.0;
    double minacc_rand, maxacc_rand;
    double minaccnonRE = 0.0, maxaccnonRE = 0.0;
        
    //set proposals 
    for(k = 0; k < npars; k++) pars_prop[k] = pars[k];
    if(nrand > 0) for(k = 0; k < nrand; k++) for(j = 0; j < nrandlevels[k]; j++) rand_prop[k][j] = rand[k][j];
    
    //initialise timer
    Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;
    
    // run chain
    Rprintf("\nStarting run:\n");
    for(i = 0; i < niter; i++)
    {  
        //check for user interruptions
        if (i % 10 == 0) R_CheckUserInterrupt();
        
        ///////////////////////////////////////////////////////
        //////////////      CENTRED UPDATES       /////////////
        ///////////////////////////////////////////////////////
        
        //cycle through update blocks
        for(m = 0; m < nblocks; m++)
        {
            if(nblock(m) > 1)
            {
                for(j = 0; j < nblock(m); j++) pars_a(m)(j) = pars[block(m)(j)];
                if(R::runif(0.0, 1.0) < scale) pars_prop_a(m) = arma::conv_to<arma::vec>::from(mvrnormArma(1, pars_a(m), propcovI(m)));
                else pars_prop_a(m) = arma::conv_to<arma::vec>::from(mvrnormArma(1, pars_a(m), propcov(m) * propsd[m]));
                for(j = 0; j < nblock(m); j++)
                {
                    nattempt[m]++;
                    pars_prop[block(m)(j)] = pars_prop_a(m)(j);
                }
            }
            else
            {
                // propose new parameter value
                if(R::runif(0.0, 1.0) < scale) pars_prop[block(m)(0)] = R::rnorm(pars[block(m)(0)], 0.1);
                else pars_prop[block(m)(0)] = R::rnorm(pars[block(m)(0)], propsd[m]);
                nattempt[m]++;
            }
            
            
            // calculate log-likelihood
            LL_prop = loglike(pars_prop, data_nrows, data_ncols, data, nsamples, nrand, rand, data_rand, logL);
            acc_prop = LL_prop;
            acc_curr = LL_curr;
            //adjust for prior
            for(j = 0; j < nblock(m); j++)
            {
                k = block(m)(j);
                acc_curr += R::dnorm(pars[k], priors(k, 0), sqrt(priors(k, 1)), 1);
                acc_prop += R::dnorm(pars_prop[k], priors(k, 0), sqrt(priors(k, 1)), 1);
            }
            acc = acc_prop - acc_curr;
            //proposals cancel
            
            //accept/reject proposal
            if(R_finite(acc) != 0)
            {
                if(log(R::runif(0.0, 1.0)) < acc)
                {
                    nacc[m]++;
                    for(j = 0; j < nblock(m); j++)
                    {
                        k = block(m)(j);
                        pars[k] = pars_prop[k];
                    }
                    LL_curr = LL_prop;
                }
                else
                {
                    for(j = 0; j < nblock(m); j++)
                    {
                        k = block(m)(j);
                        pars_prop[k] = pars[k];
                    }
                }
            }
            else
            {
                for(j = 0; j < nblock(m); j++)
                {
                    k = block(m)(j);
                    pars_prop[k] = pars[k];
                }
            }
        }
        
        ///////////////////////////////////////////////////////
        //////////////    NONCENTRED UPDATES      /////////////
        ///////////////////////////////////////////////////////
        
        if(nnoncentre > 0)
        {
            for(m = 0; m < nnoncentre; m++)
            {
                // propose new parameter value
                if(R::runif(0.0, 1.0) < scale) pars_prop[0] = R::rnorm(pars[0], 0.1);
                else pars_prop[0] = R::rnorm(pars[0], propsdnon[m]);
                nattemptnon[m]++;
                
                //update non-centred parameters
                for(j = 1; j < (npars - nrand); j++)
                {
                    pars_prop[j] = pars[j];
                    
                    if(noncentreint(m)(j) == 1)
                    {
                        double temp = pars[0] + pars[j];
                        pars_prop[j] = temp - pars_prop[0];
                    }
                }
                
                // calculate log-likelihood
                LL_prop = loglike(pars_prop, data_nrows, data_ncols, data, nsamples, nrand, rand, data_rand, logL);
                acc_prop = LL_prop;
                acc_curr = LL_curr;
                
                //adjust for priors
                acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
                acc_prop += R::dnorm(pars_prop[0], priors(0, 0), sqrt(priors(0, 1)), 1);
                for(j = 1; j < (npars - nrand); j++)
                {
                    if(noncentreint(m)(j) == 1)
                    {
                        acc_curr += R::dnorm(pars[j] + pars[0], pars[0] + priors(j, 0), sqrt(priors(j, 1)), 1);
                        acc_prop += R::dnorm(pars_prop[j] + pars_prop[0], pars_prop[0] + priors(j, 0), sqrt(priors(j, 1)), 1);
                    }
                }
                acc = acc_prop - acc_curr;
                
                //proposals cancel
                
                //accept/reject proposal
                if(R_finite(acc) != 0)
                {
                    if(log(R::runif(0.0, 1.0)) < acc)
                    {
                        naccnon[m]++;
                        for(j = 0; j < (npars - nrand); j++) pars[j] = pars_prop[j];
                        LL_curr = LL_prop;
                    }
                    else for(j = 0; j < (npars - nrand); j++) pars_prop[j] = pars[j];
                }
                else for(j = 0; j < (npars - nrand); j++) pars_prop[j] = pars[j];
            }
        }
        
        ///////////////////////////////////////////////////////
        //////////////    RANDOM INTERCEPTS       /////////////
        ///////////////////////////////////////////////////////
        
        if(nrand > 0)
        {
            //update hierarchical variance terms
            for(k = (npars - nrand); k < npars; k++)
            {
                // propose new parameter value
                if(R::runif(0.0, 1.0) < scale) pars_prop[k] = R::rnorm(pars[k], 0.1);
                else pars_prop[k] = R::rnorm(pars[k], propsd[k + nblocks - (npars - nrand)]);
                nattempt[k + nblocks - (npars - nrand)]++;
                
                if(pars_prop[k] > priors(k, 0) && pars_prop[k] < priors(k, 1))
                {
                    //add hierarchical terms
                    LL_prop = 0.0;
                    for(j = 0; j < nrandlevels[k - (npars - nrand)]; j++)
                    {
                        LL_prop += R::dnorm(rand[k - (npars - nrand)][j], 0.0, pars_prop[k], 1);
                        LL_prop -= R::dnorm(rand[k - (npars - nrand)][j], 0.0, pars[k], 1);
                    }
                    
                    //priors and proposals cancel
                    
                    acc = LL_prop;
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(R::runif(0.0, 1.0)) < acc)
                        {
                            pars[k] = pars_prop[k];
                            nacc[k + nblocks - (npars - nrand)]++;
                        }
                        else pars_prop[k] = pars[k];
                    }
                    else pars_prop[k] = pars[k];
                }
                else pars_prop[k] = pars[k];
            }
            
            ///////////////////////////////////////////////////////
            //////////////     CENTRED RE UPDATES     /////////////
            ///////////////////////////////////////////////////////
            
            //update individual hierarchical components
            for(k = 0; k < nrand; k++)
            {
                for(j = 0; j < nrandlevels[k]; j++)
                {
                    // propose new parameter value
                    if(R::runif(0.0, 1.0) < scale) rand_prop[k][j] = R::rnorm(rand[k][j], 0.1);
                    else rand_prop[k][j] = R::rnorm(rand[k][j], propsd_rand[k][j]);
                    nattempt_rand[k][j]++;
                    
                    //likelihood contributions
                    acc_prop = loglike_sub(pars, nrandindexes[k][j], data_ncols, data, nsamples, nrand, rand_prop, data_rand, randindexes, nrandindexes, k, j, logL);
                    acc_curr = loglike_sub(pars, nrandindexes[k][j], data_ncols, data, nsamples, nrand, rand, data_rand, randindexes, nrandindexes, k, j, logL);
                    
                    //hierarchical terms
                    acc_prop += R::dnorm(rand_prop[k][j], 0.0, pars[k + npars - nrand], 1);
                    acc_curr += R::dnorm(rand[k][j], 0.0, pars[k + npars - nrand], 1);
                    
                    //proposals cancel
                    
                    acc = acc_prop - acc_curr;
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(R::runif(0.0, 1.0)) < acc)
                        {
                            rand[k][j] = rand_prop[k][j];
                            nacc_rand[k][j]++;
                        }
                        else rand_prop[k][j] = rand[k][j];
                    }
                    else rand_prop[k][j] = rand[k][j];
                }
            }
            
            //calculate full log-likelihood
            LL_curr = loglike(pars, data_nrows, data_ncols, data, nsamples, nrand, rand, data_rand, logL);
            
            ///////////////////////////////////////////////////////
            //////////////   NONCENTRED RE UPDATES    /////////////
            ///////////////////////////////////////////////////////
            
            if(nnoncentreRE > 0)
            {
                for(m = 0; m < nnoncentreRE; m++)
                {
                    // propose new parameter value
                    if(R::runif(0.0, 1.0) < scale) pars_prop[0] = R::rnorm(pars[0], 0.1);
                    else pars_prop[0] = R::rnorm(pars[0], propsdnonRE[m]);
                    nattemptnonRE[m]++;
                    
                    //update non-centred parameters
                    k = noncentreintRE(m);
                    for(j = 0; j < nrandlevels[k]; j++)
                    {
                        double temp = pars[0] + rand[k][j];
                        rand_prop[k][j] = temp - pars_prop[0];
                    }
                    
                    // calculate log-likelihood
                    LL_prop = loglike(pars_prop, data_nrows, data_ncols, data, nsamples, nrand, rand_prop, data_rand, logL);
                    acc_prop = LL_prop;
                    acc_curr = LL_curr;
                    
                    //adjust for priors
                    acc_curr += R::dnorm(pars[0], priors(0, 0), sqrt(priors(0, 1)), 1);
                    acc_prop += R::dnorm(pars_prop[0], priors(0, 0), sqrt(priors(0, 1)), 1);
                    for(j = 0; j < nrandlevels[k]; j++)
                    {
                        acc_prop += R::dnorm(rand_prop[k][j] + pars_prop[0], pars_prop[0], pars[k + npars - nrand], 1);
                        acc_curr += R::dnorm(rand[k][j] + pars[0], pars[0], pars[k + npars - nrand], 1);
                    }
                        
                    //proposals cancel
                    
                    acc = acc_prop - acc_curr;
                    //accept/reject proposal
                    if(R_finite(acc) != 0)
                    {
                        if(log(R::runif(0.0, 1.0)) < acc)
                        {
                            pars[0] = pars_prop[0];
                            for(j = 0; j < nrandlevels[k]; j++) rand[k][j] = rand_prop[k][j];
                            LL_curr = LL_prop;
                            naccnonRE[m]++;
                        }
                        else for(j = 0; j < nrandlevels[k]; j++) rand_prop[k][j] = rand[k][j];
                    }
                    else for(j = 0; j < nrandlevels[k]; j++) rand_prop[k][j] = rand[k][j];
                }
            }
        }
        
        //calculate the posterior for clarity
        acc_curr = LL_curr;
        //add hierarchical terms
        if(nrand > 0)
        {
            for(j = 0; j < nrand; j++)
            {
                for(k = 0; k < nrandlevels[j]; k++) acc_curr += R::dnorm(rand[j][k], 0.0, pars[j + npars - nrand], 1);
            }
        }
        //add priors for regression parameters
        for(j = 0; j < (npars - nrand); j++) acc_curr += R::dnorm(pars[j], priors(j, 0), sqrt(priors(j, 1)), 1);
        if(nrand > 0) for(j = (npars - nrand); j < npars; j++) acc_curr += R::dunif(pars[j], priors(j, 0), priors(j, 1), 1);
        
        // save current value of chain into output matrix
        for(j = 0; j < npars; j++) output(i, j) = pars[j];
        for(j = 0; j < nrand; j++)
        {
            for(k = 0; k < nrandlevels[j]; k++) output(i, npars + cumrandlevels[j] + k) = rand[j][k];
        }
        output(i, npars + cumrandlevels[nrand - 1] + nrandlevels[nrand - 1]) = acc_curr;
        
        // record posterior mean and variances for use in proposals
        if ((i + 1) % nadapt == 0)
        {
            for(m = 0; m < nblocks; m++)
            {
                propsd[m] = adapt_scale(nacc[m] - nacc1[m], nattempt[m] - nattempt1[m], (nblock(m) > 1 ? 0.23:0.44), propsd[m], (double) i + 1, maxscale, niterdim);
            }
            if(nrand > 0)
            {
                for(m = nblocks; m < (nblocks + nrand); m++)
                {
                    propsd[m] = adapt_scale(nacc[m] - nacc1[m], nattempt[m] - nattempt1[m], 0.44, propsd[m], (double) i + 1, maxscale, niterdim);
                }
            }
            if(nnoncentre > 0)
            {
                for(m = 0; m < nnoncentre; m++) propsdnon[m] = adapt_scale(naccnon[m] - naccnon1[m], nattemptnon[m] - nattemptnon1[m], 0.23, propsdnon[m], (double) i + 1, maxscale, niterdim);
            }
            if(nrand > 0)
            {
                for(j = 0; j < nrand; j++)
                {
                    for(k = 0; k < nrandlevels[j]; k++) propsd_rand[j][k] = adapt_scale(nacc_rand[j][k] - nacc1_rand[j][k], nattempt_rand[j][k] - nattempt1_rand[j][k], 0.44, propsd_rand[j][k], (double) i + 1, maxscale, niterdim);
                }
                if(nnoncentreRE > 0)
                {
                    for(m = 0; m < nnoncentreRE; m++) propsdnonRE[m] = adapt_scale(naccnonRE[m] - naccnonRE1[m], nattemptnonRE[m] - nattemptnonRE1[m], 0.23, propsdnonRE[m], (double) i + 1, maxscale, niterdim);
                }
            }
            
            //update counts
            for(j = 0; j < (nblocks + nrand); j++)
            {
                nacc1[j] = nacc[j];
                nattempt1[j] = nattempt[j];
            }
            if(nnoncentre > 0)
            {
                for(m = 0; m < nnoncentre; m++)
                {
                    naccnon1[m] = naccnon[m];
                    nattemptnon1[m] = nattemptnon[m];
                }
            }
            if(nrand > 0)
            {
                for(j = 0; j < nrand; j++)
                {
                    for(k = 0; k < nrandlevels[j]; k++)
                    {
                        nacc1_rand[j][k] = nacc_rand[j][k];
                        nattempt1_rand[j][k] = nattempt_rand[j][k];
                    }
                }
                if(nnoncentreRE > 0)
                {
                    for(m = 0; m < nnoncentreRE; m++)
                    {
                        naccnonRE1[m] = naccnonRE[m];
                        nattemptnonRE1[m] = nattemptnonRE[m];
                    }
                }
            }
        }  
                
        //update var-cov matrix for adaptive proposal
 		if((i + 1) >= ninitial)
 		{
 		    for(m = 0; m < nblocks; m++)
 		    {
 		        if(nblock(m) > 1) 
 		        {
 		            adapt_update(i, ninitial, niter, nblock(m), &(tempmn(m)), &(meanmat(m)), &(meanmat1(m)), output, &(propcov(m)), i, &(block(m)));
//                    if((i + 1) % 100 == 0)
//                    {
//                        Rprintf("\nPropcov:\n");
//                        for(k = 0; k < nblock(m); k++)
//                        {
//                            for(j = 0; j < nblock(m); j++) Rprintf("%f ", propcov(m)(k, j));
//                            Rprintf("\n");
//                        }
    //                    Rprintf("\nPropcov_t:\n");
    //                    for(k = 0; k < npars; k++)
    //                    {
    //                        for(j = 0; j < npars; j++) Rprintf("%f ", propcov_temp(m)(k, j));
    //                        Rprintf("\n");
    //                    }
//                    }
                }
 		    }
 		}
        
        if((i + 1) % nprintsum == 0)
        {          
            // print some output to screen for book-keeping
            minacc = ((double) nacc[0]) / ((double) nattempt[0]);
            maxacc = minacc;
            for(j = 1; j < (nblocks + nrand); j++)
            {
                tempacc = ((double) nacc[j]) / ((double) nattempt[j]);
                minacc = (minacc < tempacc ? minacc:tempacc);
                maxacc = (maxacc > tempacc ? maxacc:tempacc);
            }
            if(nnoncentre > 0)
            {
                minaccnon = ((double) naccnon[0]) / ((double) nattemptnon[0]);
                maxaccnon = minaccnon;
                for(m = 1; m < nnoncentre; m++)
                {
                    tempacc = ((double) naccnon[m]) / ((double) nattemptnon[m]);
                    minaccnon = (minaccnon < tempacc ? minaccnon:tempacc);
                    maxaccnon = (maxaccnon > tempacc ? maxaccnon:tempacc);
                }
            }
            if(nrand > 0)
            {
                minacc_rand = ((double) nacc_rand[0][0]) / ((double) nattempt_rand[0][0]);
                maxacc_rand = minacc_rand;
                for(j = 0; j < nrand; j++)
                {
                    for(k = 0; k < nrandlevels[j]; k++)
                    {
                        tempacc = ((double) nacc_rand[j][k]) / ((double) nattempt_rand[j][k]);
                        minacc_rand = (minacc_rand < tempacc ? minacc_rand:tempacc);
                        maxacc_rand = (maxacc_rand > tempacc ? maxacc_rand:tempacc);
                    }
                }
                
                if(nnoncentreRE > 0)
                {
                    minaccnonRE = ((double) naccnonRE[0]) / ((double) nattemptnonRE[0]);
                    maxaccnonRE = minaccnonRE;
                    for(m = 1; m < nnoncentreRE; m++)
                    {
                        tempacc = ((double) naccnonRE[m]) / ((double) nattemptnonRE[m]);
                        minaccnonRE = (minaccnonRE < tempacc ? minaccnonRE:tempacc);
                        maxaccnonRE = (maxaccnonRE > tempacc ? maxaccnonRE:tempacc);
                    }
                }
                
                timer.step("");
                NumericVector res(timer);
                if(nnoncentre > 0)
                {
                    if(nnoncentreRE > 0)
                    {
                        Rprintf("i = %d minacc = %.2f maxacc = %.2f minaccnon = %.2f maxaccnon = %.2f\n\tminacc_rand = %.2f maxacc_rand = %.2f minacc_randnon = %.2f maxacc_randnon = %.2f\n\time = %.2f\n", i + 1, minacc, maxacc, minaccnon, maxaccnon, minacc_rand, maxacc_rand, minaccnonRE, maxaccnonRE, (res[timer_cnt] / 1e9) - prev_time);
                    }
                    else
                    {
                        Rprintf("i = %d minacc = %.2f maxacc = %.2f minaccnon = %.2f maxaccnon = %.2f\n\tminacc_rand = %.2f maxacc_rand = %.2f time = %.2f\n", i + 1, minacc, maxacc, minaccnon, maxaccnon, minacc_rand, maxacc_rand, (res[timer_cnt] / 1e9) - prev_time);
                    }
                }
                else
                {
                    if(nnoncentreRE > 0)
                    {
                        Rprintf("i = %d minacc = %.2f maxacc = %.2f minacc_rand = %.2f maxacc_rand = %.2f\n\tminacc_randnon = %.2f maxacc_randnon = %.2f time = %.2f\n", i + 1, minacc, maxacc, minaccnon, maxaccnon, minacc_rand, maxacc_rand, minaccnonRE, maxaccnonRE, (res[timer_cnt] / 1e9) - prev_time);
                    }
                    else
                    {
                        Rprintf("i = %d minacc = %.2f maxacc = %.2f minacc_rand = %.2f maxacc_rand = %.2f time = %.2f\n", i + 1, minacc, maxacc, minaccnon, maxaccnon, minacc_rand, maxacc_rand, (res[timer_cnt] / 1e9) - prev_time);
                    }
                }
                prev_time = res[timer_cnt] / 1e9;
                timer_cnt++;
            }
            else
            {
                timer.step("");
                NumericVector res(timer);
                if(nnoncentre > 0)
                {
                    Rprintf("i = %d minacc = %.2f maxacc = %.2f minaccnon = %.2f maxaccnon = %.2f time = %.2f\n", i + 1, minacc, maxacc, minaccnon, maxaccnon, (res[timer_cnt] / 1e9) - prev_time);
                }
                else
                {
                    Rprintf("i = %d minacc = %.2f maxacc = %.2f time = %.2f\n", i + 1, minacc, maxacc, minaccnon, maxaccnon, (res[timer_cnt] / 1e9) - prev_time);
                }
                prev_time = res[timer_cnt] / 1e9;
                timer_cnt++;
            }
               
            //reset counters 
            for(j = 0; j < (nblocks + nrand); j++)
            {
                nacc[j] = 0; nacc1[j] = 0;
                nattempt[j] = 0; nattempt1[j] = 0;
            }
            if(nnoncentre > 0)
            {
                for(m = 0; m < nnoncentre; m++)
                {
                    naccnon[m] = 0; naccnon1[m] = 0;
                    nattemptnon[m] = 0; nattemptnon1[m] = 0;
                }
            }
            if(nrand > 0)
            {
                for(j = 0; j < nrand; j++)
                {
                    for(k = 0; k < nrandlevels[j]; k++)
                    {
                        nacc_rand[j][k] = 0;
                        nacc1_rand[j][k] = 0;
                        nattempt_rand[j][k] = 0;
                        nattempt1_rand[j][k] = 0;
                    }
                }
                if(nnoncentreRE > 0)
                {
                    for(m = 0; m < nnoncentreRE; m++)
                    {
                        naccnonRE[m] = 0; naccnonRE1[m] = 0;
                        nattemptnonRE[m] = 0; nattemptnonRE1[m] = 0;
                    }
                }
            }
        }
    }
    
    timer.step("");
    NumericVector res(timer);
    Rprintf("Total run time = %.2f\n", res[timer_cnt] / 1e9);
    
    //return output
    return(output);
}

