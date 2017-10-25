#' @title Fits Bayesian logistic regression model
#'
#' @description Runs MCMC algorithm for fitting a logistic regression model in a Bayesian framework, includes options for performing variable selection.
#'
#' @export
#'
#' @param formula       Formula for linear regression
#' @param dat 		    Data frame containing data.
#' @param nsamples      Character specifying names of column in data to use as number of samples.
#' @param inits         List containing vectors of initial values for the regression
#' parameters. If missing then these are simulated.
#' @param priorvar      A numeric specifying the prior variance for all regression terms.
#' @param prior_rand_ub An upper bound on the prior standard deviation for the hierarchical terms.
#' @param nchains       Number of chains to run.
#' @param niter         Number of iterations to run per chain.
#' @param ninitial      Number of iterations before adaptive proposal kicks in
#' @param scale         Mixing proportion for adapt proposal.
#' @param nadapt        Number of iterations between each adaptive proposal update.
#' @param nprintsum     how often to print run time information to the screen
#' @param maxscale      the maximum scaling of the adaptive proposal variance at each update
#' @param niterdim      the iteration at which the diminishing adaptation component kicks in
#' @param blocks        a list of integer vectors defining parameters to update in blocks.
#' Missing parameters are filled in automatically as componentwise updates.
#' @param noncentreint  a vector of variable names for non-centring relative to the intercept.
#' If missing, then only centred updates performed. If \code{noncentreint == "all"} then all
#' updates are non-centred.
#'
#' @return An object of class \code{bayesLog}, which is basically a list
#' including a subset of elements:
#' \itemize{
#' \item{post:}{ an \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}} containing the posterior samples.}
#' \item{formula:}{ a \code{formula} used to define the model.}
#' \item{data:}{ a \code{data.frame} containing the data used to fit the model.}
#' \item{nregpars:}{ a scalar containing the number of regression parameters.}
#' \item{nrand:}{ a scalar containing the number of random effect variables.}
#' \item{nrandlevels:}{ a \code{numeric} containing the number of levels of each random effect.}
#' }
#'

bayesLog <- function(formula, dat, nsamples = NA, inits = NA, priorvar = 1, prior_rand_ub = 20, nchains = 2, niter = 200000, ninitial = 100, 
                     scale = 0.05, nadapt = 100, nprintsum = 1000, maxscale = 2, niterdim = 1000, blocks = NA, noncentreint = NA)
{    
    #check inputs
    stopifnot(is.data.frame(dat))
	for(j in 1:ncol(dat)) stopifnot(is.factor(dat[, j]) | is.numeric(dat[, j]) | is.logical(dat[, j]))
	
	stopifnot(is.na(nsamples) | (length(nsamples) == 1 & is.character(nsamples)))
	if(!is.na(nsamples)){
		stopifnot(length(which(colnames(dat) == nsamples)) == 1)
		form <- attributes(terms(formula))$term.labels
		form <- which(form == nsamples)
		if(length(form) != 0) {
			stop("'nsamples' term also involved in formula.")
		}
		colnames(dat)[which(colnames(dat) == nsamples)] <- "nsamples"
	} else {
		stopifnot(length(which(colnames(dat) == "nsamples")) == 0)
		form <- attributes(terms(formula))$term.labels
		form <- which(form == "nsamples")
		if(length(form) != 0) {
			stop("'nsamples' term also involved in formula.")
		}
		dat$nsamples <- rep(1, nrow(dat))
	}
    
    stopifnot(is.list(inits) | is.na(inits[1]))
    stopifnot(all(sapply(inits, is.vector)) | is.na(inits[1]))
    
    stopifnot(is.numeric(priorvar) & length(priorvar) == 1)
    stopifnot(priorvar > 0)
    
    stopifnot(is.numeric(prior_rand_ub) & length(prior_rand_ub) == 1)
    stopifnot(prior_rand_ub > 0)
    
    stopifnot(is.numeric(nchains) & length(nchains) == 1 & abs(floor(nchains) - nchains) < .Machine$double.eps ^ 0.5)
    stopifnot(is.numeric(niter) & length(niter) == 1 & abs(floor(niter) - niter) < .Machine$double.eps ^ 0.5)
    
    stopifnot(is.numeric(scale) & length(scale) == 1, scale > 0, scale < 1)
    stopifnot(is.numeric(nadapt) & length(nadapt) == 1 & abs(floor(nadapt) - nadapt) < .Machine$double.eps ^ 0.5)
    
    stopifnot(is.numeric(nprintsum) & length(nprintsum) == 1& abs(floor(nprintsum) - nprintsum) < .Machine$double.eps ^ 0.5)
    stopifnot(is.numeric(maxscale) & length(maxscale) == 1)
    stopifnot(maxscale > 0)
    stopifnot(is.numeric(niterdim) & length(niterdim) == 1 & abs(floor(niterdim) - niterdim) < .Machine$double.eps ^ 0.5)
    stopifnot(niterdim < niter)
    stopifnot(nprintsum %% nadapt == 0)
    
    stopifnot(is.list(blocks) | is.na(blocks[1]))
    stopifnot(is.character(noncentreint) | is.na(noncentreint[1]))
    
    #save formula and data set for later output
    origformula <- formula
    origdat <- dat
    
	#extract data in correct format
	temp <- extractData(formula, dat)
	mf <- temp$mf
	mf_rand <- temp$mf_rand
	formula <- temp$formula
	nrand <- temp$nrand
	nsamples <- temp$nsamples
	randnames <- temp$randnames
	origdat <- temp$origdat
	rm(temp)
	
	#extract number of parameters
	npars <- ncol(mf) - 1
	
	if(!is.na(blocks[1]))
	{
        stopifnot(all(sapply(blocks, function(x) all(!is.na(x)) & all(x > 0) & all(abs(floor(x) - x) < .Machine$double.eps ^ 0.5) & all(x <= npars))))
#        stopifnot(length(unique(do.call("c", blocks))) == length(do.call("c", blocks)))
        inds <- 1:npars
        inds <- inds[-match(unique(do.call("c", blocks)), inds)]
        if(length(inds) > 0)
        {
            inds <- as.list(inds)
            blocks <- c(blocks, inds)
        }
        nblocks <- sapply(blocks, length)
    }
    else
    {
        blocks <- as.list(1:npars)
        nblocks <- sapply(blocks, length)
    }
    blocks <- lapply(blocks, function(x) x - 1)
    
    #generate prior matrix
    priors <- matrix(rep(c(0, priorvar), npars), ncol = 2, byrow = T)
    
    #now append hierarchical variances
	if(nrand > 0)
	{
	    priors <- rbind(priors, matrix(rep(c(0, prior_rand_ub), nrand), ncol = 2, byrow = T))	    
	    #now split the data according to indexes relating to each level of factor
	    randindexes <- list(NULL)
	    for(i in 1:nrand)
	    {
	        temp <- split(1:nrow(mf_rand), factor(mf_rand[, i]))
	        temp <- lapply(temp, function(x) x - 1)
	        randindexes[[i]] <- temp
	    }
	    npars <- npars + nrand
	}
	else randindexes <- list(NULL)
	
	#generate initial values
	if(is.na(inits[1]))
	{
        inits <- list(NULL)
        for(j in 1:nchains)
        {
            inits[[j]] <- rnorm(npars, 0, 1)
            if(nrand > 0) inits[[j]][(npars - (nrand - 1)):npars] <- runif(nrand, 0, prior_rand_ub)
        }
    }
	else
	{
	    if(length(inits) != nchains) stop("List of initial values doesn't match number of chains")
	    else
	    {
	        if(!all(sapply(inits, length) == npars)) stop("Wrong number of initial values")
	    }
	}
	
	#extract variable names
	varnames <- colnames(mf)[-1]
    #set variable names
    if(nrand > 0)
    {
        sdnames <- paste0(colnames(mf_rand), "sd")
        randnames <- lapply(1:nrand, function(i, var, nam)
        {
            var <- var[i]
            randnames <- randnames[[i]]
            paste0(var, randnames)
        }, var = colnames(mf_rand), nam = randnames)
        randnames <- do.call("c", randnames)
            
        fullnames <- c(varnames, sdnames, randnames, "logposterior")
    }
    else fullnames <- c(varnames, "logposterior")
	
	# set number of levels in random intercepts terms
	nrandlevels <- sapply(randindexes, length)
	
	if(!is.na(noncentreint[1]))
	{ 
	  #extract all variables in model
	  predictors <- all.vars(origformula)
	  
	  #remove response 
	  predictors <- predictors[-match(as.character(origformula)[2], predictors)]
	  
	  if(noncentreint[1] == "all") noncentreint <- predictors
	  
	  #cross reference these against required variable names
	  noncentreint <- unique(noncentreint)
	  temp <- match(noncentreint, predictors)
	  if(!all(!is.na(temp))) stop("Variable names in required subset for non-centring do not match variables in model")
	  
	  # split j into fixed and random effects
	  randj <- noncentreint[!is.na(match(noncentreint, colnames(mf_rand)))]
	  fixj <- noncentreint[is.na(match(noncentreint, colnames(mf_rand)))]
	  
	  # if valid, then extract correct columns corresponding to columns
	  if(length(fixj) > 0)
	  {
	    fixlistj <- lapply(fixj, function(x, names, assign)
	    {
	      x <- match(x, names)
	      which(assign == x)
	    }, names = predictors, 
	    assign = attr(model.matrix(lme4::nobars(origformula), dat), "assign"))
	    
	    noncentreint <- lapply(fixlistj, function(x)
	    {
	      temp <- rep(0, npars)
	      temp[x] <- 1
	      temp
	    })
	  }
	  else noncentreint <- list(rep(0, npars))
	    
	  if(length(randj) > 0) noncentreintRE <- match(randj, colnames(mf_rand))
	  else noncentreintRE <- 0
	}
	else
	{
	  noncentreint <- list(rep(0, npars))
	  noncentreintRE <- 0
	}
	
	#run model
	model.sim <- list(NULL)
	for(j in 1:nchains)
	{
        model.sim[[j]] <- logisticMH(mf, nsamples, inits[[j]], priors, niter, ninitial, scale, nadapt, nprintsum, maxscale, niterdim, nrand, randindexes, mf_rand, nblocks, blocks, ifelse(j > 1, 0, 1), noncentreint, noncentreintRE)
        
        #set variable names
        colnames(model.sim[[j]]) <- fullnames
        
        #convert to MCMC object
        model.sim[[j]] <- as.mcmc(model.sim[[j]])
    }
    
	#return output
	if(nchains > 1) model.sim <- as.mcmc.list(model.sim)
	else model.sim <- model.sim[[1]]
	
	#append data and formula to object
	model.sim <- list(post = model.sim, formula = origformula, data = origdat, nregpars = length(varnames), nrand = nrand, nrandlevels = nrandlevels, randnames = colnames(mf_rand))
	
	#set class
	class(model.sim) <- "bayesLog"
	model.sim
}
