#' @title Fits Bayesian logistic regression model
#'
#' @description Runs MCMC algorithm for fitting a logistic regression model in a Bayesian framework, includes options for performing variable selection.
#'
#' @export
#' @import coda Rcpp
#'
#' @param formula       Formula for linear regression
#' @param dat 		    Data frame containing data.
#' @param gen_inits     Logical stating whether initial values are to be generated at random or
#'                      input by the user. If the latter then \code{inits}
#'                      must be supplied.
#' @param inits         List containing vectors of initial values for the regression parameters.
#' @param priorvar      A numeric specifying the prior variance for all regression terms.
#' @param nchains       Number of chains to run.
#' @param niter         Number of iterations to run per chain.
#' @param scale         Mixing proportion for adapt proposal.
#' @param nadapt        Number of iterations between each adaptive proposal update.
#' @param nprintsum     how often to print run time information to the screen
#' @param maxscale      the maximum scaling of the adaptive proposal variance at each update
#' @param niterdim      the iteration at which the diminishing adaptation component kicks in
#'
#' @return An object of class \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}}.
#'

bayesLog <- function(formula, dat, gen_inits = TRUE, inits = NA, priorvar = 1, nchains = 2, niter = 200000, scale = 0.05, nadapt = 100, nprintsum = 1000, maxscale = 2, niterdim = 1000)
{
    # extract formula
    form <- extractTerms(formula)[[1]]
    
    #check inputs
    stopifnot(is.data.frame(dat))
    
    stopifnot(is.logical(gen_inits) & length(gen_inits) == 1)
    
    stopifnot(is.list(inits) | is.na(inits[1]))
    stopifnot(all(sapply(inits, is.vector)) | gen_inits)
    
    stopifnot(is.numeric(priorvar) & length(priorvar) == 1)
    stopifnot(priorvar > 0)
    
    stopifnot(is.numeric(nchains) & length(nchains) == 1 & abs(floor(nchains) - nchains) < .Machine$double.eps ^ 0.5)
    stopifnot(is.numeric(niter) & length(niter) == 1 & abs(floor(niter) - niter) < .Machine$double.eps ^ 0.5)
    
    stopifnot(is.numeric(scale) & length(scale) == 1, scale > 0, scale < 1)
    stopifnot(is.numeric(nadapt) & length(nadapt) == 1 & abs(floor(nadapt) - nadapt) < .Machine$double.eps ^ 0.5)
    
    stopifnot(is.numeric(nprintsum) & length(nprintsum) == 1& abs(floor(nprintsum) - nprintsum) < .Machine$double.eps ^ 0.5)
    stopifnot(is.numeric(maxscale) & length(maxscale) == 1)
    stopifnot(maxscale > 0)
    stopifnot(is.numeric(niterdim) & length(niterdim) == 1& abs(floor(niterdim) - niterdim) < .Machine$double.eps ^ 0.5)
    stopifnot(niterdim < niter)
    stopifnot(nprintsum %% nadapt == 0)
    
	#check data
	for(j in 1:ncol(dat)) stopifnot(is.factor(dat[, j]) | is.numeric(dat[, j]) | is.logical(dat[, j]))
	
	#check data names
    mf <- model.frame(formula = form, data = dat, na.action = na.fail)
    
    #check binary response
    mf.resp <- mf[[1]]
	stopifnot(length(table(mf.resp)) == 2)
	
	#convert rest of data into design matrix
    mf <- model.matrix(model.form, mf)
    
    # check there are no columns called 'RE' or 'nsamples'
    temp <- match(c("RE", "nsamples"), attr(mf, "names"))
    if (length(temp[!is.na(temp)]) > 0) stop("Can't name variables 'RE' or 'nsamples'")
    
    # now aggregate data to speed code up
    mf <- as.data.frame(mf)
    mf <- cbind(resp = mf.resp, mf)
    mf <- aggregate(rep(1, nrow(mf)), mf, table)
    mf[, ncol(mf)] <- as.numeric(mf[, ncol(mf)])
    nsamples <- mf[, ncol(mf)]
    mf <- mf[, -ncol(mf)]
	
	#convert back to matrix for running C++ code
	mf <- as.matrix(mf)
	
	#extract number of parameters
	npars <- ncol(mf) - 1
    
    #generate prior matrix
    priors <- matrix(rep(c(0, priorvar), npars), ncol = 2, byrow = T)
	
	#generate initial values
	if(gen_inits)
	{
        inits <- list(nchains)
        for(j in 1:nchains) inits[[j]] <- rnorm(npars, 0, 10)
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
	
	#run model
	model.sim <- list(NULL)
	for(j in 1:nchains)
	{
        model.sim[[j]] <- logisticMH(mf, nsamples, inits[[j]], priors, niter, scale, nadapt, nprintsum, maxscale, niterdim)
        #set variable names
        colnames(model.sim[[j]]) <- c(varnames, "logposterior")
        #convert to MCMC object
        model.sim[[j]] <- as.mcmc(model.sim[[j]])
    }
    
	#return output
	if(nchains > 1) model.sim <- as.mcmc.list(model.sim)
	else model.sim <- model.sim[[1]]
	model.sim
}
