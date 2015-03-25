#' @title Fits Bayesian logistic regression model
#' @description Runs MCMC algorithm for fitting a logistic regression model in a Bayesian framework, includes options for performing variable selection.
#' @export
#' @import coda Rcpp
#' @param dat 		    Data frame containing data.
#' @param response      Character containing name of response variable.
#' @param gen_inits     Logical stating whether initial values are to be generated at random or
#'                      input by the user. If the latter then \code{inits} and \code{inits_sigma}
#'                      must be supplied.
#' @param inits         List containing vectors of initial values for the regression parameters.
#' @param inits_sigma   List containing vectors of initial values for the regression SD hyperparameters.
#' @param nchains       Number of chains to run.
#' @param niter         Number of iterations to run per chain.
#' @param scale         Mixing proportion for adapt proposal.
#' @param varselect     Logical stating whether variable selection is to be used.
#' @param ninitial      Number of iterations before adaptive proposal kicks in.
#' @param priorvar      Prior variance for fixed effects (currently scalar).
#' @param random        Character: taking the values "fixed", "globrand" or "locrand" to denote fixed
#'                      effects, a global random effect or local random effects.
#' @param nitertraining the number of iterations for the training run (if required)
bayesLog <- function(dat, response, gen_inits = TRUE, inits = NA, inits_sigma = NA, nchains = 2, niter = 200000, scale = 0.05, varselect = FALSE, ninitial = 10, priorvar = 10000, random = c("fixed", "globrand", "locrand"), nitertraining = NA)
{
    #check inputs
    stopifnot(is.data.frame(dat), is.character(response), length(response) == 1)
    stopifnot(is.list(inits) | is.na(inits[1]), is.list(inits_sigma) | is.na(inits[1]))
    stopifnot(all(sapply(inits, is.vector)) | gen_inits, all(sapply(inits, is.vector)) | gen_inits)
    stopifnot(is.numeric(nchains), length(nchains) == 1, is.numeric(niter), length(niter) == 1)
    stopifnot(is.numeric(scale), length(scale) == 1, scale > 0, scale < 1)
    stopifnot(is.logical(varselect), length(varselect) == 1)
    stopifnot(is.numeric(ninitial), length(ninitial) == 1, is.numeric(priorvar), length(priorvar) == 1)
    stopifnot(is.character(random), all(!is.na(match(random[1], c("fixed", "globrand", "locrand")))))
    
	#ensure response variable is in first row of data set
	stopifnot(ncol(dat) >= 2)
	respind <- which(response == colnames(dat))
	stopifnot(length(respind) > 0)
	if(respind != 1) dat <- cbind(dat[, respind], dat[, -respind])
	
	#check data
	stopifnot(is.factor(dat[, 1]))
	stopifnot(length(levels(dat[, 1])) == 2)
	for(j in 2:ncol(dat)) stopifnot(is.factor(dat[, j]) | is.numeric(dat[, j]) | is.logical(dat[, j]))
	
	#check names don't have underscores
	if(length(grep(glob2rx("*_*"), colnames(dat))) > 0) stop("Can't have variable names with underscores")
	
	#extract original number of variables
	orignpars <- ncol(dat) - 1
	
    #convert data frame into correct format for use
    #in Bayesian model
    dat <- createLinear(dat, response)
    
    #extract components necessary for MCMC
    factindex <- as.numeric(table(dat$factindex))
    cumfactindex <- c(0, cumsum(factindex)[-length(factindex)]) + 1
    varsorig <- dat$vars
    dat <- dat$data
    
    #remove missing values
    dat <- dat[!is.na(apply(dat, 1, sum)), ]
	
	#convert data to matrix
	dat <- as.matrix(dat)
	
	#extract number of parameters
	npars <- ncol(dat) - 1
	
	#set random effect indicator
	random <- ifelse(random[1] == "fixed", 0, ifelse(random[1] == "globrand", 1, 2))
	
	#generate initial values
	if(gen_inits)
	{
        gen_inits <- 1
        inits <- list(nchains)
        for(j in 1:nchains) inits[[j]] <- rep(1, ncol(dat))
        inits_sigma <- list(nchains)
        for(j in 1:nchains) inits_sigma[[j]] <- rep(1, ncol(dat))
    }
	else
	{
	    if(length(inits) != nchains | length(inits_sigma) != nchains) stop("List of initial values doesn't match number of chains")
	    else
	    {
	        if(!all(sapply(inits, length) == ncol(dat))) stop("Wrong number of initial values")
	        if(random == 2)
	        {
	            if(!all(sapply(inits_sigma, length) == ncol(dat))) stop("Wrong number of initial sigma values")
	        }
	        if(random == 1)
	        {
	            if(!all(sapply(inits_sigma, length) == 1)) stop("Wrong number of initial sigma values")
	        }
	        gen_inits <- 0
	    }
	}
	
	#set priors
	priors <- matrix(c(rep(c(0, priorvar), times = npars + 1), 0, 20), ncol = 2, byrow = T)
	
	#set training run if required
	if(is.na(nitertraining) || !varselect) nitertraining <- 0
	if(nitertraining >= niter) stop("'nitertraining' must be less than 'niter'")
	
	#run model
	model.sim <- list(NULL)
	for(j in 1:nchains)
	{
        vars <- varsorig
        model.sim[[j]] <- logisticMH(dat, factindex, cumfactindex, inits[[j]], inits_sigma[[j]], gen_inits, priors, niter, nitertraining, scale, orignpars, ifelse(varselect, 1, 0), ninitial, random)
        if(random == 0)
        {
            model.sim[[j]] <- model.sim[[j]][, -c(npars + 2)]
            colnames(model.sim[[j]]) <- c("Intercept", vars, paste0("I_", vars), "post")
        }
        if(random == 1)
        {
            model.sim[[j]] <- model.sim[[j]][, -c(npars + 3)]
            colnames(model.sim[[j]]) <- c("Intercept", vars, "sigma", paste0("I_", vars), "post")
        }
        if(random == 2)
        {
            colnames(model.sim[[j]]) <- c("Intercept", vars, paste0("sigma_", c("int", vars)), paste0("I_", c("int", vars)), "post")
            model.sim[[j]] <- model.sim[[j]][, -c(npars + 2, 2 * npars + 3)]
        }
        if(!varselect) model.sim[[j]] <- model.sim[[j]][, -grep(glob2rx("I_*"), colnames(model.sim[[j]]))]
        else
        {
            #concatenate outputs from indicator variables
            vars <- colnames(model.sim[[j]])
            varind <- grep(glob2rx("I_*"), vars)
            vars <- vars[varind]
            vars <- sapply(strsplit(vars, "_"), function(x) x[2])
            colnames(model.sim[[j]])[varind] <- paste0("I_", vars)
            model.sim[[j]] <- model.sim[[j]][, -varind[duplicated(vars)]]
        }   
        model.sim[[j]] <- as.mcmc(model.sim[[j]])
    }
    
	#return output
	model.sim <- as.mcmc.list(model.sim)
	model.sim <- list(model.sim = model.sim, varselect = varselect)
	class(model.sim) <- "bayesLog"
	model.sim
}
