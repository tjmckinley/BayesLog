library(Rcpp)
library(coda)

#source MCMC code in
sourceCpp("logistic.cpp")

#function to run model
run.mcmc <- function(dat, response, inits = NA, nchains = 2, n.iter = 50000, scale = 0.05)
{	
	#ensure response variable is in first row of data set
	respind <- which(response == colnames(dat))
	stopifnot(length(respind) > 0)
	if(respind != 1) dat <- cbind(dat[, respind], dat[, -respind])
	
	#convert to matrix
	dat <- as.matrix(dat)
	
	#extract number of parameters
	npars <- ncol(dat) - 1
	
	#generate initial values
	if(!is.list(inits))
	{
        gen_inits <- 1
        inits <- list(nchains)
        for(j in 1:nchains) inits[[j]] <- rep(1, ncol(dat) + 1)
    }
	else
	{
	    if(length(inits) != nchains) stop("List of initial values doesn't match number of chains")
	    else
	    {
	        if(!all(sapply(inits, length) == (ncol(dat) - 1))) stop("Wrong number of initial values")
	        gen_inits <- 0
	    }
	}
	
	#set priors
	priors <- matrix(c(rep(c(0, 1), times = npars + 1), 0, 20), ncol = 2, byrow = T)
	priors[1, 2] <- 10000
	
	#run model
	model.sim <- list(NULL)
	for(j in 1:nchains) model.sim[[j]] <- as.mcmc(logisticMH(dat, inits[[j]], gen_inits, priors, n.iter, scale))
	
	#return output
	model.sim <- as.mcmc.list(model.sim)
	model.sim
}

#function to output correct linear term for fast entry into WinBUGS file
createLinear <- function(data, response)
{
    y <- names(data)
    data1 <- data.frame(response = data[, y == response])
    colnames(data1) <- response
    if(is.factor(data1[, 1])) data1[, 1] <- as.numeric(data1[, 1]) - 1
    
    #remove response variable
    data <- data[, y != response]
    y <- y[y != response]
    
    for(i in 1:length(y))
    {
        if(is.numeric(data[, i]))
        {
            data1 <- cbind(data1, data[, i])
            colnames(data1)[ncol(data1)] <- y[i]
        }
        if(is.factor(data[, i]))
        {
            #expand categorical variables
            for(j in 2:length(levels(data[, i])))
            {
                data1 <- cbind(data1, ifelse(data[, i] == levels(data[, i])[j], 1, 0))
                colnames(data1)[ncol(data1)] <- paste0(y[i], j)
            }
        }
    }
    return(data1)
}
    
