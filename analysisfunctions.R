library(Rcpp)
library(coda)

#source MCMC code in
sourceCpp("logistic.cpp")

#function to run model
run.mcmc <- function(dat, response, inits = NA, inits_sigma = NA, nchains = 2, n.iter = 200000, scale = 0.05, varselect = F, ninitial = 100, priorvar = 10000, random = c("fixed", "globrand", "locrand"))
{
	#ensure response variable is in first row of data set
	respind <- which(response == colnames(dat))
	stopifnot(length(respind) > 0)
	if(respind != 1) dat <- cbind(dat[, respind], dat[, -respind])
	
	#extract original number of variables
	orignpars <- ncol(dat) - 1
	
    #convert data frame into correct format for use
    #in Bayesian model
    dat <- createLinear(dat, response)
    
    #extract components
    factindex <- as.numeric(table(dat$factindex))
    cumfactindex <- c(0, cumsum(factindex)[-length(factindex)]) + 1
    vars <- dat$vars
    dat <- dat$data
	
	#convert data to matrix
	dat <- as.matrix(dat)
	
	#extract number of parameters
	npars <- ncol(dat) - 1
	
	#set random effect indicator
	stopifnot(!is.na(match(random, c("fixed", "globrand", "locrand"))))
	random <- ifelse(random == "fixed", 0, ifelse(random == "globrand", 1, 2))
	
	#generate initial values
	if(!is.list(inits))
	{
        gen_inits <- 1
        inits <- list(nchains)
        for(j in 1:nchains) inits[[j]] <- rep(1, ncol(dat))
        inits_sigma <- list(nchains)
        for(j in 1:nchains) inits_sigma[[j]] <- rep(1, ncol(dat) - 1)
    }
	else
	{
	    if(length(inits) != nchains | length(inits_sigma) != nchains) stop("List of initial values doesn't match number of chains")
	    else
	    {
	        if(!all(sapply(inits, length) == ncol(dat))) stop("Wrong number of initial values")
	        if(random == 2)
	        {
	            if(!all(sapply(inits_sigma, length) == (ncol(dat) - 1))) stop("Wrong number of initial sigma values")
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
	
	#run model
	model.sim <- list(NULL)
	for(j in 1:nchains)
	{
        model.sim[[j]] <- logisticMH(dat, factindex, cumfactindex, inits[[j]], inits_sigma[[j]], gen_inits, priors, n.iter, scale, orignpars, ifelse(varselect, 1, 0), ninitial, random)
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
            model.sim[[j]] <- model.sim[[j]][, -c(npars + 2, 2 * npars + 3)]
            colnames(model.sim[[j]]) <- c("Intercept", vars, paste0("sigma_", vars), paste0("I_", vars), "post")
        }
        if(!varselect) model.sim[[j]] <- model.sim[[j]][, -grep(glob2rx("I_*"), colnames(model.sim[[j]]))]
        model.sim[[j]] <- as.mcmc(model.sim[[j]])
    }
    
	#return output
	model.sim <- as.mcmc.list(model.sim)
	model.sim <- list(model.sim = model.sim, varselect = varselect)
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
    data <- data[, y != response, drop = F]
    y <- y[y != response]
    
    #record factor index
    factindex <- NULL
    vars <- NULL
    for(i in 1:length(y))
    {
        if(is.numeric(data[, i]))
        {
            data1 <- cbind(data1, data[, i])
            colnames(data1)[ncol(data1)] <- y[i]
            if(is.null(factindex))
            {
                factindex <- i
                vars <- y[i]
            }
            else
            {
                factindex <- c(factindex, i)
                vars <- c(vars, y[i])
            }
        }
        if(is.factor(data[, i]))
        {
            #expand categorical variables
            for(j in 2:length(levels(data[, i])))
            {
                data1 <- cbind(data1, ifelse(data[, i] == levels(data[, i])[j], 1, 0))
                colnames(data1)[ncol(data1)] <- paste0(y[i], j)
            }
            if(is.null(factindex))
            {
                factindex <- rep(i, length(levels(data[, i])) - 1)
                for(j in 2:length(levels(data[, i]))) vars <- paste0(y[i], "_", levels(data[, i])[j])
            }
            else
            {
                factindex <- c(factindex, rep(i, length(levels(data[, i])) - 1))
                for(j in 2:length(levels(data[, i]))) vars <- c(vars, paste0(y[i], "_", levels(data[, i])[j]))
            }
        }
    }
    return(list(data = data1, factindex = factindex, vars = vars))
}

#function to summarise output
summary.varselect <- function(output, topmodels = 5, ...)
{
    #extract relevant quantities from object
    varselect <- output$varselect
    output <- output$model.sim
    
    #thin chains
    output <- window(output, ...)
    
    #now extract indicators
    if(varselect == 0) print(summary(output))
    else
    {
        nchains <- length(output)
        output <- lapply(output, as.matrix)
        indicators <- lapply(output, function(x) x[, grep(glob2rx("I_*"), colnames(x)), drop = F])
        output <- lapply(output, function(x) x[, -grep(glob2rx("I_*"), colnames(x)), drop = F])
        output <- lapply(output, as.mcmc)
        output <- as.mcmc.list(output)
        #print posterior summaries
        cat("###########    POSTERIOR SUMMARIES    ###########\n")
        print(summary(output))
        #print PPAs for variables
        cat("\n###########    PPAs for VARIABLES    ###########\n")
        PPAvar <- sapply(indicators, function(x) apply(x, 2, mean))
        if(!is.matrix(PPAvar))
        {
            PPAvar <- matrix(PPAvar, nrow = 1)
            rownames(PPAvar) <- colnames(indicators[[1]])
        }
        PPAvar <- t(apply(PPAvar, 1, function(x)
        {
            x <- x[!is.na(x)]
            x <- c(mean(x), sd(x))
            x <- signif(x, digits = 2)
            x <- c(as.character(x[1]), paste0("(", x[2], ")"))
        }))
        colnames(PPAvar) <- c("Mean", "SD")
        print(t(PPAvar), quote = F)
        #print PPAs for models
        cat("\n###########    PPAs for MODELS    ###########\n")
        indicators1 <- lapply(indicators, function(x) apply(x, 1, function(y) paste(y, collapse = ".")))
        indicators2 <- lapply(indicators1, table)
        indicators2 <- lapply(indicators2, function(x) x / sum(x))
        #match across chains
        ind <- unique(do.call("c", indicators1))
        indmatch2 <- lapply(indicators2, function(x, ind) match(ind, names(x)), ind = ind)
        indicators2 <- lapply(1:length(indicators2), function(i, x, indmatch) x[[i]][indmatch[[i]]], indmatch = indmatch2, x = indicators2)
        indicators2 <- do.call("cbind", indicators2)
        rownames(indicators2) <- ind
        indicators2 <- cbind(indicators2, t(apply(indicators2, 1, function(x)
        {
            x <- x[!is.na(x)]
            c(mean(x), sd(x))
        })))
        #sort by average
        topmodels <- ifelse(topmodels > nrow(indicators2), nrow(indicators2), topmodels)
        indicators2 <- indicators2[sort.list(indicators2[, nchains + 1], decreasing = T), ][1:topmodels, ]
        indicators1 <- rownames(indicators2)
        indicators1 <- sapply(indicators1, function(y) strsplit(y, "\\.")[[1]])
        indicators2 <- apply(indicators2, 1, function(x)
        {
            x <- signif(x, digits = 2)
            c(as.character(x[-length(x)]), paste0("(", x[length(x)], ")"))
        })
        indicators1 <- rbind(indicators1, indicators2)
        colnames(indicators1) <- paste0("M", 1:ncol(indicators1))
        var <- colnames(output[[1]])[-1]
        if(length(grep(glob2rx("sigma*"), var)) > 0) var <- var[-grep(glob2rx("sigma*"), var)]
        var <- var[-c(length(var))]
        rownames(indicators1) <- c(var, paste("Chain", 1:nchains), "Mean", "SD")
        indicators1[indicators1 == "0"] <- ""
        print(indicators1, quote = F)
    }
}
    
