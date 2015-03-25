library(Rcpp)
library(coda)

#source MCMC code in
sourceCpp("logistic.cpp")

#function to run model
run.mcmc <- function(dat, response, inits = NA, inits_sigma = NA, nchains = 2, niter = 200000, scale = 0.05, varselect = F, ninitial = 10, priorvar = 10000, random = c("fixed", "globrand", "locrand"), nitertraining = NA)
{
	#ensure response variable is in first row of data set
	respind <- which(response == colnames(dat))
	stopifnot(length(respind) > 0)
	if(respind != 1) dat <- cbind(dat[, respind], dat[, -respind])
	
	#check names don't have underscores
	if(length(grep(glob2rx("*_*"), colnames(dat))) > 0) stop("Can't have variable names with underscores")
	
	#extract original number of variables
	orignpars <- ncol(dat) - 1
	
    #convert data frame into correct format for use
    #in Bayesian model
    dat <- createLinear(dat, response)
    
    #extract components
    factindex <- as.numeric(table(dat$factindex))
    cumfactindex <- c(0, cumsum(factindex)[-length(factindex)]) + 1
    varsorig <- dat$vars
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

#function to extract conditional runs
extract.cond <- function(output)
{
    x <- lapply(output[[1]], function(x)
    {
        x <- as.matrix(x)
        vars <- colnames(x)
        
        #separate indicators
        varind <- grep(glob2rx("I_*"), vars)
        inds <- x[, varind]
        x <- x[, -varind]
        vars <- vars[-varind]
        
        #separate sigmas
        sigmaind <- grep(glob2rx("sigma_*"), vars)
        sigma <- x[, sigmaind]
        x <- x[, -sigmaind]
        vars <- vars[-sigmaind]
        
        #extract post
        postind <- which(vars == "post")
        x <- x[, -postind]
        vars <- vars[-postind]
        
        #extract intercept
        intind <- which(vars == "Intercept")
        x <- x[, -intind]
        vars <- vars[-intind]
        
        #extract conditional runs
        origvars <- vars
        vars <- sapply(strsplit(vars, "_"), function(x) x[1])
        varind <- sapply(strsplit(colnames(inds), "_"), function(x) x[2])
        x <- lapply(1:length(vars), function(i, vars, x, inds)
        {
            x <- x[, i]
            inds <- inds[, match(vars[i], varind)]
            x <- x[inds == 1]
            x
        }, vars = vars, x = x, inds = inds)
        list(x = x, origvars = origvars)
    })
    origvars <- x[[1]]$origvars
    x <- lapply(x, function(x) x$x)
    
    #now merge runs together
    x1 <- lapply(1:length(x[[1]]), function(i, x) lapply(x, function(x, i) x[[i]], i = i), x = x)
    x1 <- lapply(1:length(x1), function(i, x1)
    {
        x <- x1[[i]]
        x <- lapply(1:length(x), function(j, x) cbind(x[[j]], rep(j, length(x[[j]]))), x = x)
        x
    }, x1 = x1)
    x1 <- lapply(x1, function(x) do.call("rbind", x))
    names(x1) <- origvars
    x1
}        

#function to plot output
plot.varselect <- function(output, cond = F, trace = T, density = T, ask = T, ...)
{
    if(cond == T & output$varselect == F) cond <- F
    
    if(cond)
    {
        #plot conditional traces
        output <- extract.cond(output)
        # set up plotting parameters
        if(trace & density)
        {
            maxp <- length(output) * 2
            if (maxp >= 4) mfrow1 <- c(4, 2)
            else mfrow1 <- c(maxp, 2)
        }
        else
        {
            maxp <- length(output)
            if(maxp >= 4) mfrow1 <- c(4, 1)
            else mfrow1 <- c(maxp, 1)
        }
        cols1 <- c("black", "red", "blue", "green", "yellow", "purple")
        par(mfrow = mfrow1)
        # produce plots
        l <- 1
        for (i in 1:length(output))
        {
            # plot trace
            temp <- output[[i]]
            if(nrow(temp) > 1)
            {
                if(trace)
                {
                  plot(1:nrow(temp), temp[, 1], type = "n", main = names(output)[i], xlab = "Index", ylab = "Value")
                  segments(x0 = 1:(nrow(temp) - 1), y0 = temp[-nrow(temp), 1], x1 = 2:nrow(temp), y1 = temp[-1, 1], col = cols1[temp[, 2]])
                }
                if(density)
                {
                  # plot density
                  plot(density(temp[, 1]), main = names(output)[i], xlab = "Value", ylab = "Density")
                }
            }
            else
            {
                if(trace)
                {
                  plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                    main = names(output)[i])
                  text(0, 0, "<= 1 Sample")
                }
                if(density)
                {
                  plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                    main = names(output)[i])
                  text(0, 0, "<= 1 Sample")
                }
            }
            l <- l + 1
            if (ask == T)
            {
                if ((l - 1) %% 4 == 0 && (l - 1) != maxp) 
                readline("Press any key to continue:")
            }
        }
        par(mfrow = c(1, 1))
        return(cat(""))
    }
    else plot(output[[1]], trace = trace, density = density, ask = ask, ...)
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
        if(!is.matrix(PPAvar)) PPAvar <- matrix(PPAvar, nrow = 1)
        rownames(PPAvar) <- sapply(strsplit(colnames(indicators[[1]]), "_"), function(x) x[2])
        PPAvar <- t(apply(PPAvar, 1, function(x)
        {
            y <- x
            x <- x[!is.na(x)]
            x <- c(y, mean(x), sd(x))
            x <- signif(x, digits = 2)
            x <- c(as.character(x[-length(x)]), paste0("(", x[length(x)], ")"))
        }))
        colnames(PPAvar) <- c(paste("Chain", 1:(ncol(PPAvar) - 2)), "Mean", "SD")
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
        var <- colnames(indicators[[1]])
        var <- sapply(strsplit(var, "_"), function(x) x[2])
        rownames(indicators1) <- c(var, paste("Chain", 1:nchains), "Mean", "SD")
        indicators1[indicators1 == "0"] <- ""
        print(indicators1, quote = F)
    }
}

# window function for 'bayesord' objects
window.varselect <- function(x, start = NA, end = NA, thin = NA, chains = NA, ...) {
    y <- x
    if (!is.na(chains)) {
        x <- x[[1]][chains]
    }
    else x <- x[[1]]
    
    if(is.na(start)) {
        if(is.na(end)) {
            if(!is.na(thin)) x <- window(x, thin = thin)
        }
        else {
            if(!is.na(thin)) x <- window(x, end = end, thin = thin)
            else x <- window(x, end = end)
        }
    }
    else {
        if(is.na(end)) {
            if(!is.na(thin)) x <- window(x, start = start, thin = thin)
            else x <- window(x, start = start)
        }
        else {
            if(!is.na(thin)) x <- window(x, start = start, end = end, thin = thin)
            else x <- window(x, start = start, end = end)
        }
    }
    y[[1]] <- x
    y
} 
