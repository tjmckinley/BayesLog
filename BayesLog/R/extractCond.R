#function to extract conditional runs from MCMC output
extractCond <- function(output)
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
