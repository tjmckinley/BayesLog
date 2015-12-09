#function to extract correct variables from data set and generate
#design matrix

extractData <- function(formula, dat, agg = TRUE)
{
    #check intercept is present
    form <- terms(formula)
    stopifnot(attributes(form)$intercept == 1)
    
    # extract formula
    rand <- lme4:::findbars(formula)
    if(!is.null(rand))
    {
        rand <- sapply(rand, all.vars)
        #check that hierarchical terms are factors and
        #then convert to integers
        stopifnot(all(!is.na(match(rand, colnames(dat)))))
        temp <- dat[, match(rand, colnames(dat)), drop = F]
        dat <- dat[, -match(rand, colnames(dat)), drop = F]
        for(j in 1:ncol(temp))
        {
            stopifnot(is.factor(temp[, j]))
            temp[, j] <- as.numeric(temp[, j])
            temp[, j] <- temp[, j] - 1
        }
        dat <- cbind(dat, temp)
        nrand <- length(rand)
        
        #reformulate formula to treat extract from dataset in correct way
        formula <- terms(formula)
        formula <- drop.terms(formula, match(rand, all.vars(formula)) - 1, keep.response = T)
        formula <- reformulate(attr(formula, "term.labels"), all.vars(formula)[1])
        formula <- update.formula(formula, formula(paste("~ . +", paste(rand, collapse = "+")))) 
    }
    else nrand <- 0
	
	#check data names
    mf <- model.frame(formula = formula, data = dat, na.action = na.fail)
    
    #check binary response
    mf.resp <- mf[[1]]
	stopifnot(length(table(mf.resp)) == 2)
	
	#convert rest of data into design matrix
    mf <- model.matrix(formula, mf)
    
    # check there are no columns called 'rand' or 'nsamples'
    temp <- match("nsamples", attr(mf, "names"))
    if (length(temp[!is.na(temp)]) > 0) stop("Can't name variable 'nsamples'")
    
    if(agg)
    {
        # now aggregate data to speed code up
        mf <- as.data.frame(mf)
        mf <- cbind(resp = mf.resp, mf)
        mf <- aggregate(rep(1, nrow(mf)), mf, table)
        mf[, ncol(mf)] <- as.numeric(mf[, ncol(mf)])
        nsamples <- mf[, ncol(mf)]
        mf <- mf[, -ncol(mf)]
	
	    #convert back to matrix for running C++ code
	    mf <- as.matrix(mf)
	}
	else nsamples <- rep(1, nrow(mf))
	
	#extract hierarchical terms if necessary
	if(nrand > 0)
	{
	    temp <- (ncol(mf) - nrand + 1):ncol(mf)
	    mf_rand <- mf[, temp, drop = FALSE]
	    mf <- mf[, -temp, drop = FALSE]
	}
	else mf_rand <- matrix(NA, 1, 1)
	
	return(list(mf = mf, mf_rand = mf_rand, formula = formula, nrand = nrand, nsamples = nsamples))
}