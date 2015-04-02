#' Predict function for \code{bayesLog} objects
#' 
#' \code{predict} method for class \code{bayesLog}
#' 
#' Produce predictions for \code{\link{bayesLog}} objects.
#' 
#' @param object a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param newdata a data frame containing the new data at which to predict
#' @param \dots not used here
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{summary.mcmc}}
#'
#' @return A list of length \code{nchains}, with each element containing the 
#' posterior predictive samples
#'
#' @import coda
#' @export

predict.bayesLog <- function(object, newdata, ...)
{
    stopifnot(class(object) == "bayesLog")
    stopifnot(is.data.frame(newdata))
    
    #quick fudge to check data frames are equivalent
    test <- try(rbind(object$data[1, ], newdata[1, ]), silent = T)
    if(class(test) == "try-error") stop("Data frames not equivalent")
    
    # extract formula
    form <- extractTerms(object$formula)
    # extract random intercepts term if required
    RE <- form[[2]]
    form <- form[[1]]
    
    #check data names
    mf <- model.frame(formula = form, data = newdata, na.action = na.fail)
    
    # create vector for random intercepts if required
    if (is.null(RE)) randint <- NA 
    else
    {
        randint <- newdata[, match(RE, colnames(newdata)), drop = F]
        for(j in 1:ncol(randint)) if (!is.factor(randint[, j])) stop("Random intercepts term is not a factor")
        for(j in 1:ncol(randint)) randint[, j] <- as.numeric(randint[, j])
        randint <- randint[, 1]
        #currently only handles one random intercept term
        if(ncol(randint) > 1) print("CURRENTLY ONLY HANDLES A SINGLE RANDOM INTERCEPT VARIABLE")
        object.randint <- object$model.randint
    }
    object <- object$model.sim
    
    # condense data into succinct form for model
    newdata <- newdata[, match(attr(mf, "names"), colnames(newdata))]
        	
    #convert data frame into correct format to produce predictions
    newdata <- createLinear(newdata, colnames(newdata)[1])$data
    
    #extract posteriors in matrix form        
    nchains <- length(object)
    object <- lapply(object, as.matrix)
    #remove indicator variables if necessary
    object <- lapply(object, function(x)
    {
        cols <- grep(glob2rx("I_*"), colnames(x))
        if(length(cols) > 0) x <- x[, -cols, drop = F]
        x
    })
    #remove sigma terms if necessary
    object <- lapply(object, function(x)
    {
        cols <- grep(glob2rx("sigma_*"), colnames(x))
        if(length(cols) > 0) x <- x[, -cols, drop = F]
        x
    })
    #remove posterior term
    object <- lapply(object, function(x) x[, -ncol(x), drop = F])
    
    #extract random intercepts if necessary
    if(!is.null(RE)) object.randint <- lapply(object.randint, as.matrix)
    
    #now produce predictions
    newdata <- as.matrix(newdata)
    response <- newdata[, 1]
    newdata[, 1] <- rep(1, nrow(newdata))
    if(ncol(newdata) != ncol(object[[1]])) stop("data doesn't match outputs")
    
    #calculate predictions without random intercept
    pred <- lapply(object, function(beta, x) x %*% t(beta), x = newdata)
    
    #add random intercept if necessary
    if(!is.null(RE))
    {
        object.randint <- lapply(object.randint, function(x, randint) x[, randint], randint = randint)
        print("Need to check matrices are correct dimension")
        browser()
        pred <- lapply(1:length(pred), function(i, pred, rand)
        {
            pred <- pred[[i]]
            rand <- t(rand[[i]])
            pred + rand
        }, pred = pred, rand = object.randint)
    }
    
    #back-transform to probability scale
    pred <- lapply(pred, function(x) exp(x) / (1 + exp(x)))
    
    pred
}
