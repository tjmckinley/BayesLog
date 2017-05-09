#' Predict function for \code{bayesLog} objects
#' 
#' \code{predict} method for class \code{bayesLog}
#' 
#' Produce predictions for \code{\link{bayesLog}} objects.
#' 
#' @param object a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param newdata a data frame containing the new data at which to predict
#' @param type a character taking the value \code{"fit"} if posterior predictive
#' probabilities are required, and \code{"pred"} if binary predictions are required.
#' @param \dots not used here
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link[coda]{summary.mcmc}}
#'
#' @return A \code{bayesLog.pred} object, essentially a matrix of posterior 
#' predictive samples.
#'
#' @export

predict.bayesLog <- function(object, newdata, type = c("fit", "pred"), ...)
{
    stopifnot(class(object) == "bayesLog")
    stopifnot(is.data.frame(newdata))
    
    stopifnot(!is.na(match(type[1], c("fit", "pred"))))
    
    #extract formula, check and aggregate newdata
    origformula <- object$formula
    origdata <- object$data
    stopifnot(length(which(colnames(newdata) == "nsamples")) <= 1)
    if(length(which(colnames(newdata) == "nsamples")) == 0) {
    	newdata$nsamples <- rep(1, nrow(origdata))
    }
	temp <- extractData(object$formula, newdata, F)
	mf <- temp$mf
	mf_rand <- temp$mf_rand
	formula <- temp$formula
	nrand <- temp$nrand
	nsamples <- temp$nsamples
	newdata <- temp$origdata
	rm(temp)
    
    #extract posteriors in matrix form
    object <- as.matrix(object)
    #extract correct components if necessary
    if(nrand > 0)
    {
        #calculate number of levels of each random effect
        rand <- lme4:::findbars(origformula)
        rand <- sapply(rand, all.vars)
        origdata <- origdata[, match(rand, colnames(origdata)), drop = F]
        nlevels <- cumsum(sapply(origdata, function(x) length(levels(x))))
        #extract random effects
        rand <- object[, -(1:ncol(mf)), drop = F]
        #lose random effect variances
        rand <- rand[, -(1:nrand), drop = F]
        #lose posterior term
        rand <- rand[, -ncol(rand), drop = F]
        
        object <- object[, 1:ncol(mf), drop = F]
    }
    else object <- object[, -ncol(object), drop = F]
    
    #calculate predictions without random intercept
    pred <- mf %*% t(object)
    
    #add random intercepts if necessary
    if(nrand > 0)
    {
        j <- 1
        for(i in 1:nrand)
        {
            temp <- mf_rand[, i]
            temp <- temp + 1
            temp_rand <- rand[, j:nlevels[i], drop = F]
            temp <- temp_rand[, temp, drop = F]
            pred <- pred + t(temp)
            j <- nlevels[i] + 1
        }
    }
    
    #back-transform to probability scale
    pred <- exp(pred) / (1 + exp(pred))
    
    #perform predictive sampling
    if(type[1] == "pred") {
    	pred <- apply(pred, 1, function(x) rbinom(length(x), size = 1, prob = x))
    } else {
    	pred <- t(pred)
    }
    # return list including aggregation information
    pred <- list(pred = pred, nsamples = nsamples)
    class(pred) <- "bayesLog.pred"
    pred
}
