#' @title Generates posterior distributions for sensitivity
#' specificity, PPV and NPV
#'
#' @description Generates posterior distributions for sensitivity
#' specificity, PPV and NPV from a set of posterior predictions
#'
#' @export
#'
#' @param pred       Matrix of dimension (\code{npred} x \code{nobs}) containing predicted
#'                   values from \code{bayesLog} model.
#' @param obs 		 Vector of binary observations (of length \code{nobs})
#' @param thresh     Vector of thresholds between 0 and 1 at which to generate
#'                   classification estimates.
#'
#' @return An object of class \code{bayesLog.class}, which is a list
#' including a subset of elements: \code{sens}, \code{spec}, \code{ppv}, \code{npv}
#' which are \code{npred} x \code{nthresh} matrices, and and element \code{thresh}
#' which is a vector of thresholds.
#'
#' @alias plot.bayesLog.class
#'

classify <- function(pred, obs, thresh)
{    
    #check inputs
    stopifnot(is.matrix(pred))
    stopifnot(is.vector(obs))
    stopifnot(is.vector(thresh))
    
    stopifnot(max(pred) <= 1.0 & min(pred) >= 0)
    stopifnot(length(obs[obs != 0 & obs != 1]) == 0)
    stopifnot(max(thresh) <= 1.0 & min(thresh) >= 0)
    thresh <- unique(sort(thresh))
    
    output <- classification(pred, obs, thresh)
    
    #add thresholds
    output$thresh <- thresh
    
	#set class
	class(output) <- "bayesLog.class"
	output
}

##' @describeIn classify
##' @export
#plot.bayesLog.class <- function(x, ...)
#{
#    stopifnot(class(x) == "bayesLog.class")
#    
#    x <- lapply(x, function(x) c(mean(x), quantile(x, probs = c(0.025, 0.975))))
#    x <- lapply(x, as.data.frame)
#    x <- lapply(1:length(x), function(i, x, name)
#    {
#        x <- x[[i]]
#        x$class <- rep(name[i], nrow(x))
#        x
#    }, x = x, name = names(x))
#}


