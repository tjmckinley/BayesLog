#' Conversion methods for \code{bayesLog} objects
#' 
#' Methods for generic functions \code{as.matrix}, \code{as.array} and
#' \code{as.mcmc} for class \code{bayesLog}
#' 
#' Converts \code{\link{bayesLog}} objects to different forms.
#' Essentially a wrapper for \code{\link[coda]{mcmc.convert}} functions 
#' from the \code{coda}
#' package.
#' 
#' @param x a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param \dots additional arguments to be passed to \code{\link[coda]{mcmc.convert}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{mcmc.convert}}
#'
#' @return A \code{matrix}, \code{array} or \code{\link[coda]{mcmc}} object
#'
#' @aliases as.matrix.bayesLog as.array.bayesLog as.mcmc.bayesLog
#'
#' @export
#'

as.matrix.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    #extract 'mcmc' object
    y <- x$post
    
    #convert
    y <- as.matrix(y, ...)
    
    #return matrix object
    y
}

as.array.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    #extract 'mcmc' object
    y <- x$post
    
    #convert
    y <- as.array(y, ...)
    
    #return array object
    y
}

as.mcmc.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    #extract 'mcmc' object
    y <- x$post
    
    #return mcmc object
    y
}

