#' Subset method for \code{bayesLog} objects
#' 
#' \code{subset} method for class \code{bayesLog}
#' 
#' Extracts a subset from a \code{\link{bayesLog}} object.
#' Essentially a wrapper for \code{\link[coda]{[.mcmc}} functions from the \code{coda}
#' package.
#' 
#' @param x a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param \dots additional arguments to be passed to \code{\link[coda]{[.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{[.mcmc}}
#'
#' @return A \code{bayesLog} object
#'
#' @export
#'

"[.bayesLog" <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    #extract 'mcmc' object
    y <- x$post
    
    #extract subset
    y <- y[...]
    
    #return bayesLog object
    x$post <- y
    x
}

