#' Summary statistics for \code{bayesLog} objects
#' 
#' \code{summary} method for class \code{bayesLog}
#' 
#' Produce summary statistics for \code{\link{bayesLog}} objects.
#' Essentially a wrapper for \code{\link[coda]{summary.mcmc}} functions from the \code{coda}
#' package.
#' 
#' @param object a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param \dots additional arguments to be passed to \code{\link[coda]{summary.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{summary.mcmc}}
#'
#' @return A \code{summary.mcmc} object
#'
#' @export
#'
#' @aliases print.bayesLog

summary.bayesLog <- function(object, ...)
{
    stopifnot(class(object) == "bayesLog")
    
    #extract 'mcmc' object
    object <- object$post
    
    #summarise
    summary(object, ...)
}

#generic print function for \code{bayesLog} objects
#' @export
print.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    print(summary(x, ...))
}
