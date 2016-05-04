#' Plots \code{bayesLog} objects
#' 
#' \code{plot} method for class \code{bayesLog}
#' 
#' Produces trace and/or density plots for \code{\link{bayesLog}} objects.
#' These are based on the \code{\link[coda]{plot.mcmc}} functions from the \code{coda}
#' package
#' 
#' @param x a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param \dots additional arguments to be passed to \code{\link[coda]{plot.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{window.bayesLog}} 
#' \code{\link[coda]{plot.mcmc}}
#'
#' @export

plot.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    #extract 'mcmc' object
    y <- x$post
    
    plot(y, ...)
}
