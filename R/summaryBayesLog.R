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
#' @param vars a character specifying whether just the regression terms are to be
#' plotted (\code{"reg"}), or all of the terms (\code{"all"}).
#' @param \dots additional arguments to be passed to \code{\link[coda]{summary.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{summary.mcmc}}
#'
#' @return A \code{summary.mcmc} object
#'
#' @export
#'
#' @aliases print.bayesLog

summary.bayesLog <- function(object, vars = c("reg", "all"), ...)
{
    stopifnot(class(object) == "bayesLog")
    
    stopifnot(is.character(vars))
    stopifnot(!is.na(match(vars[1], c("reg", "all"))))
    
    #extract 'mcmc' x
    y <- object$post
    
    if(vars[1] == "reg") y <- y[, 1:(object$nregpars + object$nrand)]
    
    #summarise
    summary(y, ...)
}

#generic print function for \code{bayesLog} objects
#' @export
print.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    print(x, ...)
}
