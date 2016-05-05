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
#' @param vars a character specifying whether just the regression terms are to be
#' plotted (\code{"reg"}), or all of the terms (\code{"all"}).
#' @param \dots additional arguments to be passed to \code{\link[coda]{plot.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{window.bayesLog}} 
#' \code{\link[coda]{plot.mcmc}}
#'
#' @export

plot.bayesLog <- function(x, vars = c("reg", "all"), ...)
{
    stopifnot(class(x) == "bayesLog")
  
    stopifnot(is.character(vars))
    stopifnot(!is.na(match(vars[1], c("reg", "all"))))
    
    #extract 'mcmc' x
    y <- x$post
    
    if(vars[1] == "reg") y <- y[, 1:(x$nregpars + x$nrand)]
    
    plot(y, ...)
}
