#' Subset method for \code{bayesLog} objects
#' 
#' \code{subset} method for class \code{bayesLog}
#' 
#' Extracts a subset from a \code{\link{bayesLog}} object.
#' Differs slightly from \code{\link[coda]{mcmc.subset}} in that it
#' always returns a \code{bayesLog} object, but with the iterations
#' reset (if \code{i} is not missing).
#' 
#' @param x a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param i the iterations to extract
#' @param j the variables to extract
#' @param \dots not used here
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}}
#'
#' @return A \code{bayesLog} object
#'
#' @export
#'

"[.bayesLog" <- function(x, i, j)
{
    stopifnot(class(x) == "bayesLog")
    
    #extract 'mcmc' object
    y <- as.matrix(x$post)
    
    #extract subset
    y <- y[i, j, drop = F]
    
    #return bayesLog object
    x$post <- as.mcmc(y)
    x
}

