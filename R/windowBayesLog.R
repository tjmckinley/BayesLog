#' Time windows for \code{bayesLog} objects
#' 
#' \code{summary} method for class \code{bayesLog}
#' 
#' Acts as a wrapper function for \code{\link[coda]{window.mcmc}} from the \code{coda}
#' package
#' 
#' @param x a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param start the first iteration of interest.
#' @param end the last iteration of interest.
#' @param thin the required interval between successive samples.
#' @param chains the chains to extract (if \code{NA} then extracts all chains).
#' @param \dots futher arguments for future methods.
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{window.mcmc}}
#'
#' @import coda
#' @export window.bayesLog

window.bayesLog <- function(x, start = NA, end = NA, thin = NA, chains = NA, ...)
{
    stopifnot(class(x) == "bayesLog")
    y <- x
    if (!is.na(chains)) {
        x <- x[[1]][chains]
    }
    else x <- x[[1]]
    
    if(is.na(start)) {
        if(is.na(end)) {
            if(!is.na(thin)) x <- window(x, thin = thin)
        }
        else {
            if(!is.na(thin)) x <- window(x, end = end, thin = thin)
            else x <- window(x, end = end)
        }
    }
    else {
        if(is.na(end)) {
            if(!is.na(thin)) x <- window(x, start = start, thin = thin)
            else x <- window(x, start = start)
        }
        else {
            if(!is.na(thin)) x <- window(x, start = start, end = end, thin = thin)
            else x <- window(x, start = start, end = end)
        }
    }
    y[[1]] <- x
    class(y) <- "bayesLog"
    y
} 
