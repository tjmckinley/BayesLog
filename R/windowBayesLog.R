#' Time windows for \code{bayesLog} xs
#' 
#' \code{summary} method for class \code{bayesLog}
#' 
#' Acts as a wrapper function for \code{\link[coda]{window.mcmc}} from the \code{coda}
#' package
#' 
#' @param x a \code{bayesLog} x, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param \dots arguments to pass to \code{\link{window.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{window.mcmc}}
#'
#' @export

window.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    #extract 'mcmc' x
    y <- x$post
    
    #extract subset
    y <- window(y, ...)
    
    #generate new 'bayesLog' x
    x$post <- y
    x
}
