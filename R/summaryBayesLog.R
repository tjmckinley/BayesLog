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
#' @param vars a character specifying whether \code{"all"} parameters are to be plotted,
#' or just the regression terms (\code{"reg"}), or one of the random effect terms (\code{"rand"}).
#' If \code{vars = "rand"}, then the user must also specify a \code{rand} argument denoting
#' which term to extract.
#' @param rand an integer denoting which random effect to extract.
#' @param \dots additional arguments to be passed to \code{\link[coda]{summary.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{summary.mcmc}}
#'
#' @return A \code{summary.mcmc} object
#'
#' @export
#'
#' @aliases print.bayesLog

summary.bayesLog <- function(object, vars = c("all", "reg", "rand"), rand, ...)
{
    stopifnot(class(object) == "bayesLog")
    
    stopifnot(is.character(vars))
    stopifnot(!is.na(match(vars[1], c("reg", "rand", "all"))))
    
    if(vars[1] != "all") object <- window(object, vars = vars[1], rand = rand)
    
    #extract 'mcmc' object
    y <- object$post
    
    
    #summarise
    summary(y, ...)
}

#generic print function for \code{bayesLog} objects
#' @export
print.bayesLog <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog")
    
    print(summary(x, ...))
}
