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
    
    #extract 'mcmc' object
    y <- object$post
    
    if(vars[1] == "reg") y <- y[, 1:(object$nregpars + object$nrand)]
    else
    {
        if(vars[1] == "rand")
        {
            stopifnot(object$nrand > 0 & !missing(rand))
            stopifnot(is.numeric(rand) & length(rand) == 1 & abs(floor(rand) - rand) < .Machine$double.eps ^ 0.5)
            stopifnot(!is.na(match(rand, 1:object$nrand)))
            
            z <- rep(1:object$nrand, object$nrandlevels)
            randnames <- object$randnames[rand]
            inds <- which(!is.na(match(z, rand)))
            inds <- inds + (object$nregpars + object$nrand)
            y <- y[, inds]
        }
    }
    
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
