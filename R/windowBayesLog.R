#' Time windows for \code{bayesLog} xs
#' 
#' \code{summary} method for class \code{bayesLog}
#' 
#' Acts as a wrapper function for \code{\link[coda]{window.mcmc}} from the \code{coda}
#' package
#' 
#' @param x a \code{bayesLog} x, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param vars a character specifying whether \code{"all"} parameters are to be plotted,
#' or just the regression terms (\code{"reg"}), or one of the random effect terms (\code{"rand"}).
#' If \code{vars = "rand"}, then the user must also specify a \code{rand} argument denoting
#' which term to extract.
#' @param rand an integer denoting which random effect to extract.
#' @param \dots arguments to pass to \code{\link{window.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{window.mcmc}}
#'
#' @export

window.bayesLog <- function(x, vars = c("all", "reg", "rand"), rand, ...)
{
    stopifnot(class(x) == "bayesLog")
    stopifnot(is.character(vars))
    stopifnot(!is.na(match(vars[1], c("reg", "rand", "all"))))
    
    if(vars[1] == "rand" & missing(rand)) stop("Need to specify which random effects to extract")
    
    #extract 'mcmc' x
    y <- x$post
    
    if(vars[1] == "reg") y <- y[, 1:(x$nregpars + x$nrand)]
    else
    {
        if(vars[1] == "rand")
        {
            stopifnot(x$nrand > 0 & !missing(rand))
            stopifnot(is.numeric(rand) & length(rand) == 1 & abs(floor(rand) - rand) < .Machine$double.eps ^ 0.5)
            stopifnot(!is.na(match(rand, 1:x$nrand)))
            
            z <- rep(1:x$nrand, x$nrandlevels)
            inds <- which(!is.na(match(z, rand)))
            inds <- inds + (x$nregpars + x$nrand)
            y <- y[, inds]
        }
    }
    
    #extract subset
    y <- window(y, ...)
    
    #generate new 'bayesLog' x
    x$post <- y
    x
}
