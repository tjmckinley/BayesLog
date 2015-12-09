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
#' @param vars a character specifying whether \code{"all"} parameters are to be plotted,
#' or just the regression terms (\code{"reg"}), or one of the random effect terms (\code{"rand"}).
#' If \code{vars = "rand"}, then the user must also specify a \code{rand} argument denoting
#' which term to extract.
#' @param rand an integer denoting which random effect to extract.
#' @param type a character defining the type of plot to be drawn: \code{"mcmc"} denotes
#' trace/density plots, whereas \code{"cat"} denotes a caterpillar plot, useful for
#' comparing the random effects terms. (This latter plot can only be drawn for RE term currently.)
#' @param \dots additional arguments to be passed to \code{\link[coda]{plot.mcmc}}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{plot.mcmc}}
#'
#' @export

plot.bayesLog <- function(x, vars = c("all", "reg", "rand"), rand, type = c("mcmc", "cat"), ...)
{
    stopifnot(class(x) == "bayesLog")
    
    stopifnot(is.character(vars))
    stopifnot(!is.na(match(vars[1], c("reg", "rand", "all"))))
    
    stopifnot(is.character(type))
    stopifnot(!is.na(match(type[1], c("mcmc", "cat"))))
    
    #extract 'mcmc' object
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
            randnames <- x$randnames[rand]
            inds <- which(!is.na(match(z, rand)))
            inds <- inds + (x$nregpars + x$nrand)
            y <- y[, inds]
        }
    }
    
    #plot
    if(type[1] == "mcmc") plot(y, ...)
    else
    {
        stopifnot(vars[1] == "rand")
        
        #plot posteriors for random intercepts 
        y <- as.matrix(y)
        #turn into data frame to use in ggplot
        y <- as.data.frame(y)
        y <- melt(y)
        colnames(y) <- c("var", "nu")
        #reset names
        y$var <- as.character(y$var)
        y$var <- sapply(strsplit(y$var, randnames), function(x) x[[2]])
        y$var <- as.factor(y$var)
        #sort according to posterior mean
        y$var <- with(y, reorder(var, nu, mean))
        ggplot(y, aes(var, nu)) + geom_boxplot() + xlab(randnames)
    }
}
