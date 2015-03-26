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
#' @param cond should the xs should be subsetted conditional on
#' inclusion in the model at each iteration.
#' @param trace plot trace of each variable.
#' @param density plot density of each variable.
#' @param ask Prompt user before each page of plots
#' @param \dots additional arguments to be passed to \code{\link[coda]{plot.mcmc}} if \code{cond = FALSE}
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{summary.bayesLog}} \code{\link{window.bayesLog}} \code{\link[coda]{plot.mcmc}}
#'
#' @import coda
#' @export

plot.bayesLog <- function(x, cond = F, trace = T, density = T, ask = T, ...)
{
    stopifnot(class(x) == "bayesLog")
    if(cond == T & length(grep(glob2rx("I_*"), colnames(x$model.sim[[1]]))) > 0) cond <- F
    
    if(cond)
    {
        #plot conditional traces
        x <- extractCond(x)
        # set up plotting parameters
        if(trace & density)
        {
            maxp <- length(x) * 2
            if (maxp >= 4) mfrow1 <- c(4, 2)
            else mfrow1 <- c(maxp, 2)
        }
        else
        {
            maxp <- length(x)
            if(maxp >= 4) mfrow1 <- c(4, 1)
            else mfrow1 <- c(maxp, 1)
        }
        cols1 <- c("black", "red", "blue", "green", "yellow", "purple")
        par(mfrow = mfrow1)
        # produce plots
        l <- 1
        for (i in 1:length(x))
        {
            # plot trace
            temp <- x[[i]]
            if(nrow(temp) > 1)
            {
                if(trace)
                {
                  plot(1:nrow(temp), temp[, 1], type = "n", main = names(x)[i], xlab = "Index", ylab = "Value")
                  segments(x0 = 1:(nrow(temp) - 1), y0 = temp[-nrow(temp), 1], x1 = 2:nrow(temp), y1 = temp[-1, 1], col = cols1[temp[, 2]])
                }
                if(density)
                {
                  # plot density
                  plot(density(temp[, 1]), main = names(x)[i], xlab = "Value", ylab = "Density")
                }
            }
            else
            {
                if(trace)
                {
                  plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                    main = names(x)[i])
                  text(0, 0, "<= 1 Sample")
                }
                if(density)
                {
                  plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                    main = names(x)[i])
                  text(0, 0, "<= 1 Sample")
                }
            }
            l <- l + 1
            if (ask == T)
            {
                if ((l - 1) %% 4 == 0 && (l - 1) != maxp) 
                readline("Press any key to continue:")
            }
        }
        par(mfrow = c(1, 1))
        return(cat(""))
    }
    else plot(x[[1]], trace = trace, density = density, ask = ask, ...)
}
