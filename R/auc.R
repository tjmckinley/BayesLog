#' @title Generates AUC for \code{bayesLog.pred} objects
#' 
#' @description Produce posterior for AUC from \code{\link{bayesLog.pred}} objects.
#' 
#' @export
#' 
#' @param x    An object of class \code{bayesLog.pred}.
#' @param obs  A vector of binary observations.
#' @param \dots not used here.
#' @author TJ McKinley
#'
#' @return A \code{bayesLog.auc} object, essentially a vector of posterior 
#' samples.
#'

AUC <- function(x, ...) UseMethod("AUC")

AUC.bayesLog.pred <- function(x, obs, ...)
{
    if (!requireNamespace("caTools", quietly = TRUE))
    {
        stop("'caTools' package needed for function to work. Please install it.", 
             call. = FALSE)
    }
    require(caTools)
    
    stopifnot(class(x) == "bayesLog.pred")
    stopifnot(is.vector(obs))
    stopifnot(length(obs[obs != 0 & obs != 1]) == 0)
    
    class(x) <- "matrix"
    
    stopifnot(ncol(x) == length(obs))
    
    auc <- colAUC(t(x), obs)
    auc <- as.numeric(auc)
    class(auc) <- "bayesLog.auc"
    auc
}

#' @describeIn AUC Plot method for \code{bayesLog.auc} objects.
#' @param x    An object of class \code{bayesLog.auc}
#' @export

plot.bayesLog.auc <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog.auc")
    class(x) <- "numeric"
    ggplot(data.frame(AUC = x), aes(x = AUC)) + 
        geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
        ggtitle("Posterior predictive distribution") + ylab("Density")
}

#' @describeIn AUC Print method for \code{bayesLog.auc} objects.
#' @param x    An object of class \code{bayesLog.auc}
#' @export

print.bayesLog.auc <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog.auc")
    class(x) <- "numeric"
    print(summary(x))
}