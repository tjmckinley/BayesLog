#' @title Generates AUC for \code{bayesLog.pred} objects
#' 
#' @description Produce posterior for AUC from \code{bayesLog.pred} objects.
#' 
#' @param object    An object of class \code{bayesLog.pred}.
#' @param obs       A binary vector of observations.
#' @param \dots not used here.
#' @author TJ McKinley
#'
#' @return A \code{bayesLog.auc} object, essentially a vector of posterior 
#' samples.
#'
#' @export
AUC <- function(object, obs, ...) UseMethod("AUC")

#' @describeIn AUC AUC method for \code{bayesLog.pred} objects
#' @export
AUC.bayesLog.pred <- function(object, obs, ...)
{
    if (!requireNamespace("caTools", quietly = TRUE))
    {
        stop("'caTools' package needed for function to work. Please install it.", 
             call. = FALSE)
    }
    require(caTools)
    
    stopifnot(class(object) == "bayesLog.pred")
    stopifnot(is.vector(obs))
    stopifnot(length(obs[obs != 0 & obs != 1]) == 0)
    
    class(object) <- "matrix"
    
    stopifnot(ncol(object) == length(obs))
    
    auc <- colAUC(t(object), obs)
    auc <- as.numeric(auc)
    class(auc) <- "bayesLog.auc"
    auc
}

#' @describeIn AUC Plot method for \code{bayesLog.auc} objects
#' @param x A \code{bayesLog.auc} object
#' @export
plot.bayesLog.auc <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog.auc")
    class(x) <- "numeric"
    ggplot(data.frame(AUC = x), aes(x = AUC)) + 
        geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
        ggtitle("Posterior predictive distribution") + ylab("Density")
}

#' @describeIn AUC Print method for \code{bayesLog.auc} objects
#' @export
print.bayesLog.auc <- function(x, ...)
{
    stopifnot(class(x) == "bayesLog.auc")
    class(x) <- "numeric"
    print(summary(x))
}
