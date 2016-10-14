#' @title Generates posterior distributions for sensitivity
#' specificity, PPV and NPV
#'
#' @description Generates posterior distributions for sensitivity
#' specificity, PPV and NPV from a set of posterior predictions.
#'
#' @export
#'
#' @param pred       A \code{bayesLog.pred} object.
#' @param obs 		 Vector of binary observations (of length \code{nobs})
#' @param thresh     Vector of thresholds between 0 and 1 at which to generate
#'                   classification estimates.
#' @param \dots        Not currently used.
#'
#' @return An object of class \code{bayesLog.class}, which is a list
#' including a subset of elements: \code{sens}, \code{spec}, \code{ppv}, \code{npv}
#' which are \code{npred} x \code{nthresh} matrices, and and element \code{thresh}
#' which is a vector of thresholds.
#'

classify <- function(pred, obs, thresh, ...) UseMethod("classify")

#' @describeIn classify Classify method for \code{bayesLog.pred} objects
#' @export
classify.bayesLog.pred <- function(pred, obs, thresh, ...)
{    
    #check inputs
    stopifnot(class(pred) == "bayesLog.pred")
    stopifnot(is.vector(obs))
    stopifnot(is.vector(thresh))
    
    stopifnot(max(pred) <= 1.0 & min(pred) >= 0)
    stopifnot(length(obs[obs != 0 & obs != 1]) == 0)
    stopifnot(max(thresh) <= 1.0 & min(thresh) >= 0)
    thresh <- unique(sort(thresh))
    
    output <- classification(pred, obs, thresh)
    
    #add thresholds
    output$thresh <- thresh
    
	#set class
	class(output) <- "bayesLog.class"
	output
}

#' @describeIn classify Summary method for \code{bayesLog.class} objects
#' @param object    An object of class \code{bayesLog.class}
#' @export
summary.bayesLog.class <- function(object, ...)
{
    stopifnot(class(object) == "bayesLog.class")
    
    thresh <- object$thresh
    object <- object[-which(names(object) == "thresh")]
    
    mn <- t(sapply(object, function(x) apply(x, 2, mean)))
    colnames(mn) <- thresh
    
    quantile <- lapply(object, function(x) apply(x, 2, function(x)
    {
        x <- x[is.finite(x)]
        quantile(x, probs = c(0.025, 0.975))
    }))
    quantile <- lapply(quantile, function(x, thresh){ colnames(x) <- thresh; x}, thresh = thresh)
    
    x <- list(mean = mn, quantiles = quantile)
    class(x) <- "summary.bayesLog.class"
    x
}

#' @describeIn classify Print method for \code{bayesLog.class} objects
#' @param x    An object of class \code{bayesLog.class}
#' @export
print.bayesLog.class <- function(x, ...)
{
    stopifnot(class(x) == "summary.bayesLog.class")
    
    x <- summary(x)
    print(x)
}

#' @describeIn classify Plot method for \code{bayesLog.class} objects
#' @export
#' @param type      Character vector defining whether to plot
#'                  individual or comparative plots.
#' @param comp_type Character containing subplot to extract when doing comparative plotting.
#'                  If \code{"all"} then plots both ROC and Predictive Values plots, else one
#'                  can choose to produce either.
plot.bayesLog.class <- function(x, type = c("ind", "comp"), comp_type = c("all", "ROC", "PV"), ...)
{
   stopifnot(class(x) == "bayesLog.class")
   stopifnot(length(which(is.na(match(type[1], c("ind", "comp"))))) == 0)
   stopifnot(length(which(is.na(match(comp_type[1], c("all", "ROC", "PV"))))) == 0)
   
   #extract thresholds
   thresh <- x$thresh
   x <- x[-which(names(x) == "thresh")]
   
   x <- lapply(x, function(x) apply(x, 2, function(x)
   {
        x <- x[is.finite(x)]
        c(mean(x), quantile(x, probs = c(0.025, 0.975)))
   }))
   x <- lapply(x, function(x)
   {
        x <- t(x)
        x <- as.data.frame(x)
        colnames(x) <- c("mean", "LCI", "UCI")
        x
   })
   x <- lapply(1:length(x), function(i, x, name, thresh)
   {
       x <- x[[i]]
       x$class <- rep(name[i], nrow(x))
       x$thresh <- thresh
       x
   }, x = x, name = names(x), thresh = thresh)
   
   x <- na.omit(x)
   
   if(type[1] == "ind")
   {
       x <- do.call("rbind", x)
       x$class <- factor(x$class)
       
       p <- ggplot(x, aes(x = thresh, y = mean)) + geom_point() +
           geom_errorbar(aes(x = thresh, ymax = UCI, ymin = LCI)) +
           facet_wrap(~class)
       return(p)
    }
    else
    {
        #extract sens and spec
        x1 <- data.frame(sens = x[[1]]$mean, spec = x[[2]]$mean, thresh = x[[1]]$thresh)
        x2 <- data.frame(ppv = x[[3]]$mean, npv = x[[4]]$mean, thresh = x[[1]]$thresh)
        
        #setup plots
        p1 <- ggplot(x1, aes(x = 1 - spec, y = sens, label = thresh)) + geom_point() + geom_line() +
            xlab("1 - Specificity") + ylab("Sensitivity") + coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
            geom_text(hjust = 0, nudge_x = 0.05, check_overlap = T)
        p2 <- ggplot(x2, aes(x = 1 - npv, y = ppv, label = thresh), ylim = c(0, 1), xlim = c(0, 1)) + geom_point() + geom_line() +
            xlab("1 - NPV") + ylab("PPV") + coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
            geom_text(hjust = 0, nudge_x = 0.05, check_overlap = T)
        
        #in order to plot the object on output, we add a
        #print function below as a small hack
        if(comp_type[1] == "all")
        {
            p <- arrangeGrob(grobs = list(p1, p2), layout_matrix = matrix(c(1, 2), nrow = 1))
            class(p) <- c("grArrange", class(p))
            return(p)
        }
        else
        {
            if(comp_type[1] == "ROC") return(p1)
            else return(p2)
        }
    }
}

#' @export
print.grArrange <- function(x, ...)
{
    stopifnot(class(x)[1] == "grArrange")
    class(x) <- class(x)[-1]
    
    #set width and height of plot
    if(dev.cur() == 1) dev.new(width = 7, height = 3.5)
    
    grid.draw(x)
}

#' @export
print.summary.bayesLog.class <- function(x, ...)
{
    stopifnot(class(x) == "summary.bayesLog.class")
    
    cat("Classification means:\n\n")
    print(x$mean)
    cat("\nQuantiles:\n\n")
    print(x$quantiles)
}


