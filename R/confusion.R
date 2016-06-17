#' @title Generates confusion matrix from \code{bayesLog.pred} objects
#' 
#' @description Produces confusion matrix from \code{bayesLog.pred} objects.
#' 
#' @param x    An object of class \code{bayesLog.pred}.
#' @param obs  A vector of binary observations.
#' @param thresh A scalar between 0 and 1 defining threshold for classification.
#' @param plot   A logical denoting whether to plot confusion matrix or not.
#' @param \dots not used here.
#' @author TJ McKinley (based on code from Mario Recker)
#'
#' @return A confusion matrix, and a plot if required.
#'
#' @export

confusion <- function(x, obs, thresh, plot = T, ...) UseMethod("confusion")

#' @describeIn confusion Produces confusion matrix for \code{bayesLog.pred} objects
#' @export
confusion.bayesLog.pred <- function(x, obs, thresh, plot = T, ...)
{
    stopifnot(class(x) == "bayesLog.pred")
    stopifnot(is.vector(obs))
    stopifnot(length(obs[obs != 0 & obs != 1]) == 0)
    
    stopifnot(length(thresh) == 1)
    stopifnot(is.numeric(thresh) & thresh >= 0 & thresh <= 1)
    stopifnot(is.logical(plot))
    
    class(x) <- "matrix"
    
    stopifnot(ncol(x) == length(obs))
    
    #extract posterior mean
    x <- apply(x, 2, mean)
    Pred <- ifelse(x > thresh, 1, 0)
    Obs <- obs
    x <- t(table(Pred, Obs))
    if(plot) plot.confusion.matrix(x, c(0, 1), thresh)
    x
}

# @param x          A 2 x 2 contingency table
plot.confusion.matrix <- function(x, labs, thresh, ...)
{
    if (!requireNamespace("circlize", quietly = TRUE))
    {
        stop("'circlize' package needed for plot function to work. Please install it.", 
             call. = FALSE)
    }
    
    #load required libraries
    require(circlize)
    
    #create proportion
    r <- as.vector(prop.table(x, margin = 1))
    
    #set plot options
    par(mar = c(2, 5, 6, 1), pty = "s")
    
    #produce plot
    plot(NULL, type = "n", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', 
         xlim = c(0, 4), ylim = c(0, 4), xaxs = "i", yaxs = "i", asp = 1, cex.lab = 1.7, font.lab = 2)
    axis(3, at = c(1, 3), labels = labs, cex.axis = 1.5)
    axis(2, at = c(1, 3), labels = rev(labs), cex.axis = 1.5)
    mtext(paste0("Predicted (p = ", thresh, ")"), side = 3, line = 3, cex = 1.5)
    mtext("Observed", side = 2, line = 3, cex = 1.5)
    abline(h = 2, lwd = 0.8, col = 'blue')
    abline(v = 2, lwd = 0.8, col = 'blue')
    
    segments(x0 = seq(0.05, 1.95, by = 0.38), y0 = 1.98, y1 = 2.02, col = 'red', lwd = 0.8)
    segments(x0 = seq(2.05, 3.95, by = 0.38), y0 = 1.98, y1 = 2.02, col = 'red', lwd = 0.8)
    segments(y0 = seq(0.05, 1.95, by = 0.38), x0 = 1.98, x1 = 2.02, col = 'red', lwd = 0.8)
    segments(y0 = seq(2.05, 3.95, by = 0.38), x0 = 1.98, x1 = 2.02, col = 'red', lwd = 0.8)
    
    draw.sector(start.degree = 90, end.degree = 180, center = c(1.95, 2.05), rou1 = 1.9 * r[1], col = 'darkblue', clock.wise = F)
    draw.sector(start.degree = 180, end.degree = 270, center = c(1.95, 1.95), rou1 = 1.9 * r[2], col = 'lightblue', clock.wise = F)
    draw.sector(start.degree = 270, end.degree = 360, center = c(2.05, 1.95), rou1 = 1.9 * r[4], col = 'darkblue', clock.wise = F)
    draw.sector(start.degree = 0, end.degree = 90, center = c(2.05, 2.05), rou1 = 1.9 * r[3], col = 'lightblue', clock.wise = F)
    
    text(0.1, 3.9, paste(format(100 * r[1], digits = 1), "% (N = ", x[1, 1], ")", sep = ''), adj = c(0, 1), font = 2)
    text(0.1, 0.1, paste(format(100 * r[2], digits = 1), "% (N = ", x[2, 1], ")", sep = ''), adj = c(0, 0), font = 2)
    text(3.9, 3.9, paste(format(100 * r[3], digits = 1), "% (N = ", x[1, 2], ")", sep = ''), adj = c(1, 1), font = 2)
    text(3.9, 0.1, paste(format(100 * r[4], digits = 1), "% (N = ", x[2, 2], ")", sep = ''), adj = c(1, 0), font = 2)
    
    par(mar = c(5.1, 4.1, 4.1, 2.1), pty = "m")
}
