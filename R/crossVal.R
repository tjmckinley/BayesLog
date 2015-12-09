#' @title Function to perform cross-validation
#'
#' @description Fits logistic regression model using K-fold cross-validation
#' in order to explore robustness of predictions
#'
#' @export
#'
#' @param data      A data matrix.
#' @param model     A \code{glm} model object fitted to \code{data}.
#' @param K         Defines how many subsets to sue for cross-validation
#' @param plot      Logical defining whether to produce plot of ROC curves.
#'
#' @return An vector containing the AUC values for the full model and the
#' cross-validation fits.
#'

crossVal <- function(data, model, K = 5, plot = TRUE)
{
    #check inputs
    stopifnot(!missing(data) & !missing(model))
    stopifnot(is.data.frame(data))
#    stopifnot(class(model) == "glm" | class(model) == "lm")
    stopifnot(is.numeric(K) & length(K) == 1)
    stopifnot(K > 0 & abs(floor(K) - K) < .Machine$double.eps^0.5)
    stopifnot(is.logical(plot) & length(plot) == 1)
    
    #calculate numbers of observations required in each group
    nsubset <- floor(nrow(data) / K)
    subsets <- rep(nsubset, K)
    subsets[K] <- nrow(data) - (K - 1) * nsubset
    stopifnot(sum(subsets) == nrow(data))
    
    #now randomly sort data set
    data <- data[sample(1:nrow(data), nrow(data), replace = F), ]
    
    #now split data into K subsets
    subsets <- split(data, rep(1:K, times = subsets))
    
    #now fit model using cross-validation and plot AUCs
    model <- update(model, data = data)
    temp.orig <- roc(data$binres, fitted(model))
    if(plot) plot(temp.orig, legacy.axes = T)
    AUC <- numeric(K + 1)
    AUC[1] <- auc(temp.orig)
    
    #now perform cross-validation
    for(i in 1:K)
    {
        temp.data <- do.call("rbind", subsets[-i])
        temp.model <- update(model, data = temp.data)
        temp <- predict(temp.model, subsets[[i]], type = "response")
        temp <- roc(subsets[[i]]$binres, temp)
        if(plot) plot(temp, add = T, col = "red")
        AUC[i + 1] <- auc(temp)
    }
    #just replot original line for clarity
    if(plot) plot(temp.orig, add = T)
    names(AUC) <- c("full", 1:K)
    AUC
}

