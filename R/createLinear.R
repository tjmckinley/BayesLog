#function to output correct linear term for use in MCMC chain
createLinear <- function(data, response)
{
    y <- names(data)
    data1 <- data.frame(response = data[, y == response])
    colnames(data1) <- response
    if(is.factor(data1[, 1])) data1[, 1] <- as.numeric(data1[, 1]) - 1
    
    #remove response variable
    data <- data[, y != response, drop = F]
    y <- y[y != response]
    
    #record factor index
    factindex <- NULL
    vars <- NULL
    for(i in 1:length(y))
    {
        if(is.numeric(data[, i]))
        {
            data1 <- cbind(data1, data[, i])
            colnames(data1)[ncol(data1)] <- y[i]
            if(is.null(factindex))
            {
                factindex <- i
                vars <- y[i]
            }
            else
            {
                factindex <- c(factindex, i)
                vars <- c(vars, y[i])
            }
        }
        if(is.factor(data[, i]))
        {
            #expand categorical variables
            for(j in 2:length(levels(data[, i])))
            {
                data1 <- cbind(data1, ifelse(data[, i] == levels(data[, i])[j], 1, 0))
                colnames(data1)[ncol(data1)] <- paste0(y[i], j)
            }
            if(is.null(factindex))
            {
                factindex <- rep(i, length(levels(data[, i])) - 1)
                for(j in 2:length(levels(data[, i]))) vars <- paste0(y[i], "_", levels(data[, i])[j])
            }
            else
            {
                factindex <- c(factindex, rep(i, length(levels(data[, i])) - 1))
                for(j in 2:length(levels(data[, i]))) vars <- c(vars, paste0(y[i], "_", levels(data[, i])[j]))
            }
        }
    }
    return(list(data = data1, factindex = factindex, vars = vars))
}
