#function to summarise MCMC output
summary.varselect <- function(output, topmodels = 5, ...)
{
    #extract relevant quantities from object
    varselect <- output$varselect
    output <- output$model.sim
    
    #thin chains
    output <- window(output, ...)
    
    #now extract indicators
    if(varselect == 0) print(summary(output))
    else
    {
        nchains <- length(output)
        output <- lapply(output, as.matrix)
        indicators <- lapply(output, function(x) x[, grep(glob2rx("I_*"), colnames(x)), drop = F])
        output <- lapply(output, function(x) x[, -grep(glob2rx("I_*"), colnames(x)), drop = F])
        output <- lapply(output, as.mcmc)
        output <- as.mcmc.list(output)
        #print posterior summaries
        cat("###########    POSTERIOR SUMMARIES    ###########\n")
        print(summary(output))
        #print PPAs for variables
        cat("\n###########    PPAs for VARIABLES    ###########\n")
        PPAvar <- sapply(indicators, function(x) apply(x, 2, mean))
        if(!is.matrix(PPAvar)) PPAvar <- matrix(PPAvar, nrow = 1)
        rownames(PPAvar) <- sapply(strsplit(colnames(indicators[[1]]), "_"), function(x) x[2])
        PPAvar <- t(apply(PPAvar, 1, function(x)
        {
            y <- x
            x <- x[!is.na(x)]
            x <- c(y, mean(x), sd(x))
            x <- signif(x, digits = 2)
            x <- c(as.character(x[-length(x)]), paste0("(", x[length(x)], ")"))
        }))
        colnames(PPAvar) <- c(paste("Chain", 1:(ncol(PPAvar) - 2)), "Mean", "SD")
        print(t(PPAvar), quote = F)
        #print PPAs for models
        cat("\n###########    PPAs for MODELS    ###########\n")
        indicators1 <- lapply(indicators, function(x) apply(x, 1, function(y) paste(y, collapse = ".")))
        indicators2 <- lapply(indicators1, table)
        indicators2 <- lapply(indicators2, function(x) x / sum(x))
        #match across chains
        ind <- unique(do.call("c", indicators1))
        indmatch2 <- lapply(indicators2, function(x, ind) match(ind, names(x)), ind = ind)
        indicators2 <- lapply(1:length(indicators2), function(i, x, indmatch) x[[i]][indmatch[[i]]], indmatch = indmatch2, x = indicators2)
        indicators2 <- do.call("cbind", indicators2)
        rownames(indicators2) <- ind
        indicators2 <- cbind(indicators2, t(apply(indicators2, 1, function(x)
        {
            x <- x[!is.na(x)]
            c(mean(x), sd(x))
        })))
        #sort by average
        topmodels <- ifelse(topmodels > nrow(indicators2), nrow(indicators2), topmodels)
        indicators2 <- indicators2[sort.list(indicators2[, nchains + 1], decreasing = T), ][1:topmodels, ]
        indicators1 <- rownames(indicators2)
        indicators1 <- sapply(indicators1, function(y) strsplit(y, "\\.")[[1]])
        indicators2 <- apply(indicators2, 1, function(x)
        {
            x <- signif(x, digits = 2)
            c(as.character(x[-length(x)]), paste0("(", x[length(x)], ")"))
        })
        indicators1 <- rbind(indicators1, indicators2)
        colnames(indicators1) <- paste0("M", 1:ncol(indicators1))
        var <- colnames(indicators[[1]])
        var <- sapply(strsplit(var, "_"), function(x) x[2])
        rownames(indicators1) <- c(var, paste("Chain", 1:nchains), "Mean", "SD")
        indicators1[indicators1 == "0"] <- ""
        print(indicators1, quote = F)
    }
}
