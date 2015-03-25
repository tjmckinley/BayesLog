#function to plot MCMC output
plot.varselect <- function(output, cond = F, trace = T, density = T, ask = T, ...)
{
    if(cond == T & output$varselect == F) cond <- F
    
    if(cond)
    {
        #plot conditional traces
        output <- extractCond(output)
        # set up plotting parameters
        if(trace & density)
        {
            maxp <- length(output) * 2
            if (maxp >= 4) mfrow1 <- c(4, 2)
            else mfrow1 <- c(maxp, 2)
        }
        else
        {
            maxp <- length(output)
            if(maxp >= 4) mfrow1 <- c(4, 1)
            else mfrow1 <- c(maxp, 1)
        }
        cols1 <- c("black", "red", "blue", "green", "yellow", "purple")
        par(mfrow = mfrow1)
        # produce plots
        l <- 1
        for (i in 1:length(output))
        {
            # plot trace
            temp <- output[[i]]
            if(nrow(temp) > 1)
            {
                if(trace)
                {
                  plot(1:nrow(temp), temp[, 1], type = "n", main = names(output)[i], xlab = "Index", ylab = "Value")
                  segments(x0 = 1:(nrow(temp) - 1), y0 = temp[-nrow(temp), 1], x1 = 2:nrow(temp), y1 = temp[-1, 1], col = cols1[temp[, 2]])
                }
                if(density)
                {
                  # plot density
                  plot(density(temp[, 1]), main = names(output)[i], xlab = "Value", ylab = "Density")
                }
            }
            else
            {
                if(trace)
                {
                  plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                    main = names(output)[i])
                  text(0, 0, "<= 1 Sample")
                }
                if(density)
                {
                  plot(0, 0, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", 
                    main = names(output)[i])
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
    else plot(output[[1]], trace = trace, density = density, ask = ask, ...)
}
