# window function for MCMC output
window.varselect <- function(x, start = NA, end = NA, thin = NA, chains = NA, ...) {
    y <- x
    if (!is.na(chains)) {
        x <- x[[1]][chains]
    }
    else x <- x[[1]]
    
    if(is.na(start)) {
        if(is.na(end)) {
            if(!is.na(thin)) x <- window(x, thin = thin)
        }
        else {
            if(!is.na(thin)) x <- window(x, end = end, thin = thin)
            else x <- window(x, end = end)
        }
    }
    else {
        if(is.na(end)) {
            if(!is.na(thin)) x <- window(x, start = start, thin = thin)
            else x <- window(x, start = start)
        }
        else {
            if(!is.na(thin)) x <- window(x, start = start, end = end, thin = thin)
            else x <- window(x, start = start, end = end)
        }
    }
    y[[1]] <- x
    y
} 
