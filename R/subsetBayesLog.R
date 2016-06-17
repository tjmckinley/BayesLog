#' Subset method for \code{bayesLog} objects
#' 
#' \code{subset} method for class \code{bayesLog}
#' 
#' Extracts a subset from a \code{\link{bayesLog}} object.
#' If extracting a single iteration, then function returns either a 
#' matrix or a list of matrices. Otherwise function returns a \code{bayesLog} 
#' object (is not subsetting by column), or an \code{mcmc} or \code{mcmc.list} 
#' object otherwise. If \code{i} is not missing then the iteration counts are reset.
#' 
#' @param x a \code{bayesLog} object, usually as a result of a call to
#' \code{\link{bayesLog}}.
#' @param i the iterations to extract
#' @param j the variables to extract (can be a numerical vector, or a 
#' character vector of variables to extract)
#' @param dropcols  logical denoting whether to drop selected columns or keep them
#' @param \dots not used here
#' @author TJ McKinley
#' @seealso \code{\link{bayesLog}} \code{\link{plot.bayesLog}} \code{\link{window.bayesLog}}
#'
#' @return A \code{bayesLog} object
#'
#' @export
#'

"[.bayesLog" <- function(x, i, j, dropcols = FALSE)
{
    stopifnot(class(x) == "bayesLog")
  
    # check if the "j" subset is a character vector of names
    # and if so then check they're viable and subset accordingly
    if(!missing(j))
    {
      if(is.character(j))
      {
        #extract all variables in model
        predictors <- all.vars(x$formula)
        
        #remove response 
        predictors <- predictors[-match(as.character(x$formula)[2], predictors)]
        
        #extract intercept from names if entered
        int <- na.omit(match(c("Intercept", "(Intercept)"), j))
        if(length(int) > 0)
        {
          stopifnot(length(int) == 1)
          j <- j[-int]
        }
        else int <- NA
        
        #cross reference these against required variable names
        j <- unique(j)
        temp <- match(j, predictors)
        if(!all(!is.na(temp))) stop("Variable names in required subset do not match variables in model")
        
        #reset to get indexing correct
        predictors <- all.vars(x$formula)
        
        # split j into fixed and random effects
        randj <- j[!is.na(match(j, x$randnames))]
        fixj <- j[is.na(match(j, x$randnames))]
        
        # if valid, then extract correct columns corresponding to columns
        if(length(fixj) > 0)
        {
          fixlistj <- lapply(fixj, function(x, names, assign)
          {
            x <- match(x, names)
            which(assign == x)
          }, names = predictors, 
          assign = attr(model.matrix(lme4::nobars(x$formula), x$data), "assign") + 1)
          listj <- fixlistj
        }
        if(length(randj) > 0)
        {
          randlistj <- lapply(randj, function(rand, x)
          {
            z <- rep(1:x$nrand, x$nrandlevels)
            inds <- which(!is.na(match(z, which(rand == x$randnames))))
            inds <- inds + (x$nregpars + x$nrand)
            inds
          }, x = x)
          if(length(fixj) > 0) listj <- c(listj, randlistj)
          else listj <- randlistj
        }
        
        # reorder list of indices
        listj <- listj[match(c(fixj, randj), j)]
        
        # add intercept back in if required
        if(!is.na(int))
        {
          if(int == 1) listj <- c(1, listj)
          else listj <- c(listj[1:(int - 1)], 1, listj[int:length(listj)])
        }
        
        # extract indices as vector
        j <- do.call("c", listj)
      }
    }
    
    if(dropcols) j <- -j
    
    #extract 'mcmc' object
    y <- x$post
    
    #extract subset
    if(is.mcmc.list(y))
    {
        y <- lapply(y, as.matrix)
        
        if(missing(i)) i <- 1:nrow(y[[1]])
        
        for(k in 1:length(y))
        {
            y[[k]] <- y[[k]][i, j, drop = F]
        }
        if(length(i) > 1)
        {
          y <- lapply(y, as.mcmc)
          y <- as.mcmc.list(y)
        }
    }
    else
    {
      y <- as.matrix(y)
      if(missing(i)) i <- 1:nrow(y)
      y <- y[i, j, drop = F]
      if(length(i) > 1) y <- as.mcmc(y)
    }
    
    if(length(i) > 1)
    {
      #update bayesLog object
      x$post <- y
      
      #if subsetting by column then return an mcmc object
      if(!missing(j)) x <- y
    }
    else x <- y
    x
}

