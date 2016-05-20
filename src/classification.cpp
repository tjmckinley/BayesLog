#include "functions.hpp"

//function for producing sens, spec, ppv and npv from posterior predictions
List classification (NumericMatrix pred, IntegerVector obs, NumericVector thresh)
{
    // 'pred' is npred x nobs matrix of predicted probabilities
    // 'obs' is numeric vector of binary observations of length nobs
    // 'thresh' is vector of thresholds between 0 and 1 for classification
    
    //initialise indexing variables
    int i, j, k;
    
    // calculate number of parameters
    int nobs = obs.size();
    int npred = pred.nrow();
    int nthresh = thresh.size();
    
    if(pred.ncol() != nobs) stop("'pred' and 'obs' don't match.");
    
    //set up intermediate objects
    IntegerVector classifier(nobs);
    IntegerMatrix tab(2, 2);
    
    //set up matrices for outputs
    NumericMatrix sens(nobs, nthresh);
    NumericMatrix spec(nobs, nthresh);
    NumericMatrix ppv(nobs, nthresh);
    NumericMatrix npv(nobs, nthresh);
    
    //loop over predictions
    for(i = 0; i < npred; i++)
    {
        //loop over thresholds
        for(j = 0; j < nthresh; j++)
        {
            //set classification based on current threshold
            for(k = 0; k < nobs; k++) classifier(k) = (pred(i, k) > thresh(j) ? 1:0);
            
            //calculate sens etc.
            tab(0, 0) = 0; tab(0, 1) = 0; tab(1, 0) = 0; tab(1, 1) = 0;
            
            for(k = 0; k < nobs; k++)
            {
                tab(0, 0) += (obs(k) == 0 && classifier(k) == 0 ? 1:0);
                tab(0, 1) += (obs(k) == 1 && classifier(k) == 0 ? 1:0);
                tab(1, 0) += (obs(k) == 0 && classifier(k) == 1 ? 1:0);
                tab(1, 1) += (obs(k) == 1 && classifier(k) == 1 ? 1:0);
            }
            
            sens(k, j) = ((double) tab(1, 1)) / ((double) (tab(0, 1) + tab(1, 1)));
            spec(k, j) = ((double) tab(0, 0)) / ((double) (tab(0, 0) + tab(1, 0)));
            ppv(k, j) = ((double) tab(1, 1)) / ((double) (tab(1, 0) + tab(1, 1)));
            npv(k, j) = ((double) tab(0, 0)) / ((double) (tab(0, 0) + tab(0, 1)));
        }
    }
    
    //return output
    List output;
    output["sens"] = sens;
    output["spec"] = spec;
    output["ppv"] = ppv;
    output["npv"] = npv;
    
    return(output);
}

