#include "functions.hpp"

//function for producing sens, spec, ppv and npv from posterior predictions
List classification (NumericMatrix pred, IntegerVector obs, IntegerVector nsamples, NumericVector thresh)
{
    // 'pred' is npred x nobs matrix of predicted probabilities
    // 'obs' is numeric vector of binary observations of length nobs
    // 'nsamples' contains number of samples for each obs entry
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
    NumericMatrix sens(npred, nthresh);
    NumericMatrix spec(npred, nthresh);
    NumericMatrix ppv(npred, nthresh);
    NumericMatrix npv(npred, nthresh);
    NumericMatrix miss(npred, nthresh);
    
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
                tab(0, 0) += (obs(k) == 0 && classifier(k) == 0 ? nsamples(k):0);
                tab(0, 1) += (obs(k) == 1 && classifier(k) == 0 ? nsamples(k):0);
                tab(1, 0) += (obs(k) == 0 && classifier(k) == 1 ? nsamples(k):0);
                tab(1, 1) += (obs(k) == 1 && classifier(k) == 1 ? nsamples(k):0);
            }
            
            sens(i, j) = ((double) tab(1, 1)) / ((double) (tab(0, 1) + tab(1, 1)));
            spec(i, j) = ((double) tab(0, 0)) / ((double) (tab(0, 0) + tab(1, 0)));
            ppv(i, j) = ((double) tab(1, 1)) / ((double) (tab(1, 0) + tab(1, 1)));
            npv(i, j) = ((double) tab(0, 0)) / ((double) (tab(0, 0) + tab(0, 1)));
            miss(i, j) = ((double) tab(0, 1) + tab(1, 0)) / ((double) (nobs));
        }
    }
    
    //return output
    List output;
    output["sens"] = sens;
    output["spec"] = spec;
    output["ppv"] = ppv;
    output["npv"] = npv;
    output["miss"] = miss;
    
    return(output);
}

