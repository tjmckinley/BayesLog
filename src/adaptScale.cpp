#include "functions.h"

//function to scale proposal variance
double adapt_scale(int nacc, int niter, double desacc, double propscale, int totiter, double maxscale, double niterdim)
{
    double propscale1 = propscale;
    double accrate = (double) nacc;
    accrate = accrate / ((double) niter);

    if(niter > 0)
    {
        propscale1 = (double) totiter;
        propscale1 = 1.0 / sqrt(propscale1);
        propscale1 = ((1.0 / niterdim) < propscale1 ? (1.0 / niterdim):propscale1);
        propscale1 = niterdim * propscale1;
        
        if(accrate <= desacc)
        {
            propscale1 = pow((maxscale - (accrate / desacc)), propscale1);
            propscale1 = propscale / propscale1;
        }
        else
        {
            propscale1 = pow((maxscale - ((1.0 - accrate) / (1.0 - desacc))), propscale1);
            propscale1 = propscale * propscale1;
        }
    }
    else propscale1 = propscale;
    return propscale1;
}

