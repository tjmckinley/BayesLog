#include "functions.hpp"

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

//function to calculate means and covariance matrices for adaptive MCMC
void adapt_update(int i, int ninitial, int niter, int npars, arma::vec *tempmn, arma::mat *meanmat, arma::mat *meanmat1, NumericMatrix posterior, arma::mat *propcov, int subrow, arma::ivec *elements)
{
    int j, k, m;

    //update var-covariances for adaptive proposal
    if((i + 1) == ninitial)
    {
        //calculates variance-covariance matrix
        //first: means
        for(j = 0; j < npars; j++)
        {
            (*tempmn)(j) = 0;
            for(k = 0; k <= i; k++) (*tempmn)(j) += posterior(k, (*elements)(j));
            (*tempmn)(j) = (*tempmn)(j) / ((double) i + 1);
        }
        //matrix of product of means
        for(j = 0; j < npars; j++) for(k = 0; k < npars; k++) (*meanmat)(j, k) = (*tempmn)(j) * (*tempmn)(k);
        (*propcov).zeros();
        for(j = 0; j < npars; j++)
        {
            for(k = 0; k < npars; k++)
            {
                for(m = 0; m <= i; m++) (*propcov)(j, k) += posterior(m, (*elements)(j)) * posterior(m, (*elements)(k));
                (*propcov)(j, k) -= (i + 1) * (*meanmat)(j, k);
                (*propcov)(j, k) = (*propcov)(j, k) / ((double) i);
            }
        }
    }
    else
    {
        //recursively update mean and covariance matrix
        for(j = 0; j < npars; j++) (*tempmn)(j) = ((*tempmn)(j) * i + posterior(subrow, (*elements)(j))) / ((double) i + 1);
        //new matrix of product of means
        for(j = 0; j < npars; j++) for(k = 0; k < npars; k++) (*meanmat1)(j, k) = (*tempmn)(k) * (*tempmn)(j);
        for(j = 0; j < npars; j++)
        {
            for(k = 0; k < npars; k++)
            {
                (*propcov)(j, k) = (((double) (i - 1)) / ((double) i)) * (*propcov)(j, k) + (1.0 / ((double) i)) * (i * (*meanmat)(j, k) - (i + 1) * (*meanmat1)(j, k) + posterior(subrow, (*elements)(j)) * posterior(subrow, (*elements)(k)));
            }
        }
        for(j = 0; j < npars; j++) for(k = 0; k < npars; k++) (*meanmat)(j, k) = (*meanmat1)(j, k);
    }
//    if(npars > 1 && (i + 1) % 100 == 0)
//    {
//        Rprintf("\nPropcov:\n");
//        for(k = 0; k < npars; k++)
//        {
//            for(j = 0; j < npars; j++) Rprintf("%f ", (*propcov)(k, j) / adaptscale);
//            Rprintf("\n");
//        }
//    }
    return;
}

