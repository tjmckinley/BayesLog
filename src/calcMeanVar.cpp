#include "functions.h"

//function to calculate posterior mean and variance recursively
void calcMeanVar(int i, int ninitial, NumericVector *tempmn, NumericVector *tempvar, IntegerVector *tempcounts, NumericMatrix posterior, int postelement, int parelement, int indelement)
{
    int j, k, m;

    //update means and variances for adaptive proposal
    if((i + 1) == ninitial)
    {
        //calculate mean
        for(k = 0; k <= i; k++)
        {
            if(posterior(k, indelement) == 1)
            {
                (*tempcounts)[parelement]++;
                (*tempmn)[parelement] += posterior(k, postelement);
            }
        }
        if((*tempcounts)[parelement] > 1) (*tempmn)[parelement] = (*tempmn)[parelement] / ((double) (*tempcounts)[parelement]);
        else (*tempmn)[parelement] = 0.0;
        
        //calculate variance
        (*tempvar)[parelement] = 0.0;
        for(k = 0; k <= i; k++) if(posterior(k, indelement) == 1) (*tempvar)[parelement] += pow(posterior(k, postelement), 2.0);
        if((*tempcounts)[parelement] > 1)
        {
            (*tempvar)[parelement] -= (*tempcounts)[parelement] * pow((*tempmn)[parelement], 2.0);
            (*tempvar)[parelement] = (*tempvar)[parelement] / ((double) ((*tempcounts)[parelement] - 1.0));
        }
        else (*tempvar)[parelement] = pow(0.1, 2.0);
    }
    else
    {
        //start recursively updating variance
        if(posterior(i, indelement) == 1)
        {
            (*tempcounts)[parelement]++;
            if((*tempcounts)[parelement] > 2) (*tempvar)[parelement] = ((*tempvar)[parelement] * ((*tempcounts)[parelement] - 2.0)) + (((*tempcounts)[parelement] - 1.0) * pow((*tempmn)[parelement], 2.0));
        }
        //now update mean
        if((*tempcounts)[parelement] > 2) (*tempmn)[parelement] = (((*tempcounts)[parelement] - 1.0) * (*tempmn)[parelement] + posterior(i, postelement)) / ((double) (*tempcounts)[parelement]);
        else
        {
            if((*tempcounts)[parelement] == 2)
            {
                //calculate mean
                for(k = 0; k <= i; k++)
                {
                    if(posterior(k, indelement) == 1) (*tempmn)[parelement] += posterior(k, postelement);
                }
                (*tempmn)[parelement] = (*tempmn)[parelement] / ((double) (*tempcounts)[parelement]);
            }
            else (*tempmn)[parelement] = 0.0;
        }
        //now update variance
        if(posterior(i, indelement) == 1)
        {
            if((*tempcounts)[parelement] > 2)
            {
                (*tempvar)[parelement] += pow(posterior(i, postelement), 2.0) - (*tempcounts)[parelement] * pow((*tempmn)[parelement], 2.0);
                (*tempvar)[parelement] = (*tempvar)[parelement] / ((double) (*tempcounts)[parelement] - 1);
            }
            else
            {
                if((*tempcounts)[parelement] == 2)
                {
                    (*tempvar)[parelement] = 0.0;
                    for(k = 0; k <= i; k++) if(posterior(k, indelement) == 1) (*tempvar)[parelement] += pow(posterior(k, postelement), 2.0);
                    (*tempvar)[parelement] = (*tempvar)[parelement] - (*tempcounts)[parelement] * pow((*tempmn)[parelement], 2.0);
                    (*tempvar)[parelement] = (*tempvar)[parelement] / ((double) ((*tempcounts)[parelement] - 1));
                }
                else (*tempvar)[parelement] = pow(0.1, 2.0);
            }
        }
    }
    return;
}
