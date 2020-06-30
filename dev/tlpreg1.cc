#include <R.h>
#include <math.h>
#include "tlpreg.h"

// constrained version

void tlpreg1(double *y, double *X, double *b0, double *b, double *r, double *xtx, int *n, int *p, 
             double *tau, double *gamma, int *ngamma, int *pen_fac, double *tol, int *dc_maxit, int *cd_maxit)
{
    double b_curr[*p];
    for (int j = 0; j < *p; ++j) // may use std::copy instead
        b_curr[j] = b[j];

    double b_warm_init[*p];
    for (int j = 0; j < *p; ++j)
        b_warm_init[j] = b[j];

    double r_warm_init[*n];
    for (int i = 0; i < *n; ++i)
        r_warm_init[i] = r[i];
    
    double lambda[*ngamma];
    for (int k = 0; k < *ngamma; ++k)
        lambda[k] = gamma[k] * (*tau);

    for (int k = 0; k < *ngamma; ++k) // change loop ?
    {
        if (k > 0) 
        {
            for (int j = 0; j < *p; ++j)
                b[k * (*p) + j] = b_warm_init[j];
            for (int i = 0; i < *n; ++i)
                r[i] = r_warm_init[i];
        }
        
        int pen[*p];
        for (int j = 0; j < *p; ++j)
            pen[j] = (fabs(b[k * (*p) + j]) < *tau ? pen_fac[j] : 0);
        
        int nlambda = 1;
        lasso0(y, X, b0, &b[k * (*p)], r, xtx, n, p, &lambda[k], &nlambda, pen, tol, cd_maxit); 

        for (int j = 0; j < *p; ++j)
        {
            b_curr[j] = b[k * (*p) + j];
            b_warm_init[j] = b[k * (*p) + j];
        }
        for (int i = 0; i < *n; ++i)
            r_warm_init[i] = r[i];

        int it;
        for (it = 1; it < *dc_maxit; ++it)
        {
            for (int j = 0; j < *p; ++j)
                pen[j] = (fabs(b[k * (*p) + j]) < *tau ? pen_fac[j] : 0);

            lasso0(y, X, b0, &b[k * (*p)], r, xtx, n, p, &lambda[k], &nlambda, pen, tol, cd_maxit); 

            // termination
            double diff = 0.0;
            for (int j = 0; j < *p; ++j)
                diff += fabs(b[k * (*p) + j] - b_curr[j]);
            if (diff <= (*p) * (*tol))
                break;

            // update current b
            for (int j = 0; j < *p; ++j)
                b_curr[j] = b[k * (*p) + j];
        }

        if (it == *dc_maxit)
            Rprintf("Warning: the difference of convex functions algorithm does not converge (gamma = %f).\n", gamma[k]);
    }
}
