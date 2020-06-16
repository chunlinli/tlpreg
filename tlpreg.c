#include <R.h>
#include <math.h>

void ctlpreg0(double *y, double *X, double *b0, double *b, double *r, 
              double *xtx, int *n, int *p, double *tau, double *gamma, int *pen_fac, 
              double *tol, int *dc_maxit, int *cd_maxit)
{
    double lambda = (*gamma) * (*tau);

    double b_curr[*p];
    for(int j = 0; j < *p; ++j) 
        b_curr[j] = b[j];

    int it;
    for(it = 0; it < *dc_maxit; ++it) 
    {
        int pen[*p];
        for(int j = 0; j < *p; ++j)
            pen[j] = (fabs(b[j]) < *tau ? pen_fac[j] : 0);
        
        lasso(y, X, b0, b, r, xtx, n, p, &lambda, pen, tol, cd_maxit);

        // termination
        double diff = 0.0;
        for(int j = 0; j < *p; ++j) 
            diff += fabs(b[j] - b_curr[j]);
        if (diff <= (*p) * (*tol))
            break;
        
        // update current b
        for(int j = 0; j < *p; ++j)
            b_curr[j] = b[j];
    }

    if(it == *dc_maxit) 
        Rprintf("Warning: the difference of convex functions algorithm does not converge.");

}

