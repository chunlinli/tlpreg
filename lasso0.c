#include <R.h>
#include <math.h>

void lasso0(double *y, double *X, double *b0, double *b, double *r, double *xtx, int *n, int *p,
            double *lambda, int *nlambda, int *pen_fac, double *tol, int *cd_maxit)
{
    double b_curr[*p];
    double r_init[*n];
    for (int i = 0; i < *n; ++i)
        r_init[i] = r[i];

    for (int k = 0; k < *nlambda; ++k)
    {
        int idx = k * (*p); 

        for (int j = 0; j < *p; ++j)
            b_curr[j] = b[idx + j];

        // coordinate descent iteration
        int it;
        double thresh = lambda[k] * (*n);
        for (it = 0; it < *cd_maxit; ++it)
        {
            // b0 = mean(y - X * b_next)
            *b0 = 0.0;
            for (int i = 0; i < *n; ++i)
                *b0 += r[i];
            *b0 = *b0 / *n;

            for (int j = 0; j < *p; ++j)
            {
                //b.next[j] <- sum(X[,j] * (y - b0 - X[,-j]%*%b.next[-j]))
                b[idx + j] = b[idx + j] * xtx[j];
                for (int i = 0; i < *n; ++i)
                    b[idx + j] += X[j * (*n) + i] * (r[i] - *b0);

                //b.next[j] <- soft(b[j],lambda *n / sum(X[,j]*X[,j]))
                if (pen_fac[j] == 1)
                {
                    if (b[idx + j] > thresh)
                        b[idx + j] -= thresh;
                    else if (b[idx + j] < -thresh)
                        b[idx + j] += thresh;
                    else
                        b[idx + j] = 0.0;
                }
                b[idx + j] = b[idx + j] / xtx[j];

                for (int i = 0; i < *n; ++i)
                    r[i] += X[j * (*n) + i] * (b_curr[j] - b[idx + j]);
            }

            // termination
            double diff = 0.0;
            for (int j = 0; j < *p; ++j)
                diff += fabs(b[idx + j] - b_curr[j]);
            if (diff <= (*p) * (*tol))
                break;

            // update current b
            for (int j = 0; j < *p; ++j)
                b_curr[j] = b[idx + j];
        }

        if (it < *cd_maxit)
        {
            // warm start
            if (k < *nlambda - 1)
            {
                for (int j = 0; j < *p; ++j)
                    b[(k + 1) * (*p) + j] = b[idx + j];
            }
        }
        else 
        {
            Rprintf("Warning: the coordinate descent algorithm does not converge (lambda = %f).\n", lambda[k]);
            // restore r
            if (k < *nlambda - 1)
            {
                for (int i = 0; i < *n; ++i)
                    r[i] = r_init[i];
            }
        }
    }
}
