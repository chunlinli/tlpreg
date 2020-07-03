#include <cmath>
#include <iostream>
#include "tlpreg.h"

// default regularized version: difference of convex functions algorithm + active set termination

void tlpreg0(const double *__restrict X, double *__restrict b0, double *__restrict b, double *__restrict r, 
             const double *__restrict xtx, const int *__restrict n, const int *__restrict p, 
             const double *__restrict tau, const double *__restrict gamma, const int *__restrict ngamma, 
             const int *__restrict pen_fac, const double *__restrict tol, const int *__restrict dc_maxit, const int *__restrict cd_maxit)
{
    const int p_ = *p;
    const int n_ = *n;
    const int ngamma_ = *ngamma;
    const double tau_ = *tau;
    //const double tol_ = *tol;
    const int dc_maxit_ = *dc_maxit;
    const int nlambda = 1;

    //double b_curr[p_];
    double b0_warm_init;
    double b_warm_init[p_];
    //std::copy(b, b + p_, b_curr);
    //for (int j = 0; j != p_; ++j)
    //    b_curr[j] = b[j];
    std::copy(b, b + p_, b_warm_init);
    //for (int j = 0; j < *p; ++j)
    //    b_warm_init[j] = b[j];

    double r_warm_init[n_];
    std::copy(r, r + n_, r_warm_init);
    //for (int i = 0; i < *n; ++i)
    //    r_warm_init[i] = r[i];
    
    double lambda[ngamma_];
    for (int k = 0; k != ngamma_; ++k)
        lambda[k] = gamma[k] * tau_;

    int pen[p_];
    int pen0[p_];

    for (int k = 0; k != ngamma_; ++k) 
    {
        int kp = k * p_;
        
        if (k != 0) 
        {
            *b0 = b0_warm_init;
            std::copy(b_warm_init, b_warm_init + p_, b + kp);
            //for (int j = 0; j < *p; ++j)
            //    b[k * (*p) + j] = b_warm_init[j];
            std::copy(r_warm_init, r_warm_init + n_, r);
            //for (int i = 0; i < *n; ++i)
            //    r[i] = r_warm_init[i];
        }
        
        
        for (int j = 0; j != p_; ++j)
            pen[j] = (fabs(b[kp + j]) < tau_ ? pen_fac[j] : 0);

        std::copy(pen, pen + p_, pen0);
        
        lasso0(X, b0, b + kp, r, xtx, n, p, &lambda[k], &nlambda, pen, tol, cd_maxit); 
        
        if (k != ngamma_ - 1) {
            // store b0, b, r
            b0_warm_init = *b0;
            std::copy(b + kp, b + kp + p_, b_warm_init);
            std::copy(r, r + n_, r_warm_init);
        }
        
        //std::copy(b + kp, b + kp + p_, b_curr);
        
        //for (int j = 0; j != p_; ++j)
        //{
        //    b_curr[j] = b[kp + j];
        //    b_warm_init[j] = b[kp + j];
        //}
        
        //for (int i = 0; i < *n; ++i)
        //    r_warm_init[i] = r[i];

        for (int it = 1; it != dc_maxit_; ++it)
        {
            for (int j = 0; j != p_; ++j)
                pen[j] = (fabs(b[kp + j]) < tau_ ? pen_fac[j] : 0);
            
            // termination: active set
            if (std::equal(pen, pen + p_, pen0)) 
                break;

            //std::copy(b + kp, b + kp + p_, b_curr);
            std::copy(pen, pen + p_, pen0);

            lasso0(X, b0, b + kp, r, xtx, n, p, &lambda[k], &nlambda, pen, tol, cd_maxit); 

            // termination: L1 norm difference
            //double diff_max = 0.0;
            //double b_change = 0.0;
            //for (int j = 0; j != p_; ++j) 
            //{
            //    b_change = fabs(b[kp + j] - b_curr[j]);
            //    diff_max = (b_change > diff_max ? b_change : diff_max);
            //}
            //if (diff_max <= tol_)
            //    break;

            // update current b
            
            //for (int j = 0; j != p_; ++j)
            //    b_curr[j] = b[kp + j];
        }

        //if (it == dc_maxit_ || k % 10 == 0)
            //printf("Warning: the difference of convex functions algorithm does not converge (gamma = %f).\n", gamma[k]);
    }
}



