
#include <algorithm> // transform
#include <cmath>     // fabs
#include <numeric>   // inner_product, accumulate (need replaced with reduce for concurrency)
#include <iostream>  // printf
#include <vector>    // vector

/// default version: coordinate descent algorithm + active set update + SAFE screening rule

/// probably need to declare heap arrays for large-scale problems

/// use __restrict 

void lasso0(const double *__restrict X, double *__restrict b0, double *__restrict b, double *__restrict r, 
            const double *__restrict xtx, const int *__restrict n, const int *__restrict p,
            const double *__restrict lambda, const int *__restrict nlambda, const int *__restrict pen_fac, 
            const double *__restrict tol, const int *__restrict cd_maxit)
{
    // declare local variables
    const int p_ = *p;
    const int n_ = *n;
    const int nlambda_ = *nlambda;
    const double tol_ = *tol;
    const int cd_maxit_ = *cd_maxit;
    const double one_over_n = 1.0 / n_;

    double r0[n_];
    std::copy(r, r + n_, r0);

    int act_set[p_];
    int len_act_set = 0;
    int act_set_curr[p_];
    int len_act_set_curr = 0;

    double b_j;      // working coordinate value
    double b_j_curr; // past coordinate value
    double b_change;

    for (int k = 0; k != nlambda_; ++k) // lambda iteration
    {
        int kp = k * p_;

        // coordinate descent iteration
        double thresh = lambda[k] * n_;
        int it = 0;
        for (; it != cd_maxit_; ++it)
        {
            len_act_set = 0; // empty active set

            // update coordinate
            for (int j = 0; j != p_; ++j)
            {
                // update jth coordinate
                int jn = j * n_;
                b_j_curr = b[kp + j];
                b_j = std::inner_product(r, r + n_, X + jn, b_j_curr * xtx[j]);

                // soft thresholding and record active set
                if (pen_fac[j] != 0)
                {
                    if (b_j > thresh)
                    {
                        b_j = (b_j - thresh) / xtx[j];
                        act_set[len_act_set++] = j;
                    }
                    else if (b_j < -thresh)
                    {
                        b_j = (b_j + thresh) / xtx[j];
                        act_set[len_act_set++] = j;
                    }
                    else
                        b_j = 0.0;
                }
                else
                {
                    b_j = b_j / xtx[j];
                    act_set[len_act_set++] = j; // act_set
                }

                // update residuals
                b_change = b_j - b_j_curr; // change in b
                if (b_change != 0.0)
                {
                    std::transform(r, r + n_, X + jn, r,
                                   [b_change](double x1, double x2) { return x1 - x2 * b_change; }); // replaced with fma

                    b[kp + j] = b_j;
                }
            }

            // update intercept and residuals
            b_change = std::accumulate(r, r + n_, 0.0) * one_over_n;
            std::transform(r, r + n_, r, [b_change](double x1) { return x1 - b_change; });
            *b0 += b_change;

            // active set update: start from the 2nd iteration, seems a good strategy
            if (it != 0)
            {

                for (int it_act = 0; it_act != cd_maxit_; ++it_act)
                {
                    double diff = 0.0;
                    // update active coordinates
                    for (int m = 0; m != len_act_set; ++m)
                    {
                        int j = act_set[m];
                        int jn = j * n_;
                        b_j_curr = b[kp + j];
                        b_j = std::inner_product(r, r + n_, X + jn, b_j_curr * xtx[j]);

                        // soft thresholding and record active set
                        if (pen_fac[j] != 0)
                        {
                            if (b_j > thresh)
                                b_j = (b_j - thresh) / xtx[j];
                            else if (b_j < -thresh)
                                b_j = (b_j + thresh) / xtx[j];
                            else
                                b_j = 0.0;
                        } else {
                            b_j = b_j / xtx[j];
                        }

                        b_change = b_j - b_j_curr; // change in b
                        if (b_change != 0.0)
                        {
                            diff += fabs(b_change);
                            std::transform(r, r + n_, X + jn, r,
                                           [b_change](double x1, double x2) { return x1 - x2 * b_change; });
                            b[kp + j] = b_j;
                        }
                    }

                    // update intercept and residuals
                    b_change = std::accumulate(r, r + n_, 0.0) * one_over_n;
                    std::transform(r, r + n_, r, [b_change](double x1) { return x1 - b_change; });
                    *b0 += b_change;

                    if (diff <= len_act_set * tol_)
                        break;
                }
            }

            // terminate if active set does not change
            if (len_act_set == len_act_set_curr && std::equal(act_set, act_set + len_act_set, act_set_curr))
                break;

            // update active set
            len_act_set_curr = len_act_set;
            std::copy(act_set, act_set + len_act_set, act_set_curr);
        }

        if (it != cd_maxit_)
        {
            // previous estimate converges
            // warm start
            if (k != nlambda_ - 1)
                std::copy(b + kp, b + kp + p_, b + kp + p_);
        }
        else
        {
            printf("Warning: the coordinate descent algorithm does not converge (lambda = %f).\n", lambda[k]);

            // previous estimate does not converge
            // cold start: restore residuals, restore active set
            if (k != nlambda_ - 1)
            {
                std::copy(r0, r0 + n_, r);
                len_act_set_curr = 0;
            }
        }
    }
}
