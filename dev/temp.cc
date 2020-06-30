
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <iostream> // printf

/// default version: coordinate descent algorithm + active set update

void lasso0(double *y, double *X, double *b0, double *b, double *r_init, double *xtx, int *n, int *p,
            double *lambda, int *nlambda, int *pen_fac, double *tol, int *cd_maxit)
{
    // declare local variables
    double b_j; // working coordinate value
    double b_j_curr; // past coordinate value
    size_t p_ = *p;
    size_t n_ = *n;
    size_t nlambda_ = *nlambda;
    double tol_ = *tol;
    size_t it;
    double one_over_n = 1.0 / n_;
 
    std::vector<double> r(r_init, r_init + n_); // residual

    std::vector<size_t> act_set_curr{}; // active set 
    std::vector<size_t> act_set_0{};

    for (size_t k = 0; k < nlambda_; ++k) // lambda iteration
    {
        size_t kp = k * p_;

        // coordinate descent iteration
        double thresh = lambda[k] * n_;
        for (it = *cd_maxit; it; --it)
        {
            act_set_0.clear(); // empty active set

            // update coordinate
            double b_change;
            for (size_t j = p_; j; --j)
            {
                // update jth coordinate
                size_t jn = j * n_;
                b_j_curr = b[kp + j];
                b_j = std::inner_product(&X[jn], &X[jn + n_], r, b_j_curr * xtx[j]);

                // soft thresholding and record active set
                if (pen_fac[j] == 1)
                {
                    if (b_j > thresh)
                    {
                        b_j = (b_j - thresh) / xtx[j];
                        act_set_0.push_back(j);
                    }
                    else if (b_j < -thresh)
                    {
                        b_j = (b_j + thresh) / xtx[j];
                        act_set_0.push_back(j);
                    }
                    else
                        b_j = 0.0;
                }

                // update residuals
                b_change = b_j - b_j_curr; // change in b

                std::transform(&X[jn], &X[jn + n_], r.begin(), r.begin(),
                               [b_change](double x1, double x2) { return x2 - x1 * b_change; });

                b[kp + j] = b_j;
            }

            // update intercept and residuals
            b_change = std::accumulate(r.begin(), r.end(), 0.0) * one_over_n;
            *b0 += b_change;
            std::transform(r.begin(), r.end(), r.begin(), [b_change](double x2) { return x2 - b_change; });

            // active set update: high accuracy
            // may shuffle the order
            for (size_t it_act_set = *cd_maxit; it_act_set > 0; --it_act_set)
            {
                double diff = 0.0;
                // update active coordinates
                for (auto j : act_set_0)
                {
                    size_t jn = j * n_;
                    b_j_curr = b[kp + j];
                    b_j = std::inner_product(&X[jn], &X[jn + n_], r, b_j_curr * xtx[j]);

                    // soft thresholding and record active set
                    if (pen_fac[j] == 1)
                    {
                        if (b_j > thresh)
                        {
                            b_j = (b_j - thresh) / xtx[j];
                        }
                        else if (b_j < -thresh)
                        {
                            b_j = (b_j + thresh) / xtx[j];
                        }
                        else
                            b_j = 0.0;
                    }

                    
                    b_change = b_j - b_j_curr; // change in b
                    diff += fabs(b_change);

                    // update residuals
                    std::transform(&X[jn], &X[jn + n_], r.begin(), r.begin(),
                                   [b_change](double x1, double x2) { return x2 - x1 * b_change; });

                    // update b
                    b[kp + j] = b_j;
                }

                // update intercept and residuals
                b_change = std::accumulate(r.begin(), r.end(), 0.0) * one_over_n;
                *b0 += b_change;
                std::transform(r.begin(), r.end(), r.begin(), [b_change](double x2) { return x2 - b_change; });

                if (diff <= act_set_0.size() * tol_)
                    break;
            }

            // update working active set: need optimized
            std::vector<size_t> act_set{};
            for (auto j: act_set_0) {
                if (b[kp + j] != 0) {
                    act_set.push_back(j);
                }
            }
            
            // terminate if active set does not change
            if (act_set_0 == act_set_curr)
                break;

            // update active set
            act_set_curr = act_set_0;
        }

        if (it != 0)
        {
            // warm start
            if (k < nlambda_ - 1)
                std::copy(&b[kp], &b[kp + p_], &b[kp + p_]);
        }
        else
        {
            printf("Warning: the coordinate descent algorithm does not converge (lambda = %f).\n", lambda[k]);

            // previous estimate does not converge
            // no warm start: restore residuals
            if (k < *nlambda - 1)
                std::copy(r_init, r_init + n_, r.begin());
        }
    }
}
