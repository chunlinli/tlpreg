#include <algorithm>
#include <cmath>
#include <numeric>
#include <queue>
#include "tlpreg.h"

// constrained version

void tlpreg1(const double *__restrict X, double *__restrict b0, double *__restrict b, double *__restrict r,
             const double *__restrict xtx, const int *__restrict n, const int *__restrict p,
             const double *__restrict tau, const int *__restrict K, const int *__restrict nK,
             const double *__restrict gamma, const int *__restrict ngamma, const int *__restrict pen_fac,
             const double *__restrict tol, const int *__restrict dc_maxit, const int *__restrict cd_maxit)
{
    const int p_ = *p;
    const int n_ = *n;
    const int nK_ = *nK;
    const int K_max = K[nK_ - 1];
    const int ngamma_ = *ngamma;
    const double tau_ = *tau;
    const double tol_ = *tol;
    const int cd_maxit_ = *cd_maxit;
    const int dc_maxit_ = *dc_maxit;
    const int one_over_n = 1.0 / n_;
    const int nlambda = 1;

    double lambda[ngamma_];
    for (int k = 0; k != ngamma_; ++k)
        lambda[k] = gamma[k] * tau_;
    
    int pen[p_];
    int pen0[p_];

    double b_change;
    double b_j;
    double b_j_curr;

    // for regularized DC iteration
    double b_dc[p_];
    std::copy(b, b + p_, b_dc);
    double b0_dc = *b0;

    double b0_warm_init;
    double b_warm_init[p_];
    //std::copy(b, b + p_, b_warm_init);
    double r_warm_init[n_];
    //std::copy(r, r + n_, r_warm_init);

    // for constrained projection 
    int act_set[p_];
    int len_act_set;
    double rss[nK_];
    std::fill(rss, rss + nK_, INFINITY);

    double b0_proj_init = *b0;
    double b_proj_init[p_];
    std::fill(b_proj_init, b_proj_init + p_, 0.0);
    double b0_proj;
    double b_proj[p_];

    double r_proj[n_];
    double r_proj_init[n_];
    std::copy(r, r + n_, r_proj_init);

    for (int k = 0; k != ngamma_; ++k)
    {
        if (k != 0)
        {
            b0_dc = b0_warm_init;
            std::copy(b_warm_init, b_warm_init + p_, b_dc);
            std::copy(r_warm_init, r_warm_init + n_, r);
        }

        for (int j = 0; j != p_; ++j)
            pen[j] = (fabs(b_dc[j]) < tau_ ? pen_fac[j] : 0);

        std::copy(pen, pen + p_, pen0);

        
        lasso0(X, &b0_dc, b_dc, r, xtx, n, p, &lambda[k], &nlambda, pen, tol, cd_maxit);

        if (k != ngamma_ - 1)
        {
            b0_warm_init = b0_dc;
            std::copy(b_dc, b_dc + p_, b_warm_init);
            std::copy(r, r + n_, r_warm_init);
        }

        for (int it = 1; it != dc_maxit_; ++it)
        {
            for (int j = 0; j != p_; ++j)
                pen[j] = (fabs(b_dc[j]) < tau_ ? pen_fac[j] : 0);

            // termination: active set
            if (std::equal(pen, pen + p_, pen0))
                break;

            std::copy(pen, pen + p_, pen0);
            lasso0(X, &b0_dc, b_dc, r, xtx, n, p, &lambda[k], &nlambda, pen, tol, cd_maxit);
        }

        // next, project b to lower dimensional constrained set
        // sort the coef of b: use priority queue
        std::priority_queue<std::pair<double, int>> queue;
        for (int j = 0; j != p_; ++j)
        {
            double fabs_b_j = fabs(b_dc[j]);
            if (fabs_b_j >= tol_)
            {
                queue.push(std::pair<double, int>(fabs_b_j, j));
            }
        }
        int M = K_max < queue.size() ? K_max : queue.size(); // M is upper bound of size of candidates
        for (len_act_set = 0; len_act_set != M; ++len_act_set)
        {
            act_set[len_act_set] = queue.top().second;
            queue.pop();
        }
        
        std::copy(r_proj_init, r_proj_init + n_, r_proj);
        std::copy(b_proj_init, b_proj_init + p_, b_proj);
        b0_proj = b0_proj_init;

        int s = 0;
        int K_s = K[0]; // number of nonzero coef
        while (K_s <= M) 
        {
            // coordinate descent
            for (int it_act = 0; it_act != cd_maxit_; ++it_act)
            {
                double diff = 0.0;
                // update active coordinates
                for (int m = 0; m != K_s; ++m)
                {
                    int j = act_set[m];
                    int jn = j * n_;
                    b_j_curr = b_proj[j];
                    b_j = std::inner_product(r_proj, r_proj + n_, X + jn, b_j_curr * xtx[j]) / xtx[j];
                    b_change = b_j - b_j_curr; // change in b
                    diff += fabs(b_change);
                        std::transform(r_proj, r_proj + n_, X + jn, r_proj,
                                       [b_change](double x1, double x2) { return x1 - x2 * b_change; });
                    b_proj[j] = b_j;
                }

                // update unpenalized coef

                // update intercept and residuals
                b_change = std::accumulate(r_proj, r_proj + n_, 0.0) * one_over_n;
                std::transform(r_proj, r_proj + n_, r_proj, [b_change](double x1) { return x1 - b_change; });
                b0_proj += b_change;

                if (diff <= len_act_set * tol_)
                    break;
            }

            double loss = 0.0;
            for(int i = 0; i != n_; ++i) 
            {
                loss += r_proj[i] * r_proj[i];
            }
            if (loss < rss[s]) {
                std::copy(b_proj, b_proj + p_, b + s * p_);
                rss[s] = loss;
            }

            ++s;
            if (s == nK_) {
                break;
            }
            K_s = K[s];
        }
    }
}
