
#ifndef _TLPREG_H_
#define _TLPREG_H_

void lasso0(const double *__restrict X, double *__restrict b0, double *__restrict b, double *__restrict r, 
            const double *__restrict xtx, const int *__restrict n, const int *__restrict p,
            const double *__restrict lambda, const int *__restrict nlambda, const int *__restrict pen_fac, 
            const double *__restrict tol, const int *__restrict cd_maxit);

void tlpreg0(const double *__restrict X, double *__restrict b0, double *__restrict b, double *__restrict r, 
             const double *__restrict xtx, const int *__restrict n, const int *__restrict p, 
             const double *__restrict tau, const double *__restrict gamma, const int *__restrict ngamma, 
             const int *__restrict pen_fac, const double *__restrict tol, const int *__restrict dc_maxit, const int *__restrict cd_maxit);

void tlpreg1(const double *__restrict X, double *__restrict b0, double *__restrict b, double *__restrict r,
             const double *__restrict xtx, const int *__restrict n, const int *__restrict p,
             const double *__restrict tau, const int *__restrict K, const int *__restrict nK,
             const double *__restrict gamma, const int *__restrict ngamma, const int *__restrict pen_fac,
             const double *__restrict tol, const int *__restrict dc_maxit, const int *__restrict cd_maxit);

#endif

