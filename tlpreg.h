
#ifndef _TLPREG_H_
#define _TLPREG_H_

void lasso0(const double *X, double *b0, double *b, double *r, const double *xtx, const int *n, const int *p,
            const double *lambda, const int *nlambda, const int *pen_fac, const double *tol, const int *cd_maxit);

void tlpreg0(const double *X, double *b0, double *b, double *r, const double *xtx, const int *n, const int *p, 
             const double *tau, const double *gamma, const int* ngamma, const int *pen_fac, 
             const double *tol, const int *dc_maxit, const int *cd_maxit);

//void tlpreg1(double *y, const double *X, double *b, double *r, double *xtx, int *n, int *p, 
//             double *tau, double *gamma, int *ngamma, int *pen_fac, double *tol, int *dc_maxit, int *cd_maxit);

#endif

