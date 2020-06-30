

#ifndef _TLPREG_H_
#define _TLPREG_H_

void lasso0(double *y, double *X, double *b0, double *b, double *r, double *xtx, int *n, int *p,
            double *lambda, int *nlambda, int *pen_fac, double *tol, int *cd_maxit);

void lasso1(double *y, double *X, double *b0, double *b, double *r_b, double *xtx, int *n, int *p,
            double *lambda, int *nlambda, int *pen_fac, double *tol, int *cd_maxit);

void tlpreg0(double *y, double *X, double *b0, double *b, double *r, double *xtx, int *n, int *p, 
             double *tau, double *gamma, int* ngamma, int *pen_fac, double *tol, int *dc_maxit, int *cd_maxit);

#endif

