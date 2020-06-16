#ifndef _TLPREG_H_
#define _TLPREG_H_

void lasso0(double *y, double *X, double *b0, double *b, double *r,
           double *xtx, int *n, int *p,
           double *lambda, int *nlambda,
           int *pen_fac, double *tol, int *cd_maxit);


#endif