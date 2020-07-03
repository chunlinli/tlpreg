#include <R.h>
#include <Rinternals.h>
#include "tlpreg.h"

extern "C" {

SEXP tlpreg_c(SEXP X, SEXP b0, SEXP b, SEXP r,
           SEXP xtx, SEXP n, SEXP p, SEXP tau, SEXP K, SEXP nK, 
           SEXP gamma, SEXP ngamma,
           SEXP pen_fac, SEXP tol, SEXP dc_maxit, SEXP cd_maxit)
{
    double *X_ptr;
    double *b0_ptr;
    double *b_ptr;
    double *r_ptr;
    double *xtx_ptr;
    int    *n_ptr;
    int    *p_ptr;
    double *tau_ptr;
    int    *K_ptr;
    int    *nK_ptr;
    double *gamma_ptr;
    int    *ngamma_ptr;
    int    *pen_fac_ptr;
    double *tol_ptr;
    int    *dc_maxit_ptr;
    int    *cd_maxit_ptr;

    X_ptr        = REAL(X);
    b0_ptr       = REAL(b0);
    b_ptr        = REAL(b);
    r_ptr        = REAL(r);
    xtx_ptr      = REAL(xtx);
    n_ptr        = INTEGER(n);
    p_ptr        = INTEGER(p);
    tau_ptr      = REAL(tau);
    K_ptr        = INTEGER(K);
    nK_ptr       = INTEGER(nK);
    gamma_ptr    = REAL(gamma);
    ngamma_ptr   = INTEGER(ngamma);
    pen_fac_ptr  = INTEGER(pen_fac);
    tol_ptr      = REAL(tol);
    dc_maxit_ptr = INTEGER(dc_maxit);
    cd_maxit_ptr = INTEGER(cd_maxit);

    tlpreg1(X_ptr, b0_ptr, b_ptr, r_ptr, xtx_ptr, n_ptr, p_ptr, tau_ptr, K_ptr, nK_ptr, gamma_ptr, ngamma_ptr, pen_fac_ptr, tol_ptr, dc_maxit_ptr, cd_maxit_ptr);

    return R_NilValue;
}

} // extern "C"