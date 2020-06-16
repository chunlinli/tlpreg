## reference: Li, Shen, Pan "Simultaneous inference of directed relations with interventions"
## author: Chunlin Li (li000007@umn.edu)

dyn.load("tlpreg.so")


tlpreg0 <- function(y, X, b.init=rep(0,ncol(X)), tau=0.01, gamma=0.5, pen.fac=rep(1,ncol(X)), tol=1e-5, dc.maxit=20, cd.maxit=1e+4) {

    n <- nrow(X)
    p <- ncol(X)
    xtx <- colSums(X*X)

    .C('tlpreg0', y = as.double(y),
                   X = as.double(X),
                   b0 = as.double(numeric(1)),
                   b = as.double(b.init),
                   r = as.double(y),
                   xtx = as.double(xtx),
                   n = as.integer(n),
                   p = as.integer(p),
                   tau = as.double(tau),
                   gamma = as.double(gamma),
                   pen_fac = as.integer(pen.fac),
                   tol = as.double(tol),
                   dc_maxit = as.integer(dc.maxit),
                   cd_maxit = as.integer(cd.maxit))$b
                   
}



