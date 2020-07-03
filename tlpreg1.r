
dyn.load("tlpreg")

tlpreg1 <- function(y, X, b.init=NULL, tau=0.5*sqrt(log(p)/n), K=NULL, pen.fac=rep(1,ncol(X)), 
                    tol=1e-4, dc.maxit=as.integer(max(5,1+log2(n/log(p)))), cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)

    # initialize K sequence
    if(is.null(K))
        K <- 1:min(p, as.integer(0.5*n/log(p)))
    nK <- as.integer(length(K))

    # initialize gamma sequence 
    b0 <- mean(y)
    r <- y - b0
    gamma <- exp(seq(from=log(max(abs(crossprod(X, r)))/n),
                     to=log(ifelse(p < n, .001, .08)), 
                     length.out=100))/tau
    ngamma <- as.integer(length(gamma))

    # initialize working residuals, b matrix, intercept
    if(is.null(b.init)) {
        b <- matrix(0, p, nK)
    } else {
        r <- r - X %*% b.init
        b <- matrix(b.init, p, nK)
    }

    pen.fac <- as.integer(pen.fac)

    # call regularized version
    .Call('tlpreg_c', X, b0, b, r, xtx, n, p, tau, K, nK, gamma, ngamma, pen.fac, tol, as.integer(dc.maxit), as.integer(cd.maxit))

    list(b = b, K = K)
}
