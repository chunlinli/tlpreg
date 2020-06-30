#dyn.load("./dev/tlpreg")
dyn.load("./tlpreg")

# regularized version

tlpreg0 <- function(y, X, b.init=NULL, tau=0.5*sqrt(log(p)/n), gamma=NULL, pen.fac=rep(1,ncol(X)), 
                    tol=1e-4, dc.maxit=as.integer(max(5,log2(n/log(p)))), cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)

    # initialize gamma sequence 
    b0 <- mean(y)
    r <- y - b0
    if(is.null(gamma)) {
        lambda.max <- max(abs(crossprod(X, r)))/n
        gamma <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p < n, .001, .08)), 
                          length.out=100))/tau
    }
    ngamma <- as.integer(length(gamma))

    # initialize working residuals, b matrix, intercept
    if(is.null(b.init)) {
        b <- matrix(0, p, ngamma)
    } else {
        r <- r - X %*% b.init
        b <- matrix(b.init, p, ngamma)
    }

    pen.fac <- as.integer(pen.fac)

    # call regularized version
    .Call('tlpreg_r', X, b0, b, r, xtx, n, p, tau, gamma, ngamma, pen.fac, tol, as.integer(dc.maxit), as.integer(cd.maxit))

    list(b = b, gamma = gamma)
}
