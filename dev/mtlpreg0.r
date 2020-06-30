
dyn.load("tlpreg")

# multiple regression
library(foreach)

mtlpreg0 <- function(Y, X, b.init=NULL, tau=0.5*sqrt(log(p)/n), gamma=NULL, pen.fac=rep(1,ncol(X)), tol=1e-4, dc.maxit=20, cd.maxit=1e+4) {

    # Y is n by q mat, q < n
    # X is n by p mat

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)

    B <- foreach(k = 1:ncol(Y)) %do% {
        # initialize gamma sequence 
        r <- Y[,k] - mean(Y[,k])
        if(is.null(gamma)) {
            lambda.max <- max(abs(crossprod(X, r)))/n
            gamma <- exp(seq(from=log(lambda.max),
                              to=log(ifelse(p < n, 1e-4, .01)*lambda.max), 
                              length.out=100))/tau
        }
        ngamma <- as.integer(length(gamma))

        # initialize working residuals, b matrix, intercept
        if(is.null(b.init)) {
            b <- matrix(0, p, ngamma)
        } else {
            r <- r - X %*% b.init[,k]
            b <- matrix(b.init, p, ngamma)
        }
        b0 <- 0

        pen.fac <- as.integer(pen.fac)

        # call regularized version
        .Call('tlpreg_r', Y[,k], X, b0, b, r, xtx, n, p, tau, gamma, ngamma, pen.fac, tol, as.integer(dc.maxit), as.integer(cd.maxit))

        b
    }

    covar <- crossprod(Y - X %*% B)/n
    list(B = B, covar = covar, gamma = gamma)
}

