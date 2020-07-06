
dyn.load("tlpreg")

# multiple regression
library(foreach)

mtlpreg0 <- function(Y, X, b.init=NULL, tau=0.5*sqrt(log(p)/n), gamma=NULL, pen.fac=rep(1,ncol(X)), 
                    tol=1e-4, dc.maxit=as.integer(max(5,1+log2(n/log(p)))), cd.maxit=1e+4) {

    # Y is n by q mat, q < n
    # X is n by p mat

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)

    B <- foreach(k = 1:ncol(Y), .combine="cbind") %do% {
        # initialize gamma sequence 
        r <- Y[,k] - mean(Y[,k])
        if(is.null(gamma)) {
            gamma <- exp(seq(from=log(max(abs(crossprod(X, r)))/n),
                              to=log(ifelse(p < n, 0.001, .08)), 
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
        .Call('tlpreg_r', X, b0, b, r, xtx, n, p, tau, gamma, ngamma, pen.fac, tol, as.integer(dc.maxit), as.integer(cd.maxit))

        b
    }

    covar <- crossprod(Y - X %*% B)/n
    list(B = B, covar = covar, gamma = gamma)
}


fit.mtlpreg0 <- function(Y, X, b.init = NULL, pen.fac = rep(1,ncol(X)), 
                      tau=0.5*sqrt(log(p)/n), gamma=NULL, nfold=10, 
                      tol=1e-4, dc.maxit=as.integer(max(5,1+log2(n/log(p)))), cd.maxit=1e+4) {
    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")
    
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    foreach(k = 1:ncol(Y)) %do% {
        if(is.null(gamma)) {
            gamma <- exp(seq(from=log(max(abs(crossprod(X, Y[,k] - mean(Y[,k]))))/n),
                             to=log(ifelse(p < n, .001, .08)), 
                             length.out=100))/tau
        } 
        # nfold-cross-validation
        
        cv <- foreach(fold = 1:nfold, .combine = "cbind") %do% {
            #source("tlpreg.r")
            b <- tlpreg0(y=Y[obs.fold!=fold,k], X=X[obs.fold!=fold,], b.init = b.init, 
                         tau=tau, gamma = gamma, pen.fac = pen.fac, 
                         tol = tol, dc.maxit = dc.maxit, cd.maxit = cd.maxit)$b
        
            apply(b, 2, function(beta) {
                res <- (Y[obs.fold==fold,k] - X[obs.fold==fold,] %*% beta)
                mean((res - mean(res))^2)
            })
        }

        tlpreg0(y = Y[,k], X = X, b.init = b.init, pen.fac = pen.fac, 
                tau=tau, gamma=gamma[which.min(rowMeans(cv))], tol=tol, dc.maxit=dc.maxit, cd.maxit=cd.maxit)
    }

}