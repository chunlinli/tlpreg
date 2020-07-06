dyn.load("tlpreg")

# multiple regression
library(foreach)

fit.mtlpreg1 <- function(Y, X, b.init = NULL, pen.fac = rep(1,ncol(X)), 
                      tau=0.5*sqrt(log(p)/n), K=NULL, nfold=10, tol=1e-4, dc.maxit=as.integer(max(5,1+log2(n/log(p)))), cd.maxit=1e+4) {
    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")

    if(is.null(K))
        K <- 1:min(p, as.integer(0.5*n/log(p)))
    nK <- as.integer(length(K))
    
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    foreach(k = 1:ncol(Y)) %do% {

        cv <- foreach(fold = 1:nfold, .combine = "cbind") %do% {
            b <- tlpreg1(y=Y[obs.fold!=fold,k], X=X[obs.fold!=fold,], b.init=NULL, 
                         tau=tau, K=K, pen.fac=pen.fac, 
                         tol=tol, dc.maxit=dc.maxit, cd.maxit=cd.maxit)$b
            apply(b, 2, function(beta) {
                res <- (Y[obs.fold==fold,k] - X[obs.fold==fold,] %*% beta)
                mean((res - mean(res))^2)
            })
        }
        
        tlpreg1(y=Y[,k], X=X, b.init = b.init, pen.fac = pen.fac, 
                      tau=tau, K=K[which.min(rowMeans(cv))], tol=tol, dc.maxit=dc.maxit, cd.maxit=cd.maxit)
    }

}