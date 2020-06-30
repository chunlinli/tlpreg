

library(foreach)


cv.tlpreg1 <- function(y, X, b.init = NULL, pen.fac = rep(1,ncol(X)), 
                      tau=0.4*sqrt(log(p)/n), K=NULL, nfold=10, tol=1e-4, dc.maxit=15, cd.maxit=1e+4) {
    
    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")

    if(is.null(K))
        K <- 1:min(p, as.integer(n/log(p)))

    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cvm <- foreach(fold = 1:nfold, .combine = "cbind") %do% {
        b <- tlpreg1(y=y[obs.fold!=fold], X=X[obs.fold!=fold,], b.init=NULL, 
                     tau=tau, K=K, pen.fac=pen.fac, 
                     tol=tol, dc.maxit=dc.maxit, cd.maxit=cd.maxit)$b
        apply(b, 2, function(beta) {
            res <- (y[obs.fold==fold] - X[obs.fold==fold,] %*% beta)
            mean((res - mean(res))^2)
        })
    }

    # output CV graph
    cvm <- rowMeans(cvm)
    list(K.min=K[which.min(cvm)], cvm = cvm, K=K)
}
