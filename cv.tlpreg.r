library(foreach)

# regularized version 
cv.tlpreg0 <- function(y, X, b.init = NULL, pen.fac = rep(1,ncol(X)), 
                      tau=0.4*sqrt(log(p)/n), gamma=NULL, nfold=10, tol=1e-4, dc.maxit=15, cd.maxit=1e+4) {
    
    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")

    if(is.null(gamma)) {
        lambda.max <- max(abs(crossprod(X, y - mean(y))))/n
        gamma <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p<n, 1e-4, .01)*lambda.max), 
                          length.out=100))/tau
    } 
    
    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cvm <- foreach(fold = 1:nfold, .combine = "cbind") %do% {
        b <- tlpreg0(y=y[obs.fold!=fold], X=X[obs.fold!=fold,], b.init = b.init, 
                     tau=tau, gamma = gamma, pen.fac = pen.fac, 
                     tol = tol, dc.maxit = dc.maxit, cd.maxit = cd.maxit)$b
        
        apply(b, 2, function(beta) {
            res <- (y[obs.fold==fold] - X[obs.fold==fold,] %*% beta)
            mean((res - mean(res))^2)
        })
    }

    # output CV graph
    cvm <- rowMeans(cvm)
    list(gamma.min=gamma[which.min(cvm)], cvm = cvm, gamma=gamma)
}



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
