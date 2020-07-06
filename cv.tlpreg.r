library(foreach)
#library(doParallel)

# regularized version 
cv.tlpreg0 <- function(y, X, b.init = NULL, pen.fac = rep(1,ncol(X)), 
                      tau=0.5*sqrt(log(p)/n), gamma=NULL, nfold=10, 
                      tol=1e-4, dc.maxit=as.integer(max(5,1+log2(n/log(p)))), cd.maxit=1e+4) {
    #cl <- makeCluster(5)
    #registerDoParallel(cl)

    n <- nrow(X)
    p <- ncol(X)
    if (n < nfold) stop("nfold cannot be larger than n.")

    if(is.null(gamma)) {
        gamma <- exp(seq(from=log(max(abs(crossprod(X, y - mean(y))))/n),
                          to=log(ifelse(p < n, .001, .08)), 
                          length.out=100))/tau
    } 
    
    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cv <- foreach(fold = 1:nfold, .combine = "cbind") %do% {
        #source("tlpreg.r")
        b <- tlpreg0(y=y[obs.fold!=fold], X=X[obs.fold!=fold,], b.init = b.init, 
                     tau=tau, gamma = gamma, pen.fac = pen.fac, 
                     tol = tol, dc.maxit = dc.maxit, cd.maxit = cd.maxit)$b
        
        apply(b, 2, function(beta) {
            res <- (y[obs.fold==fold] - X[obs.fold==fold,] %*% beta)
            mean((res - mean(res))^2)
        })
    }

    #stopCluster(cl)

    # output CV graph
    cvm <- rowMeans(cv)
    list(gamma.min=gamma[which.min(cvm)], cvm = cvm, gamma=gamma)
}

