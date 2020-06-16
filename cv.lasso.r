

cv.lasso <- function(y, X, b.init = rep(0,ncol(X)), pen.fac = rep(1,ncol(X)), lambda = NULL, 
                     nfold=10, tol=1e-5, cd.maxit=1e+4) {  
    # tuning lambda
    p <- ncol(X)
    n <- nrow(X)
    if (n < nfold) stop("nfold cannot be larger than n.")
    if(is.null(lambda)) {
        lambda.max <- max(abs(t(X)%*%(y-mean(y))))/n
        lambda <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p<n, 1e-4, .01)*lambda.max), 
                          length.out=100))
    } 

    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cvm <- matrix(0,length(lambda),nfold)
    for(k in 1:nfold) {
        X.tr <- X[obs.fold!=k,]
        y.tr <- y[obs.fold!=k]
        X.te <- X[obs.fold==k,]
        y.te <- y[obs.fold==k]
        
        # warm starts CV: C code 
        for(j in 1:length(lambda)) {
            b <- lasso(y.tr, X.tr, b.init=b.init, lambda=lambda[j], tol=tol, cd.maxit=cd.maxit)
            b.init <- b
            cvm[j,k] <- mean((y.te - X.te%*%b)^2)
        }
    }
    
    # output CV graph
    cvm <- rowMeans(cvm)
    list(lambda.min=lambda[which.min(cvm)], cvm = cvm, lambda=lambda)
}
