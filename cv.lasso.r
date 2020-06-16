library(foreach)

cv.lasso <- function(y, X, b.init = NULL, pen.fac = NULL, lambda = NULL, 
                     nfold=10, tol=1e-5, cd.maxit=1e+4) {  
    # tuning lambda
    p <- ncol(X)
    n <- nrow(X)
    if (n < nfold) stop("nfold cannot be larger than n.")

    if(is.null(lambda)) {
        lambda.max <- max(abs(crossprod(X, y - mean(y))))/n
        lambda <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p < n, 1e-4, .01)*lambda.max), 
                          length.out=100))
    } 
    
    # nfold-cross-validation
    obs.fold <- rep(1:nfold,n)[sample.int(n)]
    cvm <- foreach(fold = 1:nfold, .combine = "cbind") %do% {
        b <- lasso0(y[obs.fold!=fold], X[obs.fold!=fold,], 
                   b.init, pen.fac, lambda, tol, cd.maxit)$b

        apply(b, 2, function(beta) mean((y[obs.fold==fold] - X[obs.fold==fold,] %*% beta)^2))
    }

    # output CV graph
    cvm <- rowMeans(cvm)
    list(lambda.min=lambda[which.min(cvm)], cvm = cvm, lambda=lambda)
}
