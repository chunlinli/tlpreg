library(foreach)
#library(doParallel)

cv.lasso <- function(y, X, b.init = NULL, pen.fac = rep(1,p), lambda = NULL, 
                     nfold=5, tol=1e-4, cd.maxit=1e+4) {  

    #cl <- makeCluster(5)
    #registerDoParallel(cl)

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
    cv <- foreach(fold = 1:nfold, .combine = "cbind") %do% {
        #source("lasso.r")
        b <- lasso0(y=y[obs.fold!=fold], X=X[obs.fold!=fold,], b.init=b.init, lambda=lambda, 
                    pen.fac=pen.fac, tol=tol, cd.maxit=cd.maxit)$b

        apply(b, 2, function(beta) {
            res <- (y[obs.fold==fold] - X[obs.fold==fold,] %*% beta)
            mean((res - mean(res))^2)
        })
    }

    #stopCluster(cl)

    # output CV graph
    cvm <- rowMeans(cv)



    list(lambda.min=lambda[which.min(cvm)], cvm = cvm, lambda=lambda)
}
