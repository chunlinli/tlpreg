dyn.load("tlpreg")

lasso0 <- function(y, X, b.init=NULL, lambda=NULL, pen.fac=rep(1,p), tol=1e-4, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)
    
    r <- y - mean(y)
    if(is.null(lambda)) {
        lambda.max <- max(abs(crossprod(X, r)))/n
        lambda <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p < n, 1e-4, .01)*lambda.max), 
                          length.out=100))
    } 
    nlambda <- as.integer(length(lambda))

    if(is.null(b.init)) {
        b <- matrix(0, p, nlambda)
    } else {
        r <- r - X %*% b.init
        b <- matrix(b.init, p, nlambda)
    }
    b0 <- 0
    
    pen.fac <- as.integer(pen.fac)

    .Call('lasso', y, X, b0, b, r, xtx, n, p, lambda, nlambda, pen.fac, tol, as.integer(cd.maxit))

    list(b = b, lambda = lambda)
}
