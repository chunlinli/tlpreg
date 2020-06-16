dyn.load("lasso")

lasso0 <- function(y, X, b.init=NULL, lambda=NULL, pen.fac=NULL, tol=1e-5, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)
    
    if(is.null(lambda)) {
        lambda.max <- max(abs(crossprod(X, y - mean(y))))/n
        lambda <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p < n, 1e-4, .01)*lambda.max), 
                          length.out=100))
    } 
    nlambda <- as.integer(length(lambda))

    #if(is.null(b.init)) {
        r <- y
        b0 <- 0
        b <- matrix(0, p, nlambda)
    #} else {
    #    r <- y - X %*% b.init
    #}
    pen_fac <- rep(TRUE, p)

    .Call('lasso', y, X, b0, b, r, xtx, n, p, lambda, nlambda, pen_fac, tol, as.integer(cd.maxit))

    list(b = b, lambda=lambda)
}
