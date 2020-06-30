#dyn.load("./dev/tlpreg")
dyn.load("./tlpreg")

lasso0 <- function(y, X, b.init=NULL, lambda=NULL, pen.fac=rep(1,p), tol=1e-4, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)
    
    b0 <- mean(y)
    r <- y - b0
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
    
    pen.fac <- as.integer(pen.fac)

    .Call('lasso', X, b0, b, r, xtx, n, p, lambda, nlambda, pen.fac, tol, as.integer(cd.maxit))

    list(b = b, b0 = b0, lambda = lambda)
}

# b0 should be a vec: length(b0) == length(lambda)