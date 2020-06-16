dyn.load("lasso")

lasso0 <- function(y, X, b.init=NULL, lambda=NULL, pen.fac=NULL, tol=1e-5, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)
    
    if(is.null(lambda)) {
        lambda.max <- max(abs(crossprod(X, y - mean(y))))/n
        lambda <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p<n, 1e-4, .01)*lambda.max), 
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
    cd_maxit <- 1e+4

    .Call('lasso', y, X, b0, b, r, xtx, n, p, lambda, nlambda, pen_fac, tol, as.integer(cd_maxit))

    b
}


# lasso0.old <- function(y, X, b.init=NULL, lambda=NULL, pen.fac=rep(1,ncol(X)), tol=1e-5, cd.maxit=1e+4) {

#     n <- nrow(X)
#     p <- ncol(X)
#     xtx <- colSums(X*X)
    
#     if(is.null(lambda)) {
#         lambda.max <- max(abs(crossprod(X, y - mean(y))))/n
#         lambda <- exp(seq(from=log(lambda.max),
#                           to=log(ifelse(p<n, 1e-3, .01)*lambda.max), 
#                           length.out=100))
#     } 
#     nlambda <- length(lambda)

#     #if(is.null(b.init)) {
#     #    r <- y
#         b.init <- rep(0, p)
#     #} else {
#     #    r <- y - X %*% b.init
#     #}

#     .C('lasso0', y = as.double(y),
#                 X = as.double(X),
#                 b0 = as.double(numeric(1)),
#                 b = as.double(rep(b.init,p*nlambda)),
#                 r = as.double(y),
#                 xtx = as.double(xtx),
#                 n = as.integer(n),
#                 p = as.integer(p),
#                 lambda = as.double(lambda),
#                 nlambda = as.integer(nlambda),
#                 pen_fac = as.integer(pen.fac),
#                 tol = as.double(tol),
#                 cd_maxit = as.integer(cd.maxit))$b

# }


