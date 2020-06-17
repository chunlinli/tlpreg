## reference: Li, Shen, Pan "Simultaneous inference of directed relations with interventions"
## author: Chunlin Li (li000007@umn.edu)

dyn.load("tlpreg")


tlpreg0 <- function(y, X, b.init=rep(0,ncol(X)), tau=0.01, gamma=NULL, pen.fac=rep(1,ncol(X)), tol=1e-4, dc.maxit=20, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)

    # initialize gamma sequence 
    if(is.null(gamma)) {
        lambda.max <- max(abs(crossprod(X, y - mean(y))))/n
        gamma <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p < n, 1e-4, .01)*lambda.max), 
                          length.out=100))/tau
    }
    ngamma <- as.integer(length(gamma))

    # initialize working residuals, b matrix, intercept
    if(is.null(b.init)) {
        r <- y
        b <- matrix(0, p, ngamma)
    } else {
        r <- y - X %*% b.init
        b <- matrix(b.init, p, ngamma)
    }
    b0 <- 0

    pen.fac <- as.integer(pen.fac)

    # call regularized version
    .Call('tlpreg_r', y, X, b0, b, r, xtx, n, p, tau, gamma, ngamma, pen.fac, tol, as.integer(dc.maxit), as.integer(cd.maxit))

    list(b = b, gamma = gamma)
}





# tlpreg0.old <- function(y, X, b.init=rep(0,ncol(X)), tau=0.01, gamma=0.5, pen.fac=rep(1,ncol(X)), tol=1e-5, dc.maxit=20, cd.maxit=1e+4) {

#     n <- nrow(X)
#     p <- ncol(X)
#     xtx <- colSums(X*X)

#     .C('tlpreg0', y = as.double(y),
#                    X = as.double(X),
#                    b0 = as.double(numeric(1)),
#                    b = as.double(b.init),
#                    r = as.double(y),
#                    xtx = as.double(xtx),
#                    n = as.integer(n),
#                    p = as.integer(p),
#                    tau = as.double(tau),
#                    gamma = as.double(gamma),
#                    pen_fac = as.integer(pen.fac),
#                    tol = as.double(tol),
#                    dc_maxit = as.integer(dc.maxit),
#                    cd_maxit = as.integer(cd.maxit))$b
                   
# }



