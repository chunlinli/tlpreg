

dyn.load("tlpreg")

# regularized version

tlpreg0 <- function(y, X, b.init=NULL, tau=0.4*sqrt(log(p)/n), gamma=NULL, pen.fac=rep(1,ncol(X)), tol=1e-4, dc.maxit=15, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)

    # initialize gamma sequence 
    r <- y - mean(y)
    if(is.null(gamma)) {
        lambda.max <- max(abs(crossprod(X, r)))/n
        gamma <- exp(seq(from=log(lambda.max),
                          to=log(ifelse(p < n, 1e-4, .01)*lambda.max), 
                          length.out=100))/tau
    }
    ngamma <- as.integer(length(gamma))

    # initialize working residuals, b matrix, intercept
    if(is.null(b.init)) {
        b <- matrix(0, p, ngamma)
    } else {
        r <- r - X %*% b.init
        b <- matrix(b.init, p, ngamma)
    }
    b0 <- 0

    pen.fac <- as.integer(pen.fac)

    # call regularized version
    .Call('tlpreg_r', y, X, b0, b, r, xtx, n, p, tau, gamma, ngamma, pen.fac, tol, as.integer(dc.maxit), as.integer(cd.maxit))

    list(b = b, gamma = gamma)
}





# constrained version

tlpreg1 <- function(y, X, b.init=NULL, tau=0.4*sqrt(log(p)/n), K=NULL, pen.fac=rep(1,ncol(X)), tol=1e-4, dc.maxit=15, cd.maxit=1e+4) {

    n <- as.integer(nrow(X))
    p <- as.integer(ncol(X))
    xtx <- colSums(X*X)

    # initialize K sequence
    if(is.null(K))
        K <- 1:min(p, as.integer(n/log(p)))

    # initialize gamma sequence 
    r <- y - mean(y)
    lambda.max <- max(abs(crossprod(X, r)))/n
    gamma <- exp(seq(from=log(lambda.max),
                     to=log(ifelse(p < n, 1e-4, .01)*lambda.max), 
                     length.out=100))/tau
    ngamma <- as.integer(length(gamma))

    # initialize working residuals, br matrix, intercept
    if(is.null(b.init)) {
        br <- matrix(0, p, ngamma)
    } else {
        r <- r - X %*% b.init
        br <- matrix(b.init, p, ngamma)
    }
    b0 <- 0

    pen.fac <- as.integer(pen.fac)

    # call regularized version
    .Call('tlpreg_r', y, X, b0, br, r, xtx, n, p, tau, gamma, ngamma, pen.fac, tol, as.integer(dc.maxit), as.integer(cd.maxit))

    # replace with C/C++: qsort or C++ stdlib
    k.br <- colSums(br!=0)
    idx <- apply(abs(br), 2, function(a) order(a, decreasing=T))
    b <- matrix(0, p, length(K))
    for(t in 1:length(K)) {
        loss <- rep(mean((y - mean(y))^2), ncol(br))
        for(l in 1:ncol(br)) {
            act.set <- idx[1:K[t],l]
            Z <- cbind(rep(1,n),X[,act.set])
            a <- solve(crossprod(Z), crossprod(Z, y))
            loss[l] <- mean((y - Z %*% a)^2)
        }
        act.set <- idx[1:K[t],which.min(loss)]
        Z <- cbind(rep(1,n),X[,act.set])
        b[act.set,t] <- solve(crossprod(Z), crossprod(Z, y))[-1]
    }

    list(b = b, K = K, br=br)
}




