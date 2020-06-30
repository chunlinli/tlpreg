
# README: Compare TLP with SCAD and MCP (ncvreg)
library(ncvreg)
source("tlpreg.r")
source("cv.tlpreg.r")

n <- 500
p <- 1000
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- Z + rnorm(n)
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- 1 + 0.5*(X[,1] + X[,2] + X[,10]) + rnorm(n) 


# SPEED TEST
t0 <- Sys.time()
m1 <- cv.ncvreg(X=X, y=y, penalty="SCAD")
mm1 <- ncvreg(X=X, y=y, lambda=m1$lambda.min, penalty="SCAD")
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
m2 <- cv.ncvreg(X=X, y=y, penalty="MCP")
mm2 <- ncvreg(X=X, y=y, lambda=m2$lambda.min, penalty="MCP")
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
m3 <- cv.tlpreg0(X=X, y=y, nfold=10)
mm3 <- tlpreg0(X=X, y=y, gamma=m3$gamma.min)
t1 <- Sys.time()
t1-t0

# ACCURACY TEST
which(mm1$beta[-1]!=0)
which(mm2$beta[-1]!=0)
which(as.numeric(mm3$b)!=0)

# Note: informative variables are X1, X2, X10










######################################################## DEV TEST

# example: test tlpreg
library(mvtnorm)
source("lasso.r") 
source("cv.lasso.r")
source("tlpreg.r")
source("cv.tlpreg.r")

n <- 400
p <- 500
S <- matrix(0.5,p,p) 
diag(S) <- rep(1,p)
X <- rmvnorm(n, sigma=S)
y <- 1 + 0.5*(X[,1] + X[,2]) + rnorm(n)
X <- apply(X, 2, function(z) (z - mean(z))/sd(z))


mm3 <- cv.tlpreg1(y=y, X=X)
m3 <- tlpreg1(y=y,X=X,K=mm3$K.min)
as.numeric(m3$b)

