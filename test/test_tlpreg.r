# README: Compare TLP with SCAD and MCP (ncvreg)
library(tictoc)
library(ncvreg)
source("tlpreg.r")
source("cv.tlpreg.r")
source("tlpreg1.r")
source("cv.tlpreg1.r")

set.seed(2020)
n <- 500
p <- 50000
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j] - mean(X[,j])) / sd(X[,j])
}
y <- 1 + 0.5*(X[,1] - X[,2] + X[,10] - X[,50] + X[,200]) + rnorm(n) 

# SPEED TEST
tic("SCAD")
m1 <- cv.ncvreg(X=X, y=y, penalty="SCAD")
mm1 <- ncvreg(X=X, y=y, lambda=m1$lambda.min, penalty="SCAD")
toc()

tic("MCP")
m2 <- cv.ncvreg(X=X, y=y, penalty="MCP")
mm2 <- ncvreg(X=X, y=y, lambda=m2$lambda.min, penalty="MCP")
toc()

tic("TLP-Regularized")
m3 <- cv.tlpreg0(X=X, y=y, nfold=10)
mm3 <- tlpreg0(X=X, y=y, gamma=m3$gamma.min)
toc()

tic("L0-Constrained")
m4 <- cv.tlpreg1(X=X, y=y, nfold=10)
mm4 <- tlpreg1(X=X, y=y, K=m4$K.min)
toc()

# ACCURACY TEST
which(mm1$beta[-1]!=0)
which(mm2$beta[-1]!=0)
which(as.numeric(mm3$b)!=0)
which(as.numeric(mm4$b)!=0)

# Note: informative variables are X1, X2, X10, X50, X200