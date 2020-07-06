# README: Compare TLP with SCAD and MCP (ncvreg)
library(tictoc)
source("tlpreg.r")
source("tlpreg1.r")
source("cv.tlpreg.r")
source("cv.tlpreg1.r")
source("./dev/mtlpreg0.r")
source("./dev/mtlpreg1.r")

#set.seed(2020)
n <- 500
p <- 500
q <- 50
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- X[,j] + Z
  X[,j] <- (X[,j] - mean(X[,j])) / sd(X[,j])
}
Y <- matrix(rnorm(n*q), n, q)
for (k in 1:q) {
    Y[,k] <- Y[,k] + 1 + 0.5*(X[,1 + k] - X[,2 + k] + X[,10 + k] - X[,15 + k]) 
}

# SPEED TEST
tic("TLP-Regularized")
m1 <- fit.mtlpreg0(X=X, Y=Y, nfold=5)
toc()

tic("L0-Constrained")
m2 <- fit.mtlpreg1(X=X, Y=Y, K=1:10, nfold=5)
toc()


# ACCURACY TEST
m11 <- m1[[1]]
m12 <- m1[[2]]
m13 <- m1[[3]]
m14 <- m1[[4]]
m15 <- m1[[5]]

which(m11$b!=0)
which(m12$b!=0)
which(m13$b!=0)
which(m14$b!=0)
which(m15$b!=0)

m21 <- m2[[1]]
m22 <- m2[[2]]
m23 <- m2[[3]]
m24 <- m2[[4]]
m25 <- m2[[5]]

which(m21$b!=0)
which(m22$b!=0)
which(m23$b!=0)
which(m24$b!=0)
which(m25$b!=0)

# Note: informative variables are X(1+k), X(2+k), X(10+k), X(50+k); k = 1,2,3,4,5.