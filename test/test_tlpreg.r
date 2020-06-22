library(mvtnorm)
library(ncvreg)

# example: test tlpreg
source("lasso.r") 
source("cv.lasso.r")
source("tlpreg.r")
source("cv.tlpreg.r")

n <- 400
p <- 2000
S <- matrix(0.5,p,p) 
diag(S) <- rep(1,p)
X <- rmvnorm(n, sigma=S)
#X <- matrix(rnorm(n*p),n,p)
y <- 1 + 0.5*(X[,1] + X[,2]) + rnorm(n)
X <- apply(X, 2, function(z) (z - mean(z))/sd(z))

t0 <- Sys.time()
m1 <- cv.ncvreg(X=X, y=y, penalty="SCAD")
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
m2 <- cv.ncvreg(X=X, y=y, penalty="MCP")
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
m3 <- cv.tlpreg0(X=X, y=y, nfold=10)
t1 <- Sys.time()
t1-t0

t0 <- Sys.time()
mm1 <- ncvreg(X=X, y=y, lambda=m1$lambda.min, penalty="SCAD")
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
mm2 <- ncvreg(X=X, y=y, lambda=m2$lambda.min, penalty="MCP")
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
mm3 <- tlpreg0(X=X, y=y, gamma=m3$gamma.min)
t1 <- Sys.time()
t1-t0

which(mm1$beta[-1]!=0)
which(mm2$beta[-1]!=0)
which(as.numeric(mm3$b)!=0)







########################################################

# example: test tlpreg
library(mvtnorm)
source("lasso.r") 
source("cv.lasso.r")
source("tlpreg.r")
source("tlpreg1.r")
source("cv.tlpreg.r")

n <- 400
p <- 30
S <- matrix(0.5,p,p) 
diag(S) <- rep(1,p)
X <- rmvnorm(n, sigma=S)
y <- 1 + X[,1] + X[,2] + rnorm(n)
X <- apply(X, 2, function(z) (z - mean(z))/sd(z))

z <- y
A <- X[,1:2]

summary(lm(y ~ X[,1:2]))
summary(lm(z ~ A[,1:2]))

mm3 <- tlpreg1(y=y, X=X)
mm3$br[,20]
mm3$b[,2]

summary(lm(y ~ X[,1:2]))
summary(lm(z ~ A[,1:2]))


