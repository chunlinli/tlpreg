library(glmnet)
library(ncvreg)
source("lasso.r")
source("cv.lasso.r")

# SPEED TEST: CROSS-VALIDATION

n <- 500
p <- 1000
sparse <- as.integer(0.005 * p)
idx <- sample(1:p, sparse, replace = FALSE)
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- Z + rnorm(n)
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
bt <- rnorm(length(idx))
y <- X[,idx] %*% bt + rnorm(n)

t0 <- Sys.time()
mm1 <- cv.glmnet(x=X, y=y, standardize=FALSE)
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
mm2 <- cv.ncvreg(X=X, y=y, penalty="lasso", lambda = 2*mm1$lambda)
t1 <- Sys.time()
t1-t0
t0 <- Sys.time()
mm3 <- cv.lasso(X=X, y=y, nfold =10)
t1 <- Sys.time()
t1-t0

mm1$lambda.min
mm2$lambda.min
mm3$lambda.min



# SPEED TEST: SINGLE PATH

t0 <- Sys.time()
m1 <- glmnet(x=X, y=y, standardize=FALSE)
t1 <- Sys.time()
t1-t0

t0 <- Sys.time()
m2 <- ncvreg(X=X, y=y)
t1 <- Sys.time()
t1-t0

t0 <- Sys.time()
m3 <- lasso0(X=X, y=y)
t1 <- Sys.time()
t1-t0



# SPEED TEST: SINGLE VALUE

n <- 500
p <- 10000
sparse <- as.integer(0.005 * p)
idx <- sample(1:p, sparse, replace = FALSE)
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- Z + rnorm(n)
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
bt <- rnorm(length(idx))
y <- X[,idx] %*% bt + rnorm(n)

lambda = .07

t0 <- Sys.time()
m1 <- glmnet(x=X, y=y, lambda=lambda, standardize=FALSE)
t1 <- Sys.time()
t1-t0

t0 <- Sys.time()
m2 <- ncvreg(X=X, y=y, lambda=2*lambda)
t1 <- Sys.time()
t1-t0

t0 <- Sys.time()
m3 <- lasso0(X=X, y=y, lambda=lambda)
t1 <- Sys.time()
t1-t0

b1 <- as.numeric(m1$beta)
a1 <- as.numeric(m1$a)
b2 <- as.numeric(m2$beta[-1])
a2 <- as.numeric(m2$beta[1])
b3 <- as.numeric(m3$b)
a3 <- as.numeric(m3$b0)

mean((y - X%*%b1 - a1)^2/2) + lambda*sum(abs(b1))
mean((y - X%*%b2 - a2)^2/2) + lambda*sum(abs(b2))
mean((y - X%*%b3 - a3)^2/2) + lambda*sum(abs(b3))

(mean((y - X%*%b3 - a3)^2/2) + lambda*sum(abs(b3))) - (mean((y - X%*%b1 - a1)^2/2) + lambda*sum(abs(b1)))
(mean((y - X%*%b3 - a3)^2/2) + lambda*sum(abs(b3))) - (mean((y - X%*%b2 - a2)^2/2) + lambda*sum(abs(b2)))



# ACCURACY TEST

n <- 400
p <- 500
sparse <- as.integer(0.005 * p)
idx <- sample(1:p, sparse, replace = FALSE)
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- Z + rnorm(n)
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
bt <- rnorm(length(idx))
y <- X[,idx] %*% bt + rnorm(n)


lambda = 0.07
m1 <- glmnet(x=X, y=y, lambda=lambda, standardize=FALSE)
m2 <- ncvreg(X=X, y=y, lambda=2*lambda)
m3 <- lasso0(X=X, y=y, lambda=lambda)

b1 <- as.numeric(m1$beta)
a1 <- as.numeric(m1$a)
b2 <- as.numeric(m2$beta[-1])
a2 <- as.numeric(m2$beta[1])
b3 <- as.numeric(m3$b)
a3 <- as.numeric(m3$b0)

mean((y - X%*%b1 - a1)^2/2) + lambda*sum(abs(b1))
mean((y - X%*%b2 - a2)^2/2) + lambda*sum(abs(b2))
mean((y - X%*%b3 - a3)^2/2) + lambda*sum(abs(b3))

(mean((y - X%*%b3 - a3)^2/2) + lambda*sum(abs(b3))) - (mean((y - X%*%b1 - a1)^2/2) + lambda*sum(abs(b1)))
(mean((y - X%*%b3 - a3)^2/2) + lambda*sum(abs(b3))) - (mean((y - X%*%b2 - a2)^2/2) + lambda*sum(abs(b2)))
