# `tlpreg`
Functions for regressions with truncated lasso penalty

# Usage

Compilation:
```sh
# In terminal, type 
R CMD SHLIB lasso.cc rlasso.cc rtlpreg.cc tlpreg.cc rtlpreg1.cc tlpreg1.cc -o tlpreg
```

Example: 
```r
# Generate Data
n <- 1000
p <- 3000
X <- matrix(rnorm(n*p),n,p)
Z <- rnorm(n)
for(j in 1:p) {
  X[,j] <- Z + rnorm(n)
  X[,j] <- (X[,j]-mean(X[,j]))/sd(X[,j])
}
y <- 1 + 0.5*(X[,1] + X[,2] + X[,10] + X[,50]) + rnorm(n) 

# Estimation by TLP
source("tlpreg.r")
source("cv.tlpreg.r")
m0.cv <- cv.tlpreg0(X=X, y=y) # default is 10-fold CV
m0 <- tlpreg0(X=X, y=y, gamma=m0.cv$gamma.min)

# Estimation by LASSO
source("lasso.r")
source("cv.lasso.r")
m1.cv <- cv.lasso(X=X, y=y, nfold=5) # Run 5-fold CV
m1 <- lasso0(X=X, y=y, lambda=m1.cv$lambda)
```

# To-dos
- Strong set update
- Multiple regression
- Logit regression
- Precision matrix estimation
- Inference