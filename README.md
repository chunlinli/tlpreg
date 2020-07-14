# `tlpreg`

Functions for regressions with truncated lasso penalty (TLP, also known as capped $\ell_1$). 
The implementation is partially based on Li et al. (2020+).

This project originated from my personal research need for solving high dimensional constrained regression.

*The code requires a modern C++ compiler supporting C++11 or later.* 

## Models

#### Gaussian regression

- Lasso regularized regression: `lasso0`
  
$$\frac{1}{2n} \sum_{i=1}^n(y_i -\beta_0 - \bm x_{i}^\top \bm\beta)^2 + \lambda \sum_{j=1}^{p} w_j |\beta_{j}|.$$

- Truncated Lasso penalized regression: `tlpreg0`

$$\frac{1}{2n} \sum_{i=1}^n(y_i -\beta_0 - \bm x_{i}^\top \bm\beta)^2 + \gamma \sum_{j=1}^{p} w_j \min(|\beta_{j}|,\tau).$$

- $\ell_0$-constrained regression: `tlpreg1`

$$\frac{1}{2n} \sum_{i=1}^n(y_i -\beta_0 - \bm x_{i}^\top \bm\beta)^2 \quad \text{s.t.} \quad \|\bm\beta\|_0 \leq K.$$

*By default, the penalty weight $w_j = 1$.*

## Usage

Compilation: (should replace this by a makefile)
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

# Estimation by TLP regularized regression
source("tlpreg.r")
source("cv.tlpreg.r")
m0.cv <- cv.tlpreg0(X=X, y=y) # Run 10-fold CV by default
m0 <- tlpreg0(X=X, y=y, gamma=m0.cv$gamma.min)

# Estimation by L0 constrained regression
source("tlpreg1.r")
source("cv.tlpreg1.r")
m1.cv <- cv.tlpreg1(X=X, y=y) # Run 10-fold CV by default
m1 <- tlpreg1(X=X, y=y, K=m1.cv$K.min)

# Estimation by LASSO
source("lasso.r")
source("cv.lasso.r")
m2.cv <- cv.lasso(X=X, y=y, nfold=5) # Run 5-fold CV
m2 <- lasso0(X=X, y=y, lambda=m2.cv$lambda)
```


## To-dos

- Inference
- Strong set update
- Elastic-net
- Accelerate multi-response regression
- Logistic/Poisson/Cox regressions
