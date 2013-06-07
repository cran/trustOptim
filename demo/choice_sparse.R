## Runs the hierarchical binary choice example

library(Matrix)
library(mvtnorm)
library(trustOptim)
library(sparseHessianFD)

set.seed(123)

N <- 250  ## number of heterogeneous units
k <- 8  ## number of covariates
T <- 20  ## number of "purchase opportunities per unit


## Simulate data and set priors

x.mean <- rep(0,k)
x.cov <- diag(k)
mu <- rnorm(k,0,10)
Omega <- diag(k)
inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]
inv.Omega <- solve(Omega)
X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
XB <- colSums(X * B)
log.p <- XB - log1p(exp(XB))
Y <- apply(as.matrix(log.p), 1,function(q) return(rbinom(1,T,exp(q))))

nvars <- N*k + k
start <- rnorm(nvars) ## random starting values
hess.struct <- demo.get.hess.struct(N, k)

## Setting up function to compute Hessian using sparseHessianFD package.
obj <- new.sparse.hessian.obj(start, fn=demo.get.f, gr=demo.get.grad,
                              hs=hess.struct, Y=Y, X=X,
                              inv.Omega=inv.Omega,
                              inv.Sigma=inv.Sigma, T=T)


cat("Using Sparse method in trust.optim\n")

td <- system.time(opt <- trust.optim(start, fn=obj$fn,
                                     gr = obj$gr,
                                     hs = obj$hessian,
                                     method = "Sparse",
                                     control = list(
                                       start.trust.radius=5,
                                       stop.trust.radius = 1e-7,
                                       prec=1e-7,
                                       report.freq=1L,
                                       report.level=4L,
                                       report.precision=1L,
                                       maxit=500L,
                                       preconditioner=1L
                                       ) 
                                     )
                  )
print(td)
