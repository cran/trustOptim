## Runs the binary choice with homogeneous coefficeints example

library(Matrix)
library(mvtnorm)
library(trustOptim)
library(sparseHessianFD)

set.seed(123)

N <- 250  ## number of heterogeneous units
k <- 100   ## number of covariates
T <- 20  ## number of "purchase opportunities per unit


## Simulate data and set priors

x.mean <- rep(0,k)
x.cov <- 10*diag(k)
beta <- seq(-3,3,length=k)
Omega <- 100*diag(k)
inv.Omega <- solve(Omega)
X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
XB <- t(X) %*% beta
log.p <- XB - log1p(exp(XB))
Y <- apply(as.matrix(log.p), 1,function(q) return(rbinom(1,T,exp(q))))

nvars <- k
start <- rnorm(nvars) ## random starting values

hess.struct <- Matrix.to.Coord(matrix(1,k,k))

## Setting up function to compute Hessian using sparseHessianFD package.
obj <- new.sparse.hessian.obj(start, fn=demo.get.f.dense, gr=demo.get.grad.dense,
                              hs=hess.struct, Y=Y, X=X,
                              inv.Omega=inv.Omega, T=T)


cat("Using Sparse method in trust.optim\n")
td <- system.time(opt <- trust.optim(start, fn=obj$fn,
                                     gr = obj$gr,
                                     hs = obj$hessian,
                                     method = "Sparse",                             
                                     control = list(
                                       start.trust.radius=5,
                                       stop.trust.radius = 1e-5,
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



