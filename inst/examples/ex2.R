## ex.2.R -- this file is part of the trustOptim contributed package for R
##


## example of using trust.optim function when the Hessian is known analytically

rm(list=ls())
gc()

library(plyr)
library(Matrix)
library(mvtnorm)
##library(numDeriv)
library(trustOptim)

source("ex_funcs.R")

set.seed(123)

N <- 200
k <- 5
T <- 10

method <- "Sparse"
control.list <- list(start.trust.radius=5,
                     stop.trust.radius = 1e-5,
                     prec=1e-7,
                     report.freq=10L,
                     report.level=4L,
                     report.precision=2L,
                     maxit=5000L,
                     function.scale.factor = as.numeric(-1),                           
                     preconditioner=0L
                     ) 

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
Y <- laply(log.p, function(q) return(rbinom(1,T,exp(q))))

XX.list <- vector("list",length=N)
for (i in 1:N) {
  XX.list[[i]] <- tcrossprod(X[,i])
}


nvars <- N*k + k
start <- rnorm(nvars) ## random starting values
 
cat("running ",method, "\n")
t1 <- Sys.time()
opt <- trust.optim(start, fn=get.f,
                   gr = get.grad,
                   hs = get.hess,  ## used only for Sparse method
                   method = method,
                   control = control.list,
                   Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma, XX.list=XX.list
                   )
t2 <- Sys.time()
td <- difftime(t2,t1)
print(td,units="secs")

