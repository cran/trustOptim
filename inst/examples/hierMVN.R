rm(list=ls())
library(mvtnorm)
library(plyr)
library(trustOptim)

source("hierMVN_funcs.R")

set.seed(12345)

N <- 200
T <- 10
k <- 6 # number of covariates, including intercept
s1 <- 5

B.mean <- seq(-5,5,length=k)
X <- rbind(1,matrix(rnorm(N*(k-1),0,1),k-1,N))
B <- laply(B.mean,function(mu) return(rnorm(N,mu,s1)))
y.mean <- colSums(B*X)
Y <- t(laply(y.mean, function(mu) return(rnorm(T,mu,1))))

data.list <- list(X=X, Y=Y)


nvars <- N*k + k + (k*(k+1)/2)

mu.prior.mean <- rep(0,k)
mu.prior.chol.prec <- t(chol(diag(k)))
nu <- as.numeric(k+6)
G.mean <- diag(k)
chol.inv.mean.G <- t(chol(solve(G.mean)))

prior.list <- list(mu.prior.mean=mu.prior.mean,
                   mu.prior.chol.prec=mu.prior.chol.prec,
                   nu=nu,
                   chol.inv.mean.G=chol.inv.mean.G
                   )

method <- "SR1"

if (method=="SparseFD") {
  hs.list <- get.hess.struct(nvars, N, k)  ## indexing starts at 1
  iRow <- as.integer(hs.list$iRow) 
  jCol <- as.integer(hs.list$jCol)
  hs.chk <- new("ntTMatrix",i=as.integer(iRow-1),j=as.integer(jCol-1),
                Dim=as.integer(c(nvars, nvars)),uplo="L")
} else {
  hs.list <- NULL
}

params.list <- list(data=data.list, priors=prior.list)

## find posterior mode

mu.start <- rep(0,k)
B.start <- rep(0,N*k)
G.start <- vech(chol.inv.mean.G)

startX <- c(B.start, mu.start, G.start)

control.list <- list(start.trust.radius=5.0,
                     stop.trust.radius=sqrt(.Machine$double.eps),
                     cg.tol=sqrt(.Machine$double.eps),
                     prec=sqrt(.Machine$double.eps)*10,
                     report.freq=10L,
                     report.level=6L,
                     report.precision=5L,
                     maxit=5000L,
                     contract.factor = 0.5,
                     expand.factor = 3,
                     contract.threshold = 0.25,
                     expand.threshold.ap = 0.8,
                     expand.threshold.radius = 0.8,
                     function.scale.factor = as.numeric(-1),                
                     precond.refresh.freq=1L,
                     preconditioner=1L,
                     fd.method=0L,
                     trust.iter=2000L
                     )


opt <- trust.optim(startX, fn=get.f, gr=get.df,
                   method=method,
                   hess.struct=hs.list,
                   control=control.list,
                   data.list=params.list)

post.mode <- opt$solution
gr <- opt$gradient
