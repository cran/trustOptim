rm(list=ls())

library(plyr)
library(Matrix)
library(mvtnorm)
library(numDeriv)
library(trustOptim)
library(foreach)
library(doParallel)
library(trust)

source("testfuncs.R")

set.seed(123)
registerDoParallel(cores=10)

save.file <- "~/Documents/R_packages/docs/trustOptim/test1.Rdata"

run.optim.test <- function(k, N) {

  out <- array(dim=c(n.meth,n.stats))
  dimnames(out) <- list(method=meth.vec, stat=stats)
  
  x.mean <- rep(0,k)
  x.cov <- 0.1*diag(k)

  mu <- rnorm(k,0,10)
  Omega <- diag(k)
  inv.Sigma <- rWishart(1,k+5,diag(k))[,,1]
  inv.Omega <- solve(Omega)

  X <- t(rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
  B <- t(rmvnorm(N, mean=mu, sigma=Omega)) ## k x N

  XB <- colSums(X * B)
  log.p <- XB - log1p(exp(XB))
  Y <- laply(log.p, function(q) return(rbinom(1,T,exp(q))))
  
  nvars <- N*k + k

  start <- rnorm(nvars)
    
  hess.struct <- get.hess.struct(N, k)

##  s1 <- system.time(HH <- get.hess(start, Y, X, inv.Omega, inv.Sigma))
##  print(s1)
##  s2 <- system.time(get.hess2(start, Y, X, inv.Omega, inv.Sigma))
##  print(s2)
##  return(HH)
  
  do.next.method <- TRUE
  if (do.next.method) {
 
    cat("Running CG\n")
    t1 <- Sys.time()
    opt <- optim(start, fn=log.f, gr= get.grad,
                 method="CG", hessian=FALSE,
                 control = list(fnscale= -1, maxit=100000, REPORT=0, trace=0),
                 Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                 )
    t2 <- Sys.time()
    tm <- difftime(t2, t1, units="secs")
    sol <- opt$par
    gr <- get.grad(sol, Y, X, inv.Omega, inv.Sigma)
    n.iter <- opt$counts[2]
    nrm.gr <- sqrt(sum(gr^2))
    max.abs.gr <- max(abs(gr))
    conv <- opt$convergence
    if (opt$convergence==0) out["CG",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  }

  do.next.method <- TRUE
  if (do.next.method) {
    cat("running BFGS.optim")
    t1 <- Sys.time()
    opt3 <- optim(start, fn=log.f, gr= get.grad,
                  method="BFGS", hessian=FALSE,
                  control = list(fnscale= -1, maxit=100000, REPORT=100, trace=1),
                  Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                  )
    t2 <- Sys.time()
    tm <- difftime(t2, t1, units="secs")
    sol <- opt3$par
    gr <- get.grad(sol, Y, X, inv.Omega, inv.Sigma)
    n.iter <- opt3$counts[2]
    nrm.gr <- sqrt(sum(gr^2))
    max.abs.gr <- max(abs(gr))
    conv <- opt3$convergence
    if (opt3$convergence==0) out["BFGS.optim",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  }

  do.next.method <- TRUE
  if (do.next.method) {
 
    cat("running Sparse\n")
    t1 <- Sys.time()
    opt <- trust.optim(start, fn=log.f,
                       gr = get.grad,
                       hs = get.hess,
                       hess.struct = hess.struct,
                       method = "Sparse",
                       control = control.list,
                       Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                       )
    t2 <- Sys.time()
    tm <- difftime(t2, t1, units="secs")
    sol <- opt$solution
    gr <- opt$gradient
    n.iter <- opt$iterations
    nrm.gr <- sqrt(sum(gr^2))
    max.abs.gr <- max(abs(gr))
    conv <- opt$convergence
    if (opt$status == "Success") out["Sparse",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
   }

  do.next.method <- TRUE
  if (do.next.method) {
    
    cat("running SparseFD\n")
    t1 <- Sys.time()
    opt <- trust.optim(start, fn=log.f,
                       gr = get.grad,
                       hs = get.hess,
                       hess.struct = hess.struct,
                       method = "SparseFD",
                     control = control.list,
                     Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                     )
    t2 <- Sys.time()
    tm <- difftime(t2, t1, units="secs")
    sol <- opt$solution
    gr <- opt$gradient
    n.iter <- opt$iterations
    nrm.gr <- sqrt(sum(gr^2))
    max.abs.gr <- max(abs(gr))
    conv <- opt$convergence
    if (opt$status == "Success") out["SparseFD",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  }
  
  
  do.next.method <- FALSE
  if (do.next.method) {
    cat("running BFGS\n")
    t1 <- Sys.time()
    opt <- trust.optim(start, fn=log.f,
                       gr = get.grad,
                       hs = get.hess,
                       hess.struct = hess.struct,
                       method = "BFGS",
                       control = control.list,
                       Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                       )
    t2 <- Sys.time()
    tm <- difftime(t2, t1, units="secs")
    sol <- opt$solution
    gr <- opt$gradient
    n.iter <- opt$iterations
    nrm.gr <- sqrt(sum(gr^2))
    max.abs.gr <- max(abs(gr))
    conv <- opt$convergence
    if (opt$status == "Success") out["BFGS.trust",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  }


  do.next.method <- FALSE
  if (do.next.method) {
    cat("running SR1\n")
    t1 <- Sys.time()
    opt2 <- trust.optim(start, fn=log.f,
                       gr = get.grad,
                       hs = get.hess,
                       hess.struct = hess.struct,
                       method = "SR1",
                       control = control.list,
                       Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                       )
    t2 <- Sys.time()
    tm <- difftime(t2, t1, units="secs")
    sol <- opt2$solution
    gr <- opt2$gradient
    n.iter <- opt2$iterations
    nrm.gr <- sqrt(sum(gr^2))
    max.abs.gr <- max(abs(gr))
    conv <- opt2$convergence
    if (opt2$status == "Success") out["SR1.trust",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  }


  do.next.method <- TRUE
  if (do.next.method) {
    cat("running trust\n")    
    t1 <- Sys.time()
    opt2 <- trust(objfun=trust.func,
                  parinit=start,
                  rinit=control.list$start.trust.radius,
                  rmax=10000,
                  parscale=rep(1,length(start)),
                  iterlim=50000,
                  minimize=FALSE,
                  Y=Y, X=X, inv.Omega=inv.Omega, inv.Sigma=inv.Sigma
                  )
    t2 <- Sys.time()
    tm <- difftime(t2, t1, units="secs")
    sol <- opt2$argument
    gr <- opt2$gradient
    n.iter <- opt2$iterations
    nrm.gr <- sqrt(sum(gr^2))
    max.abs.gr <- max(abs(gr))
    conv <- opt2$converged
    if (opt2$converged == TRUE) out["trust",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  }
  

  return(out)
  
}
  
T <- 25
N.vec <- c(100, 1000)
k.vec <- c(2,5,10)

meth.vec <- c("Sparse","SparseFD","BFGS.optim","trust","CG")
n.meth <- length(meth.vec)

control.list <- list(start.trust.radius=5,
                     stop.trust.radius = 1e-5,
                     prec=1e-5,
                     cg.tol=1e-5,
                     report.freq=100L,
                     report.level=6L,
                     report.precision=5L,
                     maxit=2000000L,
                     function.scale.factor = as.numeric(-1),                           
                     preconditioner=0L,
                     trust.iter=5000L
                     )


stats <- c("time","iters","nrm.gr","max.abs.gr")
n.stats <- length(stats)
n.trials <- 1

res <- array(dim=c(length(N.vec), length(k.vec), n.meth, n.stats,n.trials))
dimnames(res) <- list(N=N.vec,k=k.vec, method=meth.vec, stat=stats, trial=1:n.trials)

for (ik in 1:length(k.vec)) {
  for (ii in 1:length(N.vec)) {

    k <- k.vec[ik]
    N <- N.vec[ii]
    
    cat("\tk = ",k,"  N = ",N,"\n")  

##    HH <- run.optim.test(k,N)
  
    
    zz <- times(n.trials) %dopar% run.optim.test(k, N)
 
    res[ii, ik,,,] <- array(zz,dim=c(n.meth,n.stats, n.trials))
    
    save(res, file=save.file)
  }
}

save(res, file=save.file)
