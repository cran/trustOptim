rm(list=ls())
library(mvtnorm)
library(plyr)
library(trustOptim)
library(foreach)
library(doParallel)

registerDoParallel(cores=10)

source("hierMVN_funcs.R")

set.seed(12345)


run.optim.test <- function(k, N) {

  out <- array(dim=c(n.meth,n.stats))
  dimnames(out) <- list(method=meth.vec, stat=stats)
  
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
  
  params.list <- list(data=data.list, priors=prior.list)

  hs.list <- get.hess.struct(nvars, N, k)  ## indexing starts at 1
   
  mu.start <- rmvnorm(1,mu.prior.mean,10*diag(k))
  B.start <- as.vector(t(rmvnorm(N,mu.start,10*diag(k))))
  G.start <- vech(chol.inv.mean.G)
  
  startX <- c(B.start, mu.start, G.start)   
  

  t1 <- Sys.time()
  opt <- trust.optim(startX, fn=get.f, gr=get.df,
                     method="SparseFD",
                     hess.struct=hs.list,
                     control=control.list,
                     data.list=params.list)
  t2 <- Sys.time()
  post.mode <- opt$solution
  gr <- opt$gradient
  cat("convergence: ",opt$status,"\n")
  n.iter <- opt$iterations
  nrm.gr <- sqrt(sum(gr^2))
  tm <- difftime(t2,t1)
  max.abs.gr <- max(abs(gr))
  
  if (opt$status == "Success") out["SparseFD",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
    
  
  t1 <- Sys.time()
  opt <- optim(startX,  fn=get.f,
               gr=get.df,
               method="BFGS",
               control=list(fnscale=-1,trace=1,REPORT=50,maxit=100000),
               hessian=FALSE,
               data.list=params.list)
  t2 <- Sys.time()
  post.mode <- opt$par
  gr <- get.df(post.mode, data.list=params.list)
  cat("convergence: ",opt$convergence,"\n")
  n.iter <- opt$counts[2]
  nrm.gr <- sqrt(sum(gr^2))
  tm <- difftime(t2,t1)
  max.abs.gr <- max(abs(gr))
  if (opt$convergence==0) out["BFGS.optim",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  
 
  t1 <- Sys.time()
  opt <- optim(startX,  fn=get.f,
               gr=get.df,
               method="CG",
               control=list(fnscale=-1,trace=0,REPORT=0,maxit=100000),
               hessian=FALSE,
               data.list=params.list)
  t2 <- Sys.time()
  post.mode <- opt$par
  gr <- get.df(post.mode, data.list=params.list)
  cat("convergence: ",opt$convergence,"\n")
  n.iter <- opt$counts[2]
  nrm.gr <- sqrt(sum(gr^2))
  tm <- difftime(t2,t1)
  max.abs.gr <- max(abs(gr))
  if (opt$convergence==0) out["CG",] <- c(tm, n.iter, nrm.gr, max.abs.gr)
  
  return(out)
}


N.vec <- c(50, 200, 500)
##N.vec <- c(10,20,30)

T <- 20
k.vec <- c(2, 5, 20) # number of covariates, including intercept
##k.vec <- c(2,3,4)
s1 <- 5
meth.vec <- c("SparseFD","BFGS.optim","CG")
n.meth <- length(meth.vec)

control.list <- list(report.freq=50L,
                     report.level=6L,
                     report.precision=5L,
                     maxit=1000000L,
                     function.scale.factor = as.numeric(-1),                           
                     preconditioner=1L,
                     trust.iter=5000L
                     )


stats <- c("time","iters","nrm.gr","max.abs.gr")
n.stats <- length(stats)
n.trials <- 3

res <- array(dim=c(length(N.vec), length(k.vec), n.meth, n.stats,n.trials))
dimnames(res) <- list(N=N.vec,k=k.vec, method=meth.vec, stat=stats, trial=1:n.trials)

for (ik in 1:length(k.vec)) {
  for (ii in 1:length(N.vec)) {

    k <- k.vec[ik]
    N <- N.vec[ii]
    
    cat("\tk = ",k,"  N = ",N,"\n")  

    zz <- times(n.trials) %dopar% run.optim.test(k, N)
    res[ik, ii,,,] <- array(zz,dim=c(n.meth,n.stats, n.trials))
    
    
  }
}

save(res, file="~/Documents/R_packages/docs/trustOptim/test1.Rdata")
                            
        
