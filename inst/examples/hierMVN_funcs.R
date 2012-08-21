


get.log.wish <- function(vech.chol.G, nu, chol.inv.mean.G) {

  k <- dim(chol.inv.mean.G)[1]
  chol.G <- xpnd(vech.chol.G,k)*lower.tri(diag(k),diag=TRUE)
  
  C <- nu*sum(log(diag(chol.inv.mean.G))) + (k*nu*log(nu)/2) - nu*k*log(2)/2 - k*(k-1)*log(pi)/4 - sum(lgamma((nu+1-c(1:k))/2))
  
  
  log.det <- sum(log(diag(chol.G)))
  Z <- t(chol.inv.mean.G) %*% chol.G
  tr <- sum(diag(crossprod(Z)))
  res <- C+(nu-k-1)*log.det - nu*tr/2
  
  return(res)
  
}



get.f <- function(par, data.list){

  res <- .Call("HierMVN_get_f",par,data.list)

}



get.df <- function(par, data.list){

  res <- .Call("HierMVN_get_fdf",par,data.list)$grad

}


get.log.post.R <- function(P, Y, X, mu.prior.mean, mu.prior.chol.prec, nu, chol.inv.mean.G) {


  B <- matrix(P[1:(N*k)],k,N)
  mu <- P[(N*k+1):(N*k+k)]
  vech.chol.G <- P[(N*k+k+1):length(P)]
  chol.G <- xpnd(vech.chol.G,k)
  chol.G <- chol.G*lower.tri(chol.G,diag=TRUE)
  G <-  tcrossprod(chol.G)
  S <- solve(G)

  y.mean <- matrix(colSums(X*B),1,N)

  data.LL <- 0
  for (i in 1:N) {
    data.LL <- data.LL + sum(dnorm(Y[,i],y.mean[i],1,log=TRUE))
  }

  prior.B <- sum(dmvnorm(t(B), mu, S, log=TRUE))
  prior.mu <- dmvnorm(mu, mu.prior.mean, mu.prior.chol.prec, log=TRUE)
  prior.G <- get.log.wish(vech.chol.G, nu, chol.inv.mean.G)

  f <- data.LL + prior.B + prior.mu + prior.G
  

}



get.hess.struct <- function(nvars, N, k){

  ## construct hessian structure -- block-arrow structure
  
  ## rang
  
  b.range <- 1:(k*N)
  mu.range <- (k*N+1):(k*(N+1))
  G.range <- (k*(N+1)+1):(k*(N+1)+k*(k+1)/2)
 
  
  HS <- Matrix(0,nrow=nvars, ncol=nvars)
  HS <- as(HS,"nMatrix")
  
  HS[b.range,b.range] <- as(kronecker(Diagonal(N,TRUE),Matrix(TRUE,k,k)),"lsparseMatrix") ## check conditional dependence
  HS[mu.range,b.range] <- TRUE
  HS[G.range,b.range] <- TRUE
  HS[mu.range,mu.range] <- TRUE
  HS[G.range,mu.range] <- TRUE
  HS[G.range,G.range] <- TRUE

  HS <- tril(HS)
  HS <- as(HS,"TsparseMatrix")
  iRow <- as.integer(HS@i+1)
  jCol <- as.integer(HS@j+1)

  return(list(iRow=iRow, jCol=jCol))
}


