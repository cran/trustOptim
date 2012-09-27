## callTrust.R --   Part of the trustOptim package for the R programming language.
## This file is part of trustOptim, a nonlinear optimization package
## for the R statistical programming platform.
##
## Copyright (C) 2012 Michael Braun
##
## This Source Code Form is subject to the license terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, you can obtain one at http://mozilla.org/MPL/2.0/.
## See the trustOptim LICENSE file for more information.

get.fdfh <- function(x, fn, gr, hess.struct, fd.method=0L, fd.eps=sqrt(.Machine$double.eps), ...){

  if (!is.function(fn)) stop ("Error in get.fdfh:  fn must be a function")
  if (!is.function(gr)) stop ("Error in get.fdfh:  gr must be a function")
  if (is.null(hess.struct)) stop("Error in trust.optim:  for SparseFD method, you must supply structure of the Hessian.")
  if (!is.list(hess.struct)) stop ("Error in get.fdfh:  hs struct must be a list")
  if (names(hess.struct)!=c("iRow","jCol")) stop ("Error in get.fdfh.  Names of hess.struct must be iRow and jCol")

  if (!is.integer(fd.method) || fd.method<0 || fd.method>1 || !is.finite(fd.method)) {
    stop("Error in trust.optim:  fd.method must be an integer, and either 0 or 1.")
  }
  
  nars <- length(x)
  
  fn1 <- function(x) fn(x,...)  ## currying the data in
  gr1 <- if (!is.null(gr)) function(x) gr(x,...)

  ## test fn and gr

  r1 <- fn1(x)
  if (!is.finite(r1)) stop("Error in get.fdfh:  fn at starting values is not finite.")
  r2 <- gr1(x)
  if (any(!is.finite(r2))) stop("Error in get.fdfh:  at least one element of gr is not finite.")
  if (length(r2)!=length(x)) stop("Error in get.fdfh:  length of gradient does not match length of starting values.")
    
  res <- .Call("get_fdfh", x, fn1, gr1, hess.struct, fd.method, fd.eps)
  res$hessian <- Matrix:::t(as(res$hessian,"symmetricMatrix")) 
  

  if (res$fval != r1) stop("Error in get.fdfh:  routine returned function value that does not equal fn(x).  Contact package maintainer.")
  if (!all.equal(res$gradient, r2)) stop("Error in get.fdfh:  routine returned gradient that does not equal gr(x).  Contact package maintainer.")

  return(res)
}

