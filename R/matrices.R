## matrices.R --   Part of the trustOptim package for the R programming language.
## This file is part of trustOptim, a nonlinear optimization package
## for the R statistical programming platform.
##
## Copyright (C) 2012 Michael Braun
##
## This Source Code Form is subject to the license terms of the Mozilla
## Public License v. 2.0. If a copy of the MPL was not distributed
## with this file, you can obtain one at http://mozilla.org/MPL/2.0/.
## See the trustOptim LICENSE file for more information.


Sym.CSC.to.Matrix <- function(H,nvars) {

## H is a list of data stored in compressed sparse column (CSC) format
## returns a sparse Matrix object

  require(Matrix)
  M <- new("dsCMatrix", i = H$indrow, p = H$jpntr, x = H$vals, Dim=c(nvars, nvars),uplo="L")
  return(M)

}


Matrix.to.Coord <- function(M) {

  res <- vector("list",length=2)
  names(res) <- c("iRow","jCol")
  M <- tril(as(M,"TsparseMatrix"))
  res$iRow <- as.integer(M@i)+1 ## return to 1-based indexing
  res$jCol <- as.integer(M@j)+1
  return(res)

}

Coord.to.Pattern.Matrix <- function(H,nrows, ncols=nrows) {

  ## H is a list with two integer vectors:  iRow and jCol
  
  res <- sparseMatrix(i=H$iRow,j=H$jCol, dims=c(as.integer(nrows), as.integer(ncols)))
  return(res)

}

Coord.to.Sym.Pattern.Matrix <- function(H,nvars) {

## coords are for lower triangle, but coerces to symmetric pattern matrix
## H is a list with two integer vectors:  iRow and jCol

  
  res <- new("nsTMatrix",i=as.integer(H$iRow-1), j=as.integer(H$jCol-1),
             Dim=c(as.integer(nvars), as.integer(nvars)),uplo="L")
  return(res)

}
