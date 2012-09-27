\name{get.fdfh}
\alias{get.fdfh}
\title{Return function, gradient, and finite-differenced sparse Hessian.}
\description{
	Returns function, gradient, and finite-differenced sparse Hessian through the trust.optim routine.
	}
\usage{
get.fdfh(x, fn, gr, hess.struct, fd.method=0L, fd.eps=sqrt(.Machine$double.eps), ...)
}

\arguments{
	\item{x}{a numeric vector of starting values for the optimizer.}
	\item{fn}{an R function that takes \code{x} as its first argument.
            Returns the value of the objective function at \code{x}.  Note
            that the optimizer will \emph{minimize} \code{fn} (see
            function.scale.factor under control)}
        \item{gr}{an R function that takes x as its first argument.  Returns a
            numeric vector that is the gradient of \code{fn} at \code{x}.  Naturally, the
            length of the gradient must be the same as the length of \code{x}.  The
            user must supply this function.  If an analytic gradient is not
            available, and the method is SR1 or BFGS, the user should consider a
            numerical approximation using finite differencing (see the
            numDeriv package).} 
	\item{hess.struct}{A list with exactly two entries: integer vectors iRow and jCol. For details, see documentation for trust.optim function.}
	\item{fd.method}{Direct (1) or indirect method (0) for computing the sparse matrix.  Default is 0, which is faster, but might be less accurate.}
	\item{fd.eps}{Perturbation amount for finite differencing.  Defaults to sqrt(.Machine$double.eps).}
	\item{\dots}{Additonal arguments passed to \code{fn} and \code{gr}.  All arguments must be named.}
}


\value{
	\item{fval}{Value of the objective function}
	\item{gradient}{Gradient at the optimum} 
	\item{hessian}{Estimate of the Hessian at the optimum (returned only for Sparse and SparseFD methods).}
}


\section{Details}{

When the exact Hessian is not available, but the structure of the Hessian is known, this function will return a finite-differenced estimate of the Hessian (along with the objective function value and gradient).  Computation is performed by the same routines that are used by the SparseFD method in the trust.optim function. This function is useful for when you are using SparseFD to optimize (which requires the Hessian structure anyway, and returns the sparse Hessian at the optimum), but might need to compute the Hessian for a different x.


See documentation for the trust.optim function, and specifically the SparseFD matrix, for details about how to provide Hessian structure, and the meanings of the fd.method and fd.eps parameters.
}



