\name{trust.optim}
\alias{trust.optim}
\title{Nonlinear optimizers using trust regions}
\description{
	Run nonlinear minimizer using trust region algorithm with conjugate gradient search directions and quasi-Hessian updates
	}
\usage{
trust.optim(x, fn, gr, hs=NULL, method=c("SR1","BFGS","Sparse","SparseFD"),
	        hess.struct=NULL, control=list(), \dots)
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
            numDeriv package).  Do not use a finite-differenced gradient with 
	    the SparseFD method.  That will cause a world of hurt.}
	 \item{hs}{an R function that takes x as its first argument.
	   Returns a Hessian matrix object of class "dgCMatrix" (see the Matrix package).  This function is called only if the selected method is "Sparse."}      
\item{method}{Valid arguments are "SR1","BFGS","Sparse" or "SparseFD"}
\item{hess.struct}{For  method SparseFD, a list with exactly two entries: integer vectors iRow and jCol. See details.  Will be ignored for methods SR1, BFGS and Sparse. 
    }
\item{control}{A list containing control parameters for the
    optimizer. See
    details.}
\item{\dots}{Additonal arguments passed to \code{fn}, \code{gr} and \code{hs}.  All arguments must be
  named.}
}


\value{
	\item{fval}{Value of the objective function}
	\item{solution}{Parameter vector at the optimum}
	\item{gradient}{Gradient at the optimum} 
	\item{hessian}{Estimate of the Hessian at the optimum (returned only for Sparse and SparseFD methods).}
	\item{iterations}{Number of iterations before stopping}
	\item{status}{A message describing the last state of the iterator}
}

\section{Control parameters}{
The control list should include the following parameters.
  \describe{
\item{start.trust.radius}{Initial radius of the trust region. Default is 5.  If the algorithm returns non-finite values of the objective function early in the process, try a lower number.}
\item{stop.trust.radius}{Minimum radius of trust region.  Algorithm will terminate if radius is below this value.  This is because it may not be possible to get the norm of the gradient smaller than prec, and this is another way to get the algorithm to stop.}
\item{cg.tol}{tolerance for the conjugate gradient algorithm that is used for the trust region subproblem.  Set it to something very small.  Default is sqrt(.Machine$double.eps)}
\item{prec}{Precision for how close the norm of the gradient at the
    solution should be to zero, before the algorithm halts.  It is possible that the algorithm
    will not get that far, so it will also stop when the radius of the
    trust region is smaller thanstop.trust.radius.  If the trust
    region radius collapses, but the norm of the gradient really isn't
    close to zero, then something terrible has happened.}
\item{report.freq}{An integer. The frequency at which the algorithm
  will display the current iteration number or function value, among
  other things (see \code{report.level}).  Defaults to 1.}
\item{report.level}{The amount of detail in each report.  Defaults to 2.}
\item{report.precision}{The number of significant digits used in each
    report. Defaults to 5.}
\item{maxit}{Maximum number of iterations.  Defaults to 100.}
\item{contract.factor}{When the algorithm decides to shrink the trust
    region, it will multiply the trust radius by this factor. Defaults to 0.5.}
\item{expand.factor}{When the algorithm decides to expand the trust
    region, it will multiply the algorithm by this factor. Defaults to 3.}
\item{contract.threshold}{The algorithm with accept a proposed move if the ratio of the actual improvement in the objective function, to the predicted improvement from the trust region subproblem, is greater than this amount.  Otherwise, the trust region will contract.  Default is 0.25.}
\item{expand.threshold.ap}{First criterion to determine if the trust region should expand.  If the ratio of the actual and proposed improvements in the objective function is less than this factor, the algorithm will consider expanding the trust region.  See \code{expand.threshold.radius}. DEfault is 0.8.}
\item{expand.threshold.radius}{If the ratio of the actual and proposed improvement in the objective function is less than \code{expand.threshold.ap}, then, if the normed distance of the proposed move is greater than  \code{expand.threshold.radius}, times the current trust region radius, the trust region will expand.  Default is 0.8.}
\item{function.scale.factor}{The algorithm will minimize \code{fn} times this
    factor. If you want to maximize \code{fn}, this value should be negative
    (usually -1).  Default is 1.}
\item{precond.refresh.freq}{Frequency at which the preconditioner
    for the conjugate gradiate estimation of the trust region
    subproblem is reestimated.  Preconditioners can help the convergence properties of the algoriothm.  Default is 1.}
\item{preconditioner}{ID for choice of preconditoner.  0 is the
    identity matrix (default), and 1 is the diagonal of the Hessian. For the Sparse and SparseFD method, 2 is a modified Cholesky preconditioner. For the BFGS method, 2 is the full Cholesky decomposition.  If you select 2 for the SR1 method, the algorithm will use the diagonal preconditioner instead.}
\item{fd.method}{This flag is 1 if using the direct method for computing a finite-differenced Hessian.  Otherwise, the indirect method is used. Default value is 0. See details.  This argument is ignored for methods SR1, BFGS and Sparse.}
\item{trust.iter}{Maximum number of conjugate gradient iterations to run when solving the trust region subproblem.  A higher number will lead to more accurate solutions to the subproblem, but may also lead to longer run times. Defaults to 2000.} 
}
}

\section{Hessian information}{
For the SparseFD method, you must provide information about the structure of the Hessian. There are two elements to the hess.struct list:
\describe{
	\item{iRow}{row indices of non-zero entries in the \emph{lower triangle} of the Hessian}
	\item{jCol}{column indices of non-zero entries in the \emph{lower triangle} of the Hessian}
}
Indexing starts at 1, and should include the diagonal elements, even though it is obvious that the diagonal elements are non-zero.  Entries must be ordered by column, and by row within columns.  Do not include any entries for the upper triangle.

The algorithms used for estimating sparse Hessians using finite differencing are described in 

Coleman, Thomas F, Burton S Garbow, and Jorge J More. 1985. Software
for Estimating Sparse Hessian Matrices. ACM Transaction on
Mathematical Software 11 (4) (December): 363-377.

This method involves a partitioning and permutation of the Hessian matrix to reduce the number of distinct finite differencing operations of the gradient.  There are two methods for computing the partition and permutation:  direct and indirect.  The direct method tends to require more computation, but may be more accurate.  We recommend the indirect method to start, so the default value is 0.


}
\section{Report levels}{
The \code{report.level} control parameter determines how much information is displayed each time the algorithm reports the current state.  Possible values are 

\describe{
\item{<=0}{No information (a quiet run)}
\item{1}{Current iteration number, and current value of the objective function.}
\item{2}{Information from level 1, plus the current norm of the gradient and a status message.}
\item{3}{Information from levels 1 and 2, plus the current normed radius of the trust region.}
\item{4}{Information from levels 1, 2, and 3, plus information from each estimate of the trust region subproblem (number of conjugate gradient iterations and how/why the CG algorithm terminated).}
}

Default level is 2.  Levels 3 and 4 are available primarily for debugging purposes.

}

\section{Stopping criteria}{
The algorithm will stop when one of the following three conditions are met:
\itemize{
\item{The norm of the gradient, divided by the square root of the number of parameters, is less than \code{prec}.}
\item{The trust region collapse to a radius smaller than machine precision}
\item{The algorithm proposes zero or negative improvement in the objective function (should never happen)}
\item{The number of iterations reaches the control parameter \code{maxit}}
}

If the algorithm appears to have stopped prematurely (i.e., the norm of the gradient is still too large), then one might just restart the algorithm.  For the quasi-Newton algorithms (SR1 and BFGS), this will refresh the Hessian, and might allow more progress to be made.
}




 


