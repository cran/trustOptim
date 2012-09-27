// Rinterface.cpp.  Part of the trustOptim package for the R programming language.
// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2012 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.
// See the trustOptim LICENSE file for more information.


#ifndef __TRUST_OPTIM_RINTERFACE
#define __TRUST_OPTIM_RINTERFACE

#define IS_R

#include <RcppEigen.h>
#include <common_R.hpp>
#include <include_eigen.hpp>

#include <CG-sparse.h>
#include <CG-quasi.h>
#include <Rfunc.cpp>
#include <RfuncHess.cpp>


RcppExport SEXP sparseTR(SEXP start_, SEXP fn_, SEXP gr_, SEXP hs_,
			 SEXP control_) {
  
  BEGIN_R_INTERFACE
    
    using Rcpp::NumericVector;
  using Rcpp::IntegerVector;
  using Rcpp::Function;
  using Rcpp::List;
  using Rcpp::as;
  
  using Eigen::VectorXi;
  using Eigen::Map;

  typedef double Scalar;

  typedef SparseMatrix<Scalar> optHessType;
  typedef SimplicialLLT<optHessType> optPrecondType; 

  NumericVector start(start_);
  int nvars = start.size();
  if (nvars<=0) throw MyException("Number of variables (starting values) must be positive\n",__FILE__,__LINE__);

  Function fn(fn_);
  Function gr(gr_);

 // Control parameters for optimizer

  List control(control_);
  const double rad = as<double>(control["start.trust.radius"]);
  const double min_rad = as<double>(control["stop.trust.radius"]);
  const double tol = as<double>(control["cg.tol"]);
  const double prec = as<double>(control["prec"]);
  const int report_freq = as<int>(control["report.freq"]);
  const int report_level = as<int>(control["report.level"]);
  const int report_precision = as<int>(control["report.precision"]);
  const int maxit = as<int>(control["maxit"]);
  const double contract_factor = as<double>(control["contract.factor"]);
  const double expand_factor = as<double>(control["expand.factor"]);
  const double contract_threshold = as<double>(control["contract.threshold"]);
  const double expand_threshold_rad = as<double>(control["expand.threshold.radius"]);
  const double expand_threshold_ap = as<double>(control["expand.threshold.ap"]);
  const double function_scale_factor = as<double>(control["function.scale.factor"]);
  const int precond_refresh_freq = as<int>(control["precond.refresh.freq"]);
  const int precond_ID = as<int>(control["preconditioner"]);
  int fd_method = as<int>(control["fd.method"]);
  const double fd_eps = as<double>(control["fd.eps"]);
  const int trust_iter = as<int>(control["trust.iter"]);

  // get Hessian structure

  List hs(hs_);

  IntegerVector iRow_((SEXP)hs["iRow"]);
  IntegerVector jCol_((SEXP)hs["jCol"]);
  int nnz = iRow_.size();

  if (nnz<=0) throw MyException("number of nonzeros in Hessian must be positive\n",__FILE__,__LINE__);
  if (jCol_.size()!=nnz) throw MyException("index vectors for Hessian must be of same length\n",__FILE__,__LINE__);


  Map<VectorXi> iRow(iRow_.begin(),nnz);
  Map<VectorXi> jCol(jCol_.begin(),nnz);

  Rfunc func(nvars, fn, gr);
  func.hessian_init(iRow, jCol, fd_method, fd_eps);

  Map<VectorXd> startX(start.begin(),nvars); 

  Trust_CG_Sparse<Map<VectorXd>, Rfunc, optHessType, optPrecondType> opt(func,
									 startX, rad, min_rad, tol, prec,
									 report_freq,
									 report_level, report_precision,
									 maxit, contract_factor, expand_factor,
									 contract_threshold,
									 expand_threshold_rad,
									 expand_threshold_ap,
									 function_scale_factor,
									 precond_refresh_freq,
									 precond_ID,
									 trust_iter);
  
  opt.run();

 // collect results and return

  VectorXd P(nvars);
  VectorXd grad(nvars);
  SparseMatrix<double> hess(nvars, nvars);
  hess.reserve(nnz);

  double fval, radius;
  int iterations;
  MB_Status status;

  status = opt.get_current_state(P, fval, grad, hess,
				 iterations, radius);

  List res;
  res = List::create(Rcpp::Named("fval") = Rcpp::wrap(fval),
		     Rcpp::Named("solution") = Rcpp::wrap(P),
		     Rcpp::Named("gradient") = Rcpp::wrap(grad),	
		     Rcpp::Named("hessian") = Rcpp::wrap(hess),
		     Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		     Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		     Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		     Rcpp::Named("nnz") = Rcpp::wrap(nnz),
		     Rcpp::Named("method") = Rcpp::wrap("SparseFD")
		     );
   
  return(res);

  END_R_INTERFACE
  
}

RcppExport SEXP quasiTR(SEXP start_, SEXP fn_, SEXP gr_,
			SEXP control_) {
  
    
  BEGIN_R_INTERFACE

  using Rcpp::NumericVector;
  using Rcpp::IntegerVector;
  using Rcpp::Function;
  using Rcpp::List;
  using Rcpp::as;
  
  using Eigen::VectorXi;
  using Eigen::Map;

  typedef MatrixXd optHessType;
  typedef LLT<optHessType> optPrecondType;

  NumericVector start(start_);
  int nvars = start.size();

  List control(control_);
  double rad = as<double>(control["start.trust.radius"]);
  const double min_rad = as<double>(control["stop.trust.radius"]);
  const double tol = as<double>(control["cg.tol"]);
  const double prec = as<double>(control["prec"]);
  const int report_freq = as<int>(control["report.freq"]);
  const int report_level = as<int>(control["report.level"]);
  const int report_precision = as<int>(control["report.precision"]);
  const int maxit = as<int>(control["maxit"]);
  const double contract_factor = as<double>(control["contract.factor"]);
  const double expand_factor = as<double>(control["expand.factor"]);
  const double contract_threshold = as<double>(control["contract.threshold"]);
  const double expand_threshold_rad = as<double>(control["expand.threshold.radius"]);
  const double expand_threshold_ap = as<double>(control["expand.threshold.ap"]);
  const double function_scale_factor = as<double>(control["function.scale.factor"]);
  const int precond_refresh_freq = as<int>(control["precond.refresh.freq"]);
  const int precond_ID = as<int>(control["preconditioner"]);
  const int quasi_newton_method = as<int>(control["quasi.newton.method"]);
  const int trust_iter = as<int>(control["trust.iter"]);

  Function fn(fn_);
  Function gr(gr_);

  Rfunc func(nvars, fn, gr);
  
  Map<VectorXd> startX(start.begin(),nvars); 
  
  // Control parameters for optimizer
  
  Trust_CG_Optimizer<Map<VectorXd>, Rfunc, optHessType, optPrecondType> opt(func,
									    startX, rad, min_rad, tol,
									    prec, report_freq,
									    report_level,
									    report_precision,
									    maxit, contract_factor,
									    expand_factor,
									    contract_threshold,
									    expand_threshold_rad,
									    expand_threshold_ap,
									    function_scale_factor,
									    precond_refresh_freq,
									    precond_ID,
									    quasi_newton_method,
									    trust_iter);

  opt.run();
  
 // collect results and return

  VectorXd P(nvars);
  VectorXd grad(nvars);

  // return sparse hessian information in CSC format

  double fval, radius;
  int iterations;
  MB_Status status;

  status = opt.get_current_state(P, fval, grad,
  				 iterations, radius);
  
  List res;
  res = List::create(Rcpp::Named("fval") = Rcpp::wrap(fval),
		     Rcpp::Named("solution") = Rcpp::wrap(P),
		     Rcpp::Named("gradient") = Rcpp::wrap(grad),	
		     Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		     Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		     Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		     Rcpp::Named("method") = Rcpp::wrap("quasi-newton"),
		     Rcpp::Named("hessian.update.method") = Rcpp::wrap(quasi_newton_method)
		     );
   
 return(res);

 END_R_INTERFACE 

    }


RcppExport SEXP get_fdfh(SEXP x_, SEXP fn_, SEXP gr_, SEXP hs_struct_, SEXP fd_method_, SEXP fd_eps_)
{
  BEGIN_R_INTERFACE

  using Rcpp::List;
  using Rcpp::NumericVector;
  using Rcpp::IntegerVector;
  using Rcpp::NumericMatrix;
  using Rcpp::Function;
  using Rcpp::as;
  using Eigen::Map;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::VectorXi;

  NumericVector x2(x_);
  int nvars = x2.size();

  Function fn(fn_);
  Function gr(gr_);

  List hs_struct(hs_struct_);
  int fd_method = Rcpp::as<int>(fd_method_);
  double fd_eps = Rcpp::as<double>(fd_eps_);

  // get hessian structure

  IntegerVector iRow_((SEXP)hs_struct["iRow"]);
  IntegerVector jCol_((SEXP)hs_struct["jCol"]);
  int nnz = iRow_.size();
  if (jCol_.size()!=nnz) throw MyException("Error:  Hessian index vectors must have same length\n",__FILE__, __LINE__);
 
  Map<VectorXi> iRow(iRow_.begin(),nnz);
  Map<VectorXi> jCol(jCol_.begin(),nnz);

  Rfunc func(nvars, fn, gr);
  func.hessian_init(iRow, jCol, fd_method, fd_eps);

  Map<VectorXd> x(x2.begin(),nvars);

  double val;
  VectorXd grad(nvars);
  SparseMatrix<double> hess(nvars, nvars);
  hess.reserve(nnz);

  func.get_fdf(x, val, grad);
  func.get_hessian(x, hess);
  
  List res = List::create(Rcpp::Named("fval") = val,
			  Rcpp::Named("gradient") = Rcpp::wrap(grad),
			  Rcpp::Named("hessian") = Rcpp::wrap(hess)
			  );

  return(res);

END_R_INTERFACE
}



RcppExport SEXP sparseTR2(SEXP start_, SEXP fn_, SEXP gr_, SEXP hs_,
			 SEXP control_) {
  
  BEGIN_R_INTERFACE
    
  // use this version when the user supplies his own Hessian function.
  //  Hessian function must return a Matrix of class dgCMatrix


    using Rcpp::NumericVector;
  using Rcpp::IntegerVector;
  using Rcpp::Function;
  using Rcpp::List;
  using Rcpp::as;
  
  using Eigen::VectorXi;
  using Eigen::Map;

  typedef double Scalar;

  typedef SparseMatrix<Scalar> optHessType;
  typedef SimplicialLLT<optHessType> optPrecondType; 

  NumericVector start(start_);
  int nvars = start.size();
  if (nvars<=0) throw MyException("Number of variables (starting values) must be positive\n",__FILE__,__LINE__);
 
  Function fn(fn_);
  Function gr(gr_);
  Function hs(hs_);

 // Control parameters for optimizer

  List control(control_);
  double rad = as<double>(control["start.trust.radius"]);
  const double min_rad = as<double>(control["stop.trust.radius"]);
  const double tol = as<double>(control["cg.tol"]);
  const double prec = as<double>(control["prec"]);
  const int report_freq = as<int>(control["report.freq"]);
  const int report_level = as<int>(control["report.level"]);
  const int report_precision = as<int>(control["report.precision"]);
  const int maxit = as<int>(control["maxit"]);
  const double contract_factor = as<double>(control["contract.factor"]);
  const double expand_factor = as<double>(control["expand.factor"]);
  const double contract_threshold = as<double>(control["contract.threshold"]);
  const double expand_threshold_rad = as<double>(control["expand.threshold.radius"]);
  const double expand_threshold_ap = as<double>(control["expand.threshold.ap"]);
  const double function_scale_factor = as<double>(control["function.scale.factor"]);
  const int precond_refresh_freq = as<int>(control["precond.refresh.freq"]);
  const int precond_ID = as<int>(control["preconditioner"]);
  const int trust_iter = as<int>(control["trust.iter"]);

  Map<VectorXd> startX(start.begin(),nvars); 
  Rcpp::S4 sh_  = hs(startX);
  MappedSparseMatrix<double> sh = Rcpp::as<MappedSparseMatrix<double> >(sh_);
  int nnz = (sh.nonZeros() + nvars)/2;

  RfuncHess func(nvars, nnz, fn, gr, hs);
  
  Trust_CG_Sparse<Map<VectorXd>, RfuncHess, optHessType, optPrecondType> opt(func,
									     startX, rad, min_rad, tol, prec,
									     report_freq,
									     report_level, report_precision,
									     maxit, contract_factor,
									     expand_factor,
									     contract_threshold,
									     expand_threshold_rad,
									     expand_threshold_ap,
									     function_scale_factor,
									     precond_refresh_freq,
									     precond_ID,
									     trust_iter);
  
  opt.run();

 // collect results and return

  VectorXd P(nvars);
  VectorXd grad(nvars);
  SparseMatrix<double> hess(nvars,nvars);
  hess.reserve(nnz);

  double fval, radius;
  int iterations;
  MB_Status status;

  status = opt.get_current_state(P, fval, grad, hess,
				 iterations, radius);

  List res;
  res = List::create(Rcpp::Named("fval") = Rcpp::wrap(fval),
		     Rcpp::Named("solution") = Rcpp::wrap(P),
		     Rcpp::Named("gradient") = Rcpp::wrap(grad),	
		     Rcpp::Named("hessian") = Rcpp::wrap(hess),
		     Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		     Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		     Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		     Rcpp::Named("nnz") = Rcpp::wrap(nnz),
		     Rcpp::Named("method") = Rcpp::wrap("Sparse")
		     );
   
  return(res);

  END_R_INTERFACE
  
}



#endif
