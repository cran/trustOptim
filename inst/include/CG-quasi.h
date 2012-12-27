// CG-quasi.cpp.  Part of the trustOptim package for the R programming language.
//
// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2012 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef __TRUST_OPTIM_CGQUASI
#define __TRUST_OPTIM_CGQUASI

#include <iostream>
#include <iomanip>

#include <include_eigen.hpp>

using namespace Eigen;
using namespace std;

template<typename TP, typename TFunc, typename THess, typename TPreLLt>  // TP is type for the parameter vector
  class Trust_CG_Optimizer : public Trust_CG_Base<TP, TFunc, THess, TPreLLt> {

  typedef typename TP::Scalar Scalar;
  typedef typename TP::Index Index;
  

  typedef Matrix<Scalar, Dynamic, 1> Vec;
  typedef Matrix<Scalar, Dynamic, Dynamic> Mat;

  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::Bk;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::nvars;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::func;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::function_scale_factor;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::xk;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::precond_ID;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::PrecondLLt;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::status;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::yk;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::sk;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::rad;
  using Trust_CG_Base<TP, TFunc, THess, TPreLLt>::iter;

private:

  const int & quasi_newton_method;
  
  Mat Precond; // current preconditioning matrix for trust sub-problem
  void init_precond();
  void init_precond_identity();
  void init_precond_Cholesky();

  void update_precond();
  void update_precond_identity();
  void update_precond_Cholesky();

  void updateHessian_SR1();
  void updateHessian_BFGS();

  Vec work; // some workspace

  void update_hessian();

public:



  // control values for the optimizer

  Trust_CG_Optimizer(TFunc&, const MatrixBase<TP>&,	
		     const double &, const double &, const double &,
		     const double &, const int &, const int , const int &,
		     const int &, const double &,
		     const double &, const double &,
		     const double &, const double &, const double &,
		     const int &, const int &, const int &, const int &);


 
  template<typename Tvec>
  MB_Status get_current_state(const MatrixBase<Tvec> &, double &, const MatrixBase<Tvec> &,
			      int &, double &);


};


template<typename TP, typename TFunc, typename THess, typename TPreLLt>
  Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::Trust_CG_Optimizer(TFunc & func_,
								    const MatrixBase<TP>& startX_,
								    const double & rad_,
								    const double & min_rad_,
								    const double & tol_,
								    const double & prec_,
								    const int & report_freq_,
								    const int  report_level_,
								    const int & report_precision_,
								    const int & maxit_,
								    const double & contract_factor_,
								    const double & expand_factor_,
								    const double & contract_threshold_,
								    const double & expand_threshold_rad_,
								    const double & expand_threshold_ap_,
								    const double & function_scale_factor_,
								    const int & precond_refresh_freq_,
								    const int & precond_ID_,
								    const int & quasi_newton_method_,
								    const int & trust_iter_) :
  Trust_CG_Base<TP, TFunc, THess, TPreLLt>(func_, startX_, rad_, min_rad_, tol_, prec_, report_freq_,
					   report_level_,
					   report_precision_, maxit_, contract_factor_,
					   expand_factor_, contract_threshold_,
					   expand_threshold_rad_, expand_threshold_ap_,
					   function_scale_factor_, precond_refresh_freq_, precond_ID_, trust_iter_),
    quasi_newton_method(quasi_newton_method_)    
  {


  Bk.resize(nvars, nvars);  
  Bk.setIdentity();
  Precond.resize(nvars, nvars);
  init_precond();
  work.resize(nvars);

  }

template<typename TP, typename TFunc, typename THess, typename TPreLLt>
  template<typename Tvec>
  MB_Status Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::get_current_state(const MatrixBase<Tvec> & pars_,
							     double & fval,
							     const MatrixBase<Tvec> & grad_,
							     int & iterations,
							     double & radius)
 {

   MatrixBase<Tvec> & pars = const_cast<MatrixBase<Tvec>& >(pars_);
   MatrixBase<Tvec> & grad = const_cast<MatrixBase<Tvec>& >(grad_);


   // note:  results are NOT multiplied by the function scale factor

   pars = xk;

   func.get_fdf(xk, fval, grad);

   iterations = iter;
   radius = rad;

   return status;
}

template<typename TP, typename TFunc, typename THess, typename TPreLLt>
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::init_precond() {

  switch (precond_ID) {

  case 0:
    init_precond_identity();
    break;
  case 1:
    init_precond_Cholesky();
    break;
  default:
    init_precond_identity();
  }

  return;
}

template<typename TP, typename TFunc, typename THess, typename TPreLLt>
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::update_precond() {

  switch (precond_ID) {

  case 0:
    update_precond_identity();
    break;
  case 1:
    update_precond_Cholesky();
    break;
  default:
    update_precond_identity();
  }
    
  return;
}


template<typename TP, typename TFunc, typename THess, typename TPreLLt> 
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::init_precond_Cholesky() {

  // this is a special case of a diagonal preconditioner
  // in general, preconditioner must be lower triangular and positive definite

  update_precond_Cholesky();

  return;
}


template<typename TP, typename TFunc, typename THess, typename TPreLLt> 
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::update_precond_Cholesky() {
 
  PrecondLLt.compute(Bk);

  return;
}


template<typename TP, typename TFunc, typename THess, typename TPreLLt> 
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::init_precond_identity() {

  // this is a special case of a diagonal preconditioner
  // in general, preconditioner must be upper triangular and positive definite

  Precond.setIdentity();
  PrecondLLt.compute(Precond);
  return;

}


template<typename TP, typename TFunc, typename THess, typename TPreLLt> 
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::update_precond_identity() {

  // this is a special case of a diagonal preconditioner
  
  // Do nothing.  Preconditioner will always be identity matrix
 
  return;
}


template<typename TP, typename TFunc, typename THess, typename TPreLLt>
    void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::update_hessian() {
   
  if (quasi_newton_method==1) {
    updateHessian_SR1();
  } else { 
    if (quasi_newton_method==2) {
      updateHessian_BFGS();
    }
  } 
}



template<typename TP, typename TFunc, typename THess, typename TPreLLt>
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::updateHessian_SR1()
{

 // // SR1 update of Hessian for next iteration
  double denom, crit;
 
  // using wd workspace for y'Bs
  work = yk - Bk.template selfadjointView<Lower>()*sk;
  denom = work.dot(sk);

  crit = 1e-7 * sk.norm()  * work.norm();
    
  if (abs(denom) > crit) {   // do the SR1 update only if denominator is large enough
  
    Bk.template selfadjointView<Lower>().rankUpdate(work,1/denom);
  }

  return;
}


template<typename TP, typename TFunc, typename THess, typename TPreLLt>
void Trust_CG_Optimizer<TP, TFunc, THess, TPreLLt>::updateHessian_BFGS()
{

 // // BFGS update of Hessian for next iteration
 // using wd workspace for Bs
  Scalar ys = yk.dot(sk);
  work = Bk.template selfadjointView<Lower>() * sk;
  Scalar sBs = sk.dot(work);


  Bk.template selfadjointView<Lower>().rankUpdate(work,-1.0/sBs);
  Bk.template selfadjointView<Lower>().rankUpdate(yk,1.0/ys);

  return;
}


#endif
