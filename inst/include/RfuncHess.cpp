// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2012 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef __TRUST_OPTIM_RFUNC_HESS__
#define __TRUST_OPTIM_RFUNC_HESS__

// uses a user-supplied function for the Hessian

#include<RcppEigen.h>
#include <common_R.hpp>
#include <include_eigen.hpp>

using Eigen::Matrix;
using Eigen::MatrixBase;
using Eigen::MappedSparseMatrix;
using Eigen::Dynamic;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class RfuncHess {

  int nvars; 
  const Rcpp::Function & fn;
  const Rcpp::Function & gr;
  const Rcpp::Function & hs;
  
 protected:
  
  VectorXi iRow; // row indices of nonzero elements
  VectorXi jCol; // col indices of nonzero elements
  VectorXi jpntr; // for each col, pointer to first element
  VectorXd fhes; // values for nonzero elements
  
  int nnz;
  
  VectorXd tmp1;
  VectorXd tmp2;

 public:

  RfuncHess(const int, const int, const Rcpp::Function&, const Rcpp::Function&, const Rcpp::Function&);

  ~RfuncHess();

  template <typename Tpars>
    void get_f(const Eigen::MatrixBase<Tpars>&, const double&);
  
  template <typename Tpars, typename Tgrad>
    void get_df(const Eigen::MatrixBase<Tpars>&, const Eigen::MatrixBase<Tgrad>&);
  
  template <typename Tpars, typename Tgrad>
    void get_fdf(const Eigen::MatrixBase<Tpars>&, const double&, const Eigen::MatrixBase<Tgrad>&);


  template<typename Tpars, typename Tout>
    void get_hessian(const Eigen::MatrixBase<Tpars>&, const Eigen::SparseMatrixBase<Tout>&);

  int get_nnz();  
};

RfuncHess::RfuncHess(const int nvars_, const int nnz_,
		 const Rcpp::Function & fn_,
		 const Rcpp::Function & gr_,
		 const Rcpp::Function & hs_) :
  nvars(nvars_), fn(fn_), gr(gr_), hs(hs_), nnz(nnz_)
{
}


int RfuncHess::get_nnz() {
  return(nnz);
}


template<typename Tpars>
void RfuncHess::get_f(const MatrixBase<Tpars>& P_, const double& f_) {
  
  Eigen::MatrixBase<Tpars>& P = const_cast<Eigen::MatrixBase<Tpars>& >(P_);
  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n", __FILE__, __LINE__);
 
  double &f = const_cast<double&>(f_);

  Rcpp::NumericVector pars(P.derived().data(), P.derived().data() + P.derived().size());
  
  double res = Rcpp::as<double>(fn(pars));
  f = res;
  return;

}


template<typename Tpars, typename Tgrad>
void RfuncHess::get_df(const MatrixBase<Tpars>& P_, const MatrixBase<Tgrad>& df_) {
  
  using Rcpp::NumericVector;
  using Eigen::VectorXd;

  Eigen::MatrixBase<Tpars>& P = const_cast<Eigen::MatrixBase<Tpars>& >(P_);
  Eigen::MatrixBase<Tgrad> & df = const_cast<Eigen::MatrixBase<Tgrad>& >(df_);
  
  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n", __FILE__, __LINE__);
  if (df.size()!=nvars) throw MyException("Incorrect gradient length\n", __FILE__, __LINE__);
  
  NumericVector pars(P.derived().data(), P.derived().data() + P.size());
  
  NumericVector grad_  = gr(pars);
  
  VectorXd grad = VectorXd::Map(grad_.begin(), nvars);
  
  df = grad;
  
  return;  
}


template<typename Tpars, typename Tgrad>
void RfuncHess::get_fdf(const Eigen::MatrixBase<Tpars>& P_, const double& f_,
		    const Eigen::MatrixBase<Tgrad>& df_)
{

  get_f(P_, f_);
  get_df(P_, df_);
  return;    
  
}

template<typename TP, typename Tout>
void RfuncHess::get_hessian(const MatrixBase<TP>& P,
			    const SparseMatrixBase<Tout>& out_
			    )
{
  
  SparseMatrixBase<Tout>& out = const_cast<SparseMatrixBase<Tout>&>(out_);
  
  using Rcpp::NumericVector;
  using Eigen::VectorXd;

  if (P.size()!=nvars) throw MyException("Incorrect number of parameters\n", __FILE__, __LINE__);
  
  NumericVector pars(P.derived().data(), P.derived().data() + P.size());
  
  Rcpp::S4 sh_ = hs(pars);
  MappedSparseMatrix<double> sh(Rcpp::as<MappedSparseMatrix<double> >(sh_));
  out = sh.selfadjointView<Lower>();
}



RfuncHess::~RfuncHess(){
}

#endif




