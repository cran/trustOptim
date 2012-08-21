//
// mvnorm.cpp.  Part of the trustOptim package for the R programming language.
//
// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2012 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// See the trustOptim LICENSE file for more information.

#ifndef MVNORM
#define MVNORM

#include <iostream>
#include <Eigen/Core>
#include "common_R.hpp"

#define LOG_2_PI 1.837877066409345339082

using namespace std;
using namespace Eigen;

template<typename yBase, typename muBase>
class MVN_Distr {

  typedef typename yBase::Scalar Scalar;
  typedef typename yBase::Index Index;

  typedef Matrix<Scalar, yBase::RowsAtCompileTime, yBase::ColsAtCompileTime> yType;
  typedef Matrix<Scalar, muBase::RowsAtCompileTime, muBase::ColsAtCompileTime> muType;
  typedef Matrix<Scalar, Dynamic, Dynamic> precType;


private:
  
  const MatrixBase<yBase> & y;
  const MatrixBase<muBase> & mu;
  Index N, k;

  precType chol_prec; // Cholesky of precision matrix, lower triangle only
  yType Z;
  Scalar sum_log_diag_chol_prec, normConst;
  int isStandardized;
  void standardize_MVN();

public:

  template<typename cholType>
  MVN_Distr(const MatrixBase<yBase>&, const MatrixBase<muBase>&, const MatrixBase<cholType>&);

  Scalar get_log_pdf();

  template<typename outType>
  void get_dlogpdf_dmu(const MatrixBase<outType>&); // not really const

  template<typename outType>
  void get_dlogpdf_dy(const MatrixBase<outType>&); // not really const

    template<typename outType>
   void get_dlogpdf_dvechcholprec(const MatrixBase<outType>&);
};

template<typename yBase, typename muBase>
void
MVN_Distr<yBase, muBase>::standardize_MVN(){
  using namespace Eigen;
  if (mu.cols()==N) {
    Z = chol_prec.template triangularView<Lower>().transpose() * (y-mu); // do I need the noalias here?
    
  } else {
    if (mu.cols()==1) {
      Z = chol_prec.template triangularView<Lower>().transpose() * (y.colwise()-mu.col(0)); // do I need the noalias here?
     
 
    } else {
      throw MyException("error in standardize_MVN",__FILE__, __LINE__);
    }
  }
}

template<typename yBase, typename muBase>
template<typename cholType>
MVN_Distr<yBase, muBase>::MVN_Distr(const MatrixBase<yBase> &y_, const MatrixBase<muBase> &mu_,
				    const MatrixBase<cholType> &S_)

  : y(y_), mu(mu_)
    
{
  k = y.rows();
  N = y.cols();

  chol_prec = precType::Zero(k, k);
  
  if (S_.cols()==1 && S_.rows()==k*(k+1)/2) {
    // S_ is vech(chol(prec))
    Index i, j, idx=0;
    for (j=0; j<k; j++) {
      for (i=j; i<k; i++) {
	chol_prec(i,j)=S_(idx);
	idx++;
      }
    }
  } else {
    if (S_.cols()==k && S_.rows()==k) {
      // S_ is chol(prec) (lower triangle)
      chol_prec = S_; // make a copy?
    } else { 
      // throw exception
    }
  }
 
  sum_log_diag_chol_prec = chol_prec.diagonal().array().log().sum();
 
  Z = yType::Zero(k,N);
  isStandardized = 0;
  normConst = -N*k*LOG_2_PI/2;
 
}


template<typename yBase, typename muBase>
typename MVN_Distr<yBase, muBase>::Scalar MVN_Distr<yBase, muBase>::get_log_pdf() {
 
  if (isStandardized==0) {
    standardize_MVN();
    isStandardized = 1;
  }
  
  Scalar res =  normConst + N*sum_log_diag_chol_prec - Z.squaredNorm()/2;  
 
  return(res); 
    
}

template<typename yBase, typename muBase>
template<typename outType>
void MVN_Distr<yBase, muBase>::get_dlogpdf_dmu(const MatrixBase<outType>  & out){

  MatrixBase<outType> & dmu = const_cast< MatrixBase<outType>& >(out);

  // summing across multiple y observations

  if (isStandardized==0) {
    standardize_MVN();
    isStandardized = 1;
  }
 
  if (mu.cols()==1) {
  
    yType tmp =  chol_prec.template triangularView<Lower>() * Z;
   
    dmu = tmp.rowwise().sum();   // sum across N people
  } else {
    if (mu.cols()==N) {
      dmu = chol_prec.template triangularView<Lower>() * Z;  // result for each person     
    } else {
      // throw exception
    }
  }
}

template<typename yBase, typename muBase>
template<typename outType>
void MVN_Distr<yBase, muBase>::get_dlogpdf_dy(const MatrixBase<outType>  & out){

  MatrixBase<outType> & dy = const_cast< MatrixBase<outType>& >(out);
  // summing across multiple y observations

  if (isStandardized==0) {
    standardize_MVN();
    isStandardized = 1;
  }
  dy = chol_prec.template triangularView<Lower>() * Z;  // result for each person  
  dy.array() *= -1.0;
}

template<typename yBase, typename muBase>
template<typename outType>
void MVN_Distr<yBase, muBase>::get_dlogpdf_dvechcholprec(const MatrixBase<outType>  & out_) {

    MatrixBase<outType> & out = const_cast< MatrixBase<outType>& >(out_);

  if (isStandardized==0) {
    standardize_MVN();
    isStandardized = 1;
  }
    
  Matrix<Scalar, Dynamic, Dynamic> tmp(k,k);
  tmp.setIdentity(k,k);
  tmp *= N;
  tmp -=  Z*Z.transpose();

  chol_prec.template triangularView<Lower>().transpose().solveInPlace(tmp);

  int idx = 0;
  for (int j=0; j<k; j++) {
    for (int i=j; i<k; i++) {
      out(idx) = tmp(i,j);
      idx++;
    }
  }

}

#endif
