//
// wishart.cpp.   this file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2012 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.
// See the trustOptim LICENSE file for more information.


#ifndef __WISHART
#define __WISHART

#include <iostream>
#include <Eigen/Core>

#define LOG_2_PI 1.837877066409345339082 // log(2*pi)
#define LOG_2 0.69314718056
#define LOG_PI_4 0.28618247146 // log(pi)/4

using namespace std;
using namespace Eigen;

template<typename GBase, typename meanBase>
class Wishart_Distr {


/* Computes   Wishart
Pass in the cholesky of the Sigma (variable of interest)

G = variable of interest.  Typically the precision matrix
L = Cholesky of G
B  = Cholesky of the inverse of mean of G
nu = degrees of freedom

If y is a vector, it is the vech of L.  If it is a matrix, it is L (lower triangle)

*/



  typedef typename GBase::Scalar Scalar;
  typedef typename GBase::Index Index;

  typedef Matrix<Scalar, Dynamic, Dynamic> GType;

private:

  GType chol_G; // Cholesky of precision matrix, lower triangle only
  Index k;
  const MatrixBase<meanBase> & chol_inv_mean;
  const Scalar & nu; 
  Scalar sum_log_diag_chol_G, normConst;
  
  int isStandardized;
  GType Z; // LZ = chol_mean
  
public:


  Wishart_Distr(const MatrixBase<GBase> &, const Scalar &, const MatrixBase<meanBase> &);

  void standardize_Wishart();
  typename Wishart_Distr<GBase, meanBase>::Scalar get_log_pdf();

  template<typename DType>
  void get_dlogpdf_dvechcholG(const MatrixBase<DType> &);

};

template<typename GBase, typename meanBase>
Wishart_Distr<GBase, meanBase>::Wishart_Distr(const MatrixBase<GBase> & chol_G_,
						    const Scalar & nu_,
						    const MatrixBase<meanBase> & chol_mean_) :
  chol_inv_mean(chol_mean_), nu(nu_) {
  
 k = chol_inv_mean.cols();
  
  if (chol_G_.cols()==1) {

    // input is vech(chol(S))
    chol_G = GType::Zero(k,k);
    Index idx=0;
    for (Index j=0; j<k; j++) {
      for (Index i=j; i<k; i++) {
	chol_G(i,j) = chol_G_(idx);
	idx++;
      }
    }

  } else {
    if ( (chol_G_.cols()==k) & (chol_G_.rows()==k) ) {

          chol_G = chol_G_; // create a copy?

    } else {
      throw MyException("Error constructing Wishart distribution",__FILE__,__LINE__);
      // throw exception
    }
  }

  sum_log_diag_chol_G = chol_G.diagonal().array().log().sum();

  normConst = nu*k*((log(nu)-LOG_2)/2) + nu*chol_inv_mean.diagonal().array().log().sum() - k*(k-1)*LOG_PI_4;
 
  for (int i=1; i<=k; i++) {
    normConst -= lgamma((nu+1-i)/2);
  }
 
  Z = GType::Zero(k,k);
  isStandardized = 0;

}


template<class GBase, typename meanBase>
void  Wishart_Distr<GBase, meanBase>::standardize_Wishart(){

   Z = chol_inv_mean.template triangularView<Lower>().transpose() * chol_G;  




}


template<typename GBase, typename meanBase>
typename Wishart_Distr<GBase, meanBase>::Scalar Wishart_Distr<GBase, meanBase>::get_log_pdf() {

  if (isStandardized==0) {
    standardize_Wishart();
    isStandardized = 1;
  }

  Scalar res = normConst + (nu-k-1)*sum_log_diag_chol_G;
  res -= nu * Z.squaredNorm()/2;



  return(res);

}

template<typename GBase, typename meanBase>
template<typename DType>
void Wishart_Distr<GBase, meanBase>::get_dlogpdf_dvechcholG(const MatrixBase<DType> & out_){
  
  MatrixBase<DType> & out = const_cast< MatrixBase<DType>& >(out_);
  
    if (isStandardized==0) {
      standardize_Wishart();
      isStandardized = 1;
    }

    GType tmp;
    tmp.setIdentity(k,k);
    tmp.array() *= (k+1-nu)/nu;
    tmp += Z.transpose() * Z;
    chol_G.template triangularView<Lower>().transpose().template solveInPlace<OnTheLeft>(tmp);
    tmp.array() *= -nu;

    Index idx=0;
    for (Index j=0; j<k; j++) {
      for (Index i=j; i<k; i++) {
	out(idx) = tmp(i,j);
	idx++;
      }
    }
 
}


#endif



