//
// norm.cpp. --  This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2012 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// See the trustOptim LICENSE file for more information. Part of the trustOptim package for the R programming language.

#ifndef __NORM
#define __NORM

#include <iostream>
#include <Eigen/Core>

#define LOG_2_PI 1.837877066409345339082 // log(2*pi)

using namespace std;
using namespace Eigen;

template<typename yBase, typename meanBase>
class norm_Distr {

  typedef typename yBase::Scalar Scalar; 
  typedef typename yBase::Index Index;
  typedef Matrix<Scalar, yBase::RowsAtCompileTime, yBase::ColsAtCompileTime> yType;
  
private:

  const MatrixBase<yBase> & y;  
  const MatrixBase<meanBase> & mu;  
  const Scalar & log_sig;  

  yType ymu;

  Scalar sig, sig2, normConst;
  Index N, k;

  int is_swept;
  void sweep_mean();

public:

  norm_Distr(const MatrixBase<yBase> &, const MatrixBase<meanBase> &, const Scalar &);
  Scalar get_log_pdf();

  template<typename outType>
  void get_dlogpdf_dmu(const MatrixBase<outType> &);

  template<typename outType>
  void get_dlogpdf_dy(const MatrixBase<outType> &);

  void get_dlogpdf_dsig(Scalar &);
  void get_dlogpdf_dlogsig(Scalar &);
  void get_dlogpdf_dlogsig2(Scalar &);

};

template<typename yBase, typename meanBase>
norm_Distr<yBase, meanBase>::norm_Distr(const MatrixBase<yBase> & y_,
					const MatrixBase<meanBase> & mu_,
					const Scalar & log_sig_) :
  y(y_), mu(mu_), log_sig(log_sig_)

{

  N = y.cols();
  k = y.rows();
  sig = exp(log_sig);
  sig2 = exp(2*log_sig);
  normConst = -k*N*LOG_2_PI/2;
  
  ymu = yType::Zero(k,N);
  is_swept = 0;


}

template<typename yBase, typename meanBase>
void norm_Distr<yBase, meanBase>::sweep_mean() {

  if (mu.rows()==k) {
    if (mu.cols()==N) {
      ymu = y-mu;
    } else {
      if (mu.cols()==1) {
	ymu = y.colwise() - mu.col(0);
      } else {
	// throw exception
      }
    }
  } else {
    if (mu.rows()==1) {
      if (mu.cols()==N) {
	ymu = y.rowwise() - mu.row(0);	
      } else {
	if (mu.cols()==1) {
	  ymu = y-mu;
	} else {
	  // throw exception
	}
      }
    }
  }
}
template<typename yBase, typename meanBase>
typename norm_Distr<yBase, meanBase>::Scalar norm_Distr<yBase, meanBase>::get_log_pdf() {

  if (is_swept==0) {
    sweep_mean();
    is_swept = 1;
  }

  Scalar res = normConst - k*N*log_sig - ymu.squaredNorm()/(2*sig2);
  return(res);

}
 
template<typename yBase, typename meanBase>
template<typename outType> 
void norm_Distr<yBase, meanBase>::get_dlogpdf_dy(const MatrixBase<outType> & out_) {

  MatrixBase<outType> & out = const_cast< MatrixBase<outType>& >(out_);

  if (is_swept==0) {
    sweep_mean();
    is_swept = 1;
  }
  
  out = -ymu / sig2;
}


template<typename yBase, typename meanBase>
template<typename outType> 
void norm_Distr<yBase, meanBase>::get_dlogpdf_dmu(const MatrixBase<outType> & out_) {

  MatrixBase<outType> & out = const_cast< MatrixBase<outType>& >(out_);

  if (is_swept==0) {
    sweep_mean();
    is_swept = 1;
  }

  if (mu.rows()==k) {
    if (mu.cols()==N) {
      out = ymu / sig2;
    } else {
      if (mu.cols()==1) {
	out = ymu.rowwise().sum() / sig2;
      } else {
	// throw exception
      }
    }
  } else {
    if (mu.rows()==1) {
      if (mu.cols()==N) {
	out = ymu.colwise().sum() / sig2;
      } else {
	if (mu.cols()==1) {
	  out = ymu / sig2;
	} else {
	  // throw exception
	}}
    } else {
      // throw exception
    }
  } 
} 

template<typename yBase, typename meanBase>
void norm_Distr<yBase, meanBase>::get_dlogpdf_dsig(Scalar & out) {
  
  if (is_swept==0) {
    sweep_mean();
    is_swept = 1;
  }

  out = (ymu.squaredNorm()/sig2 - k*N) / sig;

}

template<typename yBase, typename meanBase>
void norm_Distr<yBase, meanBase>::get_dlogpdf_dlogsig(Scalar & out) {
  
  if (is_swept==0) {
    sweep_mean();
    is_swept = 1;
  }

  get_dlogpdf_dsig(out);
  out = out*sig;
  
}

template<typename yBase, typename meanBase>
void norm_Distr<yBase, meanBase>::get_dlogpdf_dlogsig2(Scalar & out) {
  

  if (is_swept==0) {
    sweep_mean();
    is_swept = 1;
  }

  get_dlogpdf_dsig(out);
  out = out * sig / 2;
  
}


#endif


