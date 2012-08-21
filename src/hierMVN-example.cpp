// hierMVN-example.cpp -- This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2012 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// See the trustOptim LICENSE file for more information.



#ifndef __TRUST_OPTIM_HIERMVN
#define __TRUST_OPTIM_HIERMVN

#include <Rcpp.h>
#include <Eigen/Core>
#include <norm.cpp>
#include <mvnorm.cpp>
#include <wishart.cpp>

template<typename TD>
class HierMVN  {


  typedef double Scalar; // Tpars is the type of the parameter vector
  typedef int Index; // Tpars is the type of the parameter vector
  
private:
  
  // data and priors
  
  const int & nvars;
  MatrixXd Y;
  MatrixXd X;
  MatrixXd mu_prior_mean;
  MatrixXd mu_prior_chol_prec;
  Scalar nu;
  MatrixXd chol_G_inv_prior_mean;

  
  Index N, T, k; // number of observations per unit, number of units, number of regressors
  
  
public:

  HierMVN(const int&, const TD&);

  template <typename Tpars>
  void get_f(const Eigen::MatrixBase<Tpars>&, const Scalar&);
  
  template <typename Tpars, typename Tgrad>
  void get_df(const Eigen::MatrixBase<Tpars>&, const Eigen::MatrixBase<Tgrad>&);
  
  template <typename Tpars, typename Tgrad>
  void get_fdf(const Eigen::MatrixBase<Tpars>&, const Scalar&, const Eigen::MatrixBase<Tgrad>&, const int&);

  template <typename Tpars>
  void get_data_LL(const Eigen::MatrixBase<Tpars>&, const Scalar&);

  int get_nvars();
  
};

template<typename TD>
HierMVN<TD>::HierMVN(const int& nvars_, const TD& params) :
  nvars(nvars_)
{

  using Rcpp::List;
  using Rcpp::NumericMatrix;
  using Rcpp::NumericVector;
  using Rcpp::as;

  List & pars = static_cast<List &>(const_cast<List &>(params));

  List data = as<List>(pars["data"]);
  List priors = as<List>(pars["priors"]);

  NumericMatrix Y_((SEXP)data["Y"]);
  NumericMatrix X_((SEXP)data["X"]);
  NumericVector mu_prior_mean_((SEXP)priors["mu.prior.mean"]);
  NumericMatrix mu_prior_chol_prec_((SEXP)priors["mu.prior.chol.prec"]);
  nu = as<double>((SEXP)priors["nu"]);
  NumericMatrix chol_G_inv_prior_mean_((SEXP)priors["chol.inv.mean.G"]);
  
  N = Y_.ncol();
  T = Y_.nrow();
  k = X_.nrow();

  X = MatrixXd::Map(X_.begin(),k,N);
  Y = MatrixXd::Map(Y_.begin(),T,N);
  mu_prior_mean = MatrixXd::Map(mu_prior_mean_.begin(),k,1);
  mu_prior_chol_prec = MatrixXd::Map(mu_prior_chol_prec_.begin(),k,k);
  chol_G_inv_prior_mean = MatrixXd::Map(chol_G_inv_prior_mean_.begin(),k,k);

}


template<typename TD>
int HierMVN<TD>::get_nvars()
{

  return(nvars);

}

template<typename TD>
template<typename Tpars>
void HierMVN<TD>::get_f(const MatrixBase<Tpars>& P, const Scalar& f_) {
  

  // return value of objective function in f
  // still computes the gradient.  need to find a way around that.
  
  Scalar & f = const_cast<Scalar&>(f_);
  VectorXd tmp(nvars);
  get_fdf(P,f,tmp,0);
  
}

template<typename TD>
template<typename Tpars, typename Tgrad>
void HierMVN<TD>::get_df(const MatrixBase<Tpars>& P, const MatrixBase<Tgrad>& df_) {
  
 
  // return gradient in df_
  
  MatrixBase<Tgrad> & df = const_cast<MatrixBase<Tgrad>& >(df_);
  
  Scalar f;  // get f anyway, since it's cheap relative to df.
  get_fdf(P, f, df,1);
  
}

template<typename TD>
template <typename Tpars, typename Tgrad>
void HierMVN<TD>::get_fdf(const MatrixBase<Tpars>& P_, const Scalar& f_, const MatrixBase<Tgrad>& df_, const int& gradFlag) {
  
  // return both value in f and gradient in df_
  // note:  do not call with an expression for P_.  It must be a column vector, alone.

  MatrixBase<Tpars>& P = const_cast<MatrixBase<Tpars>& >(P_);
  
  const Index kk12 = k*(k+1)/2;
  
 
  Scalar & f = const_cast<Scalar&>(f_);
  MatrixBase<Tgrad> & df = const_cast<MatrixBase<Tgrad>& >(df_);

  // Maps to parameter vectors

  Map< MatrixXd>B(&P(0), k, N); // map for k*N elements of P to B matrix

  Block<Tpars> mu = P.derived().block(k*N,0,k,1);
  Block<Tpars> vech_chol_G = P.derived().block(k*N+k, 0, kk12,1);


  RowVectorXd y_mean = (B.array() * X.array()).colwise().sum().matrix();

  norm_Distr< MatrixXd, RowVectorXd >  obj_data(Y, y_mean, 0.0);
  MVN_Distr< Map<MatrixXd >, Block<Tpars> >obj_B_prior(B, mu, vech_chol_G);
  MVN_Distr< Block<Tpars >, MatrixXd > obj_mu_prior(mu, mu_prior_mean, mu_prior_chol_prec);
  Wishart_Distr< Block<Tpars>, MatrixXd > obj_G_prior(vech_chol_G, nu, chol_G_inv_prior_mean);
  
  Scalar data_LL = obj_data.get_log_pdf();
  Scalar B_prior = obj_B_prior.get_log_pdf();
  Scalar mu_prior = obj_mu_prior.get_log_pdf();
  Scalar G_prior = obj_G_prior.get_log_pdf();
  
  f = data_LL + B_prior + mu_prior + G_prior;
  
  // computing gradients
  if (gradFlag>0) { 
    
    RowVectorXd dLL_dmean(N);
    MatrixXd dLL_dB(k,N);
    MatrixXd dB_prior_dB(k,N);
    VectorXd dB_prior_dmu(k);
    VectorXd dB_prior_dvech_chol_G(kk12);
    VectorXd dmu_prior_dmu(k);
    VectorXd dG_prior_dvech_chol_G(kk12);
 
    obj_data.get_dlogpdf_dmu(dLL_dmean);
 
    dLL_dB = X * dLL_dmean.row(0).asDiagonal(); // chain rule

    obj_B_prior.get_dlogpdf_dy(dB_prior_dB);
    obj_B_prior.get_dlogpdf_dmu(dB_prior_dmu);
    obj_B_prior.get_dlogpdf_dvechcholprec(dB_prior_dvech_chol_G);
    obj_mu_prior.get_dlogpdf_dy(dmu_prior_dmu);
    obj_G_prior.get_dlogpdf_dvechcholG(dG_prior_dvech_chol_G);

    // maps to gradient vectors

    Map< MatrixXd >df_B(&df(0), k, N);
    Map< MatrixXd >df_mu(&df(k*N), k, 1);
    Map< MatrixXd >df_dG(&df(k*(N+1)), kk12, 1);

    df_B = dLL_dB + dB_prior_dB;
    df_mu = dB_prior_dmu + dmu_prior_dmu;
    df_dG = dB_prior_dvech_chol_G + dG_prior_dvech_chol_G;

  }
				    
}


template<typename TD>
template <typename Tpars>
void HierMVN<TD>::get_data_LL(const MatrixBase<Tpars>& P_, const Scalar& f_) {
  
  // return both value in f and gradient in df_
  // note:  do not call with an expression for P_.  It must be a column vector, alone.

  MatrixBase<Tpars>& P = const_cast<MatrixBase<Tpars>& >(P_);
  Scalar & f = const_cast<Scalar&>(f_);
 
  // Maps to parameter vectors

  Map< MatrixXd>B(&P(0), k, N); // map for k*N elements of P to B matrix
  
  RowVectorXd y_mean = (B.array() * X.array()).colwise().sum().matrix();

  norm_Distr< MatrixXd, RowVectorXd >  obj_data(Y, y_mean, 0.0);
  Scalar data_LL = obj_data.get_log_pdf();
 
  f = data_LL;
				    
}

RcppExport SEXP HierMVN_get_fdf(SEXP PR, SEXP params_)
{
  
  using Rcpp::List;
  using Rcpp::NumericVector;
  using Eigen::Map;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
 
  typedef HierMVN<List> FuncType;

  NumericVector P_(PR);
  int nvars = P_.size();
  List params(params_);

  Map<VectorXd> P(P_.begin(),nvars);
  
  FuncType mod(nvars, params);
   
  double val;
  VectorXd grad(nvars);
  NumericVector grad_(nvars);
 
   
  mod.get_fdf(P, val, grad, 1);
 
  std::copy(grad.data(),grad.data()+grad.size(),grad_.begin());  
  List res = List::create(Rcpp::Named("val") = val,
			  Rcpp::Named("grad") = Rcpp::wrap(grad_)
			  );
  return(res);

}


RcppExport SEXP HierMVN_get_f(SEXP PR, SEXP params_)
{
  
  using Rcpp::List;
  using Rcpp::NumericVector;
  using Eigen::Map;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
 
  typedef HierMVN<List> FuncType;

  NumericVector P_(PR);
  int nvars = P_.size();
  List params(params_);

  Map<VectorXd> P(P_.begin(),nvars);
  
  FuncType mod(nvars, params);
   
  double val;   
  mod.get_f(P, val);
   
  return(Rcpp::wrap(val));

}




#endif



