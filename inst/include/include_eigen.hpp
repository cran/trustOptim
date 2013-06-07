// include_eigen.hpp - part of the trustOptim package for the R programming language.
//
// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2013 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef __TRUST_OPTIM_INCLUDE_EIGEN
#define __TRUST_OPTIM_INCLUDE_EIGEN

#include "exceptions.hpp"

inline bool my_ret_bool(bool x) {return(x);}

#define my_assert(x) do { \
    if(!my_ret_bool(x) ) throw MyException(EIGEN_MAKESTRING(x), __FILE__, __LINE__); \
} \
while (false);

#ifdef eigen_assert
#undef eigen_assert
#endif
#define eigen_assert(x) my_assert(x)



#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>




#endif
