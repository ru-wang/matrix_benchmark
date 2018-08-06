#pragma once

#include <cstdlib>
#include <iostream>

#define NONIUS_RUNNER
#include <nonius/nonius_single.h++>



/*****************************************************
 *                     Armadillo                     *
 *****************************************************/

#ifndef ARMA_DONT_USE_WRAPPER
#   define ARMA_DONT_USE_WRAPPER
#endif

#ifndef ARMA_DONT_USE_OPENMP
#   define ARMA_DONT_USE_OPENMP
#endif

#ifdef ARMA_DONT_USE_LAPACK
#   undef ARMA_DONT_USE_LAPACK
#endif
#ifndef ARMA_USE_LAPACK
#   define ARMA_USE_LAPACK
#endif

#ifdef ARMA_DONT_USE_BLAS
#   undef ARMA_DONT_USE_BLAS
#endif
#ifndef ARMA_USE_BLAS
#   define ARMA_USE_BLAS
#endif

#ifdef NDEBUG
#   ifndef ARMA_NO_DEBUG
#   define ARMA_NO_DEBUG
#   endif
#endif

/*****************************************************
 *                       Eigen                       *
 *****************************************************/

#ifndef EIGEN_DONT_PARALLELIZE
#   define EIGEN_DONT_PARALLELIZE
#endif

#include <armadillo>
#include <Eigen/Eigen>

#ifndef STRING
#   define STRING(str) #str
#endif
