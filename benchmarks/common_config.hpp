#pragma once

#include <cstdlib>
#include <iostream>

#define NONIUS_RUNNER
#include <nonius/nonius_single.h++>

#ifdef ARMA_DONT_USE_BLAS
#undef ARMA_DONT_USE_BLAS
#endif

#ifndef ARMA_USE_BLAS
#define ARMA_USE_BLAS
#endif

#ifndef ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_OPENMP
#endif

#ifndef ARMA_NO_DEBUG
#define ARMA_NO_DEBUG
#endif

#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <armadillo>
#include <Eigen/Eigen>

#ifndef STRING
#define STRING(str) #str
#endif

