#include <cstdlib>

#define NONIUS_RUNNER
#include <nonius/nonius_single.h++>

#ifdef ARMA_DONT_USE_BLAS
#undef ARMA_DONT_USE_BLAS
#endif
#ifndef ARMA_USE_BLAS
#define ARMA_USE_BLAS
#endif
#include <armadillo>

#include <Eigen/Eigen>

#define STRING(str) #str

template<int dims>
struct bm_chol_arma {
    bm_chol_arma(
            const Eigen::Matrix<double, dims, dims>& A)
        : A(arma::mat::fixed<dims, dims>(A.data())) { }
    const double* impl(const arma::mat& A, arma::mat& ATA) {
        ATA = A.t() * A;
        return ATA.memptr();
    }
    arma::mat A;
    arma::mat ATA = arma::mat::fixed<dims, dims>();
};

template<int dims>
struct bm_chol_eigen {
    bm_chol_eigen(
            const Eigen::Matrix<double, dims, dims>& A)
        : A(A) { }
    const double* impl(
            const Eigen::Matrix<double, dims, dims>& A,
                  Eigen::Matrix<double, dims, dims>& ATA) {
        ATA.noalias() = A.transpose().lazyProduct(A);
        return ATA.data();
    }
    Eigen::Matrix<double, dims, dims> A;
    Eigen::Matrix<double, dims, dims> ATA;
};
