#pragma once

#include "../common_config.hpp"

#ifndef GEMM_ARMA_MAT
#define GEMM_ARMA_MAT(N)                                                \
NONIUS_BENCHMARK(STRING(ARMA_GEMM_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, -1, -1> A =                                   \
            Eigen::Matrix<double, -1, -1>::Random(N, N);                \
    gemm_arma_mat<N> bm(A);                                             \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_, bm._a_t_a_); });                   \
})
#endif

template<int N>
struct gemm_arma_mat {
    gemm_arma_mat(const Eigen::Matrix<double, -1, -1>& A)
            : _a_(arma::mat::fixed<N, N>(A.data())) { }
    double impl(const arma::mat& A, arma::mat& ATA) {
        ATA = A.t() * A;
        return arma::accu(ATA);
    }
    arma::mat _a_;
    arma::mat _a_t_a_ = arma::mat::fixed<N, N>();
};
