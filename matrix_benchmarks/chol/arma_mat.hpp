#pragma once

#include "../common_config.hpp"

#ifndef CHOL_ARMA_MAT
#define CHOL_ARMA_MAT(N)                                                \
NONIUS_BENCHMARK(STRING(ARMA_CHOL_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, -1, -1> A =                                   \
            Eigen::Matrix<double, -1, -1>::Random(N, N);                \
    chol_arma_mat<N> bm(A);                                             \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_t_a_, bm._l_t_); });                 \
    std::cerr << bm._l_t_(0, 0);                                        \
})
#endif

template<int N>
struct chol_arma_mat {
    chol_arma_mat(const Eigen::Matrix<double, -1, -1>& A)
            : _a_(arma::mat::fixed<N, N>(A.data())) {
        _a_t_a_.eye();
        _a_t_a_ += _a_.t() * _a_;
    }
    void impl(const arma::mat& ATA, arma::mat& LT) {
        LT = arma::chol(ATA);
    }
    arma::mat _a_;
    arma::mat _a_t_a_ = arma::mat::fixed<N, N>();
    arma::mat _l_t_   = arma::mat::fixed<N, N>();
};
