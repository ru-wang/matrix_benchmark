#pragma once

#include "../common_config.hpp"

#ifndef CHOL_EIGEN_SMAT
#define CHOL_EIGEN_SMAT(N)                                              \
NONIUS_BENCHMARK(STRING(SMAT_CHOL_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, N, N> A =                                     \
            Eigen::Matrix<double, N, N>::Random();                      \
    chol_eigen_mat<N> bm(A);                                            \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_t_a_, bm._l_t_); });                 \
    std::cerr << bm._l_t_(0, 0);                                        \
})
#endif

#ifndef CHOL_EIGEN_DMAT
#define CHOL_EIGEN_DMAT(N)                                              \
NONIUS_BENCHMARK(STRING(DMAT_CHOL_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, -1, -1> A =                                   \
            Eigen::Matrix<double, -1, -1>::Random(N, N);                \
    chol_eigen_mat<-1> bm(A);                                           \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_t_a_, bm._l_t_); });                 \
    std::cerr << bm._l_t_(0, 0);                                        \
})
#endif

template<int N>
struct chol_eigen_mat {
    chol_eigen_mat(const Eigen::Matrix<double, N, N>& A) : _a_(A) {
        _a_t_a_.setIdentity();
        _a_t_a_.noalias() += _a_.transpose().lazyProduct(_a_);
    }
    void impl(
            const Eigen::Matrix<double, N, N>& ATA,
                  Eigen::Matrix<double, N, N>& LT) {
        LT = ATA.llt().matrixU();
    }
    Eigen::Matrix<double, N, N> _a_;
    Eigen::Matrix<double, N, N> _a_t_a_;
    Eigen::Matrix<double, N, N> _l_t_;
};

template<>
struct chol_eigen_mat<-1> {
    chol_eigen_mat(const Eigen::Matrix<double, -1, -1>& A)
            : _a_(A), _l_t_(A.rows(), A.cols()) {
        _a_t_a_.setIdentity(A.rows(), A.cols());
        _a_t_a_.noalias() += _a_.transpose().lazyProduct(_a_);
    }
    void impl(
            const Eigen::Matrix<double, -1, -1>& ATA,
                  Eigen::Matrix<double, -1, -1>& LT) {
        LT = ATA.llt().matrixU();
    }
    Eigen::Matrix<double, -1, -1> _a_;
    Eigen::Matrix<double, -1, -1> _a_t_a_;
    Eigen::Matrix<double, -1, -1> _l_t_;
};
