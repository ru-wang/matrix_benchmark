#pragma once

#include "../common_config.hpp"

#ifndef GEMM_EIGEN_SMAT
#define GEMM_EIGEN_SMAT(N)                                              \
NONIUS_BENCHMARK(STRING(SMAT_GEMM_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, N, N> A =                                     \
            Eigen::Matrix<double, N, N>::Random();                      \
    gemm_eigen_mat<N> bm(A);                                            \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_, bm._a_t_a_); });                   \
    std::cerr << bm._a_t_a_(0, 0);                                      \
})
#endif

#ifndef GEMM_EIGEN_DMAT
#define GEMM_EIGEN_DMAT(N)                                              \
NONIUS_BENCHMARK(STRING(DMAT_GEMM_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, -1, -1> A =                                   \
            Eigen::Matrix<double, -1, -1>::Random(N, N);                \
    gemm_eigen_mat<-1> bm(A);                                           \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_, bm._a_t_a_); });                   \
    std::cerr << bm._a_t_a_(0, 0);                                      \
})
#endif

template<int N>
struct gemm_eigen_mat {
    gemm_eigen_mat(const Eigen::Matrix<double, N, N>& A) : _a_(A) { }
    void impl(
            const Eigen::Matrix<double, N, N>& A,
                  Eigen::Matrix<double, N, N>& ATA) {
        ATA.noalias() = A.transpose() * A;
    }
    Eigen::Matrix<double, N, N> _a_;
    Eigen::Matrix<double, N, N> _a_t_a_;
};

template<>
struct gemm_eigen_mat<-1> {
    gemm_eigen_mat(const Eigen::Matrix<double, -1, -1>& A)
            : _a_(A), _a_t_a_(A.rows(), A.cols()) { }
    void impl(
            const Eigen::Matrix<double, -1, -1>& A,
                  Eigen::Matrix<double, -1, -1>& ATA) {
        ATA.noalias() = A.transpose() * A;
    }
    Eigen::Matrix<double, -1, -1> _a_;
    Eigen::Matrix<double, -1, -1> _a_t_a_;
};
