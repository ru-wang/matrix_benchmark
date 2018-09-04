#pragma once

#include "../common_config.hpp"

#ifndef GEMM_EIGEN_SMAP
#define GEMM_EIGEN_SMAP(N)                                              \
NONIUS_BENCHMARK(STRING(SMAP_GEMM_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, N, N> A =                                     \
            Eigen::Matrix<double, N, N>::Random();                      \
    gemm_eigen_map<N> bm(A);                                            \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_.data(), bm._a_t_a_.data()); });     \
    std::cerr << bm._a_t_a_(0, 0);                                      \
})
#endif

#ifndef GEMM_EIGEN_DMAP
#define GEMM_EIGEN_DMAP(N)                                              \
NONIUS_BENCHMARK(STRING(DMAP_GEMM_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, -1, -1> A =                                   \
            Eigen::Matrix<double, -1, -1>::Random(N, N);                \
    gemm_eigen_map<-1> bm(A);                                           \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_.data(), bm._a_t_a_.data(), N); });  \
    std::cerr << bm._a_t_a_(0, 0);                                      \
})
#endif

template<int N>
struct gemm_eigen_map {
    gemm_eigen_map(const Eigen::Matrix<double, N, N>& A) : _a_(A) { }
    void impl(const double* A, double* ATA) {
        Eigen::Map<const Eigen::Matrix<double, N, N>> A_wrap(A);
        Eigen::Map<Eigen::Matrix<double, N, N>> ATA_wrap(ATA);
        ATA_wrap.noalias() = A_wrap.transpose() * A_wrap;
    }
    Eigen::Matrix<double, N, N> _a_;
    Eigen::Matrix<double, N, N> _a_t_a_;
};

template<>
struct gemm_eigen_map<-1> {
    gemm_eigen_map(const Eigen::Matrix<double, -1, -1>& A)
            : _a_(A), _a_t_a_(A.rows(), A.cols()) { }
    void impl(const double* A, double* ATA, int dims) {
        Eigen::Map<const Eigen::Matrix<double, -1, -1>> A_wrap(A, dims, dims);
        Eigen::Map<Eigen::Matrix<double, -1, -1>> ATA_wrap(ATA, dims, dims);
        ATA_wrap.noalias() = A_wrap.transpose() * A_wrap;
    }
    Eigen::Matrix<double, -1, -1> _a_;
    Eigen::Matrix<double, -1, -1> _a_t_a_;
};
