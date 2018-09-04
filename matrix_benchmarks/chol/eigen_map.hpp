#pragma once

#include "../common_config.hpp"

#ifndef CHOL_EIGEN_SMAP
#define CHOL_EIGEN_SMAP(N)                                              \
NONIUS_BENCHMARK(STRING(SMAP_CHOL_##N), [](nonius::chronometer meter) { \
    srand(0);                                                           \
    Eigen::Matrix<double, N, N> A =                                     \
            Eigen::Matrix<double, N, N>::Random();                      \
    chol_eigen_map<N> bm(A);                                            \
    meter.measure([&bm]                                                 \
            { return bm.impl(bm._a_t_a_.data(), bm._l_t_.data()); });   \
    std::cerr << bm._l_t_(0, 0);                                        \
})
#endif

#ifndef CHOL_EIGEN_DMAP
#define CHOL_EIGEN_DMAP(N)                                               \
NONIUS_BENCHMARK(STRING(DMAP_CHOL_##N), [](nonius::chronometer meter) {  \
    srand(0);                                                            \
    Eigen::Matrix<double, -1, -1> A =                                    \
            Eigen::Matrix<double, -1, -1>::Random(N, N);                 \
    chol_eigen_map<-1> bm(A);                                            \
    meter.measure([&bm]                                                  \
            { return bm.impl(bm._a_t_a_.data(), bm._l_t_.data(), N); }); \
    std::cerr << bm._l_t_(0, 0);                                         \
})
#endif

template<int N>
struct chol_eigen_map {
    chol_eigen_map(const Eigen::Matrix<double, N, N>& A) : _a_(A) {
        _a_t_a_.setIdentity();
        _a_t_a_.noalias() += _a_.transpose().lazyProduct(_a_);
    }
    void impl(const double* ATA, double* LT) {
        Eigen::Map<const Eigen::Matrix<double, N, N>> ATA_wrap(ATA);
        Eigen::Map<Eigen::Matrix<double, N, N>> LT_wrap(LT);
        LT_wrap = ATA_wrap.llt().matrixU();
    }
    Eigen::Matrix<double, N, N> _a_;
    Eigen::Matrix<double, N, N> _a_t_a_;
    Eigen::Matrix<double, N, N> _l_t_;
};

template<>
struct chol_eigen_map<-1> {
    chol_eigen_map(const Eigen::Matrix<double, -1, -1>& A)
            : _a_(A), _l_t_(A.rows(), A.cols()) {
        _a_t_a_.setIdentity(A.rows(), A.cols());
        _a_t_a_.noalias() += _a_.transpose().lazyProduct(_a_);
    }
    void impl(
            const double* ATA, double* LT, int dims) {
        Eigen::Map<const Eigen::Matrix<double, -1, -1>> ATA_wrap(ATA, dims, dims);
        Eigen::Map<Eigen::Matrix<double, -1, -1>> LT_wrap(LT, dims, dims);
        LT_wrap = ATA_wrap.llt().matrixU();
    }
    Eigen::Matrix<double, -1, -1> _a_;
    Eigen::Matrix<double, -1, -1> _a_t_a_;
    Eigen::Matrix<double, -1, -1> _l_t_;
};

