#pragma once

#include "./common_config.hpp"

BENCH_STRONG_INLINE Eigen::MatrixXd
pseudo_inverse_cholesky(
        const Eigen::MatrixXd& A, double epsilon) {
    const int n = A.rows();
    const auto& ldlt = A.ldlt();
    const auto& P = ldlt.transpositionsP();
    const auto& D = ldlt.vectorD();
    const auto& L = ldlt.matrixL();
    Eigen::MatrixXd A_inv = Eigen::MatrixXd::Identity(n, n);
    for (int j = 0; j < n - 1; ++j)
        for (int i = j + 1; i < n; ++i)
            A_inv.block(i, 0, 1, j+1) -=
                    L(i, j) * A_inv.block(j, 0, 1, j+1);
    A_inv = P.transpose() * (A_inv.transpose() *
            (D.array()>epsilon).select(D.cwiseInverse(), 0.0).asDiagonal() *
                    A_inv) * P.transpose();
    return A_inv;
}

template<int N> class cholesky_ldlt {
public:
    explicit cholesky_ldlt(
            const sparse_symmetric_block_matrix<N>& ssbm)
            : _ssbm_(ssbm), _ldlt_(ssbm) { run(); }

    sparse_symmetric_block_matrix<N>&
    matrix_ldlt() { return _ldlt_; }

    const sparse_symmetric_block_matrix<N>&
    matrix_ldlt() const { return _ldlt_; }

    BENCH_STRONG_INLINE sparse_symmetric_block_vector<N>
    solve(const sparse_symmetric_block_vector<N>& ssbv) {
        sparse_symmetric_block_vector<N> ret = ssbv;
        for (int j = 1; j < N; ++j)
            for (int i = j; i < N; ++i)
                ret[i].noalias() -=
                        _ldlt_.at(i, j-1).lazyProduct(ret[j-1]);
        for (int i = 0; i < N; ++i)
            ret[i].applyOnTheLeft(_ldlt_.at(i, i));
        for (int j = N - 2; j >= 0; --j)
            for (int i = 0; i <= j; ++i)
                ret[i].noalias() -= _ldlt_.at(j+1, i).transpose()
                        .lazyProduct(ret[j+1]);
        return ret;
    }

    BENCH_STRONG_INLINE void
    solve_inplace(sparse_symmetric_block_vector<N>& ssbv) {
        for (int j = 1; j < N; ++j)
            for (int i = j; i < N; ++i)
                ssbv[i].noalias() -=
                        _ldlt_.at(i, j-1).lazyProduct(ssbv[j-1]);
        for (int i = 0; i < N; ++i)
            ssbv[i].applyOnTheLeft(_ldlt_.at(i, i));
        for (int j = N - 2; j >= 0; --j)
            for (int i = 0; i <= j; ++i)
                ssbv[i].noalias() -= _ldlt_.at(j+1, i).transpose()
                        .lazyProduct(ssbv[j+1]);
    }

private:
    BENCH_STRONG_INLINE void run() {
        std::array<Eigen::MatrixXd, N> col_vector;
        for (int k = 0; k < N; ++k) {
            auto& diag = _ldlt_(k, k);
            diag = pseudo_inverse_cholesky(diag, EPSILON);
            for (int i = k + 1; i < N; ++i) {
                col_vector[i].noalias() = _ldlt_(i, k).transpose();
                _ldlt_(i, k).applyOnTheRight(diag);
            }
            for (int i = k + 1; i < N; ++i)
                for (int j = k + 1; j <= i; ++j)
                    _ldlt_(i, j).noalias() -=
                            _ldlt_(i, k).lazyProduct(col_vector[j]);
        }
    }

    static constexpr double EPSILON = 1.0e-10;

    const sparse_symmetric_block_matrix<N>& _ssbm_;
    sparse_symmetric_block_matrix<N> _ldlt_;
};

#ifndef CHOL_MAT
#define CHOL_MAT(N)                                                            \
NONIUS_BENCHMARK(STRING(CHOL_MAT_##N), [](nonius::chronometer meter) {         \
    using namespace Eigen;                                                     \
    srand(0);                                                                  \
    constexpr int c_rows = (1 + N) * N / 2;                                    \
    const MatrixXd J = MatrixXd::Random(c_rows, c_rows);                       \
    const MatrixXd A = J.transpose() * J + MatrixXd::Identity(c_rows, c_rows); \
    const VectorXd b = VectorXd::Random(c_rows);                               \
    VectorXd d(c_rows);                                                        \
    meter.measure([&A, &b, &d] { d = A.ldlt().solve(b); });                    \
    std::cerr << d[0];                                                         \
})
#endif

#ifndef CHOL_SBM
#define CHOL_SBM(N)                                                            \
NONIUS_BENCHMARK(STRING(CHOL_SBM_##N), [](nonius::chronometer meter) {         \
    using namespace Eigen;                                                     \
    srand(0);                                                                  \
    constexpr int c_rows = (1 + N) * N / 2;                                    \
    const MatrixXd J = MatrixXd::Random(c_rows, c_rows);                       \
    const MatrixXd A = J.transpose() * J + MatrixXd::Identity(c_rows, c_rows); \
    const VectorXd b = VectorXd::Random(c_rows);                               \
    sparse_symmetric_block_matrix<N> S;                                        \
    sparse_symmetric_block_vector<N> t;                                        \
    sparse_symmetric_block_vector<N> r;                                        \
    for (int i = 0, istart = 0; i < N; ++i, istart += i) {                     \
        for (int j = 0, jstart = 0; j < N; ++j, jstart += j)                   \
            S(i, j) = A.block(istart, jstart, i+1, j+1);                       \
        t[i] = b.segment(istart, i+1);                                         \
    }                                                                          \
    meter.measure([&S, &t, &r] { r = cholesky_ldlt<N>(S).solve(t); });         \
    std::cerr << r[0][0];                                                      \
})
#endif
