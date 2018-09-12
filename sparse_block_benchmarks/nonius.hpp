#pragma once

#include <Eigen/Eigen>

#define NONIUS_RUNNER
#include <nonius/nonius_single.h++>

#ifndef CHOL_ALT
#define CHOL_ALT(N)                                                            \
NONIUS_BENCHMARK(STRING(CHOL_ALT_##N), [](nonius::chronometer meter) {         \
    using namespace Eigen;                                                     \
    srand(0);                                                                  \
    constexpr int c_rows = (1 + N) * N / 2;                                    \
    const MatrixXd J = MatrixXd::Random(c_rows, c_rows);                       \
    const MatrixXd A = J.transpose() * J + MatrixXd::Identity(c_rows, c_rows); \
    const VectorXd b = VectorXd::Random(c_rows);                               \
    VectorXd d = b;                                                            \
    meter.measure([&A, &d] { cholesky_solver(A).solve_inplace(d); });          \
    std::cerr << d[0];                                                         \
/*    std::cerr << d.transpose() << std::endl;                                   */\
})
#endif

#ifndef CHOL_MAT
#define CHOL_MAT(N)                                                            \
NONIUS_BENCHMARK(STRING(CHOL_MAT_##N), [](nonius::chronometer meter) {         \
    using namespace Eigen;                                                     \
    srand(0);                                                                  \
    constexpr int c_rows = (1 + N) * N / 2;                                    \
    const MatrixXd J = MatrixXd::Random(c_rows, c_rows);                       \
    const MatrixXd A = J.transpose() * J + MatrixXd::Identity(c_rows, c_rows); \
    const VectorXd b = VectorXd::Random(c_rows);                               \
    VectorXd d = b;                                                            \
    meter.measure([&A, &d] { A.ldlt().solveInPlace(d); });                     \
    std::cerr << d[0];                                                         \
/*    std::cerr << d.transpose() << std::endl;                                   */\
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
    for (int i = 0, istart = 0; i < N; ++i, istart += i) {                     \
        for (int j = 0, jstart = 0; j < N; ++j, jstart += j)                   \
            S(i, j) = A.block(istart, jstart, i+1, j+1);                       \
        t[i] = b.segment(istart, i+1);                                         \
    }                                                                          \
    meter.measure([&S, &t] { cholesky_ldlt<N>(S).solve_inplace(t); });         \
    std::cerr << t[0][0];                                                      \
/*    VectorXd d(c_rows);                                                        */\
/*    for (int i = 0, istart = 0; i < N; ++i, istart += i)                       */\
/*        d.segment(istart, i+1) = t[i];                                         */\
/*    std::cerr << d.transpose() << std::endl;                                   */\
})
#endif

struct cholesky_solver {
    cholesky_solver(const Eigen::MatrixXd& A) : ldlt(A) {
        const int n = ldlt.rows();
        for (int k = 0; k < n; ++k) {
            const int rs_size = n - k - 1;
            ldlt(k, k) = 1.0 / ldlt(k, k);
            ldlt.row(k).tail(rs_size).noalias() =
                    ldlt.col(k).tail(rs_size).transpose();
            ldlt.col(k).tail(rs_size) *= ldlt(k, k);
            for (int j = k + 1; j < n; ++j)
                ldlt.col(j).tail(n-j).noalias() -=
                        ldlt(k, j) * ldlt.col(k).tail(n-j);
        }
    }
    void solve_inplace(Eigen::VectorXd& b) {
        const int n = ldlt.rows();
        for (int j = 1; j < n; ++j) {
            b.tail(n-j).noalias() -= b[j-1] * ldlt.col(j-1).tail(n-j);
            b[j-1] *= ldlt(j-1, j-1);
        }
        b[n-1] *= ldlt(n-1, n-1);
        for (int j = n - 2; j >= 0; --j)
            b.head(j+1).noalias() -= b[j+1] *
                    ldlt.row(j+1).head(j+1).transpose();
    }
    Eigen::MatrixXd ldlt;
};
