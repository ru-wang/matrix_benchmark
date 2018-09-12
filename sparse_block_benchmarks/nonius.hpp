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
    cholesky_solver(const Eigen::MatrixXd& A)
            : ldlt(A) { run(); }
    void run() {
        using namespace Eigen;
        const int n = ldlt.rows();
        if (n <= 0) return;
        ldlt.col(0).tail(n-1) /= ldlt(0, 0);
        VectorXd temp(n-1);
        for (int k = 1; k < n; ++k) {
            const int rs = n - k - 1;
            Block<MatrixXd, Dynamic, 1> A21(ldlt, k+1, k, rs, 1);
            Block<MatrixXd, 1, Dynamic> A10(ldlt, k, 0, 1, k);
            Block<MatrixXd, Dynamic, Dynamic> A20(ldlt, k+1, 0, rs, k);
            temp.head(k).noalias() = ldlt.diagonal().head(k).asDiagonal() * A10.transpose();
            ldlt(k, k) -= (A10 * temp.head(k)).value();
            if (rs > 0) {
                A21.noalias() -= A20.lazyProduct(temp.head(k));
                A21 /= ldlt(k, k);
            }
        }
    }
    void solve_inplace(Eigen::VectorXd& b)
    {
        const int n = ldlt.rows();
        for (int j = 1; j < n; ++j) {
            b.tail(n-j).noalias() -= b[j-1] * ldlt.col(j-1).tail(n-j);
            b[j-1] /= ldlt(j-1, j-1);
        }
        b[n-1] /= ldlt(n-1, n-1);
        for (int j = n - 2; j >= 0; --j)
            b.head(j+1).noalias() -= b[j+1] *
                    ldlt.row(j+1).head(j+1).transpose();
    }
    Eigen::MatrixXd ldlt;
};
