#pragma once

#include "./common_config.hpp"

template<int N> class cholesky_ldlt {
public:
    explicit cholesky_ldlt(
            const sparse_symmetric_block_matrix<N>& ssbm)
            : _ssbm_(ssbm), _ldlt_(ssbm) { run(); }

    sparse_symmetric_block_matrix<N>&
    matrix_ldlt() { return _ldlt_; }

    const sparse_symmetric_block_matrix<N>&
    matrix_ldlt() const { return _ldlt_; }

    void pseudo_inverse_cholesky_inplace(Eigen::MatrixXd& A) {
        const int n = A.rows();
        const auto& ldlt = A.ldlt();
        const auto& P = ldlt.transpositionsP();
        const auto& D = ldlt.vectorD();
        const auto& L = ldlt.matrixL();
        A.setIdentity();
        for (int j = 0; j < n - 1; ++j)
            for (int i = j + 1; i < n; ++i)
                A.block(i, 0, 1, j+1) -=
                    L(i, j) * A.block(j, 0, 1, j+1);
        A = P.transpose() *
                (A.transpose() * D.cwiseInverse().asDiagonal() * A)
                        * P.transpose();
    }

    sparse_symmetric_block_vector<N> solve(
            const sparse_symmetric_block_vector<N>& ssbv) {
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

    void solve_inplace(
            sparse_symmetric_block_vector<N>& ssbv) {
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
    void run() {
        std::array<Eigen::MatrixXd, N> col_vector;
        for (int k = 0; k < N; ++k) {
            auto& diag = _ldlt_(k, k);
            pseudo_inverse_cholesky_inplace(diag);
            for (int i = k + 1; i < N; ++i) {
                col_vector[i].noalias() = _ldlt_(i, k).transpose();
                _ldlt_(i, k).applyOnTheRight(diag);
                for (int j = k + 1; j <= i; ++j)
                    _ldlt_(i, j).noalias() -=
                            _ldlt_(i, k).lazyProduct(col_vector[j]);
            }
        }
    }

    const sparse_symmetric_block_matrix<N>& _ssbm_;
    sparse_symmetric_block_matrix<N> _ldlt_;
};
