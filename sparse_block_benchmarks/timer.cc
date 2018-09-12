#include <algorithm>
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifndef EIGEN_DONT_PARALLELIZE
#   define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Eigen>

#define TIMER_START(name) timer::start(name)
#define TIMER_STOP(name)  timer::stop (name)
#define TIMER_RUN(expr)   expr
#define TIMER_REPORT      timer::report()

#define TIMER_FUNC_START TIMER_START("[func] " + std::string(__FUNCTION__) + "()")
#define TIMER_FUNC_STOP  TIMER_STOP ("[func] " + std::string(__FUNCTION__) + "()")

#define TIMER_EXPR(expr)                         \
    TIMER_START("[expr] " + std::string(#expr)); \
    TIMER_RUN  (expr);                           \
    TIMER_STOP ("[expr] " + std::string(#expr))

#ifndef DEFINE_SPARSE_BLOCK_MATRIX_ENTRY_ACCESSOR
#   define DEFINE_SPARSE_BLOCK_MATRIX_ENTRY_ACCESSOR                                                    \
        const Eigen::MatrixXd& at(int row, int col)         const { return blocks.at(col*b_rows+row); } \
        const Eigen::MatrixXd& at(int id)                   const { return blocks.at(id); }             \
        const Eigen::MatrixXd& operator[](int id)           const { return blocks[id]; }                \
        const Eigen::MatrixXd& operator()(int row, int col) const { return blocks[col*b_rows+row]; }    \
              Eigen::MatrixXd& at(int row, int col)               { return blocks.at(col*b_rows+row); } \
              Eigen::MatrixXd& at(int id)                         { return blocks.at(id); }             \
              Eigen::MatrixXd& operator[](int id)                 { return blocks[id]; }                \
              Eigen::MatrixXd& operator()(int row, int col)       { return blocks[col*b_rows+row]; }
#endif

#ifndef DEFINE_SPARSE_BLOCK_VECTOR_ENTRY_ACCESSOR
#   define DEFINE_SPARSE_BLOCK_VECTOR_ENTRY_ACCESSOR                              \
        const Eigen::VectorXd& at(int id)         const { return blocks.at(id); } \
        const Eigen::VectorXd& operator[](int id) const { return blocks[id]; }    \
        const Eigen::VectorXd& operator()(int id) const { return blocks[id]; }    \
              Eigen::VectorXd& at(int id)               { return blocks.at(id); } \
              Eigen::VectorXd& operator[](int id)       { return blocks[id]; }    \
              Eigen::VectorXd& operator()(int id)       { return blocks[id]; }
#endif

class timer {
public:
    using clk = std::chrono::steady_clock;
    using dur = std::chrono::duration<double, std::milli>;
    using tp  = std::chrono::time_point<clk>;

    static void start(const std::string& name)
    { start_timestamps().emplace_back(name, clk::now()); }
    static void stop(const std::string& name)
    { stop_timestamps().emplace_back(name, clk::now()); }

    static std::string report() {
        std::map<std::string,
                std::pair<double,
                        std::vector<std::reference_wrapper<const tp>>>
        > timers;

        auto& start_ts = start_timestamps();
        auto& stop_ts = stop_timestamps();
        for (auto it = start_ts.begin(); it != start_ts.end(); ++it) {
            timers[it->first].first = 0;
            timers[it->first].second.emplace_back(it->second);
        }
        for (auto it = stop_ts.begin(); it != stop_ts.end(); ++it) {
            timers[it->first].first +=
                dur(it->second - timers[it->first].second.back().get()).count();
            timers[it->first].second.pop_back();
        }

        int max_len = 0;
        std::vector<std::pair<std::string, double>> timer_list;
        for (auto& entry : timers) {
            timer_list.emplace_back(entry.first, entry.second.first);
            max_len = std::max(max_len, (int)entry.first.length());
        }

        auto comp = [](
                const std::pair<std::string, double>& a,
                const std::pair<std::string, double>& b) {
            return a.second > b.second;
        };
        std::sort(timer_list.begin(), timer_list.end(), comp);

        std::stringstream ss;
        ss << std::fixed << std::setprecision(4) << "[timer report]\n";
        for (auto& entry : timer_list)
            ss << ">| " << std::setw(max_len) << std::left << entry.first
               << " : " << entry.second << "ms\n";
        start_timestamps().clear();
        stop_timestamps().clear();
        return ss.str();
    }

private:
    static std::vector<std::pair<std::string, tp>>& start_timestamps() {
        static std::vector<std::pair<std::string, tp>> timesteamps;
        return timesteamps;
    }
    static std::vector<std::pair<std::string, tp>>& stop_timestamps() {
        static std::vector<std::pair<std::string, tp>> timesteamps;
        return timesteamps;
    }
};

template<int b_rows> struct sparse_symmetric_block_matrix {
    static_assert(b_rows > 0, "block size must positive");
    DEFINE_SPARSE_BLOCK_MATRIX_ENTRY_ACCESSOR
    bool nonzero(int id)           const { return at(id).size(); }
    bool nonzero(int row, int col) const { return at(row, col).size(); }
    std::array<Eigen::MatrixXd, b_rows*b_rows> blocks;
};

template<int b_rows> struct sparse_symmetric_block_vector {
    static_assert(b_rows > 0, "block size must positive");
    DEFINE_SPARSE_BLOCK_VECTOR_ENTRY_ACCESSOR
    bool nonzero(int id) const { return at(id).size(); }
    std::array<Eigen::VectorXd, b_rows> blocks;
};

template<int N> class cholesky_ldlt {
public:
    explicit cholesky_ldlt(
            const sparse_symmetric_block_matrix<N>& ssbm)
            : _ssbm_(ssbm), _ldlt_(ssbm) { run(); }

     void pseudo_inverse_cholesky_inplace(
            Eigen::MatrixXd& A) {
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
        TIMER_FUNC_START;
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
        TIMER_FUNC_STOP;
    }

private:
     void run() {
        TIMER_FUNC_START;
        std::array<Eigen::MatrixXd, N> col_vector;
        for (int k = 0; k < N; ++k) {
            auto& diag = _ldlt_(k, k);
            pseudo_inverse_cholesky_inplace(diag);
            for (int i = k + 1; i < N; ++i) {
                col_vector[i].noalias() = _ldlt_(i, k);
                _ldlt_(i, k).applyOnTheRight(diag);
                for (int j = k + 1; j <= i; ++j)
                    _ldlt_(i, j).noalias() -=
                            _ldlt_(i, k).lazyProduct(col_vector[j].transpose());
            }
        }
        TIMER_FUNC_STOP;
    }

    const sparse_symmetric_block_matrix<N>& _ssbm_;
    sparse_symmetric_block_matrix<N> _ldlt_;
};

template<int N> auto test_case() {
    constexpr int c_rows = (1 + N) * N / 2;
    const Eigen::MatrixXd J = Eigen::MatrixXd::Random(c_rows, c_rows);
    const Eigen::MatrixXd A = J.transpose() * J + Eigen::MatrixXd::Identity(c_rows, c_rows);
    const Eigen::VectorXd b = Eigen::VectorXd::Random(c_rows);

    sparse_symmetric_block_matrix<N> S;
    sparse_symmetric_block_vector<N> t;
    for (int i = 0, istart = 0; i < N; ++i, istart += i) {
        for (int j = 0, jstart = 0; j < N; ++j, jstart += j)
            S(i, j) = A.block(istart, jstart, i+1, j+1);
        t[i] = b.segment(istart, i+1);
    }

    cholesky_ldlt<N>(S).solve_inplace(t);

    Eigen::VectorXd d1(c_rows);
    for (int i = 0, istart = 0; i < N; ++i, istart += i)
        d1.segment(istart, i+1) = t[i];

    return d1;
}

int main() {
    constexpr int K = 1000;
    for (int k = 0; k < K; ++k) {
        auto r1 = test_case<20>();
        std::cerr << r1.sum() << std::endl;
    }
    std::cout << TIMER_REPORT << std::endl;
    return 0;
}
