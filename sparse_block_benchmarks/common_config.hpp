#pragma once

#include <array>
#include <cstdlib>
#include <iostream>

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Eigen>

#include "./inline.hpp"

#ifndef STRING
#define STRING(str) #str
#endif

#ifndef DEFINE_SPARSE_BLOCK_MATRIX_ENTRY_ACCESSOR
#define DEFINE_SPARSE_BLOCK_MATRIX_ENTRY_ACCESSOR                                                   \
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
#define DEFINE_SPARSE_BLOCK_VECTOR_ENTRY_ACCESSOR                             \
    const Eigen::VectorXd& at(int id)         const { return blocks.at(id); } \
    const Eigen::VectorXd& operator[](int id) const { return blocks[id]; }    \
    const Eigen::VectorXd& operator()(int id) const { return blocks[id]; }    \
          Eigen::VectorXd& at(int id)               { return blocks.at(id); } \
          Eigen::VectorXd& operator[](int id)       { return blocks[id]; }    \
          Eigen::VectorXd& operator()(int id)       { return blocks[id]; }
#endif

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
