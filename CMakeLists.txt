project(MatrixBenchmark)
cmake_minimum_required(VERSION 3.5.1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always -Wfatal-errors -Wall -Wno-sign-compare")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -march=native")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

include_directories(.)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/nonius/include)

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

find_package(BLAS REQUIRED)
include_directories(SYSTEM ${BLAS_INCLUDE_DIR})

find_package(LAPACK REQUIRED)
include_directories(SYSTEM ${LAPACK_INCLUDE_DIR})

set(LIBRARIES ${ARMADILLO_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} pthread)

add_executable(gemm_arma_mat_vs_eigen_smat benchmarks/gemm_arma_mat_vs_eigen_smat.cc)
add_executable(gemm_arma_mat_vs_eigen_smap benchmarks/gemm_arma_mat_vs_eigen_smap.cc)
add_executable(gemm_arma_mat_vs_eigen_dmat benchmarks/gemm_arma_mat_vs_eigen_dmat.cc)
add_executable(gemm_arma_mat_vs_eigen_dmap benchmarks/gemm_arma_mat_vs_eigen_dmap.cc)
add_executable(gemm_arma_vs_eigen_small    benchmarks/gemm_arma_vs_eigen_small.cc)
add_executable(gemm_arma_vs_eigen_big      benchmarks/gemm_arma_vs_eigen_big.cc)

add_executable(chol_arma_mat_vs_eigen_smat benchmarks/chol_arma_mat_vs_eigen_smat.cc)
add_executable(chol_arma_mat_vs_eigen_smap benchmarks/chol_arma_mat_vs_eigen_smap.cc)
add_executable(chol_arma_mat_vs_eigen_dmat benchmarks/chol_arma_mat_vs_eigen_dmat.cc)
add_executable(chol_arma_mat_vs_eigen_dmap benchmarks/chol_arma_mat_vs_eigen_dmap.cc)

target_link_libraries(gemm_arma_mat_vs_eigen_smat ${LIBRARIES})
target_link_libraries(gemm_arma_mat_vs_eigen_smap ${LIBRARIES})
target_link_libraries(gemm_arma_mat_vs_eigen_dmat ${LIBRARIES})
target_link_libraries(gemm_arma_mat_vs_eigen_dmap ${LIBRARIES})
target_link_libraries(gemm_arma_vs_eigen_small    ${LIBRARIES})
target_link_libraries(gemm_arma_vs_eigen_big      ${LIBRARIES})

target_link_libraries(chol_arma_mat_vs_eigen_smat ${LIBRARIES})
target_link_libraries(chol_arma_mat_vs_eigen_smap ${LIBRARIES})
target_link_libraries(chol_arma_mat_vs_eigen_dmat ${LIBRARIES})
target_link_libraries(chol_arma_mat_vs_eigen_dmap ${LIBRARIES})

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
