#pragma once

#define CVODE_USE_EIGEN_SOLVER 1
#define CVODE_USE_SPARSE_SOLVER 0

#if CVODE_USE_EIGEN_SOLVER
#include "nvector_serial_eigen.h"
inline N_Vector MakeCVodeVector(size_t N) { return N_VNew_Eigen((long)N); }
#if CVODE_USE_SPARSE_SOLVER
#include "sunmatrix_sparse_eigen.h"
#include "sunlinsol_sparse_eigen.h"
inline SUNMatrix MakeCVodeMatrix(size_t N, size_t M) { return SUNSparseEigenMatrix(M, N); }
inline SUNLinearSolver MakeCVodeLinearSolver(N_Vector y, SUNMatrix J) { return SUNLinSol_Sparse_Eigen(y, J); }
#else
#include "sunmatrix_dense_eigen.h"
#include "sunlinsol_dense_eigen.h"
inline SUNMatrix MakeCVodeMatrix(size_t N, size_t M) { return SUNDenseEigenMatrix(M, N); }
inline SUNLinearSolver MakeCVodeLinearSolver(N_Vector y, SUNMatrix J) { return SUNLinSol_Dense_Eigen(y, J); }
#endif
#else
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
inline N_Vector MakeCVodeVector(size_t N) { return N_VNew_Serial((long)N); }
inline SUNMatrix MakeCVodeMatrix(size_t N, size_t M) { return SUNDenseMatrix(M, N); }
inline SUNLinearSolver MakeCVodeLinearSolver(N_Vector y, SUNMatrix J) { return SUNLinSol_Dense(y, J); }
#endif
