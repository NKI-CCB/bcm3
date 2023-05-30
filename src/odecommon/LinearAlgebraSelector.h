#pragma once

#define ODE_SINGLE_PRECISION 0 // If this is changed, make sure it's also changed in sundials_config.h
#define CVODE_USE_EIGEN_SOLVER 0
#define CVODE_USE_SPARSE_SOLVER 0

#if ODE_SINGLE_PRECISION
typedef float OdeReal;
typedef Eigen::VectorXf OdeVectorReal;
typedef Eigen::MatrixXf OdeMatrixReal;
#else
typedef double OdeReal;
typedef Eigen::VectorXd OdeVectorReal;
typedef Eigen::MatrixXd OdeMatrixReal;
#endif

#if CVODE_USE_EIGEN_SOLVER
#include "nvector_serial_eigen.h"
#if CVODE_USE_SPARSE_SOLVER
#include "sunmatrix_sparse_eigen.h"
#include "sunlinsol_sparse_eigen.h"
#else
#include "sunmatrix_dense_eigen.h"
#include "sunlinsol_dense_eigen.h"
#endif
#else
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#endif

N_Vector MakeCVodeVector(size_t N);
SUNMatrix MakeCVodeMatrix(size_t N, size_t M);
SUNLinearSolver MakeCVodeLinearSolver(N_Vector y, SUNMatrix J);
