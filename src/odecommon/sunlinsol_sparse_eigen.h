#pragma once

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include "sunmatrix_sparse_eigen.h"

struct _SUNLinearSolverContent_Sparse_Eigen {
	Eigen::SparseLU< SparseMat > lu;
	int analyzed_nonzero_count;
};

typedef struct _SUNLinearSolverContent_Sparse_Eigen* SUNLinearSolverContent_Sparse_Eigen;

SUNLinearSolver SUNLinSol_Sparse_Eigen(N_Vector y, SUNMatrix A);

SUNLinearSolver_Type SUNLinSolGetType_Sparse_Eigen(SUNLinearSolver S);
SUNLinearSolver_ID SUNLinSolGetID_Sparse_Eigen(SUNLinearSolver S);
int SUNLinSolInitialize_Sparse_Eigen(SUNLinearSolver S);
int SUNLinSolSetup_Sparse_Eigen(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSolve_Sparse_Eigen(SUNLinearSolver S, SUNMatrix A,
                         N_Vector x, N_Vector b, realtype tol);
sunindextype SUNLinSolLastFlag_Sparse_Eigen(SUNLinearSolver S);
int SUNLinSolSpace_Sparse_Eigen(SUNLinearSolver S,
                         long int *lenrwLS,
                         long int *leniwLS);
int SUNLinSolFree_Sparse_Eigen(SUNLinearSolver S);
