#pragma once

#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include "EigenPartialPivLUSomewhatSparse.h"

struct _SUNLinearSolverContent_Dense_Eigen {
	PartialPivLUExtended<OdeMatrixReal> lu;
    Eigen::Matrix3f inverse;
};

typedef struct _SUNLinearSolverContent_Dense_Eigen* SUNLinearSolverContent_Dense_Eigen;

SUNLinearSolver SUNLinSol_Dense_Eigen(N_Vector y, SUNMatrix A);

SUNLinearSolver_Type SUNLinSolGetType_Dense_Eigen(SUNLinearSolver S);
SUNLinearSolver_ID SUNLinSolGetID_Dense_Eigen(SUNLinearSolver S);
int SUNLinSolInitialize_Dense_Eigen(SUNLinearSolver S);
int SUNLinSolSetup_Dense_Eigen(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSetup_Dense_Eigen2x2(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSetup_Dense_Eigen3x3(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSolve_Dense_Eigen(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
int SUNLinSolSolve_Dense_Eigen2x2(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
int SUNLinSolSolve_Dense_Eigen3x3(SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);
sunindextype SUNLinSolLastFlag_Dense_Eigen(SUNLinearSolver S);
int SUNLinSolSpace_Dense_Eigen(SUNLinearSolver S, long int *lenrwLS, long int *leniwLS);
int SUNLinSolFree_Dense_Eigen(SUNLinearSolver S);
