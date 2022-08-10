#pragma once

#include <sundials/sundials_matrix.h>

/* ------------------------------------
 * Macros for access to SUNMATRIX_DENSE
 * ------------------------------------ */

#define EIGMAT(A) (*(MatrixReal*)(A->content))
#define SM_ELEMENT_D(A, i, j) ((*(MatrixReal*)(A->content))(i,j))

/* ---------------------------------------
 * Exported Functions for SUNMATRIX_DENSE
 * --------------------------------------- */

SUNMatrix SUNDenseEigenMatrix(sunindextype M, sunindextype N);

void SUNDenseEigenMatrix_Print(SUNMatrix A, FILE* outfile);

SUNMatrix_ID SUNMatGetID_DenseEigen(SUNMatrix A);
SUNMatrix SUNMatClone_DenseEigen(SUNMatrix A);
void SUNMatDestroy_DenseEigen(SUNMatrix A);
int SUNMatZero_DenseEigen(SUNMatrix A);
int SUNMatCopy_DenseEigen(SUNMatrix A, SUNMatrix B);
int SUNMatScaleAdd_DenseEigen(realtype c, SUNMatrix A, SUNMatrix B);
int SUNMatScaleAddI_DenseEigen(realtype c, SUNMatrix A);
int SUNMatMatvec_DenseEigen(SUNMatrix A, N_Vector x, N_Vector y);
int SUNMatSpace_DenseEigen(SUNMatrix A, long int *lenrw, long int *leniw);
