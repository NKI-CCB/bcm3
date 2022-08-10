#pragma once

#include <sundials/sundials_matrix.h>
#include <Eigen/Sparse>

/* ------------------------------------
 * Macros for access to SUNMATRIX_DENSE
 * ------------------------------------ */

typedef Eigen::SparseMatrix<double> SparseMat;
#define EIGMAT(A) (*(SparseMat*)(A->content))
#define SM_ELEMENT_D(A, i, j) (((SparseMat*)(A->content))->coeffRef(i,j))

/* ---------------------------------------
 * Exported Functions for SUNMATRIX_DENSE
 * --------------------------------------- */

SUNMatrix SUNSparseEigenMatrix(sunindextype M, sunindextype N);

void SUNSparseEigenMatrix_Print(SUNMatrix A, FILE* outfile);

SUNMatrix_ID SUNMatGetID_SparseEigen(SUNMatrix A);
SUNMatrix SUNMatClone_SparseEigen(SUNMatrix A);
void SUNMatDestroy_SparseEigen(SUNMatrix A);
int SUNMatZero_SparseEigen(SUNMatrix A);
int SUNMatCopy_SparseEigen(SUNMatrix A, SUNMatrix B);
int SUNMatScaleAdd_SparseEigen(realtype c, SUNMatrix A, SUNMatrix B);
int SUNMatScaleAddI_SparseEigen(realtype c, SUNMatrix A);
int SUNMatMatvec_SparseEigen(SUNMatrix A, N_Vector x, N_Vector y);
int SUNMatSpace_SparseEigen(SUNMatrix A, long int *lenrw, long int *leniw);
