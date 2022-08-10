#include "Utils.h"
#include "sunmatrix_sparse_eigen.h"
#include "nvector_serial_eigen.h"
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* ----------------------------------------------------------------------------
 * Function to create a new sparse matrix
 */

SUNMatrix SUNSparseEigenMatrix(sunindextype M, sunindextype N)
{
	SUNMatrix A;
	SparseMat* content;
	sunindextype j;

	/* return with NULL matrix on illegal dimension input */
	if ( (M <= 0) || (N <= 0) ) return(NULL);

	/* Create an empty matrix object */
	A = NULL;
	A = SUNMatNewEmpty();
	if (A == NULL) return(NULL);

	/* Attach operations */
	A->ops->getid     = SUNMatGetID_SparseEigen;
	A->ops->clone     = SUNMatClone_SparseEigen;
	A->ops->destroy   = SUNMatDestroy_SparseEigen;
	A->ops->zero      = SUNMatZero_SparseEigen;
	A->ops->copy      = SUNMatCopy_SparseEigen;
	A->ops->scaleadd  = SUNMatScaleAdd_SparseEigen;
	A->ops->scaleaddi = SUNMatScaleAddI_SparseEigen;
	A->ops->matvec    = SUNMatMatvec_SparseEigen;
	A->ops->space     = SUNMatSpace_SparseEigen;

	/* Create content */
	content = NULL;
	content = new SparseMat(M, N);
	if (content == NULL) { SUNMatDestroy(A); return(NULL); }

	/* Attach content */
	A->content = content;

	return(A);
}

/* ----------------------------------------------------------------------------
 * Function to print the Sparse matrix 
 */
 
void SUNSparseMatrix_Print(SUNMatrix A, FILE* outfile)
{
  sunindextype i, j;
  
  /* should not be called unless A is a Sparse matrix; 
     otherwise return immediately */
  if (SUNMatGetID(A) != SUNMATRIX_EIGEN_SPARSE)
    return;

  /* perform operation */
  fprintf(outfile,"\n");
  for (i=0; i<EIGMAT(A).rows(); i++) {
    for (j=0; j<EIGMAT(A).cols(); j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(outfile,"%12Lg  ", EIGMAT(A)(i, j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(outfile,"%12g  ", EIGMAT(A).coeffRef(i,j));
#else
      fprintf(outfile,"%12g  ", EIGMAT(A)(i, j));
#endif
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  return;
  SparseMat m;
  
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_SparseEigen(SUNMatrix A)
{
	return SUNMATRIX_EIGEN_SPARSE;
}

SUNMatrix SUNMatClone_SparseEigen(SUNMatrix A)
{
	SUNMatrix B = SUNSparseEigenMatrix(EIGMAT(A).rows(), EIGMAT(A).cols());
	return(B);
}

void SUNMatDestroy_SparseEigen(SUNMatrix A)
{
	if (A == NULL) return;

	/* free content */
	if (A->content != NULL) {
		delete (SparseMat*)(A->content);
		A->content = NULL;
	}

	/* free ops and matrix */
	if (A->ops) { free(A->ops); A->ops = NULL; }
	free(A); A = NULL;

	return;
}

int SUNMatZero_SparseEigen(SUNMatrix A)
{
	EIGMAT(A).setZero();
	return SUNMAT_SUCCESS;
}

int SUNMatCopy_SparseEigen(SUNMatrix A, SUNMatrix B)
{
	EIGMAT(B) = EIGMAT(A);
	return SUNMAT_SUCCESS;
}

int SUNMatScaleAddI_SparseEigen(realtype c, SUNMatrix A)
{
	EIGMAT(A) *= c;
	for (int i = 0; i < EIGMAT(A).cols(); i++) {
		SM_ELEMENT_D(A, i, i) += ONE;
	}
	return SUNMAT_SUCCESS;
}

int SUNMatScaleAdd_SparseEigen(realtype c, SUNMatrix A, SUNMatrix B)
{
	EIGMAT(A) = c * EIGMAT(A) + EIGMAT(B);
	return SUNMAT_SUCCESS;
}

int SUNMatMatvec_SparseEigen(SUNMatrix A, N_Vector x, N_Vector y)
{
	EIGV(y) = EIGMAT(A) * EIGV(x);
	return SUNMAT_SUCCESS;
}

int SUNMatSpace_SparseEigen(SUNMatrix A, long int *lenrw, long int *leniw)
{
	*lenrw = 0;
	*leniw = 0;
	return SUNMAT_SUCCESS;
}
