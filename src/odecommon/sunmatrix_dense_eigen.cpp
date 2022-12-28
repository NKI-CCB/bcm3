#include "Utils.h"
#include "LinearAlgebraSelector.h"

#if CVODE_USE_EIGEN_SOLVER
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)

/* ----------------------------------------------------------------------------
 * Function to create a new dense matrix
 */

SUNMatrix SUNDenseEigenMatrix(sunindextype M, sunindextype N)
{
	SUNMatrix A;
	OdeMatrixReal* content;
	sunindextype j;

	/* return with NULL matrix on illegal dimension input */
	if ( (M <= 0) || (N <= 0) ) return(NULL);

	/* Create an empty matrix object */
	A = NULL;
	A = SUNMatNewEmpty();
	if (A == NULL) return(NULL);

	/* Attach operations */
	A->ops->getid     = SUNMatGetID_DenseEigen;
	A->ops->clone     = SUNMatClone_DenseEigen;
	A->ops->destroy   = SUNMatDestroy_DenseEigen;
	A->ops->zero      = SUNMatZero_DenseEigen;
	A->ops->copy      = SUNMatCopy_DenseEigen;
	A->ops->scaleadd  = SUNMatScaleAdd_DenseEigen;
	A->ops->scaleaddi = SUNMatScaleAddI_DenseEigen;
	A->ops->matvec    = SUNMatMatvec_DenseEigen;
	A->ops->space     = SUNMatSpace_DenseEigen;

	/* Create content */
	content = NULL;
	content = new OdeMatrixReal(M, N);
	if (content == NULL) { SUNMatDestroy(A); return(NULL); }

	/* Attach content */
	A->content = content;

	return(A);
}

/* ----------------------------------------------------------------------------
 * Function to print the dense matrix 
 */
 
void SUNDenseMatrix_Print(SUNMatrix A, FILE* outfile)
{
  sunindextype i, j;
  
  /* should not be called unless A is a dense matrix; 
     otherwise return immediately */
  if (SUNMatGetID(A) != SUNMATRIX_EIGEN_DENSE)
    return;

  /* perform operation */
  fprintf(outfile,"\n");
  for (i=0; i<EIGMAT(A).rows(); i++) {
    for (j=0; j<EIGMAT(A).cols(); j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
      fprintf(outfile,"%12Lg  ", EIGMAT(A)(i, j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
      fprintf(outfile,"%12g  ", EIGMAT(A)(i,j));
#else
      fprintf(outfile,"%12g  ", EIGMAT(A)(i, j));
#endif
    }
    fprintf(outfile,"\n");
  }
  fprintf(outfile,"\n");
  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of matrix operations
 * -----------------------------------------------------------------
 */

SUNMatrix_ID SUNMatGetID_DenseEigen(SUNMatrix A)
{
	return SUNMATRIX_EIGEN_DENSE;
}

SUNMatrix SUNMatClone_DenseEigen(SUNMatrix A)
{
	SUNMatrix B = SUNDenseEigenMatrix(EIGMAT(A).rows(), EIGMAT(A).cols());
	return(B);
}

void SUNMatDestroy_DenseEigen(SUNMatrix A)
{
	if (A == NULL) return;

	/* free content */
	if (A->content != NULL) {
		delete (OdeMatrixReal*)(A->content);
		A->content = NULL;
	}

	/* free ops and matrix */
	if (A->ops) { free(A->ops); A->ops = NULL; }
	free(A); A = NULL;

	return;
}

int SUNMatZero_DenseEigen(SUNMatrix A)
{
	EIGMAT(A).setZero();
	return SUNMAT_SUCCESS;
}

int SUNMatCopy_DenseEigen(SUNMatrix A, SUNMatrix B)
{
	EIGMAT(B) = EIGMAT(A);
	return SUNMAT_SUCCESS;
}

int SUNMatScaleAddI_DenseEigen(realtype c, SUNMatrix A)
{
	EIGMAT(A) *= c;
	EIGMAT(A).diagonal().array() += ONE;
	return SUNMAT_SUCCESS;
}

int SUNMatScaleAdd_DenseEigen(realtype c, SUNMatrix A, SUNMatrix B)
{
	EIGMAT(A) = c * EIGMAT(A) + EIGMAT(B);
	return SUNMAT_SUCCESS;
}

int SUNMatMatvec_DenseEigen(SUNMatrix A, N_Vector x, N_Vector y)
{
	EIGV(y) = EIGMAT(A) * EIGV(x);
	return SUNMAT_SUCCESS;
}

int SUNMatSpace_DenseEigen(SUNMatrix A, long int *lenrw, long int *leniw)
{
	*lenrw = 0;
	*leniw = 0;
	return SUNMAT_SUCCESS;
}

#endif
