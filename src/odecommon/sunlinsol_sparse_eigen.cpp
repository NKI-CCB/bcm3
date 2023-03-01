#include "Utils.h"
#include "LinearAlgebraSelector.h"

#if CVODE_USE_EIGEN_SOLVER && CVODE_USE_SPARSE_SOLVER

#include <sundials/sundials_math.h>

#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Sparse solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

inline Eigen::SparseLU< SparseMat >& EIGSOL(SUNLinearSolver S) {
	return reinterpret_cast<SUNLinearSolverContent_Sparse_Eigen>(S->content)->lu;
}

/* ----------------------------------------------------------------------------
 * Function to create a new sparse linear solver
 */

SUNLinearSolver SUNLinSol_Sparse_Eigen(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_Sparse_Eigen content;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_EIGEN_SPARSE) return(NULL);

  if ((N_VGetVectorID(y) != SUNDIALS_NVEC_EIGEN))
    return(NULL);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_Sparse_Eigen;
  S->ops->getid      = SUNLinSolGetID_Sparse_Eigen;
  S->ops->initialize = SUNLinSolInitialize_Sparse_Eigen;
  S->ops->setup      = SUNLinSolSetup_Sparse_Eigen;
  S->ops->solve      = SUNLinSolSolve_Sparse_Eigen;
  S->ops->lastflag   = SUNLinSolLastFlag_Sparse_Eigen;
  S->ops->space      = SUNLinSolSpace_Sparse_Eigen;
  S->ops->free       = SUNLinSolFree_Sparse_Eigen;

  /* Create content */
  content = NULL;
  content = new _SUNLinearSolverContent_Sparse_Eigen;
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }
  content->analyzed_nonzero_count = 0;

  /* Attach content */
  S->content = content;

  /* Fill content */
  return(S);
}

/*
 * -----------------------------------------------------------------
 * implementation of linear solver operations
 * -----------------------------------------------------------------
 */

SUNLinearSolver_Type SUNLinSolGetType_Sparse_Eigen(SUNLinearSolver S)
{
	return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID SUNLinSolGetID_Sparse_Eigen(SUNLinearSolver S)
{
	return SUNLINEARSOLVER_CUSTOM;
}

int SUNLinSolInitialize_Sparse_Eigen(SUNLinearSolver S)
{
	/* all solver-specific memory has already been allocated */
	return(SUNLS_SUCCESS);
}

int SUNLinSolSetup_Sparse_Eigen(SUNLinearSolver S, SUNMatrix A)
{
	/* check for valid inputs */
	if ( (A == NULL) || (S == NULL) )
		return(SUNLS_MEM_NULL);

	/* Ensure that A is a Sparse eigen matrix */
	if (SUNMatGetID(A) != SUNMATRIX_EIGEN_SPARSE) {
		return(SUNLS_ILL_INPUT);
	}

	if (EIGMAT(A).nonZeros() != ((SUNLinearSolverContent_Sparse_Eigen)S->content)->analyzed_nonzero_count) {
		EIGSOL(S).analyzePattern((*(SparseMat*)(A->content)));
		((SUNLinearSolverContent_Sparse_Eigen)S->content)->analyzed_nonzero_count = EIGMAT(A).nonZeros();
	}
	EIGSOL(S).factorize(EIGMAT(A));
	if (EIGSOL(S).info() == Eigen::Success) {
		return(SUNLS_SUCCESS);
	} else {
		return(SUNLS_LUFACT_FAIL);
	}
}

int SUNLinSolSolve_Sparse_Eigen(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                         N_Vector b, realtype tol)
{
	if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) )
		return(SUNLS_MEM_NULL);

	EIGV(x) = EIGSOL(S).solve(EIGV(b));

	return(SUNLS_SUCCESS);
}

sunindextype SUNLinSolLastFlag_Sparse_Eigen(SUNLinearSolver S)
{
	/* return the stored 'last_flag' value */
	if (S == NULL) return(-1);
	return(SUNLS_SUCCESS);
}

int SUNLinSolSpace_Sparse_Eigen(SUNLinearSolver S,
                         long int *lenrwLS,
                         long int *leniwLS)
{
	*leniwLS = 0;
	*lenrwLS = 0;
	return(SUNLS_SUCCESS);
}

int SUNLinSolFree_Sparse_Eigen(SUNLinearSolver S)
{
	/* return if S is already free */
	if (S == NULL) return(SUNLS_SUCCESS);

	/* delete items from contents, then delete generic structure */
	if (S->content) {
		delete (_SUNLinearSolverContent_Sparse_Eigen*)(S->content);
		S->content = NULL;
	}
	if (S->ops) {
		free(S->ops);
		S->ops = NULL;
	}
	free(S);
	return(SUNLS_SUCCESS);
}

#endif
