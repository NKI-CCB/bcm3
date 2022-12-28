#include "Utils.h"
#include "LinearAlgebraSelector.h"

#if CVODE_USE_EIGEN_SOLVER && !CVODE_USE_SPARSE_SOLVER

#include <sundials/sundials_math.h>

#define ONE  RCONST(1.0)

/*
 * -----------------------------------------------------------------
 * Dense solver structure accessibility macros:
 * -----------------------------------------------------------------
 */

inline PartialPivLUExtended< OdeMatrixReal >& EIGSOL(SUNLinearSolver S) {
	return reinterpret_cast<SUNLinearSolverContent_Dense_Eigen>(S->content)->lu;
}



/* ----------------------------------------------------------------------------
 * Function to create a new dense linear solver
 */

SUNLinearSolver SUNLinSol_Dense_Eigen(N_Vector y, SUNMatrix A)
{
  SUNLinearSolver S;
  SUNLinearSolverContent_Dense_Eigen content;

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (SUNMatGetID(A) != SUNMATRIX_EIGEN_DENSE) return(NULL);

  if ((N_VGetVectorID(y) != SUNDIALS_NVEC_EIGEN))
    return(NULL);

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty();
  if (S == NULL) return(NULL);

  /* Attach operations */
  S->ops->gettype    = SUNLinSolGetType_Dense_Eigen;
  S->ops->getid      = SUNLinSolGetID_Dense_Eigen;
  S->ops->initialize = SUNLinSolInitialize_Dense_Eigen;
  S->ops->setup      = SUNLinSolSetup_Dense_Eigen;
  S->ops->solve      = SUNLinSolSolve_Dense_Eigen;
  S->ops->lastflag   = SUNLinSolLastFlag_Dense_Eigen;
  S->ops->space      = SUNLinSolSpace_Dense_Eigen;
  S->ops->free       = SUNLinSolFree_Dense_Eigen;

  /* Create content */
  content = NULL;
  content = new _SUNLinearSolverContent_Dense_Eigen;
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

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

SUNLinearSolver_Type SUNLinSolGetType_Dense_Eigen(SUNLinearSolver S)
{
	return SUNLINEARSOLVER_DIRECT;
}

SUNLinearSolver_ID SUNLinSolGetID_Dense_Eigen(SUNLinearSolver S)
{
	return SUNLINEARSOLVER_CUSTOM;
}

int SUNLinSolInitialize_Dense_Eigen(SUNLinearSolver S)
{
	/* all solver-specific memory has already been allocated */
	return(SUNLS_SUCCESS);
}

int SUNLinSolSetup_Dense_Eigen(SUNLinearSolver S, SUNMatrix A)
{
	/* check for valid inputs */
	if ( (A == NULL) || (S == NULL) )
		return(SUNLS_MEM_NULL);

	/* Ensure that A is a dense eigen matrix */
	if (SUNMatGetID(A) != SUNMATRIX_EIGEN_DENSE) {
		return(SUNLS_ILL_INPUT);
	}

	EIGSOL(S).compute_optimized(EIGMAT(A));

#if 0
	std::cout << EIGMAT(A) << std::endl << std::endl;
	std::cout << EIGSOL(S).matrixLU() << std::endl << std::endl;
	Eigen::PartialPivLU< MatrixReal > test;
	test.compute(EIGMAT(A));
	std::cout << test.matrixLU() << std::endl;
#endif

	return(SUNLS_SUCCESS);
}

int SUNLinSolSolve_Dense_Eigen(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                         N_Vector b, realtype tol)
{
	if ( (A == NULL) || (S == NULL) || (x == NULL) || (b == NULL) )
		return(SUNLS_MEM_NULL);

	EIGV(x) = EIGSOL(S).solve(EIGV(b));

	return(SUNLS_SUCCESS);
}

sunindextype SUNLinSolLastFlag_Dense_Eigen(SUNLinearSolver S)
{
	/* return the stored 'last_flag' value */
	if (S == NULL) return(-1);
	return(SUNLS_SUCCESS);
}

int SUNLinSolSpace_Dense_Eigen(SUNLinearSolver S,
                         long int *lenrwLS,
                         long int *leniwLS)
{
	*leniwLS = 0;
	*lenrwLS = 0;
	return(SUNLS_SUCCESS);
}

int SUNLinSolFree_Dense_Eigen(SUNLinearSolver S)
{
	/* return if S is already free */
	if (S == NULL) return(SUNLS_SUCCESS);

	/* delete items from contents, then delete generic structure */
	if (S->content) {
		delete (_SUNLinearSolverContent_Dense_Eigen*)(S->content);
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
