#include "Utils.h"
#include "LinearAlgebraSelector.h"

#if CVODE_USE_EIGEN_SOLVER

N_Vector MakeCVodeVector(size_t N)
{
	return N_VNew_Eigen((long)N);
}

#if CVODE_USE_SPARSE_SOLVER

SUNMatrix MakeCVodeMatrix(size_t N, size_t M)
{
	return SUNSparseEigenMatrix(M, N);
}
SUNLinearSolver MakeCVodeLinearSolver(N_Vector y, SUNMatrix J)
{
	return SUNLinSol_Sparse_Eigen(y, J);
}
#else

SUNMatrix MakeCVodeMatrix(size_t N, size_t M)
{
	return SUNDenseEigenMatrix(M, N);
}
SUNLinearSolver MakeCVodeLinearSolver(N_Vector y, SUNMatrix J)
{
	return SUNLinSol_Dense_Eigen(y, J);
}
#endif

#else

N_Vector MakeCVodeVector(size_t N)
{
	return N_VNew_Serial((long)N);
}
SUNMatrix MakeCVodeMatrix(size_t N, size_t M)
{
	return SUNDenseMatrix(M, N);
}
SUNLinearSolver MakeCVodeLinearSolver(N_Vector y, SUNMatrix J)
{
	return SUNLinSol_Dense(y, J);
}

#endif
