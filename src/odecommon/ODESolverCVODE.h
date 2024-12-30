#pragma once

#include <cvode/cvode.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include "ODESolver.h"

class ODESolverCVODE : public ODESolver
{
public:
	ODESolverCVODE();
	virtual ~ODESolverCVODE();

	virtual bool Initialize(size_t N, void* user);
	virtual bool Simulate(const OdeReal* initial_conditions, const OdeVectorReal& timepoints, OdeMatrixReal& output, bool verbose = false);

	void DumpStatistics(const char* filename);

	virtual OdeReal get_y(size_t i);
	virtual void set_y(size_t i, OdeReal y);

protected:
	int cvode_rhs_fn(OdeReal t, OdeReal* y, OdeReal* ydot);
	int cvode_jac_fn(OdeReal t, OdeReal* y, OdeReal* ydot, OdeMatrixReal& jac);

	void* cvode_mem;
	SUNLinearSolver LS;
	SUNNonlinearSolver NLS;
	N_Vector y;
	N_Vector tmpvector;
	SUNMatrix J;

	std::vector<size_t> step_counts;
	std::vector<size_t> failed_step_counts;

	friend int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data);
	friend int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3);
};
