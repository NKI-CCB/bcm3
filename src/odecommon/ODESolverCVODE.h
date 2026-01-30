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
	virtual bool SetSolverParameter(const std::string& parameter, int int_value, OdeReal real_value);
	virtual OdeReal GetInterpolatedY(OdeReal t, size_t i);
	virtual OdeReal get_current_y(size_t i);
	virtual void set_current_y(size_t i, OdeReal y);

protected:
	virtual bool Solve(const OdeVectorReal& initial_conditions, OdeReal end_time, bool do_interpolation, bool store_integration_points, bool verbose);

	// Settings
	int max_steps;

	// Runtime variables
	void* cvode_mem;
	SUNLinearSolver LS;
	SUNNonlinearSolver NLS;
	N_Vector y;
	N_Vector tmpvector;
	SUNMatrix J;

	int cvode_rhs_fn(OdeReal t, const OdeReal* y, OdeReal* ydot);
	friend int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data);

#if CVODE_USE_EIGEN_SOLVER
	OdeVectorReal y_copy;
	OdeVectorReal work_vector;
	bool DifferenceQuotientJacobian(OdeReal t, const OdeVectorReal& y, const OdeVectorReal& ydot, OdeMatrixReal& jac);
	int cvode_jac_fn(OdeReal t, const OdeVectorReal& y, const OdeVectorReal& ydot, OdeMatrixReal& jac);
	friend int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3);
#endif
};
