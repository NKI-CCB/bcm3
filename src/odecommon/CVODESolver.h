#pragma once

#include <cvode/cvode.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include "LinearAlgebraSelector.h"

class CVODESolver
{
public:
	typedef std::function<bool (OdeReal, const OdeReal*, OdeReal*, void*)> TDeriviativeFunction;
	typedef std::function<bool (OdeReal, const OdeReal*, const OdeReal*, OdeReal**, void*)> TJacobianFunction;
	typedef std::function<Real (OdeReal, void*)> TDiscontinuityCallback;

	CVODESolver();
	~CVODESolver();

	bool Initialize(size_t N, void* user);
	bool SetTolerance(Real relative, Real absolute);
	bool SetTolerance(Real relative, VectorReal absolute);
	void SetUserData(void* user);
	void SetDiscontinuity(Real time, TDiscontinuityCallback cb, void* user);
	void ResetDiscontinuity();

	void SetDerivativeFunction(TDeriviativeFunction f);
	void SetJacobianFunction(TJacobianFunction f);

	bool Simulate(const OdeReal* initial_conditions, const OdeVectorReal& timepoints, OdeMatrixReal& output);

	void DumpStatistics(const char* filename);

	int cvode_rhs_fn(OdeReal t, OdeReal* y, OdeReal* ydot);
	int cvode_jac_fn(OdeReal t, N_Vector y, N_Vector ydot, SUNMatrix jac);
	OdeReal get_y(size_t i);
	void set_y(size_t i, OdeReal y);

private:
	int cvode_jac_dq(OdeReal t, N_Vector y, N_Vector ydot, SUNMatrix jac);
	void* cvode_mem;
	SUNLinearSolver LS;
	SUNNonlinearSolver NLS;
	size_t N;
	N_Vector y;
	N_Vector interpolate_y;
	SUNMatrix J;
	void* user_data;
	Real discontinuity_time;
	TDiscontinuityCallback discontinuity_cb;
	void* discontinuity_user;

	TDeriviativeFunction derivative;
	TJacobianFunction jacobian;

	std::vector<size_t> step_counts;
	std::vector<size_t> failed_step_counts;
};
