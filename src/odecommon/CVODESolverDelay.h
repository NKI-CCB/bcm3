#pragma once

#include <cvode/cvode.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_nonlinearsolver.h>
#include "LinearAlgebraSelector.h"

#if ODE_SINGLE_PRECISION
	typedef float OdeReal;
	typedef Eigen::VectorXf OdeVectorReal;
	typedef Eigen::MatrixXf OdeMatrixReal;
#else
	typedef double OdeReal;
	typedef Eigen::VectorXd OdeVectorReal;
	typedef Eigen::MatrixXd OdeMatrixReal;
#endif

class CVODESolverDelay
{
public:
	typedef std::function<bool (OdeReal, const OdeReal*, const std::vector< OdeReal >&, const std::vector< OdeVectorReal >&, size_t, OdeReal*, void*)> TDeriviativeFunction;

	CVODESolverDelay();
	~CVODESolverDelay();

	bool Initialize(size_t N, void* user);
	bool SetTolerance(Real relative, Real absolute);
	bool SetTolerance(Real relative, VectorReal absolute);
	void SetUserData(void* user);
	void SetKeepHistory(Real duration);

	void SetDerivativeFunction(TDeriviativeFunction f);
	void SetDebugLogging(bool log);

	bool Simulate(const Real* initial_conditions, const VectorReal& timepoints, const VectorReal& discontinuities_t, MatrixReal& output);
	static bool InterpolateHistory(const OdeReal t, const std::vector< OdeReal >& history_t, const std::vector< OdeVectorReal >& history_y, OdeVectorReal& out);
	bool InterpolateHistory(const OdeReal t, OdeVectorReal& out);
	
	int cvode_rhs_fn(OdeReal t, void* y_nvector, void* ydot_nvector);

private:
	void* cvode_mem;
	SUNLinearSolver LS;
	SUNNonlinearSolver NLS;
	size_t N;
	N_Vector y;
	N_Vector interpolate_y;
	SUNMatrix J;
	void* user_data;
	Real history_duration;
	size_t max_steps;
	bool debug_log;

	std::vector< OdeReal > history_time;
	std::vector< OdeVectorReal > history_y;
	size_t current_dci;

	TDeriviativeFunction derivative;
};
