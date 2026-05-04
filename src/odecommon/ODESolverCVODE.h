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

	virtual bool Initialize(size_t N, void* user, size_t store_integration_points_buffer);
	virtual bool Restart();
	virtual bool SetSolverParameter(const std::string& parameter, int int_value, OdeReal real_value);
	
	virtual void RestartInterpolationIteration();
	virtual const OdeVectorReal& GetInterpolatedY(OdeReal t);
	virtual OdeReal GetInterpolatedY(OdeReal t, size_t i);

	virtual OdeVectorReal get_current_y() const;
	virtual OdeReal get_current_y(size_t i) const;
	virtual void set_current_y(size_t i, OdeReal y) const;
	virtual Real get_threshold_crossing_time(size_t species_ix, Real threshold, bool above, Real prev_time) const;

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
	OdeReal t;
	size_t current_step;

	struct CVodeTimepoint
	{
		CVodeTimepoint();
		OdeReal cvode_time;
		OdeReal cv_uround;
		OdeReal cv_tn;
		OdeReal cv_h;
		OdeReal cv_hu;
		int cv_q;
	};
	std::vector<CVodeTimepoint> cvode_timepoints;
	OdeMatrixReal cvode_timepoints_zn[6];
	size_t cvode_timepoint_iter;
	OdeVectorReal interpolation_y;
	OdeReal interpolation_time;

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
