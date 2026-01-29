#include "Utils.h"
#include "ODESolverCVODE.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_math.h>

#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_impl.h"
#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_ls_impl.h"

static const size_t MAX_CVODE_STEPS = 10000;

void static_cvode_err_fn(int error_code, const char *module, const char *function, char *msg, void *user_data)
{
	LOGERROR("CVODE error %d in module %s, function %s: %s", error_code, module, function, msg);
}

int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data)
{
	ODESolverCVODE* solver = reinterpret_cast<ODESolverCVODE*>(user_data);
	return solver->cvode_rhs_fn(t, NV_DATA_S(y), NV_DATA_S(ydot));
}

int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3)
{
	ODESolverCVODE* solver = reinterpret_cast<ODESolverCVODE*>(user_data);
#if CVODE_USE_EIGEN_SOLVER
	return solver->cvode_jac_fn(t, NV_DATA_S(y), NV_DATA_S(fy), EIGMAT(Jac));
#else
	ASSERT(false);
	return -1;
#endif
}

ODESolverCVODE::ODESolverCVODE()
	: max_steps(2000)
	, cvode_mem(NULL)
	, LS(NULL)
	, NLS(NULL)
	, y(NULL)
	, tmpvector(NULL)
	, J(NULL)
{
}

ODESolverCVODE::~ODESolverCVODE()
{
	if (y) {
		N_VDestroy((N_Vector)y);
	}
	if (tmpvector) {
		N_VDestroy((N_Vector)tmpvector);
	}
	if (J) {
		SUNMatDestroy(J);
	}
	if (NLS) {
		SUNNonlinSolFree(NLS);
	}
	if (LS) {
		SUNLinSolFree(LS);
	}
	if (cvode_mem) {
		CVodeFree(&cvode_mem);
	}
}

bool ODESolverCVODE::Initialize(size_t N, void* user)
{
	if (!ODESolver::Initialize(N, user)) {
		return false;
	}

	if (N > std::numeric_limits<long>::max()) {
		LOGERROR("System too large to simulate with CVode");
		return false;
	}

	y = MakeCVodeVector((long)N);
	tmpvector = MakeCVodeVector((long)N);

	cvode_mem = CVodeCreate(CV_BDF);
	if (!cvode_mem) {
		LOGERROR("CVodeCreate failed");
		return false;
	}

	user_data = user;

	J = MakeCVodeMatrix(N, N);
	LS = MakeCVodeLinearSolver(y, J);
	NLS = SUNNonlinSol_Newton(y);

	CVodeInit(cvode_mem, &static_cvode_rhs_fn, 0.0, y);
	CVodeSetUserData(cvode_mem, this);
	CVodeSetLinearSolver(cvode_mem, LS, J);
	CVodeSetNonlinearSolver(cvode_mem, NLS);
	CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);

	return true;
}

bool ODESolverCVODE::SetSolverParameter(const std::string& parameter, int int_value, OdeReal real_value)
{
	if (cvode_mem == nullptr) {
		LOGERROR("Solver should be initialized before specifying parameters.");
		return false;
	}

	if (parameter == "min_dt") {
		if (real_value <= 0.0) {
			LOGERROR("min_dt should be strictly positive, but %g was provided", real_value);
			return false;
		}
		CVodeSetMinStep(cvode_mem, (OdeReal)real_value);
		return true;
	} if (parameter == "max_dt") {
		if (real_value < 0.0) {
			LOGERROR("max_dt should be non-negative, but %g was provided", real_value);
			return false;
		}
		CVodeSetMaxStep(cvode_mem, (OdeReal)real_value);
		return true;
	} else if (parameter == "max_steps") {
		if (int_value <= 0) {
			LOGERROR("max_steps should be strictly positive, but %d was provided", int_value);
			return false;
		}
		max_steps = int_value;
		return true;
	} else {
		LOGERROR("Unknown solver parameter \"%s\"", parameter.c_str());
		return false;
	}
}

OdeReal ODESolverCVODE::GetInterpolatedY(OdeReal t, size_t i)
{
	// Not yet implemented
	return std::numeric_limits<Real>::quiet_NaN();
}

OdeReal ODESolverCVODE::get_current_y(size_t i)
{
	return NV_Ith_S(this->y, i);
}

void ODESolverCVODE::set_current_y(size_t i, OdeReal y)
{
	NV_Ith_S(this->y, i) = y;
}

bool ODESolverCVODE::Solve(const OdeVectorReal& initial_conditions, OdeReal end_time, bool do_interpolation, bool store_integration_points, bool verbose)
{
	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(y, i) = initial_conditions[i];
		NV_Ith_S(tmpvector, i) = absolute_tolerance(i);
	}

	if (verbose) {
		CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn, this);
	}

	CVodeSVtolerances(cvode_mem, relative_tolerance, tmpvector);
	CVodeReInit(cvode_mem, 0.0, (N_Vector)y);

	if (!std::isnan(next_discontinuity_time)) {
		CVodeSetStopTime(cvode_mem, (OdeReal)next_discontinuity_time);
	}

	// Simulate the system one step at a time
	size_t step_count = 0;
	OdeReal t = 0.0;
	size_t tpi = interpolation_timepoints_start;
	while (1) {
		// Make sure we don't take too many steps
		if (step_count >= MAX_CVODE_STEPS) {
			if (verbose) {
				LOGERROR("Maximum number of steps reached, stopping.");
			}
			return false;
		}
		step_count++;

		// Take one CVode step
		OdeReal tret;
		int result = CVode(cvode_mem, end_time, (N_Vector)y, &tret, CV_ONE_STEP);
		if (result < 0) {
			if (verbose) {
				LOGERROR("CVode failed at step %u, last y values:", tpi);
				for (size_t i = 0; i < N; i++) {
					LOGERROR("%u - %g", i, NV_Ith_S(y, i));
				}
			}
			if (do_interpolation) {
				for (size_t j = tpi; j < interpolation_timepoints->size(); j++) {
					for (size_t i = 0; i < N; i++) {
						(*interpolated_output)(i, tpi) = std::numeric_limits<Real>::quiet_NaN();
					}
				}
			}
			return false;
		} else {
			if (do_interpolation) {
				// If we pass a timepoint, interpolate back to that timepoint
				while (tret >= (*interpolation_timepoints)(tpi)) {
					result = CVodeGetDky(cvode_mem, (*interpolation_timepoints)(tpi), 0, (N_Vector)tmpvector);
					if (result != CV_SUCCESS) {
						return false;
					}

					for (size_t i = 0; i < N; i++) {
						(*interpolated_output)(i, tpi) = NV_Ith_S(tmpvector, i);
					}

					tpi++;
					if (tpi >= interpolation_timepoints->size()) {
						break;
					}
				}
			}

			if (!std::isnan(next_discontinuity_time) && (result == CV_TSTOP_RETURN || next_discontinuity_time == tret)) {
				next_discontinuity_time = discontinuity_cb(tret, discontinuity_user);
				if (!std::isnan(next_discontinuity_time) && next_discontinuity_time < std::numeric_limits<Real>::infinity()) {
					CVodeReInit(cvode_mem, tret, (N_Vector)y);
					CVodeSetStopTime(cvode_mem, next_discontinuity_time);
				} else {
					// Continue simulation but no new discontinuity
					CVodeReInit(cvode_mem, tret, (N_Vector)y);
				}
			}
		}
	}

	return true;
}

int ODESolverCVODE::cvode_rhs_fn(OdeReal t, OdeReal* y, OdeReal* ydot)
{
	if (!derivative) {
		return -2;
	}

	bool result = derivative(t, y, ydot, user_data);
	if (result) {
		return 0;
	} else {
		return -1;
	}
}

int ODESolverCVODE::cvode_jac_fn(OdeReal t, OdeReal* y, OdeReal* ydot, OdeMatrixReal& jac)
{
	if (!jacobian) {
		return -2;
	}

	bool result = jacobian(t, y, ydot, jac, user_data);
	if (result) {
		return 0;
	} else {
		return -1;
	}
}
