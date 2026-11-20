#include "Utils.h"
#include "CVODESolverDelay.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>

static void static_cvode_err_fn(int error_code, const char *module, const char *function, char *msg, void *user_data)
{
#if 0
	LOGERROR("CVODE error %d in module %s, function %s: %s", error_code, module, function, msg);
#endif
}

static int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data)
{
	CVODESolverDelay* solver = reinterpret_cast<CVODESolverDelay*>(user_data);
	return solver->cvode_rhs_fn(t, y, ydot);
}

static int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3)
{
	CVODESolverDelay* solver = reinterpret_cast<CVODESolverDelay*>(user_data);
	return solver->cvode_jac_fn(t, y, fy, Jac);
}

CVODESolverDelay::CVODESolverDelay()
	: cvode_mem(NULL)
	, LS(NULL)
	, NLS(NULL)
	, y(NULL)
	, interpolate_y(NULL)
	, J(NULL)
	, user_data(NULL)
	, max_steps(2000)
	, debug_log(false)
{
}

CVODESolverDelay::~CVODESolverDelay()
{
	if (cvode_mem) {
		CVodeFree(&cvode_mem);
	}
	if (y) {
		N_VDestroy(y);
	}
	if (interpolate_y) {
		N_VDestroy(interpolate_y);
	}
	if (NLS) {
		SUNNonlinSolFree(NLS);
	}
	if (LS) {
		SUNLinSolFree(LS);
	}
	if (J) {
		SUNMatDestroy(J);
	}
}

bool CVODESolverDelay::Initialize(size_t N, void* user)
{
	if (N == 0 || N > std::numeric_limits<long>::max() || y || interpolate_y) {
		LOGERROR("Invalid system size or delay solver already initialized.");
		return false;
	}
#if CVODE_USE_EIGEN_SOLVER
	if (!jacobian) {
		LOGERROR("The Eigen CVODE delay solver requires a Jacobian callback.");
		return false;
	}
#endif
	this->N = N;

	y = MakeCVodeVector(N);
	interpolate_y = MakeCVodeVector(N);
	if (!y || !interpolate_y) {
		LOGERROR("Could not allocate CVODE delay vectors.");
		return false;
	}

	cvode_mem = CVodeCreate(CV_BDF);
	//cvode_mem = CVodeCreate(CV_ADAMS);
	if (!cvode_mem) {
		LOGERROR("CVodeCreate failed");
		return false;
	}

	user_data = user;

	int result = CVodeInit(cvode_mem, &static_cvode_rhs_fn, 0.0, y);
	if (result != CV_SUCCESS) {
		LOGERROR("CVodeInit failed (%d).", result);
		return false;
	}

	result = CVodeSetUserData(cvode_mem, this);
	if (result != CV_SUCCESS) {
		LOGERROR("CVodeSetUserData failed (%d).", result);
		return false;
	}

	result = CVodeSStolerances(cvode_mem, (OdeReal)1e-6, (OdeReal)1e-6);
	if (result != CV_SUCCESS) {
		LOGERROR("CVodeSStolerances failed (%d).", result);
		return false;
	}

	result = CVodeSetMaxNumSteps(cvode_mem, 500);
	if (result != CV_SUCCESS) {
		LOGERROR("CVodeSetMaxNumSteps failed (%d).", result);
		return false;
	}

	result = CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn, this);
	if (result != CV_SUCCESS) {
		LOGERROR("CVodeSetErrHandlerFn failed (%d).", result);
		return false;
	}

	J = MakeCVodeMatrix(N, N);
	if (!J) {
		LOGERROR("Could not allocate the CVODE Jacobian matrix.");
		return false;
	}

	LS = MakeCVodeLinearSolver(y, J);
	if (!LS) {
		LOGERROR("Could not create the CVODE linear solver.");
		return false;
	}

	NLS = SUNNonlinSol_Newton(y);
	if (!NLS) {
		LOGERROR("SUNNonlinSol_Newton failed.");
		return false;
	}

	result = CVodeSetLinearSolver(cvode_mem, LS, J);
	if (result != CVLS_SUCCESS) {
		LOGERROR("CVodeSetLinearSolver failed (%d).", result);
		return false;
	}

	result = CVodeSetNonlinearSolver(cvode_mem, NLS);
	if (result != CV_SUCCESS) {
		LOGERROR("CVodeSetNonlinearSolver failed (%d).", result);
		return false;
	}

	if (jacobian) {
		result = CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);
		if (result != CVLS_SUCCESS) {
			LOGERROR("CVodeSetJacFn failed (%d).", result);
			return false;
		}

		result = CVodeSetMaxStepsBetweenJac(cvode_mem, 10);
		if (result != CVLS_SUCCESS) {
			LOGERROR("CVodeSetMaxStepsBetweenJac failed (%d).", result);
			return false;
		}
	}

	return true;
}

bool CVODESolverDelay::SetTolerance(OdeReal relative, OdeReal absolute)
{
	return cvode_mem && CVodeSStolerances(cvode_mem, relative, absolute) == CV_SUCCESS;
}

bool CVODESolverDelay::SetTolerance(OdeReal relative, OdeVectorReal absolute)
{
	if (!cvode_mem || absolute.size() != N) {
		return false;
	}

	N_Vector abstol = MakeCVodeVector(N);
	if (!abstol) {
		return false;
	}
	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(abstol,i) = absolute(i);
	}
	int result = CVodeSVtolerances(cvode_mem, relative, abstol);
	N_VDestroy(abstol);
	return result == CV_SUCCESS;
}

void CVODESolverDelay::SetUserData(void* user)
{
	user_data = user;
}

void CVODESolverDelay::SetDerivativeFunction(TDeriviativeFunction f)
{
	derivative = f;
}

void CVODESolverDelay::SetJacobianFunction(TJacobianFunction f)
{
	jacobian = f;
}

void CVODESolverDelay::SetDebugLogging(bool log)
{
	debug_log = log;
}

bool CVODESolverDelay::Simulate(const OdeReal* initial_conditions, const OdeVectorReal& timepoints, const OdeVectorReal& discontinuities_t, OdeReal min_delay, OdeMatrixReal& output)
{
	if (timepoints.size() <= 0) {
		LOGERROR("No time points requested");
		return false;
	}

	// Leave room for rounding when subtracting the delay from a trial time.
	const OdeReal max_step = min_delay * (1 - 4 * std::numeric_limits<OdeReal>::epsilon());
	if (!std::isfinite(max_step) || max_step <= 0.0 || CVodeSetMaxStep(cvode_mem, max_step) != CV_SUCCESS) {
		LOGERROR("Could not set the delay solver step limit from minimum delay %g.", min_delay);
		return false;
	}

	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(y,i) = (OdeReal)initial_conditions[i];
	}

	if (CVodeReInit(cvode_mem, 0.0, y) != CV_SUCCESS) {
		return false;
	}
	output.setConstant(N, timepoints.size(), std::numeric_limits<OdeReal>::quiet_NaN());

	// Copy the initial conditions to any requested time points at t=0.
	size_t tpi = 0;
	while (tpi < timepoints.size() && timepoints[tpi] == 0.0) {
		for (size_t i = 0; i < N; i++) {
			output(i, tpi) = NV_Ith_S(y, i);
		}
		tpi++;
	}

	// Store initial condition in history
	history_time.reserve(max_steps);
	history_y.setConstant(N, max_steps, std::numeric_limits<Real>::quiet_NaN());
	history_time.resize(1, 0.0);
	for (size_t i = 0; i < N; i++) {
		history_y(i,0) = NV_Ith_S(y, i);
	}

	current_dci = 0;
	if (tpi == timepoints.size()) {
		return true;
	}
	if (discontinuities_t.size() > 0 && CVodeSetStopTime(cvode_mem, discontinuities_t[current_dci]) != CV_SUCCESS) {
		return false;
	}

	// Simulate the system one step at a time
	while (1) {
		OdeReal tret;

		// CVODE estimates the initial step from tout before applying the maximum step size.
		const OdeReal tout = std::min(timepoints[tpi], history_time.back() + max_step);
		int result = CVode(cvode_mem, tout, (N_Vector)y, &tret, CV_ONE_STEP);
		if (debug_log) {
			LOG("%6g - %g - %.12g %.12g %.12g", timepoints[tpi], tret, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2));
		}
		if (result == CV_TSTOP_RETURN) {
			// We've met a discontinuity
			
			// If we pass a timepoint, interpolate back to that timepoint
			while (tret >= timepoints[tpi]) {
				result = CVodeGetDky(cvode_mem, (OdeReal)timepoints[tpi], 0, (N_Vector)interpolate_y);
				if (result != CV_SUCCESS) {
					return false;
				}
			
				for (size_t i = 0; i < N; i++) {
					output(i, tpi) = NV_Ith_S(interpolate_y, i);
				}

				tpi++;
				if (tpi >= (size_t)timepoints.size()) {
					break;
				}
			}

			if (tpi >= (size_t)timepoints.size()) {
				break;
			}

			// Restart cvode
			if (CVodeReInit(cvode_mem, tret, y) != CV_SUCCESS) {
				return false;
			}
			current_dci++;
			if (current_dci < discontinuities_t.size() && CVodeSetStopTime(cvode_mem, discontinuities_t[current_dci]) != CV_SUCCESS) {
				return false;
			}
		} else if (result != CV_SUCCESS) {
#if 0
			LOGERROR("CVode failed, last y values:");
			for (size_t i = 0; i < N; i++) {
				LOGERROR( "%u - %g", i, NV_EIGEN_Ith(y, i));
			}
#endif
			return false;
		} else {
			// If we pass a timepoint, interpolate back to that timepoint
			while (tret >= timepoints[tpi]) {
				result = CVodeGetDky(cvode_mem, (OdeReal)timepoints[tpi], 0, (N_Vector)interpolate_y);
				if (result != CV_SUCCESS) {
					return false;
				}
			
				for (size_t i = 0; i < N; i++) {
					output(i, tpi) = NV_Ith_S(interpolate_y, i);
				}

				tpi++;
				if (tpi >= (size_t)timepoints.size()) {
					break;
				}
			}
			if (tpi >= (size_t)timepoints.size()) {
				break;
			}
		}

		if (history_time.size() >= max_steps) {
			return false;
		}

		// Store history
		history_time.push_back(tret);
		for (size_t i = 0; i < N; i++) {
			history_y(i, history_time.size() - 1) = NV_Ith_S(y, i);
		}
	}

	return true;
}

bool CVODESolverDelay::InterpolateHistory(const OdeReal t, const std::vector< OdeReal >& history_t, const OdeMatrixReal& history_y, OdeVectorReal& out)
{
	if (history_t.empty()) {
		LOGERROR("No history supplied");
		return false;
	}

	if (t <= history_t.front()) {
		out = history_y.col(0);
	} else if (t >= history_t.back()) {
		// Allow rounding at the end of the history.
		const OdeReal slack = 4 * std::numeric_limits<OdeReal>::epsilon() * std::abs(t);
		if (t - history_t.back() > slack) {
			return false;
		}
		out = history_y.col(history_t.size() - 1);
	} else {
		size_t imin = 0;
		size_t imax = history_t.size();

		size_t step_count = 0;
		while (1) {
			size_t ti = imin + (imax - imin) / 2;
			step_count++;
			if (step_count == 50) {
				LOGERROR("Not found in 50 steps?!");
				return false;
			}

			if (history_t[ti] == t) {
				out = history_y.col(ti);
				break;
			} else if (history_t[ti] > t) {
				if (ti == 0) {
					out = history_y.col(0);
					break;
				} else {
					if (history_t[ti-1] < t) {
						
						OdeReal time_i = history_t[ti];
						OdeReal time_im1 = history_t[ti-1];
						OdeReal a = (t - time_im1) / (time_i - time_im1);

						const OdeVectorReal& y_i = history_y.col(ti);
						const OdeVectorReal& y_im1 = history_y.col(ti-1);
						out = y_im1;
						for (int j = 0; j < y_i.size(); j++) {
							out(j) += a * (y_i(j) - y_im1(j));
						}
						break;

					} else {
						imax = ti;
					}
				}
			} else {
				if (ti == history_t.size() - 1) {
					out = history_y.col(ti);
					break;
				} else {
					if (history_t[ti+1] > t) {
						
						OdeReal time_i = history_t[ti+1];
						OdeReal time_im1 = history_t[ti];
						OdeReal a = (t - time_im1) / (time_i - time_im1);

						const OdeVectorReal& y_i = history_y.col(ti+1);
						const OdeVectorReal& y_im1 = history_y.col(ti);
						out = y_im1;
						for (int j = 0; j < y_i.size(); j++) {
							out(j) += a * (y_i(j) - y_im1(j));
						}
						break;

					} else {
						imin = ti;
					}
				}
			}
		}
	}
	
	return true;
}

bool CVODESolverDelay::InterpolateHistory(const OdeReal t, OdeVectorReal& out)
{
	return InterpolateHistory(t, history_time, history_y, out);
}

int CVODESolverDelay::cvode_rhs_fn(OdeReal t, void* y_nvector, void* ydot_nvector)
{
	if (!derivative) {
		return -2;
	}

	bool result = derivative(t, NV_DATA_S(((N_Vector)y_nvector)), history_time, history_y, current_dci, NV_DATA_S(((N_Vector)ydot_nvector)), user_data);
	if (result) {
		return 0;
	} else {
		return -1;
	}
}

int CVODESolverDelay::cvode_jac_fn(OdeReal t, void* y_nvector, void* fy_nvector, void* Jac_matrix)
{
	if (!jacobian) {
		return -2;
	}

#if CVODE_USE_EIGEN_SOLVER
	bool result = jacobian(t, NV_DATA_S(((N_Vector)y_nvector)), NV_DATA_S(((N_Vector)fy_nvector)), history_time, history_y, current_dci, EIGMAT(((SUNMatrix)Jac_matrix)), user_data);
	if (result) {
		return 0;
} else {
		return -1;
	}
#else
	return -3;
#endif
}
