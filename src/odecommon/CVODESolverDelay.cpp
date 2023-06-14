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
	, history_duration(std::numeric_limits<Real>::infinity())
	, max_steps(2000)
	, debug_log(false)
{
}

CVODESolverDelay::~CVODESolverDelay()
{
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
	if (cvode_mem) {
		CVodeFree(&cvode_mem);
	}
}

bool CVODESolverDelay::Initialize(size_t N, void* user)
{
	if (N > std::numeric_limits<long>::max()) {
		LOGERROR("System too large");
		return false;
	}
	this->N = N;

	y = MakeCVodeVector(N);
	interpolate_y = MakeCVodeVector(N);

	cvode_mem = CVodeCreate(CV_BDF);
	//cvode_mem = CVodeCreate(CV_ADAMS);
	if (!cvode_mem) {
		LOGERROR("CVodeCreate failed");
		return false;
	}

	user_data = user;

	CVodeInit(cvode_mem, &static_cvode_rhs_fn, 0.0, (N_Vector)y);
	CVodeSetUserData(cvode_mem, this);
	CVodeSStolerances(cvode_mem, (OdeReal)1e-6, (OdeReal)1e-6);
	CVodeSetMaxNumSteps(cvode_mem, 500);
	CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn, this);

	J = MakeCVodeMatrix(N, N);
	LS = MakeCVodeLinearSolver(y, J);
	CVodeSetLinearSolver(cvode_mem, LS, J);

	//FILE* infofp = fopen("cvode_monitor.txt", "w+");

	NLS = SUNNonlinSol_Newton(y);
	//SUNNonlinSolSetPrintLevel_Newton(NLS, 1);
	//SUNNonlinSolSetInfoFile_Newton(NLS, infofp);
	CVodeSetNonlinearSolver(cvode_mem, NLS);

	if (jacobian) {
		CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);
		CVodeSetMaxStepsBetweenJac(cvode_mem, 10);
	}

	return true;
}

bool CVODESolverDelay::SetTolerance(Real relative, Real absolute)
{
	CVodeSStolerances(cvode_mem, (OdeReal)relative, (OdeReal)absolute);
	return true;
}

bool CVODESolverDelay::SetTolerance(Real relative, VectorReal absolute)
{
	if (absolute.size() != N) {
		return false;
	}

	N_Vector abstol = MakeCVodeVector(N);
	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(abstol,i) = (OdeReal)absolute(i);
	}
	CVodeSVtolerances(cvode_mem, (OdeReal)relative, abstol);
	N_VDestroy(abstol);
	return true;
}

void CVODESolverDelay::SetUserData(void* user)
{
	user_data = user;
}

void CVODESolverDelay::SetKeepHistory(Real duration)
{
	history_duration = duration;
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

bool CVODESolverDelay::Simulate(const Real* initial_conditions, const VectorReal& timepoints, const VectorReal& discontinuities_t, MatrixReal& output)
{
	if (timepoints.size() <= 0) {
		LOGERROR("No time points requested");
		return false;
	}

	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(y,i) = (OdeReal)initial_conditions[i];
	}

	CVodeReInit(cvode_mem, 0.0, (N_Vector)y);
	output.setZero(N, timepoints.size());

	// If the first time point is at t=0, copy the initial conditions to the output
	size_t tpi = 0;
	if (timepoints[tpi] == 0.0) {
		for (size_t i = 0; i < N; i++) {
			output(i, tpi) = NV_Ith_S(y, i);
		}
		tpi++;
	}

	// Store initial condition in history
	history_time.reserve(max_steps);
	history_y.reserve(max_steps);
	history_time.resize(1, 0.0);
	history_y.resize(1);
	history_y.back().resize(N);
	for (size_t i = 0; i < N; i++) {
		history_y.back()(i) = NV_Ith_S(y, i);
	}

	current_dci = 0;
	if (discontinuities_t.size() > 0) {
		CVodeSetStopTime(cvode_mem, (OdeReal)discontinuities_t[current_dci]);
	}

	// Simulate the system one step at a time
	while (1) {
		OdeReal tret;

		int result = CVode(cvode_mem, (OdeReal)timepoints[tpi], (N_Vector)y, &tret, CV_ONE_STEP);
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
			CVodeReInit(cvode_mem, tret, (N_Vector)y);
			current_dci++;
			if (current_dci < (size_t)discontinuities_t.size()) {
				CVodeSetStopTime(cvode_mem, (OdeReal)discontinuities_t[current_dci]);
			}
		} else if (result != CV_SUCCESS) {
#if 0
			LOGERROR("CVode failed, last y values:");
			for (size_t i = 0; i < N; i++) {
				LOGERROR( "%u - %g", i, NV_EIGEN_Ith(y, i));
			}
#endif
			for (size_t j = tpi; j < (size_t)timepoints.size(); j++) {
				for (size_t i = 0; i < N; i++) {
					output(i, tpi) = std::numeric_limits<Real>::quiet_NaN();
				}
			}
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
		history_y.resize(history_y.size() + 1);
		history_y.back().resize(N);
		for (size_t i = 0; i < N; i++) {
			history_y.back()(i) = NV_Ith_S(y, i);
		}
	}

	return true;
}

bool CVODESolverDelay::InterpolateHistory(const OdeReal t, const std::vector< OdeReal >& history_t, const std::vector< OdeVectorReal >& history_y, OdeVectorReal& out)
{
	if (history_t.empty()) {
		LOGERROR("No history supplied");
		return false;
	}
	if (history_t.size() != history_y.size()) {
		LOGERROR("History sizes do not match");
		return false;
	}

#if 1
	if (t <= history_t.front()) {
		out = history_y.front();
	} else if (t >= history_t.back()) {
		out = history_y.back();
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
				out = history_y[ti];
				break;
			} else if (history_t[ti] > t) {
				if (ti == 0) {
					out = history_y.front();
					break;
				} else {
					if (history_t[ti-1] < t) {
						
						OdeReal time_i = history_t[ti];
						OdeReal time_im1 = history_t[ti-1];
						OdeReal a = (t - time_im1) / (time_i - time_im1);

						const OdeVectorReal& y_i = history_y[ti];
						const OdeVectorReal& y_im1 = history_y[ti-1];
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
					out = history_y[ti];
					break;
				} else {
					if (history_t[ti+1] > t) {
						
						OdeReal time_i = history_t[ti+1];
						OdeReal time_im1 = history_t[ti];
						OdeReal a = (t - time_im1) / (time_i - time_im1);

						const OdeVectorReal& y_i = history_y[ti+1];
						const OdeVectorReal& y_im1 = history_y[ti];
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
#else
	for (size_t i = 0; i < history_t.size(); i++) {
		if (history_t[i] == t) {
			out = history_y[i];
			break;
		}

		if (t < history_t[i]) {
			if (i == 0) {
				out = history_y[i];
				return true;
			} else {
				OdeReal time_i = history_t[i];
				OdeReal time_im1 = history_t[i-1];
				OdeReal a = (t - time_im1) / (time_i - time_im1);

				const OdeVectorReal& y_i = history_y[i];
				const OdeVectorReal& y_im1 = history_y[i-1];
				out = y_im1;
				for (int j = 0; j < y_i.size(); j++) {
					out(j) += a * (y_i(j) - y_im1(j));
				}
				return true;
			}
		}
	}
	out = history_y.back();
#endif
	
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
