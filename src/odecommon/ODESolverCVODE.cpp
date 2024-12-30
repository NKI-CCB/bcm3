#include "Utils.h"
#include "ODESolverCVODE.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_math.h>

#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_impl.h"
#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_ls_impl.h"

static const size_t MAX_CVODE_STEPS = 10000;

void static_cvode_err_fn(int error_code, const char *module, const char *function, char *msg, void *user_data)
{
#if 0
	LOGERROR("CVODE error %d in module %s, function %s: %s", error_code, module, function, msg);
#endif
}

int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data)
{
	ODESolverCVODE* solver = reinterpret_cast<ODESolverCVODE*>(user_data);
	return solver->cvode_rhs_fn(t, NV_DATA_S(y), NV_DATA_S(ydot));
}

int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3)
{
	ODESolverCVODE* solver = reinterpret_cast<ODESolverCVODE*>(user_data);
	return solver->cvode_jac_fn(t, NV_DATA_S(y), NV_DATA_S(fy), EIGMAT(Jac));
}

ODESolverCVODE::ODESolverCVODE()
	: cvode_mem(NULL)
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
	CVodeSetMaxNumSteps(cvode_mem, 500); // Note that this is not actually used, since we take one step at a time.
	CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn, this);
	CVodeSetMinStep(cvode_mem, (OdeReal)1e-6);
	CVodeSetMaxStep(cvode_mem, 1.0);
	CVodeSetLinearSolver(cvode_mem, LS, J);
	CVodeSetNonlinearSolver(cvode_mem, NLS);
	CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);

	step_counts.resize(MAX_CVODE_STEPS + 1, 0);
	failed_step_counts.resize(MAX_CVODE_STEPS + 1, 0);

	return true;
}

bool ODESolverCVODE::Simulate(const OdeReal* initial_conditions, const OdeVectorReal& timepoints, OdeMatrixReal& output, bool verbose)
{
	if (timepoints.size() == 0) {
		LOGERROR("No time points requested");
		return false;
	}

	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(y,i) = (OdeReal)initial_conditions[i];
		NV_Ith_S(tmpvector, i) = abstol(i);
	}

	CVodeSVtolerances(cvode_mem, reltol, tmpvector);
	CVodeReInit(cvode_mem, 0.0, (N_Vector)y);
	output.setZero(N, timepoints.size());

	// If the first time point is at t=0, copy the initial conditions to the output
	unsigned int tpi = 0;
	if (timepoints(tpi) == 0.0) {
		for (size_t i = 0; i < N; i++) {
			output(i, tpi) = NV_Ith_S(y, i);
		}
		tpi++;
	}

	if (discontinuity_time == discontinuity_time) {
		CVodeSetStopTime(cvode_mem, (OdeReal)discontinuity_time);
	}

	// Simulate the system one step at a time
	size_t step_count = 0;
	while (tpi < timepoints.size()) {
		// Make sure we don't take too many steps
		if (step_count >= MAX_CVODE_STEPS) {
			//LOGERROR("Maximum number of steps reached, stopping - %s.");
			failed_step_counts[MAX_CVODE_STEPS]++;
			return false;
		}
		step_count++;

		// Take one CVode step
		OdeReal tret;
		int result = CVode(cvode_mem, (OdeReal)timepoints(tpi), (N_Vector)y, &tret, CV_ONE_STEP);
		if (result < 0) {
			if (verbose) {
				LOGERROR("CVode failed at step %u, last y values:", tpi);
				for (size_t i = 0; i < N; i++) {
					LOGERROR("%u - %g", i, NV_Ith_S(y, i));
				}
			}
			for (size_t j = tpi; j < timepoints.size(); j++) {
				for (size_t i = 0; i < N; i++) {
					output(i, tpi) = std::numeric_limits<Real>::quiet_NaN();
				}
			}
			failed_step_counts[step_count]++;
			return false;
		} else {
			// If we pass a timepoint, interpolate back to that timepoint
			while (tret >= timepoints(tpi)) {
				result = CVodeGetDky(cvode_mem, (OdeReal)timepoints(tpi), 0, (N_Vector)tmpvector);
				if (result != CV_SUCCESS) {
					return false;
				}

				for (size_t i = 0; i < N; i++) {
					output(i, tpi) = NV_Ith_S(tmpvector, i);
				}

				tpi++;
				if (tpi >= timepoints.size()) {
					break;
				}
			}

			if (discontinuity_time == discontinuity_time && (result == CV_TSTOP_RETURN || discontinuity_time == tret)) {
				discontinuity_time = discontinuity_cb(tret, discontinuity_user);
				if (discontinuity_time == discontinuity_time) {
					CVodeReInit(cvode_mem, tret, (N_Vector)y);
					CVodeSetStopTime(cvode_mem, (OdeReal)discontinuity_time);
				} else if (discontinuity_time == std::numeric_limits<Real>::infinity()) {
					// Continue simulation but no new discontinuity
					CVodeReInit(cvode_mem, tret, (N_Vector)y);
				}
			}
		}
	}

	step_counts[step_count]++;

	return true;
}

void ODESolverCVODE::DumpStatistics(const char* filename)
{
	FILE* file = fopen(filename, "w");
	if (file) {
		for (size_t i = 0; i <= MAX_CVODE_STEPS; i++) {
			fprintf(file, "%zu\t%zu\n", i, step_counts[i]);
		}
		for (size_t i = 0; i <= MAX_CVODE_STEPS; i++) {
			fprintf(file, "%zu\t%zu\n", i, failed_step_counts[i]);
		}

		fclose(file);
	}
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

OdeReal ODESolverCVODE::get_y(size_t i)
{
	return NV_Ith_S(this->y, i);
}

void ODESolverCVODE::set_y(size_t i, OdeReal y)
{
	NV_Ith_S(this->y, i) = y;
}
