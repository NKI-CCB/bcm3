#include "Utils.h"
#include "CVODESolver.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_math.h>

#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_impl.h"
#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_ls_impl.h"

static const size_t MAX_CVODE_STEPS = 10000;

static void static_cvode_err_fn(int error_code, const char *module, const char *function, char *msg, void *user_data)
{
#if 0
	LOGERROR("CVODE error %d in module %s, function %s: %s", error_code, module, function, msg);
#endif
}

static int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data)
{
	CVODESolver* solver = reinterpret_cast<CVODESolver*>(user_data);
	return solver->cvode_rhs_fn(t, NV_DATA_S(y), NV_DATA_S(ydot));
}

static int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3)
{
	CVODESolver* solver = reinterpret_cast<CVODESolver*>(user_data);
	return solver->cvode_jac_fn(t, y, fy, Jac);
}

CVODESolver::CVODESolver()
	: cvode_mem(NULL)
	, LS(NULL)
	, NLS(NULL)
	, y(NULL)
	, interpolate_y(NULL)
	, J(NULL)
	, user_data(NULL)
	, discontinuity_time(std::numeric_limits<Real>::quiet_NaN())
{
}

CVODESolver::~CVODESolver()
{
	if (y) {
		N_VDestroy((N_Vector)y);
	}
	if (interpolate_y) {
		N_VDestroy((N_Vector)interpolate_y);
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

bool CVODESolver::Initialize(size_t N, void* user)
{
	if (N > std::numeric_limits<long>::max()) {
		LOGERROR("System too large");
		return false;
	}
	this->N = N;

	y = MakeCVodeVector((long)N);
	interpolate_y = MakeCVodeVector((long)N);

	cvode_mem = CVodeCreate(CV_BDF);
	//cvode_mem = CVodeCreate(CV_ADAMS);
	if (!cvode_mem) {
		LOGERROR("CVodeCreate failed");
		return false;
	}

	user_data = user;

	CVodeInit(cvode_mem, &static_cvode_rhs_fn, 0.0, y);
	CVodeSetUserData(cvode_mem, this);
	CVodeSStolerances(cvode_mem, (OdeReal)1e-12, (OdeReal)1e-12);
	CVodeSetMaxNumSteps(cvode_mem, 500); // Note that this is not actually used, since we take one step at a time.
	CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn, this);
	CVodeSetMinStep(cvode_mem, (OdeReal)1e-12);

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
	} else {
#if CVODE_USE_EIGEN_SOLVER
		// The CVode reference implementation of difference quotient Jacobian approximation is not available, so use our own
		CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);
#endif
	}

	step_counts.resize(MAX_CVODE_STEPS + 1, 0);
	failed_step_counts.resize(MAX_CVODE_STEPS + 1, 0);

	return true;
}

bool CVODESolver::SetTolerance(Real relative, Real absolute)
{
	CVodeSStolerances(cvode_mem, (OdeReal)relative, (OdeReal)absolute);
	return true;
}

bool CVODESolver::SetTolerance(Real relative, VectorReal absolute)
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

void CVODESolver::SetUserData(void* user)
{
	user_data = user;
}

void CVODESolver::SetDiscontinuity(Real time, TDiscontinuityCallback cb, void* user)
{
	discontinuity_time = time;
	discontinuity_cb = cb;
	discontinuity_user = user;
}

void CVODESolver::ResetDiscontinuity()
{
	discontinuity_time = std::numeric_limits<Real>::quiet_NaN();
}

void CVODESolver::SetDerivativeFunction(TDeriviativeFunction f)
{
	derivative = f;
}

void CVODESolver::SetJacobianFunction(TJacobianFunction f)
{
	jacobian = f;
}

bool CVODESolver::Simulate(const OdeReal* initial_conditions, const OdeVectorReal& timepoints, OdeMatrixReal& output)
{
	if (timepoints.size() == 0) {
		LOGERROR("No time points requested");
		return false;
	}

	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(y,i) = (OdeReal)initial_conditions[i];
	}

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

#if 1
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
#if 0
			LOGERROR("CVode failed at step %u, last y values:", tpi);
			for (size_t i = 0; i < N; i++) {
				LOGERROR("%u - %g", i, NV_Ith_S(y, i));
			}
#endif
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
				result = CVodeGetDky(cvode_mem, (OdeReal)timepoints(tpi), 0, (N_Vector)interpolate_y);
				if (result != CV_SUCCESS) {
					return false;
				}

				for (size_t i = 0; i < N; i++) {
					output(i, tpi) = NV_Ith_S(interpolate_y, i);
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
#else
	// Simulate the system until each time point
	Real prevtime = 0.0;
	for (; tpi < timepoints.size(); tpi++) {
		if (timepoints[tpi] <= prevtime) {
			LOGERROR("Time points are not monotonically increasing");
			return false;
		}

		Real tret;
		int result = CVode(cvode_mem, timepoints[tpi], (N_Vector)y, &tret, CV_NORMAL);
		if (result != CV_SUCCESS) {
#if 0
			LOGERROR("CVode failed, last species values:");
			for (size_t i = 0; i < N; i++) {
				LOGERROR( "%u - %g", i, NV_EIGEN_Ith(y, i));
			}
#endif
			for (size_t j = tpi; j < timepoints.size(); j++) {
				for (size_t i = 0; i < N; i++) {
					output(i, tpi) = std::numeric_limits<Real>::quiet_NaN();
				}
			}
			return false;
		}
		
		for (size_t i = 0; i < N; i++) {
			output(i, tpi) = NV_EIGEN_Ith(y, i);
		}

		prevtime = timepoints[tpi];
	}
#endif

	step_counts[step_count]++;

	return true;
}

void CVODESolver::DumpStatistics(const char* filename)
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

int CVODESolver::cvode_rhs_fn(OdeReal t, OdeReal* y, OdeReal* ydot)
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

int CVODESolver::cvode_jac_fn(OdeReal t, N_Vector y, N_Vector ydot, SUNMatrix jac)
{
	if (!jacobian) {
		return cvode_jac_dq(t, y, ydot, jac);
	}
#if TODO
	bool result = jacobian(t, y, ydot, jac, user_data);
	if (result) {
		return 0;
	} else {
		return -1;
	}
#else
	return -1;
#endif
}

int CVODESolver::cvode_jac_dq(OdeReal t, N_Vector y, N_Vector ydot, SUNMatrix jac)
{
	static const realtype MIN_INC_MULT = 1000.0;
	static const realtype ZERO         = 0.0;
	static const realtype PT25         = 0.25;
	static const realtype ONE          = 1.0;

	int retval = 0;

	CVodeMem cv_mem = (CVodeMem)cvode_mem;
	CVLsMem cvls_mem = (CVLsMem)cv_mem->cv_lmem;

	/* Obtain pointers to the data for ewt, y */
	realtype* ewt_data = N_VGetArrayPointer(cv_mem->cv_ewt);
	realtype* y_data = NV_DATA_S(y);
	realtype* cns_data = nullptr;
	if (cv_mem->cv_constraintsSet) {
		cns_data = N_VGetArrayPointer(cv_mem->cv_constraints);
	}

	/* Set minimum increment based on uround and norm of f */
	realtype srur = SUNRsqrt(cv_mem->cv_uround);
	realtype fnorm = N_VWrmsNorm(ydot, cv_mem->cv_ewt);
	realtype minInc;
	if (fnorm != ZERO) {
		minInc = (MIN_INC_MULT * SUNRabs(cv_mem->cv_h) * cv_mem->cv_uround * N * fnorm);
	} else {
		minInc = ONE;
	}

	N_Vector jthCol = N_VCloneEmpty(y);

	//VectorReal jthCol;
	for (size_t j = 0; j < N; j++) {
		/* Generate the jth col of J(tn,y) */
		//N_VSetArrayPointer(SUNDenseMatrix_Column(jac, j), jthCol);

	}

#if 0
	/* access matrix dimension */
	N = SUNDenseMatrix_Rows(Jac);

	/* Rename work vector for readibility */
	ftemp = tmp1;

	/* Create an empty vector for matrix column calculations */
	jthCol = N_VCloneEmpty(tmp1);



	for (j = 0; j < N; j++) {

		/* Generate the jth col of J(tn,y) */
		N_VSetArrayPointer(SUNDenseMatrix_Column(Jac, j), jthCol);

		yjsaved = y_data[j];
		inc = SUNMAX(srur * SUNRabs(yjsaved), minInc / ewt_data[j]);

		/* Adjust sign(inc) if y_j has an inequality constraint. */
		if (cv_mem->cv_constraintsSet) {
			conj = cns_data[j];
			if (SUNRabs(conj) == ONE) { if ((yjsaved + inc) * conj < ZERO)  inc = -inc; } else if (SUNRabs(conj) == TWO) { if ((yjsaved + inc) * conj <= ZERO) inc = -inc; }
		}

		y_data[j] += inc;

		retval = cv_mem->cv_f(t, y, ftemp, cv_mem->cv_user_data);
		cvls_mem->nfeDQ++;
		if (retval != 0) break;

		y_data[j] = yjsaved;

		inc_inv = ONE / inc;
		N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

	}

	/* Destroy jthCol vector */
	N_VSetArrayPointer(NULL, jthCol);  /* SHOULDN'T BE NEEDED */
	N_VDestroy(jthCol);
#endif
	return(retval);
}

OdeReal CVODESolver::get_y(size_t i)
{
	return NV_Ith_S(this->y, i);
}

void CVODESolver::set_y(size_t i, OdeReal y)
{
	NV_Ith_S(this->y, i) = y;
}
