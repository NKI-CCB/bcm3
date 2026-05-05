#include "Utils.h"
#include "ODESolverCVODE.h"
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_math.h>

#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_impl.h"
#include "../../dependencies/cvode-5.3.0/src/cvode/cvode_ls_impl.h"

static const size_t MAX_CVODE_STEPS = 10000;

void static_cvode_err_fn(int error_code, const char* module, const char* function, char* msg, void* user_data)
{
}

void static_cvode_err_fn_log(int error_code, const char *module, const char *function, char *msg, void *user_data)
{
	LOGERROR("CVODE error %d in module %s, function %s: %s", error_code, module, function, msg);
}

int static_cvode_rhs_fn(OdeReal t, N_Vector y, N_Vector ydot, void* user_data)
{
	ODESolverCVODE* solver = reinterpret_cast<ODESolverCVODE*>(user_data);
	return solver->cvode_rhs_fn(t, NV_DATA_S(y), NV_DATA_S(ydot));
}

#if CVODE_USE_EIGEN_SOLVER
int static_cvode_jac_fn(OdeReal t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector ytmp1, N_Vector ytmp2, N_Vector ytmp3)
{
	ODESolverCVODE* solver = reinterpret_cast<ODESolverCVODE*>(user_data);
	return solver->cvode_jac_fn(t, *NV_CONTENT_S(y), *NV_CONTENT_S(fy), EIGMAT(Jac));
}
#endif

ODESolverCVODE::CVodeTimepoint::CVodeTimepoint()
	: cvode_time(0.0)
	, cv_uround(std::numeric_limits<Real>::quiet_NaN())
	, cv_tn(std::numeric_limits<Real>::quiet_NaN())
	, cv_h(std::numeric_limits<Real>::quiet_NaN())
	, cv_hu(std::numeric_limits<Real>::quiet_NaN())
	, cv_q(0)
{
}

ODESolverCVODE::ODESolverCVODE()
	: max_steps(2000)
	, cvode_mem(NULL)
	, LS(NULL)
	, NLS(NULL)
	, y(NULL)
	, tmpvector(NULL)
	, J(NULL)
	, t(0.0)
	, current_step(0)
	, cvode_timepoint_iter(0)
	, interpolation_time(std::numeric_limits<Real>::quiet_NaN())
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

bool ODESolverCVODE::Initialize(size_t N, void* user, size_t store_integration_points_buffer)
{
	if (!ODESolver::Initialize(N, user, store_integration_points_buffer)) {
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

#if CVODE_USE_EIGEN_SOLVER
	CVodeSetJacFn(cvode_mem, &static_cvode_jac_fn);
	if (!jacobian) {
		// Difference quotient approximation requires some work space
		y_copy = OdeVectorReal(N);
		work_vector = OdeVectorReal(N);
	}
#else
	if (jacobian) {
		LOGERROR("Custom jacobian is not implemented for solver other than the Eigen solver");
		return false;
	} else {
		// Standard CVODE solver has a built-in jacobian approximation, so no need to specify.
	}	
#endif

	cvode_timepoints.resize(store_integration_points_buffer);
	for (int i = 0; i < 6; i++) {
		cvode_timepoints_zn[i].resize(N, store_integration_points_buffer);
	}
	interpolation_y.resize(N);

	return true;
}

bool ODESolverCVODE::Restart()
{
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

void ODESolverCVODE::RestartInterpolationIteration()
{
	interpolation_y.setZero();
	cvode_timepoint_iter = 0;
	interpolation_time = std::numeric_limits<Real>::quiet_NaN();
}

const OdeVectorReal& ODESolverCVODE::GetInterpolatedY(OdeReal t)
{
	if (t == interpolation_time) {
		// We still have this interpolation timepoint
		return interpolation_y;
	}

	// Requests for interpolations should be in strictly increasing time, so we can continue searching where we left off
	while (cvode_timepoint_iter < current_step && cvode_timepoint_iter < cvode_timepoints.size()) {
		const CVodeTimepoint& cvt = cvode_timepoints[cvode_timepoint_iter];
		if (cvt.cvode_time > t) {
			break;
		} else {
			cvode_timepoint_iter++;
		}
	}
	if (cvode_timepoint_iter == current_step || cvode_timepoint_iter == cvode_timepoints.size()) {
		// This timepoint is beyond the solver's time
		interpolation_time = std::numeric_limits<Real>::quiet_NaN();
		interpolation_y.setConstant(std::numeric_limits<Real>::quiet_NaN());
		return interpolation_y;
	}

	interpolation_time = t;

	// Calculate the interpolated value
	// This is based on CVodeGetDky with k=0; we've stored the relevant part of cvode_mem in the CVodeTimepoint struct
	const CVodeTimepoint& cvt = cvode_timepoints[cvode_timepoint_iter];

	/* Allow for some slack */
	Real tfuzz = 100.0 * cvt.cv_uround * (fabs(cvt.cv_tn) + fabs(cvt.cv_hu));
	if (cvt.cv_hu < 0.0) {
		tfuzz = -tfuzz;
	}
	Real tp = cvt.cv_tn - cvt.cv_hu - tfuzz;
	Real tn1 = cvt.cv_tn + tfuzz;
	if ((t - tp) * (t - tn1) > 0.0) {
		LOGERROR("Time error for interpolation; %g - %zu - %g", t, cvode_timepoint_iter, cvt.cv_tn);
		interpolation_y.setConstant(std::numeric_limits<Real>::quiet_NaN());
		return interpolation_y;
	}

	/* Sum the differentiated interpolating polynomial */
	CVodeMem cv_mem = (CVodeMem)cvode_mem;
	Real s = (t - cvt.cv_tn) / cvt.cv_h;
	interpolation_y.setZero();
	for (int j = cvt.cv_q; j >= 0; j--) {
		Real cval = 1.0;
		for (int i = j; i >= j + 1; i--)
			cval *= i;
		for (int i = 0; i < j; i++)
			cval *= s;

		interpolation_y += cval * cvode_timepoints_zn[j].col(cvode_timepoint_iter);
	}

	return interpolation_y;
}

OdeReal ODESolverCVODE::GetInterpolatedY(OdeReal t, size_t i)
{
	return GetInterpolatedY(t)(i);
}

OdeVectorReal ODESolverCVODE::get_current_y() const
{
	return EIGV(this->y);
}

OdeReal ODESolverCVODE::get_current_y(size_t i) const
{
	return NV_Ith_S(this->y, i);
}

void ODESolverCVODE::set_current_y(size_t i, OdeReal y) const
{
	NV_Ith_S(this->y, i) = y;
}

Real ODESolverCVODE::get_threshold_crossing_time(size_t species_ix, Real threshold, bool above, Real prev_time) const
{
	if (current_step == 0 || current_step > cvode_timepoints.size()) {
		// This shouldn't happen.
		ASSERT(false);
		LOGERROR("Integration error.");
		return prev_time;
	}

	size_t cvode_step_ix = current_step - 1;
	const CVodeTimepoint& cvt = cvode_timepoints[cvode_step_ix];

	Real dt = (t - prev_time) * 0.5;
	Real time = prev_time + dt;

	bool crossed = true;
	for (int iter = 0; iter < 10; iter++) {
		/* Allow for some slack */
		Real tfuzz = 100.0 * cvt.cv_uround * fabs(cvt.cv_tn) + fabs(cvt.cv_hu);
		if (cvt.cv_hu < 0.0) {
			tfuzz = -tfuzz;
		}
		Real tp = cvt.cv_tn - cvt.cv_hu - tfuzz;
		Real tn1 = cvt.cv_tn + tfuzz;

		/* Sum the differentiated interpolating polynomial */
		CVodeMem cv_mem = (CVodeMem)cvode_mem;
		Real s = (time - cvt.cv_tn) / cvt.cv_h;
		Real x = 0.0;
		for (int j = cvt.cv_q; j >= 0; j--) {
			Real cval = 1.0;
			for (int i = j; i >= j + 1; i--)
				cval *= i;
			for (int i = 0; i < j; i++)
				cval *= s;

			x += cval * cvode_timepoints_zn[j](species_ix, cvode_step_ix);
		}

		dt *= 0.5;
		if (above) {
			if (x > threshold) {
				time -= dt;
			} else {
				time += dt;
			}
		} else {
			if (x < threshold) {
				time -= dt;
			} else {
				time += dt;
			}
		}
	}

	return time;
}

bool ODESolverCVODE::Solve(const OdeVectorReal& initial_conditions, OdeReal end_time, bool do_interpolation, bool store_integration_points, bool verbose)
{
	for (size_t i = 0; i < N; i++) {
		NV_Ith_S(y, i) = initial_conditions[i];
		NV_Ith_S(tmpvector, i) = absolute_tolerance(i);
	}

	if (verbose) {
		CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn_log, this);
	} else {
		CVodeSetErrHandlerFn(cvode_mem, &static_cvode_err_fn, this);
	}

	CVodeSVtolerances(cvode_mem, relative_tolerance, tmpvector);
	CVodeReInit(cvode_mem, 0.0, (N_Vector)y);

	if (!std::isnan(next_discontinuity_time)) {
		CVodeSetStopTime(cvode_mem, (OdeReal)next_discontinuity_time);
	}

	// Simulate the system one step at a time
	current_step = 0;
	min_used_step_size = std::numeric_limits<Real>::infinity();
	t = 0.0;
	size_t tpi = interpolation_timepoints_start;
	while (1) {
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
			num_steps_used = current_step;
			return false;
		} else {
			min_used_step_size = std::min(min_used_step_size, tret - t);

			t = tret;

			if (store_integration_points) {
				if (current_step >= cvode_timepoints.size()) {
					LOGERROR("CVODE integration timepoint storage buffer too small");
					return false;
				}
				CVodeTimepoint& cvt = cvode_timepoints[current_step];

				// Copy relevant variables from cvode_mem
				CVodeMem cv_mem = (CVodeMem)cvode_mem;

				cvt.cvode_time = t;
				cvt.cv_uround = cv_mem->cv_uround;
				cvt.cv_tn = cv_mem->cv_tn;
				cvt.cv_h = cv_mem->cv_h;
				cvt.cv_hu = cv_mem->cv_hu;
				cvt.cv_q = cv_mem->cv_q;

				ASSERT(cvt.cv_q <= 5);
				for (int j = cv_mem->cv_q; j >= 0; j--) {
#if CVODE_USE_EIGEN_SOLVER
					cvode_timepoints_zn[j].col(current_step) = *NV_CONTENT_S(cv_mem->cv_zn[j]);
#else
					for (int i = 0; i < NV_LENGTH_S(cv_mem->cv_zn[j]); i++) {
						cvode_timepoints_zn[j](i, cvode_steps) = NV_Ith_S(cv_mem->cv_zn[j], i);
					}
#endif
				}
			}

			if (do_interpolation) {
				// If we pass a timepoint, interpolate back to that timepoint
				while (tret >= (*interpolation_timepoints)(tpi)) {
					result = CVodeGetDky(cvode_mem, (*interpolation_timepoints)(tpi), 0, (N_Vector)tmpvector);
					if (result != CV_SUCCESS) {
						if (verbose) {
							LOG("Interpolation failed.");
						}
						num_steps_used = current_step;
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

			current_step++;

			if (integration_step_cb) {
				bool continue_integration = integration_step_cb(t, NV_DATA_S(y), user_data);
				if (!continue_integration) {
					break;
				}
			}

			if (t >= end_time) {
				break;
			}

			if (current_step == max_steps) {
				if (verbose) {
					LOG("Maximum number of time steps reached.");
				}
				num_steps_used = current_step;
				return false;
			}

			if (!std::isnan(next_discontinuity_time) && (result == CV_TSTOP_RETURN || next_discontinuity_time == t)) {
				next_discontinuity_time = discontinuity_cb(t, discontinuity_user);
				if (!std::isnan(next_discontinuity_time) && next_discontinuity_time < std::numeric_limits<Real>::infinity()) {
					CVodeReInit(cvode_mem, t, (N_Vector)y);
					CVodeSetStopTime(cvode_mem, next_discontinuity_time);
				} else {
					// Continue simulation but no new discontinuity
					CVodeReInit(cvode_mem, t, (N_Vector)y);
				}
			}
		}
	}

	num_steps_used = current_step;
	return true;
}

int ODESolverCVODE::cvode_rhs_fn(OdeReal t, const OdeReal* y, OdeReal* ydot)
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

#if CVODE_USE_EIGEN_SOLVER

int ODESolverCVODE::cvode_jac_fn(OdeReal t, const OdeVectorReal& y, const OdeVectorReal& ydot, OdeMatrixReal& jac)
{
	bool result;
	if (jacobian) {
		result = jacobian(t, y.data(), ydot.data(), jac, user_data);
	} else {
		result = DifferenceQuotientJacobian(t, y, ydot, jac);
	}
	if (result) {
		return 0;
	} else {
		return -1;
	}
}

bool ODESolverCVODE::DifferenceQuotientJacobian(OdeReal t, const OdeVectorReal& y, const OdeVectorReal& ydot, OdeMatrixReal& jac)
{
	// Adapted from cvLsDenseDQJac

	static const OdeReal MIN_INC_MULT = OdeReal(1000.0);

	CVodeMem cv_mem = reinterpret_cast<CVodeMem>(cvode_mem);
	CVLsMem cvls_mem = reinterpret_cast<CVLsMem>(cv_mem->cv_lmem);

	y_copy = y;

	// Set minimum increment based on uround and norm of f
	OdeReal srur = SUNRsqrt(cv_mem->cv_uround);
	OdeVectorReal& ewt = *NV_CONTENT_S(cv_mem->cv_ewt);
	OdeReal fnorm = sqrt((ydot.array() * ewt.array()).square().sum() / N);

	OdeReal minInc;
	if (fnorm != (OdeReal)0.0) {
		minInc = (MIN_INC_MULT * fabs(cv_mem->cv_h) * cv_mem->cv_uround * N * fnorm);
	} else {
		minInc = 1.0;
	}

	for (int j = 0; j < N; j++) {
		OdeReal inc = std::max(srur * fabs(y[j]), minInc / NV_Ith_S(cv_mem->cv_ewt, j));

		y_copy[j] += inc;

		int retval = cvode_rhs_fn(t, y_copy.data(), work_vector.data());
		cvls_mem->nfeDQ++;
		if (retval != 0) {
			return false;
		}

		y_copy[j] = y[j];

		OdeReal inc_inv = (OdeReal)1.0 / inc;
		jac.col(j) = inc_inv * (work_vector - ydot);
	}

	return true;
}

#endif
