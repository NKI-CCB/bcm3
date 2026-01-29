#include "Utils.h"
#include "ODESolver.h"

ODESolver::ODESolver()
	: user_data(nullptr)
	, N(0)
	, next_discontinuity_time(std::numeric_limits<Real>::quiet_NaN())
	, discontinuity_user(nullptr)
	, relative_tolerance(std::numeric_limits<Real>::quiet_NaN())
{
}

ODESolver::~ODESolver()
{
}

bool ODESolver::Initialize(size_t N, void* user)
{
	if (N > std::numeric_limits<size_t>::max()) {
		LOGERROR("System too large");
		return false;
	}
	this->N = N;
	user_data = user;

	return true;
}

bool ODESolver::SetSolverParameter(const std::string& parameter, int int_value, OdeReal real_value)
{
	LOGWARNING("Unhandled ODE solver parameter %s - %d / %g", parameter.c_str(), int_value, real_value);
	return false;
}

void ODESolver::SetTolerance(OdeReal relative, OdeReal absolute)
{
	ASSERT(N != 0);
	absolute_tolerance.setConstant(N, absolute);
	relative_tolerance = relative;
}

void ODESolver::SetTolerance(OdeReal relative, OdeVectorReal absolute)
{
	ASSERT(N != 0);
	ASSERT(absolute.size() == N);
	absolute_tolerance = absolute;
	relative_tolerance = relative;
}

void ODESolver::SetUserData(void* user)
{
	user_data = user;
}

void ODESolver::SetDiscontinuity(OdeReal time, TDiscontinuityCallback cb, void* user)
{
	next_discontinuity_time = time;
	discontinuity_cb = cb;
	discontinuity_user = user;
}

void ODESolver::ResetDiscontinuity()
{
	next_discontinuity_time = std::numeric_limits<OdeReal>::quiet_NaN();
}

void ODESolver::SetIntegrationStepCallback(TIntegrationStepCallback cb)
{
	integration_step_cb = cb;
}

void ODESolver::SetDerivativeFunction(TDeriviativeFunction f)
{
	derivative = f;
}

void ODESolver::SetJacobianFunction(TJacobianFunction f)
{
	jacobian = f;
}

bool ODESolver::SolveReturnSolution(const OdeVectorReal& initial_conditions, const OdeVectorReal* timepoints, OdeMatrixReal* output, bool verbose)
{
	if (timepoints == nullptr || timepoints->size() == 0) {
		LOGERROR("No timepoints requested");
		return false;
	}
	if (output == nullptr) {
		LOGERROR("No output matrix provided");
		return false;
	}

	output->resize(N, timepoints->size());

	size_t ti = 0;
	while ((*timepoints)(ti) < std::numeric_limits<OdeReal>::epsilon()) {
		for (size_t i = 0; i < N; i++) {
			(*output)(i, ti) = initial_conditions[i];
		}
		ti++;
		if (ti == timepoints->size()) {
			return true;
		}
	}

	OdeReal end_time = (*timepoints)[timepoints->size() - 1];

	interpolation_timepoints_start = ti;
	interpolation_timepoints = timepoints;
	interpolated_output = output;

	if (!Solve(initial_conditions, end_time, true, false, verbose)) {
		return false;
	}

	return true;
}

bool ODESolver::SolveStoreIntegrationPoints(const OdeVectorReal& initial_conditions, Real end_time, bool verbose)
{
	if (end_time <= std::numeric_limits<OdeReal>::epsilon()) {
		LOGERROR("end_time should be larger than zero");
		return false;
	}

	if (!Solve(initial_conditions, end_time, false, true, verbose)) {
		return false;
	}

	return true;
}
