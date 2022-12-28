#include "Utils.h"
#include "ODESolver.h"

ODESolver::ODESolver()
	: user_data(nullptr)
	, discontinuity_time(std::numeric_limits<Real>::quiet_NaN())
	, discontinuity_user(nullptr)
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

bool ODESolver::SetTolerance(OdeReal relative, OdeReal absolute)
{
	abstol = absolute;
	reltol = relative;
	return true;
}

void ODESolver::SetUserData(void* user)
{
	user_data = user;
}

void ODESolver::SetDiscontinuity(OdeReal time, TDiscontinuityCallback cb, void* user)
{
	discontinuity_time = time;
	discontinuity_cb = cb;
	discontinuity_user = user;
}

void ODESolver::ResetDiscontinuity()
{
	discontinuity_time = std::numeric_limits<OdeReal>::quiet_NaN();
}

void ODESolver::SetDerivativeFunction(TDeriviativeFunction f)
{
	derivative = f;
}

void ODESolver::SetJacobianFunction(TJacobianFunction f)
{
	jacobian = f;
}
