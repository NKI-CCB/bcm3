#pragma once

#include "LinearAlgebraSelector.h"

class ODESolver
{
public:
	typedef std::function<bool(OdeReal, const OdeReal*, OdeReal*, void*)> TDeriviativeFunction;
	typedef std::function<bool(OdeReal, const OdeReal*, const OdeReal*, OdeReal**, void*)> TJacobianFunction;
	typedef std::function<Real(OdeReal, void*)> TDiscontinuityCallback;

	ODESolver();
	virtual ~ODESolver();

	virtual bool Initialize(size_t N, void* user);
	bool SetTolerance(OdeReal relative, OdeReal absolute);
	void SetUserData(void* user);
	void SetDiscontinuity(OdeReal time, TDiscontinuityCallback cb, void* user);
	void ResetDiscontinuity();

	void SetDerivativeFunction(TDeriviativeFunction f);
	void SetJacobianFunction(TJacobianFunction f);

	virtual bool Simulate(const OdeReal* initial_conditions, const OdeVectorReal& timepoints, OdeMatrixReal& output, bool verbose = false) = 0;
	virtual OdeReal get_y(size_t i) { return std::numeric_limits<OdeReal>::quiet_NaN(); }
	virtual void set_y(size_t i, OdeReal y) {}

protected:
	size_t N;
	void* user_data;
	OdeReal discontinuity_time;
	TDiscontinuityCallback discontinuity_cb;
	void* discontinuity_user;
	TDeriviativeFunction derivative;
	TJacobianFunction jacobian;

	OdeReal abstol;
	OdeReal reltol;
};
