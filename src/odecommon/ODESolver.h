#pragma once

#include "LinearAlgebraSelector.h"

class ODESolver
{
public:
	typedef std::function<bool(OdeReal, const OdeReal*, OdeReal*, void*)> TDeriviativeFunction;
	typedef std::function<bool(OdeReal, const OdeReal*, const OdeReal*, OdeMatrixReal&, void*)> TJacobianFunction;
	typedef std::function<Real(OdeReal, void*)> TDiscontinuityCallback;
	typedef std::function<void(OdeReal, const OdeReal*, void*)> TIntegrationStepCallback;

	ODESolver();
	virtual ~ODESolver();

	virtual bool Initialize(size_t N, void* user);
	virtual bool SetSolverParameter(const std::string& parameter, int int_value, OdeReal real_value);
	void SetTolerance(OdeReal relative, OdeReal absolute);			// SetTolerance should be called after Initialize
	void SetTolerance(OdeReal relative, OdeVectorReal absolute);	// SetTolerance should be called after Initialize
	void SetUserData(void* user);
	void SetDiscontinuity(OdeReal time, TDiscontinuityCallback cb, void* user);
	void ResetDiscontinuity();
	void SetIntegrationStepCallback(TIntegrationStepCallback cb);

	void SetDerivativeFunction(TDeriviativeFunction f);
	void SetJacobianFunction(TJacobianFunction f);
	inline size_t GetNumVariables() const { return N; }

	bool SolveReturnSolution(const OdeVectorReal& initial_conditions, const OdeVectorReal* timepoints, OdeMatrixReal* output, bool verbose = false);
	bool SolveStoreIntegrationPoints(const OdeVectorReal& initial_conditions, Real end_time, bool verbose = false);

	// For use after SolveStoreIntegrationPoints
	virtual OdeReal GetInterpolatedY(OdeReal t, size_t i) = 0;

	// For use during callbacks
	virtual OdeReal get_current_y(size_t i) { return std::numeric_limits<OdeReal>::quiet_NaN(); }
	virtual void set_current_y(size_t i, OdeReal y) {}

protected:
	virtual bool Solve(const OdeVectorReal& initial_conditions, OdeReal end_time, bool do_interpolation, bool store_integration_points, bool verbose) = 0;

	// Settings
	size_t N;
	void* user_data;

	OdeReal next_discontinuity_time;
	TDiscontinuityCallback discontinuity_cb;
	void* discontinuity_user;
	TIntegrationStepCallback integration_step_cb;

	TDeriviativeFunction derivative;
	TJacobianFunction jacobian;

	OdeReal relative_tolerance;
	OdeVectorReal absolute_tolerance;

	// Runtime variables
	size_t interpolation_timepoints_start;
	const OdeVectorReal* interpolation_timepoints;
	OdeMatrixReal* interpolated_output;
};
