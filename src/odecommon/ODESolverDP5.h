#pragma once

#include "ODESolver.h"

class ODESolverDP5 : public ODESolver
{
public:
	ODESolverDP5();
	~ODESolverDP5();

	virtual bool Initialize(size_t N, void* user);
	virtual bool Restart();
	virtual bool SetSolverParameter(const std::string& parameter, int int_value, OdeReal real_value);

	virtual void RestartInterpolationIteration();
	virtual const OdeVectorReal& GetInterpolatedY(OdeReal t);
	virtual OdeReal GetInterpolatedY(OdeReal t, size_t i);

	virtual OdeVectorReal get_current_y() const;
	virtual OdeReal get_current_y(size_t i);
	virtual void set_current_y(size_t i, OdeReal y);
	virtual Real get_threshold_crossing_time(size_t species_ix, Real threshold, bool above, Real prev_time) const;

private:
	virtual bool Solve(const OdeVectorReal& initial_conditions, OdeReal end_time, bool do_interpolation, bool store_integration_points, bool verbose);
	OdeReal ApplyRK(OdeReal t, OdeReal cur_dt);

	// Settings
	OdeReal min_dt;
	OdeReal max_dt;
	int max_steps;

	// Runtime variables
	OdeReal* k[7];
	OdeReal* ytmp;
	OdeReal* yn;
	OdeVectorReal interpolated_y;

	static bool warn_discontinuity;
};
