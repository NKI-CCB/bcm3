#pragma once

#include "ODESolver.h"

class ODESolverDP5 : public ODESolver
{
public:
	ODESolverDP5();
	~ODESolverDP5();

	virtual bool Initialize(size_t N, void* user);
	virtual bool Simulate(const Real* initial_conditions, const VectorReal& timepoints, MatrixReal& output, bool verbose = false);
	virtual OdeReal get_y(size_t i);
	virtual void set_y(size_t i, OdeReal y);

private:
	OdeReal ApplyRK(OdeReal t, OdeReal cur_dt);

	OdeReal* k[7];
	OdeReal* ytmp;
	OdeReal* yn;
};
