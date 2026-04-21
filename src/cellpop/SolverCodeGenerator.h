#pragma once

#include "SBMLModel.h"

extern "C" {
    typedef void (*derivative_fn)(OdeReal* out, const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
    typedef void (*jacobian_fn)(OdeMatrixReal& out, const OdeReal* species, const OdeReal* constant_species, const OdeReal* parameters, const OdeReal* non_sampled_parameters);
}

class Experiment;

class SolverCodeGenerator {
public:
    SolverCodeGenerator();
    ~SolverCodeGenerator();

    bool GenerateAndCompileSolverCode(Experiment* experiment, SBMLModel& cell_model, const std::string& codegen_name);

	derivative_fn GetGeneratedDerivative() const { return derivative; }
	jacobian_fn GetGeneratedJacobian() const { return jacobian; }

private:
#if PLATFORM_WINDOWS
    HMODULE derivative_dll;
#else
    void* derivative_dll;
#endif

    derivative_fn derivative;
    jacobian_fn jacobian;
};
