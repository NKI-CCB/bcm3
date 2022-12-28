#pragma once

#include "VariableSet.h"

class ExperimentalCondition;
class SBMLModel;

class SBMLModelParameters
{
public:
	SBMLModelParameters(SBMLModel* model);
	~SBMLModelParameters();

	bool SetVariableSet(std::shared_ptr<const bcm3::VariableSet> variables);
	void ResetForcedVariables();
	bool SetForcedVariable(const std::string& variable_name, Real value);

	void UpdateVariableValues(Real* variables) const;
	void UpdateInitialValues(OdeReal* initial_values) const;

private:
	const SBMLModel* Model;
	std::shared_ptr<const bcm3::VariableSet> Variables;

	// Forcing variable to specific values
	std::map<size_t, Real> ForcedVariableValues;
	std::map<size_t, Real> ForcedInitialConditions;
};
