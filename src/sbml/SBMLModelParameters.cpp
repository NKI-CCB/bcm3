#include "Utils.h"
#include "SBMLModel.h"
#include "SBMLModelParameters.h"
#include "SBMLSpecies.h"

SBMLModelParameters::SBMLModelParameters(SBMLModel* model)
	: Model(model)
	, Variables(NULL)
{
}

SBMLModelParameters::~SBMLModelParameters()
{
}

bool SBMLModelParameters::SetVariableSet(std::shared_ptr<const bcm3::VariableSet> variables)
{
	Variables = variables;
	return true;
}

void SBMLModelParameters::ResetForcedVariables()
{
	ForcedVariableValues.clear();
}

bool SBMLModelParameters::SetForcedVariable(const std::string& variable_name, Real value)
{
	if (variable_name.substr(0, 14) == "initial_value_") {
		std::string varname = variable_name.substr(14);
		for (size_t i = 0; i < Model->GetNumSimulatedSpecies(); i++) {
			if (Model->GetSimulatedSpeciesFullName(i) == varname) {
				ForcedInitialConditions[i] = value;
				return true;
			}
		}
	}

	size_t ix = Variables->GetVariableIndex(variable_name);
	if (ix != std::numeric_limits<size_t>::max()) {
		ForcedVariableValues[ix] = value;
		return true;
	}

	LOGERROR("Force variable of unknown parameter \"%s\"", variable_name.c_str());
	return false;
}

void SBMLModelParameters::UpdateVariableValues(Real* variables) const
{
	for (std::map<size_t, Real>::const_iterator fvvi = ForcedVariableValues.begin(); fvvi != ForcedVariableValues.end(); ++fvvi) {
		variables[fvvi->first] = fvvi->second;
	}
}

void SBMLModelParameters::UpdateInitialValues(Real* initial_values) const
{
	for (std::map<size_t, Real>::const_iterator fici = ForcedInitialConditions.begin(); fici != ForcedInitialConditions.end(); ++fici) {
		initial_values[fici->first] = fici->second;
	}
}
