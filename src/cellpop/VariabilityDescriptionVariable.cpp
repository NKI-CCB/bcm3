#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "VariabilityDescriptionVariable.h"

VariabilityDescriptionVariable::VariabilityDescriptionVariable()
	: apply_type(EType::Invalid)
	, entry_time(false)
	, negate(false)
	, only_initial_cells(false)
{
}

VariabilityDescriptionVariable::~VariabilityDescriptionVariable()
{
}

std::unique_ptr<VariabilityDescriptionVariable> VariabilityDescriptionVariable::Create(const boost::property_tree::ptree& xml_node)
{
	std::unique_ptr<VariabilityDescriptionVariable> desc;

	try {
		desc = std::make_unique<VariabilityDescriptionVariable>();
		if (!desc->Load(xml_node)) {
			desc.reset();
		}
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		desc.reset();
	}

	return desc;
}

bool VariabilityDescriptionVariable::PostInitialize(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model)
{
	bool result = true;

	// Load the scale value reference
	result &= scale.Load(varset, non_sampled_parameter_names, model);

	// If the variability is on a parameter, get the parameter index
	if (!parameter_name.empty()) {
		// Check whether the parameter is sampled
		size_t ix = varset->GetVariableIndex(parameter_name, false);
		if (ix == std::numeric_limits<size_t>::max()) {
			// Is it non-sampled?
			auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), parameter_name);
			if (it == non_sampled_parameter_names.end()) {
				LOGERROR("Variability has been specified for parameter \"%s\", but the parameter is not find in either the sampled parameters or the non-sampled parameters", parameter_name.c_str());
				return false;
			}
		}
	}

	// If the variability is on an initial condition, get the species index
	if (!initial_condition_species_name.empty()) {
		size_t ix = model.GetODEIntegratedSpeciesByName(initial_condition_species_name);
		if (ix == std::numeric_limits<size_t>::max()) {
			return false;
		}
	}

	return true;
}

bool VariabilityDescriptionVariable::ApplyVariabilityEntryTime(Real& value, const VectorReal& pseudorandom_values, int& pseudorandom_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	if (entry_time) {
		// This variability description applies to entry time
		Real sobol_sample = pseudorandom_values[pseudorandom_ix++];
		if (!only_initial_cells || is_initial_cell) {
			Real variability = 0.0;
			variability = TransformValue(variability, transformed_values, non_sampled_parameters);
			Apply(value, variability);
		}
		return true;
	} else {
		return false;
	}
}

bool VariabilityDescriptionVariable::ApplyVariabilityParameter(const std::string& parameter, OdeReal& value, const VectorReal& pseudorandom_values, int& pseudorandom_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	if (!parameter_name.empty() && parameter == parameter_name) {
		// This variability description applies to this parameter
		Real sobol_sample = pseudorandom_values[pseudorandom_ix++];
		if (!only_initial_cells || is_initial_cell) {
			Real variability = 0.0;
			variability = TransformValue(variability, transformed_values, non_sampled_parameters);
			Real real_value = (Real)value;
			Apply(real_value, variability);
			value = (OdeReal)real_value;
		}
		return true;
	} else {
		return false;
	}
}

bool VariabilityDescriptionVariable::ApplyVariabilityInitialConditionSpecies(const std::string& species, OdeReal& value, const VectorReal& pseudorandom_values, int& pseudorandom_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	if (!initial_condition_species_name.empty() && species == initial_condition_species_name) {
		// This variability description applies to this species
		Real sobol_sample = pseudorandom_values[pseudorandom_ix++];
		if (!only_initial_cells || is_initial_cell) {
			Real variability = 0.0;
			variability = TransformValue(variability, transformed_values, non_sampled_parameters);
			Real real_value = (Real)value;
			Apply(real_value, variability);
			value = (OdeReal)real_value;
		}
		return true;
	} else {
		return false;
	}
}

bool VariabilityDescriptionVariable::Load(const boost::property_tree::ptree& xml_node)
{
	initial_condition_species_name = xml_node.get<std::string>("<xmlattr>.initial_condition_species", "");
	parameter_name = xml_node.get<std::string>("<xmlattr>.model_parameter", "");
	std::string entry_time_str = xml_node.get<std::string>("<xmlattr>.entry_time", "");
	entry_time = !entry_time_str.empty();

	int count = 0;
	if (!initial_condition_species_name.empty()) count++;
	if (!parameter_name.empty()) count++;
	if (!entry_time_str.empty()) count++;

	if (count == 0) {
		LOGERROR("Cell variability description has neither a initial_condition_species name nor a model_parameter name, nor entry_time; exactly one of these should be specified");
		return false;
	} else if (count > 1) {
		LOGERROR("Cell variability description has both a species name and a parameter name or entry_time; only one of these should be specified");
		return false;
	}

	if (entry_time) {
		only_initial_cells = xml_node.get<bool>("<xmlattr>.only_initial_cells", true);
	} else {
		only_initial_cells = xml_node.get<bool>("<xmlattr>.only_initial_cells", false);
	}

	std::string apply_str = xml_node.get<std::string>("<xmlattr>.apply");
	if (apply_str == "additive") {
		apply_type = EType::Additive;
	} else if (apply_str == "multiplicative") {
		apply_type = EType::Multiplicative;
	} else if (apply_str == "multiplicative_log2") {
		apply_type = EType::Multiplicative_Log2;
	} else if (apply_str == "replace") {
		apply_type = EType::Replace;
	} else {
		LOGERROR("Unknown variability application type \"%s\"", apply_str.c_str());
		return false;
	}
	
	scale.SetString(xml_node.get<std::string>("<xmlattr>.scale"));
	negate = xml_node.get<bool>("<xmlattr>.negate", false);

	return true;
}

Real VariabilityDescriptionVariable::TransformValue(Real value, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const
{
	Real scale_value = scale.GetValue(transformed_values, non_sampled_parameters);
	value *= scale_value;

	if (negate) {
		value = -value;
	}

	return value;
}

void VariabilityDescriptionVariable::Apply(Real& x, Real& value) const
{
	switch (apply_type) {
	case EType::Additive:
		x += value;
		break;

	case EType::Additive_Log:
		x += exp(value);
		break;

	case EType::Additive_Log2:
		x += pow(2.0, value);
		break;

	case EType::Multiplicative:
		x *= value;
		break;

	case EType::Multiplicative_Log:
		x *= exp(value);
		break;

	case EType::Multiplicative_Log2:
		x *= pow(2.0, value);
		break;

	case EType::Replace:
		x = value;
		break;

	default:
		ASSERT(false);
		break;
	}
}
