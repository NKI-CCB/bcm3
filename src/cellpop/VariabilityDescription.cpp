#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "VariabilityDescription.h"

VariabilityDescription::Parameter::Parameter()
	: fixed_value(std::numeric_limits<Real>::quiet_NaN())
	, variable_ix(std::numeric_limits<size_t>::max())
	, non_sampled_variable_ix(std::numeric_limits<size_t>::max())
{
}

VariabilityDescription::VariabilityDescription()
	: type(EType::Invalid)
	, distribution(EDistribution::Invalid)
	, negate(false)
	, only_initial_cells(false)
{

}

VariabilityDescription::~VariabilityDescription()
{

}

std::unique_ptr<VariabilityDescription> VariabilityDescription::Create(const boost::property_tree::ptree& xml_node)
{
	std::unique_ptr<VariabilityDescription> desc;

	try {
		desc = std::make_unique<VariabilityDescription>();
		if (!desc->Load(xml_node)) {
			desc.reset();
		}
	} catch (boost::property_tree::ptree_error & e) {
		LOGERROR("Error parsing likelihood file: %s", e.what());
		desc.reset();
	}

	return desc;
}

bool VariabilityDescription::PostInitialize(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model)
{
	bool result = true;

	if (distribution != EDistribution::Bernoulli) {
		result &= InitializeParameter(range, varset, non_sampled_parameter_names, model);
	}
	if (distribution == EDistribution::StudentT || distribution == EDistribution::ProportionHalfT) {
		result &= InitializeParameter(dof, varset, non_sampled_parameter_names, model);
	}
	if (distribution == EDistribution::Bernoulli || distribution == EDistribution::MultipliedBernoulli || distribution == EDistribution::ProportionExponential || distribution == EDistribution::ProportionHalfT) {
		result &= InitializeParameter(proportion, varset, non_sampled_parameter_names, model);
	}

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

	if (!species_name.empty()) {
		size_t ix = model.GetCVodeSpeciesByName(species_name);
		if (ix == std::numeric_limits<size_t>::max()) {
			return false;
		}
	}

	return true;
}

bool VariabilityDescription::ApplyVariabilityEntryTime(Real& value, const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	if (entry_time) {
		// This variability description applies to entry time
		Real sobol_sample = sobol_sequence[sobol_sequence_ix++];
		if (!only_initial_cells || is_initial_cell) {
			Real q = DistributionQuantile(sobol_sample, transformed_values, non_sampled_parameters);
			ApplyType(value, q);
		}
		return true;
	} else {
		return false;
	}
}

bool VariabilityDescription::ApplyVariabilityParameter(const std::string& parameter, OdeReal& value, const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	if (!parameter_name.empty() && parameter == parameter_name) {
		// This variability description applies to this parameter
		Real sobol_sample = sobol_sequence[sobol_sequence_ix++];
		if (!only_initial_cells || is_initial_cell) {
			Real q = DistributionQuantile(sobol_sample, transformed_values, non_sampled_parameters);
			Real real_value = (Real)value;
			ApplyType(real_value, q);
			value = (OdeReal)real_value;
		}
		return true;
	} else {
		return false;
	}
}

bool VariabilityDescription::ApplyVariabilitySpecies(const std::string& species, OdeReal& value, const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const
{
	if (!species_name.empty() && species == species_name) {
		// This variability description applies to this species
		Real sobol_sample = sobol_sequence[sobol_sequence_ix++];
		if (!only_initial_cells || is_initial_cell) {
			Real q = DistributionQuantile(sobol_sample, transformed_values, non_sampled_parameters);
			Real real_value = (Real)value;
			ApplyType(real_value, q);
			real_value = fabs(real_value);
			value = (OdeReal)real_value;
		}
		return true;
	} else {
		return false;
	}
}

bool VariabilityDescription::Load(const boost::property_tree::ptree& xml_node)
{
	species_name = xml_node.get<std::string>("<xmlattr>.species", "");
	parameter_name = xml_node.get<std::string>("<xmlattr>.parameter", "");
	std::string entry_time_str = xml_node.get<std::string>("<xmlattr>.entry_time", "");
	entry_time = !entry_time_str.empty();

	int count = 0;
	if (!species_name.empty()) count++;
	if (!parameter_name.empty()) count++;
	if (!entry_time_str.empty()) count++;

	if (count == 0) {
		LOGERROR("Cell variability description has neither a species name nor a parameter name, nor entry_time");
		return false;
	} else if (count > 1) {
		LOGERROR("Cell variability description has both a species name and a parameter name or entry_time");
		return false;
	}

	if (entry_time) {
		only_initial_cells = xml_node.get<bool>("<xmlattr>.only_initial_cells", true);
	} else {
		only_initial_cells = xml_node.get<bool>("<xmlattr>.only_initial_cells", false);
	}

	std::string type_str = xml_node.get<std::string>("<xmlattr>.type");
	if (type_str == "additive") {
		type = EType::Additive;
	} else if (type_str == "multiplicative") {
		type = EType::Multiplicative;
	} else if (type_str == "multiplicative_log2") {
		type = EType::Multiplicative_Log2;
	} else if (type_str == "replace") {
		type = EType::Replace;
	} else {
		LOGERROR("Unknown variability description type \"%s\"", type_str.c_str());
		return false;
	}

	std::string distribution_str = xml_node.get<std::string>("<xmlattr>.distribution");
	if (distribution_str == "normal") {
		distribution = EDistribution::Normal;
	} else if (distribution_str == "half_normal") {
		distribution = EDistribution::HalfNormal;
	} else if (distribution_str == "student_t") {
		distribution = EDistribution::StudentT;
	} else if (distribution_str == "bernoulli") {
		distribution = EDistribution::Bernoulli;
	} else if (distribution_str == "multiplied_bernoulli") {
		distribution = EDistribution::MultipliedBernoulli;
	} else if (distribution_str == "uniform") {
		distribution = EDistribution::Uniform;
	} else if (distribution_str == "exponential") {
		distribution = EDistribution::Exponential;
	} else if (distribution_str == "propotion_exponential") {
		distribution = EDistribution::ProportionExponential;
	} else if (distribution_str == "propotion_half_t") {
		distribution = EDistribution::ProportionHalfT;
	} else {
		LOGERROR("Unknown distribution \"%s\" in variability description", distribution_str.c_str());
		return false;
	}

	if (distribution != EDistribution::Bernoulli) {
		range.str = xml_node.get<std::string>("<xmlattr>.range");
	}
	if (distribution == EDistribution::StudentT || distribution == EDistribution::ProportionHalfT) {
		dof.str = xml_node.get<std::string>("<xmlattr>.degrees_of_freedom_minus_two");
	}
	if (distribution == EDistribution::Bernoulli || distribution == EDistribution::MultipliedBernoulli || distribution == EDistribution::ProportionExponential || distribution == EDistribution::ProportionHalfT) {
		proportion.str = xml_node.get<std::string>("<xmlattr>.proportion");
	}

	negate = xml_node.get<bool>("<xmlattr>.negate", false);

	return true;
}

bool VariabilityDescription::InitializeParameter(Parameter& p, std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model)
{
	ASSERT(!p.str.empty()); // Missing variables should have been caught in VariabilityDescription::Load already

	p.variable_ix = varset->GetVariableIndex(p.str, false);
	if (p.variable_ix == std::numeric_limits<size_t>::max()) {
		// No, it ins't - is it a non-sampled parameter?
		auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), p.str);
		if (it != non_sampled_parameter_names.end()) {
			p.non_sampled_variable_ix = it - non_sampled_parameter_names.begin();
		} else {
			// No, it isn't - try to cast it to a Real
			try {
				p.fixed_value = boost::lexical_cast<Real>(p.str);
			} catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for parameter \"%s\", and could also not cast it to a constant real value: %s", p.str.c_str(), e.what());
				return false;
			}
		}
	}

	return true;
}

Real VariabilityDescription::GetParameterValue(const Parameter& p, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const
{
	if (p.variable_ix != std::numeric_limits<size_t>::max()) {
		return transformed_values[p.variable_ix];
	} else if (p.non_sampled_variable_ix != std::numeric_limits<size_t>::max()) {
		return non_sampled_parameters[p.non_sampled_variable_ix];
	} else {
		ASSERT(!std::isnan(p.fixed_value));
		return p.fixed_value;
	}
}

Real VariabilityDescription::DistributionQuantile(Real p, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const
{
	Real q = std::numeric_limits<Real>::quiet_NaN();

	switch (distribution) {
	case EDistribution::Normal:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			q = bcm3::QuantileNormal(p, 0, range_value);
		}
		break;

	case EDistribution::HalfNormal:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			q = bcm3::QuantileNormal(0.5 + p * 0.5, 0, range_value);
		}
		break;

	case EDistribution::StudentT:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			Real dof_value = GetParameterValue(dof, transformed_values, non_sampled_parameters) + 2.0;
			q = bcm3::QuantileT(p, 0, range_value, dof_value);
		}
		break;

	case EDistribution::Bernoulli:
		{
			Real proportion_value = GetParameterValue(proportion, transformed_values, non_sampled_parameters);
			if (p < proportion_value) {
				q = 0.0;
			} else {
				q = 1.0;
			}
		}
		break;

	case EDistribution::MultipliedBernoulli:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			Real proportion_value = GetParameterValue(proportion, transformed_values, non_sampled_parameters);
			if (p < proportion_value) {
				q = 0.0;
			} else {
				q = range_value;
			}
		}
		break;

	case EDistribution::Uniform:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			q = p * range_value;
		}
		break;

	case EDistribution::Exponential:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			q = bcm3::QuantileExponential(p, 1.0 / range_value);
		}
		break;

	case EDistribution::ProportionExponential:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			Real proportion_value = GetParameterValue(proportion, transformed_values, non_sampled_parameters);
			if (p < proportion_value) {
				q = bcm3::QuantileExponential(p / proportion_value, 1.0 / range_value);
			} else {
				q = 0.0;
			}
		}
		break;

		case EDistribution::ProportionHalfT:
		{
			Real range_value = GetParameterValue(range, transformed_values, non_sampled_parameters);
			Real dof_value = GetParameterValue(dof, transformed_values, non_sampled_parameters) + 2.0;
			Real proportion_value = GetParameterValue(proportion, transformed_values, non_sampled_parameters);
			if (p < proportion_value) {
				q = bcm3::QuantileT(0.5 + 0.5 * (p / proportion_value), 0, range_value, dof_value);
			} else {
				q = 0.0;
			}
		}
		break;

	default:
		ASSERT(false);
		break;
	}

	if (negate) {
		q = -q;
	}

	return q;
}

void VariabilityDescription::ApplyType(Real& x, Real& value) const
{
	switch (type) {
	case EType::Additive:
		x += value;
		break;

	case EType::Additive_Log2:
		x += pow(2.0, value);
		break;

	case EType::Multiplicative:
		x *= value;
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
