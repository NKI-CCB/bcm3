#include "Utils.h"
#include "ProbabilityDistributions.h"
#include "VariabilityDescription.h"

VariabilityDescription::VariabilityDescription()
	: type(EType::Invalid)
	, distribution(EDistribution::Invalid)
	, fixed_range_value(std::numeric_limits<Real>::quiet_NaN())
	, range_ix(std::numeric_limits<size_t>::max())
	, non_sampled_range_ix(std::numeric_limits<size_t>::max())
	, fixed_dof_value(std::numeric_limits<Real>::quiet_NaN())
	, dof_ix(std::numeric_limits<size_t>::max())
	, non_sampled_dof_ix(std::numeric_limits<size_t>::max())
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
	if (!range_str.empty()) {
		range_ix = varset->GetVariableIndex(range_str, false);
		if (range_ix == std::numeric_limits<size_t>::max()) {
			// No, it ins't - is it a non-sampled parameter?
			auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), range_str);
			if (it != non_sampled_parameter_names.end()) {
				non_sampled_range_ix = it - non_sampled_parameter_names.begin();
			} else {
				// No, it isn't - try to cast it to a Real
				try {
					fixed_range_value = boost::lexical_cast<Real>(range_str);
				} catch (const boost::bad_lexical_cast & e) {
					LOGERROR("Could not find variable for range parameter \"%s\", and could also not cast it to a constant real value: %s", range_str.c_str(), e.what());
					return false;
				}
			}
		}
	}
	if (!dof_str.empty()) {
		dof_ix = varset->GetVariableIndex(dof_str, false);
		if (dof_ix == std::numeric_limits<size_t>::max()) {
			// No, it ins't - is it a non-sampled parameter?
			auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), dof_str);
			if (it != non_sampled_parameter_names.end()) {
				non_sampled_dof_ix = it - non_sampled_parameter_names.begin();
			} else {
				// No, it isn't - try to cast it to a Real
				try {
					fixed_dof_value = boost::lexical_cast<Real>(dof_str);
				} catch (const boost::bad_lexical_cast& e) {
					LOGERROR("Could not find variable for degrees of freedom parameter \"%s\", and could also not cast it to a constant real value: %s", dof_str.c_str(), e.what());
					return false;
				}
			}
		}
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
			Real range = GetRangeValue(transformed_values, non_sampled_parameters);
			Real dof = GetDOFValue(transformed_values, non_sampled_parameters);
			Real q = DistributionQuantile(sobol_sample, range, dof);
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
			Real range = GetRangeValue(transformed_values, non_sampled_parameters);
			Real dof = GetDOFValue(transformed_values, non_sampled_parameters);
			Real q = DistributionQuantile(sobol_sample, range, dof);
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
			Real range = GetRangeValue(transformed_values, non_sampled_parameters);
			Real dof = GetDOFValue(transformed_values, non_sampled_parameters);
			Real q = DistributionQuantile(sobol_sample, range, dof);
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
	} else if (distribution_str == "uniform") {
		distribution = EDistribution::Uniform;
	} else {
		LOGERROR("Unknown distribution \"%s\" in variability description", distribution_str.c_str());
		return false;
	}

	range_str = xml_node.get<std::string>("<xmlattr>.range");
	dof_str = xml_node.get<std::string>("<xmlattr>.degrees_of_freedom", "0.0");

	if (distribution == EDistribution::StudentT && dof_str == "0.0") {
		LOGERROR("Student T-distribution with fixed 0 degrees of freedom; likely because it was not specified - please specify degrees_of_freedom in the variability description");
		return false;
	}

	return true;
}

Real VariabilityDescription::GetRangeValue(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const
{
	if (range_ix != std::numeric_limits<size_t>::max()) {
		return transformed_values[range_ix];
	} else if (non_sampled_range_ix != std::numeric_limits<size_t>::max()) {
		return non_sampled_parameters[non_sampled_range_ix];
	} else {
		ASSERT(fixed_range_value != std::numeric_limits<Real>::quiet_NaN());
		return fixed_range_value;
	}
}

Real VariabilityDescription::GetDOFValue(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const
{
	if (dof_ix != std::numeric_limits<size_t>::max()) {
		return transformed_values[dof_ix];
	} else if (non_sampled_dof_ix != std::numeric_limits<size_t>::max()) {
		return non_sampled_parameters[non_sampled_dof_ix];
	} else {
		ASSERT(fixed_dof_value != std::numeric_limits<Real>::quiet_NaN());
		return fixed_dof_value;
	}
}

Real VariabilityDescription::DistributionQuantile(Real p, Real range, Real dof) const
{
	Real q = std::numeric_limits<Real>::quiet_NaN();

	switch (distribution) {
	case EDistribution::Normal:
		q = bcm3::QuantileNormal(p, 0, range);
		break;

	case EDistribution::HalfNormal:
		q = bcm3::QuantileNormal(0.5 + p * 0.5, 0, range);
		break;

	case EDistribution::StudentT:
		q = bcm3::QuantileT(p, 0, range, dof);
		break;

	case EDistribution::Bernoulli:
		if (p < range) {
			q = 0.0;
		} else {
			q = 1.0;
		}
		break;

	case EDistribution::Uniform:
		q = p * range;
		break;

	default:
		ASSERT(false);
		break;
	}

	return q;
}

void VariabilityDescription::ApplyType(Real& x, Real& value) const
{
	switch (type) {
	case EType::Additive:
		x += value;
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
