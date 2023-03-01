#pragma once

#include "SBMLModel.h"
#include "VariableSet.h"

#include <boost/random/sobol.hpp>

class VariabilityDescription
{
public:
	VariabilityDescription();
	~VariabilityDescription();

	static std::unique_ptr<VariabilityDescription> Create(const boost::property_tree::ptree& xml_node);
	bool PostInitialize(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model);

	bool ApplyVariabilityEntryTime(Real& value, const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;
	bool ApplyVariabilityParameter(const std::string& parameter, OdeReal& value, const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;
	bool ApplyVariabilitySpecies(const std::string& species, OdeReal& value, const VectorReal& sobol_sequence, int& sobol_sequence_ix, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;

private:
	bool Load(const boost::property_tree::ptree& xml_node);
	Real GetRangeValue(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const;
	Real DistributionQuantile(Real p, Real range) const;
	void ApplyType(Real& x, Real& value) const;

	enum class EType {
		Additive,
		Multiplicative,
		Multiplicative_Log2,
		Replace,

		Invalid,
	};
	enum class EDistribution {
		Normal,
		HalfNormal,
		Bernoulli,
		Uniform,

		Invalid
	};

	std::string species_name;
	std::string parameter_name;
	bool entry_time;

	EType type;
	EDistribution distribution;
	
	std::string range_str;
	Real fixed_range_value;
	size_t range_ix;
	size_t non_sampled_range_ix;
	bool only_initial_cells;
};
