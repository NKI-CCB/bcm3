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
	struct Parameter {
		Parameter();
		std::string str;
		Real fixed_value;
		size_t variable_ix;
		size_t non_sampled_variable_ix;
	};

	bool Load(const boost::property_tree::ptree& xml_node);
	bool InitializeParameter(Parameter& p, std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model);
	Real GetParameterValue(const Parameter& p, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const;
	Real DistributionQuantile(Real p, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const;
	void ApplyType(Real& x, Real& value) const;

	enum class EType {
		Additive,
		Additive_Log2,
		Multiplicative,
		Multiplicative_Log2,
		Replace,

		Invalid,
	};
	enum class EDistribution {
		Normal,
		HalfNormal,
		StudentT,
		Bernoulli,
		MultipliedBernoulli,
		Uniform,
		Exponential,
		ProportionExponential,
		ProportionHalfT,

		Invalid
	};

	std::string species_name;
	std::string parameter_name;
	bool entry_time;

	EType type;
	EDistribution distribution;

	Parameter range;
	Parameter dof;
	Parameter proportion;

	bool only_initial_cells;
};
