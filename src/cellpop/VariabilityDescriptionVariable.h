#pragma once

#include "SBMLModel.h"
#include "ValueReference.h"
#include "VariableSet.h"

#include <boost/random/sobol.hpp>

class VariabilityDescriptionVariable
{
public:
	VariabilityDescriptionVariable();
	~VariabilityDescriptionVariable();

	static std::unique_ptr<VariabilityDescriptionVariable> Create(const boost::property_tree::ptree& xml_node);
	bool PostInitialize(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names, const SBMLModel& model);

	bool ApplyVariabilityEntryTime(Real& value, Real pseudorandom_value, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;
	bool ApplyVariabilityParameter(const std::string& parameter, OdeReal& value, Real pseudorandom_value, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;
	bool ApplyVariabilityInitialCondition(const std::string& species, OdeReal& value, Real pseudorandom_value, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, bool is_initial_cell) const;

	const ValueReference& GetScaleReference() const { return scale; }

private:
	bool Load(const boost::property_tree::ptree& xml_node);
	Real TransformValue(Real value, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const;
	void Apply(Real& x, Real& value) const;

	enum class EType {
		Additive,
		Additive_Log,
		Additive_Log2,
		Multiplicative,
		Multiplicative_Log,
		Multiplicative_Log2,
		Replace,
		Invalid,
	};

	std::string initial_condition_species_name;
	std::string parameter_name;
	bool entry_time;

	EType apply_type;

	ValueReference scale;
	bool negate;
	bool only_initial_cells;
};
