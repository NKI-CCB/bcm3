#pragma once

#include "VariableSet.h"

class ValueReference
{
public:
	ValueReference();

	void SetString(const std::string& str);
	bool Load(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names);
	Real GetValue(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const;

private:
	std::string value_name;
	Real fixed_value;
	size_t variable_ix;
	size_t non_sampled_variable_ix;
};
