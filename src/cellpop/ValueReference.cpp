#include "Utils.h"
#include "ValueReference.h"

ValueReference::ValueReference()
	: fixed_value(std::numeric_limits<Real>::quiet_NaN())
	, variable_ix(std::numeric_limits<size_t>::max())
	, non_sampled_variable_ix(std::numeric_limits<size_t>::max())
{
}

void ValueReference::SetString(const std::string& str)
{
	value_name = str;
}

bool ValueReference::Load(std::shared_ptr<const bcm3::VariableSet> varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	if (value_name.empty()) {
		LOGERROR("Value reference is empty");
		return false;
	
	}
	// See if it references a sampled variable.
	variable_ix = varset->GetVariableIndex(value_name, false);
	if (variable_ix == std::numeric_limits<size_t>::max()) {
		// No, it isn't - is it a non-sampled parameter?
		auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), p.str);
		if (it != non_sampled_parameter_names.end()) {
			non_sampled_variable_ix = it - non_sampled_parameter_names.begin();
		} else {
			// No, it isn't - try to cast it to a Real
			try {
				p.fixed_value = boost::lexical_cast<Real>(p.str);
			}
			catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for parameter \"%s\", and also could not cast it to a constant real value: %s", value_name.c_str(), e.what());
				return false;
			}
		}
	}

	return true;
}

Real ValueReference::GetValue(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters) const
{
	if (variable_ix != std::numeric_limits<size_t>::max()) {
		return transformed_values[p.variable_ix];
	} else if (non_sampled_variable_ix != std::numeric_limits<size_t>::max()) {
		return non_sampled_parameters[p.non_sampled_variable_ix];
	} else {
		ASSERT(!std::isnan(fixed_value));
		return fixed_value;
	}
}
