#include "DataLikelihoodBase.h"
#include "DataLikelihoodDuration.h"
#include "DataLikelihoodTimeCourse.h"
#include "DataLikelihoodTimePoints.h"
#include "ProbabilityDistributions.h"

DataLikelihoodBase::DataLikelihoodBase()
	: weight(1.0)
	, error_model(ErrorModel::StudentT4)
{
}

DataLikelihoodBase::~DataLikelihoodBase()
{
}

std::unique_ptr<DataLikelihoodBase> DataLikelihoodBase::Create(const boost::property_tree::ptree& xml_node, size_t parallel_evaluations)
{
	std::string type = xml_node.get<std::string>("<xmlattr>.type", "time_course");
	std::unique_ptr<DataLikelihoodBase> dl;
	if (type == "time_course") {
		dl = std::make_unique<DataLikelihoodTimeCourse>(parallel_evaluations);
	} else if (type == "time_points") {
		dl = std::make_unique<DataLikelihoodTimePoints>(parallel_evaluations);
	} else if (type == "duration") {
		dl = std::make_unique<DataLikelihoodDuration>(parallel_evaluations);
	} else {
		LOGERROR("Unknown data likelihood type \"%s\"", type.c_str());
	}
	return dl;
}

bool DataLikelihoodBase::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm)
{
	data_name = xml_node.get<std::string>("<xmlattr>.data_name");
	
	weight = xml_node.get<Real>("<xmlattr>.weight", 1.0);

	stdev_str = xml_node.get<std::string>("<xmlattr>.stdev");
	offset_str = xml_node.get<std::string>("<xmlattr>.offset", "");
	scale_str = xml_node.get<std::string>("<xmlattr>.scale", "");

	return true;
}

bool DataLikelihoodBase::PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	bool result = true;

	std::vector<std::string> stdev_tokens;
	bcm3::tokenize(stdev_str, stdev_tokens, ";");
	stdev_ix.resize(stdev_tokens.size(), std::numeric_limits<size_t>::max());
	non_sampled_stdev_ix.resize(stdev_tokens.size(), std::numeric_limits<size_t>::max());
	fixed_stdev_value.resize(stdev_tokens.size(), std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < stdev_tokens.size(); i++) {
		result &= ParseString(stdev_tokens[i], stdev_ix[i], non_sampled_stdev_ix[i], fixed_stdev_value[i], varset, non_sampled_parameter_names);
	}

	if (!offset_str.empty()) {
		std::vector<std::string> offset_tokens;
		bcm3::tokenize(offset_str, offset_tokens, ";");
		offset_ix.resize(offset_tokens.size(), std::numeric_limits<size_t>::max());
		non_sampled_offset_ix.resize(offset_tokens.size(), std::numeric_limits<size_t>::max());
		fixed_offset_value.resize(offset_tokens.size(), std::numeric_limits<size_t>::max());
		for (size_t i = 0; i < offset_tokens.size(); i++) {
			result &= ParseString(offset_tokens[i], offset_ix[i], non_sampled_offset_ix[i], fixed_offset_value[i], varset, non_sampled_parameter_names);
		}
	}

	if (!scale_str.empty()) {
		std::vector<std::string> scale_tokens;
		bcm3::tokenize(scale_str, scale_tokens, ";");
		scale_ix.resize(scale_tokens.size(), std::numeric_limits<size_t>::max());
		non_sampled_scale_ix.resize(scale_tokens.size(), std::numeric_limits<size_t>::max());
		fixed_scale_value.resize(scale_tokens.size(), std::numeric_limits<size_t>::max());
		for (size_t i = 0; i < scale_tokens.size(); i++) {
			result &= ParseString(scale_tokens[i], scale_ix[i], non_sampled_scale_ix[i], fixed_scale_value[i], varset, non_sampled_parameter_names);
		}
	}

	return result;
}

Real DataLikelihoodBase::GetCurrentSTDev(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, size_t i)
{
	size_t ix;
	if (stdev_ix.size() == 1) {
		ix = 0;
	} else if (i < stdev_ix.size()) {
		ix = i;
	} else {
		LOGERROR("Out of bounds");
		return std::numeric_limits<Real>::quiet_NaN();
	}

	Real stdev;
	if (stdev_ix[ix] != std::numeric_limits<size_t>::max()) {
		stdev = transformed_values[stdev_ix[ix]];
	} else if (non_sampled_stdev_ix[ix] != std::numeric_limits<size_t>::max()) {
		stdev = non_sampled_parameters[non_sampled_stdev_ix[ix]];
	} else {
		stdev = fixed_stdev_value[ix];
	}
	return stdev;
}

Real DataLikelihoodBase::GetCurrentDataOffset(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, size_t i)
{
	size_t ix;
	if (offset_ix.size() == 0) {
		return 0.0;
	} else if (offset_ix.size() == 1) {
		ix = 0;
	} else if (i < offset_ix.size()) {
		ix = i;
	} else {
		LOGERROR("Out of bounds");
		return std::numeric_limits<Real>::quiet_NaN();
	}

	Real data_offset = 0.0;
	if (offset_ix[ix] != std::numeric_limits<size_t>::max()) {
		data_offset = transformed_values[offset_ix[ix]];
	} else if (non_sampled_offset_ix[ix] != std::numeric_limits<size_t>::max()) {
		data_offset = non_sampled_parameters[non_sampled_offset_ix[ix]];
	} else if (fixed_offset_value[ix] == fixed_offset_value[ix]) {
		data_offset = fixed_offset_value[ix];
	}
	return data_offset;
}

Real DataLikelihoodBase::GetCurrentDataScale(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, size_t i)
{
	size_t ix;
	if (scale_ix.size() == 0) {
		return 1.0;
	} else if (scale_ix.size() == 1) {
		ix = 0;
	} else if (i < scale_ix.size()) {
		ix = i;
	} else {
		LOGERROR("Out of bounds");
		return std::numeric_limits<Real>::quiet_NaN();
	}

	Real data_scale = 1.0;
	if (scale_ix[ix] != std::numeric_limits<size_t>::max()) {
		data_scale = transformed_values[scale_ix[ix]];
	} else if (non_sampled_scale_ix[ix] != std::numeric_limits<size_t>::max()) {
		data_scale = non_sampled_parameters[non_sampled_scale_ix[ix]];
	} else if (fixed_scale_value[ix] == fixed_scale_value[ix]) {
		data_scale = fixed_scale_value[ix];
	}
	return data_scale;
}

bool DataLikelihoodBase::ParseString(std::string str, size_t& var_ix, size_t& non_sampled_var_ix, Real& fixed_value, const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	var_ix = varset.GetVariableIndex(str, false);
	if (var_ix == std::numeric_limits<size_t>::max()) {
		// No, it ins't - is it a non-sampled parameter?
		auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), str);
		if (it != non_sampled_parameter_names.end()) {
			non_sampled_var_ix = it - non_sampled_parameter_names.begin();
		} else {
			// No, it isn't - try to cast it to a Real
			try {
				fixed_value = boost::lexical_cast<Real>(str);
			}
			catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for data parameter \"%s\", and could also not cast it to a constant real value: %s", str.c_str(), e.what());
				return false;
			}
		}
	}
	return true;
}
