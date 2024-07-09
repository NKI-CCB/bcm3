#include "DataLikelihoodBase.h"
#include "DataLikelihoodDuration.h"
#include "DataLikelihoodTimeCourse.h"
#include "DataLikelihoodTimePoints.h"
#include "ProbabilityDistributions.h"

DataLikelihoodBase::DataLikelihoodBase()
	: weight(1.0)
	, error_model(ErrorModel::StudentT4)
	, offset_ix(std::numeric_limits<size_t>::max())
	, scale_ix(std::numeric_limits<size_t>::max())
	, stdev_ix(std::numeric_limits<size_t>::max())
	, non_sampled_offset_ix(std::numeric_limits<size_t>::max())
	, non_sampled_scale_ix(std::numeric_limits<size_t>::max())
	, non_sampled_stdev_ix(std::numeric_limits<size_t>::max())
	, fixed_offset_value(std::numeric_limits<Real>::quiet_NaN())
	, fixed_scale_value(std::numeric_limits<Real>::quiet_NaN())
	, fixed_stdev_value(std::numeric_limits<Real>::quiet_NaN())
	, stdev_multiplication_factor(1.0)
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
	stdev_multiplication_factor = xml_node.get<Real>("<xmlattr>.stdev_multiplication_factor", 1.0);

	// Is the stdev a parameter in the prior?
	stdev_str = xml_node.get<std::string>("<xmlattr>.stdev");
	offset_str = xml_node.get<std::string>("<xmlattr>.offset", "");
	scale_str = xml_node.get<std::string>("<xmlattr>.scale", "");

	return true;
}

bool DataLikelihoodBase::PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	stdev_ix = varset.GetVariableIndex(stdev_str, false);
	if (stdev_ix == std::numeric_limits<size_t>::max()) {
		// No, it ins't - is it a non-sampled parameter?
		auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), stdev_str);
		if (it != non_sampled_parameter_names.end()) {
			non_sampled_stdev_ix = it - non_sampled_parameter_names.begin();
		} else {
			// No, it isn't - try to cast it to a Real
			try {
				fixed_stdev_value = boost::lexical_cast<Real>(stdev_str);
			} catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for data stdev parameter \"%s\", and could also not cast it to a constant real value: %s", stdev_str.c_str(), e.what());
				return false;
			}
		}
	}
	if (!offset_str.empty()) {
		offset_ix = varset.GetVariableIndex(offset_str, false);
		if (offset_ix == std::numeric_limits<size_t>::max()) {
			// No, it ins't - is it a non-sampled parameter?
			auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), offset_str);
			if (it != non_sampled_parameter_names.end()) {
				non_sampled_offset_ix = it - non_sampled_parameter_names.begin();
			} else {
				// No, it isn't - try to cast it to a Real
				try {
					fixed_offset_value = boost::lexical_cast<Real>(offset_str);
				} catch (const boost::bad_lexical_cast& e) {
					LOGERROR("Could not find variable for data offset parameter \"%s\", and could also not cast it to a constant real value: %s", offset_str.c_str(), e.what());
					return false;
				}
			}
		}
	}
	if (!scale_str.empty()) {
		scale_ix = varset.GetVariableIndex(scale_str, false);
		if (scale_ix == std::numeric_limits<size_t>::max()) {
			// No, it ins't - is it a non-sampled parameter?
			auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), scale_str);
			if (it != non_sampled_parameter_names.end()) {
				non_sampled_scale_ix = it - non_sampled_parameter_names.begin();
			} else {
				// No, it isn't - try to cast it to a Real
				try {
					fixed_scale_value = boost::lexical_cast<Real>(scale_str);
				} catch (const boost::bad_lexical_cast& e) {
					LOGERROR("Could not find variable for data scale parameter \"%s\", and could also not cast it to a constant real value: %s", scale_str.c_str(), e.what());
					return false;
				}
			}
		}
	}
	return true;
}

Real DataLikelihoodBase::GetCurrentSTDev(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters)
{
	Real stdev;
	if (stdev_ix != std::numeric_limits<size_t>::max()) {
		stdev = transformed_values[stdev_ix];
	} else if (non_sampled_stdev_ix != std::numeric_limits<size_t>::max()) {
		stdev = non_sampled_parameters[non_sampled_stdev_ix];
	} else {
		stdev = fixed_stdev_value;
	}
	stdev *= stdev_multiplication_factor;
	stdev += 1e-4;
	return stdev;
}

Real DataLikelihoodBase::GetCurrentDataOffset(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters)
{
	Real data_offset = 0.0;
	if (offset_ix != std::numeric_limits<size_t>::max()) {
		data_offset = transformed_values[offset_ix];
	} else if (non_sampled_offset_ix != std::numeric_limits<size_t>::max()) {
		data_offset = non_sampled_parameters[non_sampled_offset_ix];
	} else if (fixed_offset_value == fixed_offset_value) {
		data_offset = fixed_offset_value;
	}
	return data_offset;
}

Real DataLikelihoodBase::GetCurrentDataScale(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters)
{
	Real data_scale = 1.0;
	if (scale_ix != std::numeric_limits<size_t>::max()) {
		data_scale = transformed_values[scale_ix];
	} else if (non_sampled_scale_ix != std::numeric_limits<size_t>::max()) {
		data_scale = non_sampled_parameters[non_sampled_scale_ix];
	} else if (fixed_scale_value == fixed_scale_value) {
		data_scale = fixed_scale_value;
	}
	return data_scale;
}

