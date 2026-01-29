#include "Utils.h"
#include "Correlation.h"
#include "DataLikelihoodTimeCoursePopulationAverage.h"
#include "Experiment.h"
#include "ProbabilityDistributions.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/algorithm/string.hpp>

DataLikelihoodTimeCourseBase::DataLikelihoodTimeCourseBase()
{
}

DataLikelihoodTimeCourseBase::~DataLikelihoodTimeCourseBase()
{
}

bool DataLikelihoodTimeCourseBase::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm)
{
	if (!DataLikelihoodBase::Load(xml_node, experiment, varset, data_file, vm)) {
		return false;
	}

	bool result = true;

	std::string species_name = xml_node.get<std::string>("<xmlattr>.species_name");
	use_log_ratio = xml_node.get<bool>("<xmlattr>.use_log_ratio", false);
	optimize_offset_scale = xml_node.get<bool>("<xmlattr>.optimize_offset_scale", false);
	include_only_cells_that_went_through_mitosis = xml_node.get<bool>("<xmlattr>.include_only_cells_that_went_through_mitosis", false);
	missing_simulation_time_stdev_str = xml_node.get<std::string>("<xmlattr>.missing_simulation_time_stdev", "");

	saturation_scale_str = xml_node.get<std::string>("<xmlattr>.saturation_scale", "");

	if (optimize_offset_scale) {
		if (error_model != ErrorModel::Normal) {
			LOGWARNING("Using offset/scale optimization without a normal error model - this likely leads to strong convergence problems!");
		}

		optimize_offset_min = xml_node.get<Real>("<xmlattr>.optimize_offset_min", -1.0);
		optimize_offset_max = xml_node.get<Real>("<xmlattr>.optimize_offset_max", 1.0);
		optimize_scale_min = xml_node.get<Real>("<xmlattr>.optimize_scale_min", 0.1);
		optimize_scale_max = xml_node.get<Real>("<xmlattr>.optimize_scale_max", 10.0);
	}

	std::string time_dimension_name;
	size_t num_dimensions = 0;
	result &= data_file.GetDimensionName(experiment->GetName(), data_name, 0, time_dimension_name);
	result &= data_file.GetDimensionCount(experiment->GetName(), data_name, &num_dimensions);

	// Retrieve time points from the datafile
	size_t num_timepoints;
	result &= data_file.GetDimensionSize(experiment->GetName(), time_dimension_name, &num_timepoints);
	result &= data_file.GetValues(experiment->GetName(), time_dimension_name, 0, num_timepoints, timepoints);

	if (!result) {
		return false;
	}

	// See if the species reference is a sum of species
	if (species_name.find_first_of(';') != std::string::npos) {
		bcm3::tokenize(species_name, species_names, ";");
		for (size_t i = 0; i < species_names.size(); i++) {
			boost::trim(species_names[i]);
		}
	} else {
		species_names.resize(1, species_name);
	}

	return result;
}

bool DataLikelihoodTimeCourseBase::PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	if (!DataLikelihoodBase::PostInitialize(varset, non_sampled_parameter_names)) {
		return false;
	}

	if (!missing_simulation_time_stdev_str.empty()) {
		fixed_missing_simulation_time_stdev_ix = varset.GetVariableIndex(missing_simulation_time_stdev_str, false);
		if (fixed_missing_simulation_time_stdev_ix == std::numeric_limits<size_t>::max()) {
			// No, it ins't - is it a non-sampled parameter?
			auto it = std::find(non_sampled_parameter_names.begin(), non_sampled_parameter_names.end(), missing_simulation_time_stdev_str);
			if (it != non_sampled_parameter_names.end()) {
				fixed_missing_simulation_time_stdev_non_sampled_ix = it - non_sampled_parameter_names.begin();
			} else {
				// No, it isn't - try to cast it to a Real
				try {
					fixed_missing_simulation_time_stdev = boost::lexical_cast<Real>(missing_simulation_time_stdev_str);
				} catch (const boost::bad_lexical_cast& e) {
					LOGERROR("Could not find variable for data missing simulation time stdev parameter \"%s\", and could also not cast it to a constant real value: %s", missing_simulation_time_stdev_str.c_str(), e.what());
					return false;
				}
			}
		}
	}

	if (!saturation_scale_str.empty()) {
		use_signal_saturation = true;

		saturation_scale_ix = varset.GetVariableIndex(saturation_scale_str, false);
		if (saturation_scale_ix == std::numeric_limits<size_t>::max()) {
			// Try to cast it to a Real
			try {
				saturation_scale = boost::lexical_cast<Real>(saturation_scale_str);
			} catch (const boost::bad_lexical_cast& e) {
				LOGERROR("Could not find variable for saturation scale \"%s\", and could also not cast it to a constant real value: %s", saturation_scale_str.c_str(), e.what());
				return false;
			}
		}
	}

	stdevs.resize(species_names.size(), std::numeric_limits<Real>::quiet_NaN());
	proportional_stdevs.resize(species_names.size(), std::numeric_limits<Real>::quiet_NaN());
	data_offsets.resize(species_names.size(), std::numeric_limits<Real>::quiet_NaN());
	data_scales.resize(species_names.size(), std::numeric_limits<Real>::quiet_NaN());
	minus_log_sigma.resize(species_names.size(), std::numeric_limits<Real>::quiet_NaN());
	inv_two_sigma_sq.resize(species_names.size(), std::numeric_limits<Real>::quiet_NaN());

	return true;
}

bool DataLikelihoodTimeCourseBase::RequestSimulationInfo(Experiment* experiment, ESynchronizeCellTrajectory synchronize)
{
	for (size_t j = 0; j < species_names.size(); j++) {
		const std::string& species_name = species_names[j];

		if (use_log_ratio) {
			if (species_name.find_first_of('/') == std::string::npos) {
				LOGERROR("use_log_ratio is specified as true, but the species_name does not contain a division");
				return false;
			}
		}

		if (species_name.find_first_of('+') != std::string::npos) {
			std::vector<std::string> sum_names;
			bcm3::tokenize(species_name, sum_names, "+");
			for (size_t i = 0; i < sum_names.size(); i++) {
				boost::trim(sum_names[i]);
				size_t species_ix = experiment->GetCVodeSpeciesByName(sum_names[i]);
				if (species_ix == std::numeric_limits<size_t>::max()) {
					return false;
				}
				if (species_map[species_ix].empty()) {
					for (size_t i = 0; i < timepoints.size(); i++) {
						experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, synchronize);
					}
				}
				species_map[species_ix].push_back(j);
			}
		} else if (species_name.find_first_of('/') != std::string::npos) {
			if (!use_log_ratio) {
				LOGERROR("simulated species reference has a division, but use_log_ratio has not been specified; only log ratios are supported for now");
				return false;
			}

			std::vector<std::string> sum_names;
			bcm3::tokenize(species_name, sum_names, "/");
			if (sum_names.size() != 2) {
				LOGERROR("only division of exactly two species is supported");
				return false;
			}
			for (size_t i = 0; i < sum_names.size(); i++) {
				boost::trim(sum_names[i]);
				size_t species_ix = experiment->GetCVodeSpeciesByName(sum_names[i]);
				if (species_ix == std::numeric_limits<size_t>::max()) {
					return false;
				}
				if (species_map[species_ix].empty()) {
					for (size_t i = 0; i < timepoints.size(); i++) {
						experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, synchronize);
					}
				}
				species_map[species_ix].push_back(j);
				if (i == 0) {
					log_ratio_numerator_ix.resize(j + 1);
					log_ratio_numerator_ix[j] = species_ix;
				}
			}
		} else {
			size_t species_ix = experiment->GetCVodeSpeciesByName(species_names[j]);
			if (species_ix == std::numeric_limits<size_t>::max()) {
				return false;
			}
			if (species_map[species_ix].empty()) {
				for (size_t i = 0; i < timepoints.size(); i++) {
					experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, synchronize);
				}
			}
			species_map[species_ix].push_back(j);
		}
	}

	return true;
}

bool DataLikelihoodTimeCourseBase::PrepateEvaluation(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp)
{
	for (size_t i = 0; i < species_names.size(); i++) {
		stdevs[i] = GetCurrentSTDev(transformed_values, non_sampled_parameters, i);
		proportional_stdevs[i] = GetCurrentProportionalSTDev(transformed_values, non_sampled_parameters, i);
		data_offsets[i] = GetCurrentDataOffset(transformed_values, non_sampled_parameters, i);
		data_scales[i] = GetCurrentDataScale(transformed_values, non_sampled_parameters, i);
	}

	missing_simulation_time_stdev = 1.0;
	if (fixed_missing_simulation_time_stdev_ix != std::numeric_limits<size_t>::max()) {
		missing_simulation_time_stdev = transformed_values[fixed_missing_simulation_time_stdev_ix];
	}
	else if (fixed_missing_simulation_time_stdev_non_sampled_ix != std::numeric_limits<size_t>::max()) {
		missing_simulation_time_stdev = non_sampled_parameters[fixed_missing_simulation_time_stdev_non_sampled_ix];
	}
	else if (!std::isnan(fixed_missing_simulation_time_stdev)) {
		missing_simulation_time_stdev = fixed_missing_simulation_time_stdev;
	}

	saturation_scale = saturation_scale_fixed_value;
	if (saturation_scale_ix != std::numeric_limits<size_t>::max()) {
		saturation_scale = transformed_values[saturation_scale_ix];
	}

	// Precalculate
	if (error_model == ErrorModel::Normal) {
		// Can precalculate minus_log_sigma and inv_two_sigma_sq for faster evaluation
		for (size_t i = 0; i < species_names.size(); i++) {
			minus_log_sigma[i] = -log(stdevs[i]);
			inv_two_sigma_sq[i] = 1.0 / (2.0 * stdevs[i] * stdevs[i]);
		}
	} else if (error_model == ErrorModel::StudentT4) {
		// Can precalculate logC for faster evaluation
	}

	return true;
}

Real DataLikelihoodTimeCourseBase::EvaluateValue(Real simulated, Real observed, size_t stdev_ix)
{
	Real logp;

	switch (error_model) {
		case ErrorModel::Normal:
			{
				Real d = observed - simulated;
				logp = minus_log_sigma[stdev_ix] - 0.91893853320467274178032973640562 - d * d * inv_two_sigma_sq[stdev_ix];
			}
			break;

		case ErrorModel::ProportionalNormal:
			logp = bcm3::LogPdfNormal(observed, simulated, proportional_stdevs[stdev_ix] * std::max(simulated, 0.0));
			break;

		case ErrorModel::AdditiveProportionalNormal:
			logp = bcm3::LogPdfNormal(observed, simulated, stdevs[stdev_ix] + proportional_stdevs[stdev_ix] * std::max(simulated, 0.0));
			break;

		case ErrorModel::StudentT4:
			logp = bcm3::LogPdfTnu4(observed, simulated, stdevs[stdev_ix]);
			break;

		default:
			assert(false);
			logp = std::numeric_limits<Real>::quiet_NaN();
			break;
	}

	return logp;
}

Real DataLikelihoodTimeCourseBase::EvaluateMissingValue(Real time_offset)
{
	Real logp;

	if (error_model == ErrorModel::Normal || error_model == ErrorModel::ProportionalNormal || error_model == ErrorModel::AdditiveProportionalNormal) {
		logp = bcm3::LogPdfNormal(time_offset, 0, missing_simulation_time_stdev);
	} else if (error_model == ErrorModel::StudentT4) {
		logp = bcm3::LogPdfTnu4(time_offset, 0, missing_simulation_time_stdev);
	} else {
		assert(false);
		logp = std::numeric_limits<Real>::quiet_NaN();
	}

	return logp;
}

void DataLikelihoodTimeCourseBase::OptimizeOffsetScale(const MatrixReal& observed, const MatrixReal& simulated, int col_ix, Real& offset, Real& scale) const
{
	bcm3::linear_regress_columns(simulated.col(col_ix), observed.col(col_ix), offset, scale);
	scale = std::min(std::max(scale, optimize_scale_min), optimize_scale_max);
	offset = std::min(std::max(offset, optimize_offset_min), optimize_offset_max);
}
