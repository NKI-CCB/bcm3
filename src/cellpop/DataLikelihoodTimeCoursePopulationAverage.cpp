#include "Utils.h"
#include "Correlation.h"
#include "DataLikelihoodTimeCoursePopulationAverage.h"
#include "Experiment.h"
#include "ProbabilityDistributions.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/algorithm/string.hpp>

DataLikelihoodTimeCoursePopulationAverage::DataLikelihoodTimeCoursePopulationAverage(size_t parallel_evaluations)
{
	if (parallel_evaluations > 1) {
		parallel_population_averages.resize(parallel_evaluations);
	}
}

DataLikelihoodTimeCoursePopulationAverage::~DataLikelihoodTimeCoursePopulationAverage()
{
}

bool DataLikelihoodTimeCoursePopulationAverage::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm)
{
	if (!DataLikelihoodTimeCourseBase::Load(xml_node, experiment, varset, data_file, vm)) {
		return false;
	}

	bool result = true;

	relative_to_time_average = xml_node.get<bool>("<xmlattr>.relative_to_time_average", false);

	std::string use_only_cell_ix_str = vm["cellpop.use_only_cell_ix"].as<std::string>();
	std::vector<std::string> use_only_cell_ix_tokens;
	bcm3::tokenize(use_only_cell_ix_str, use_only_cell_ix_tokens, ",");

	size_t num_dimensions = 0;
	result &= data_file.GetDimensionCount(experiment->GetName(), data_name, &num_dimensions);

	size_t num_replicates;
	std::string replicate_dimension_name;
	if (num_dimensions == 1) {
		LOG("Time course population likelihood for data %s; data has 1 dimension; assuming no replicates", data_name.c_str());
		num_replicates = 1;
	} else if (num_dimensions == 2) {
		result &= data_file.GetDimensionName(experiment->GetName(), data_name, 1, replicate_dimension_name);
		result &= data_file.GetDimensionSize(experiment->GetName(), replicate_dimension_name, &num_replicates);
		LOG("Time course population likelihood for data %s; data has 2 dimensions; assuming the second dimension (\"%s\") are replicate observations (%zu replicates).", data_name.c_str(), replicate_dimension_name.c_str(), num_replicates);
	} else {
		LOGERROR("Time course population likelihood for data %s; data has %zu dimensions but can only handle 1 or 2 dimensions.", data_name.c_str(), num_dimensions);
		return false;
	}

	observed_data = MatrixReal(num_replicates, timepoints.size());
	for (size_t ri = 0; ri < num_replicates; ri++) {
		VectorReal od;
		result &= data_file.GetValues(experiment->GetName(), data_name, ri, timepoints.size(), od);
		observed_data.row(ri) = od;
	}

	population_average.setZero(timepoints.size(), species_names.size());
	for (size_t i = 0; i < parallel_population_averages.size(); i++) {
		parallel_population_averages[i].setZero(timepoints.size(), species_names.size());
	}

	result &= RequestSimulationInfo(experiment, ESynchronizeCellTrajectory::None);

	return result;
}

bool DataLikelihoodTimeCoursePopulationAverage::PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	if (!DataLikelihoodTimeCourseBase::PostInitialize(varset, non_sampled_parameter_names)) {
		return false;
	}

	return true;
}

void DataLikelihoodTimeCoursePopulationAverage::Reset()
{
	population_average.setZero();
	for (size_t i = 0; i < parallel_population_averages.size(); i++) {
		parallel_population_averages[i].setZero();
	}
}

bool DataLikelihoodTimeCoursePopulationAverage::Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp)
{
	if (!PrepateEvaluation(values, transformed_values, non_sampled_parameters, logp)) {
		return false;
	}

	// For parallel evaluations, add all the parallel evaluations together
	for (size_t i = 0; i < parallel_population_averages.size(); i++) {
		for (int j = 0; j < population_average.array().size(); j++) {
			if (!std::isnan(parallel_population_averages[i].array()(j))) {
				if (std::isnan(population_average.array()(j))) {
					population_average.array()(j) = parallel_population_averages[i].array()(j);
				} else {
					population_average.array()(j) += parallel_population_averages[i].array()(j);
				}
			}
		}
	}

	if (relative_to_time_average) {
		VectorReal time_average = population_average.colwise().mean();
		for (int i = 0; i < time_average.size(); i++) {
			//population_average.col(i) = (population_average.col(i) /= time_average(i)).array().log();
			population_average.col(i) = population_average.col(i).array() -= time_average(i);
		}
	}

	ASSERT(species_names.size() == 1);
	population_average.array() *= data_scales[0];
	population_average.array() += data_offsets[0];

	logp = 0.0;

	for (int i = 0; i < timepoints.size(); i++) {
		Real x = population_average.row(i).sum();
		if (std::isnan(x)) {
			Real cell_trajectories_firstnonan = timepoints(timepoints.size() - 1);
			for (int m = 0; m < timepoints.size(); m++) {
				if (!std::isnan(population_average.row(m).sum())) {
					cell_trajectories_firstnonan = timepoints(m);
					break;
				}
			}
			Real cell_trajectories_lastnonan = timepoints(0);
			for (int m = timepoints.size() - 1; m >= 0; m--) {
				if (!std::isnan(population_average.row(m).sum())) {
					cell_trajectories_lastnonan = timepoints(m);
					break;
				}
			}

			Real time_offset = std::min(std::abs(timepoints(i) - cell_trajectories_firstnonan), std::abs(timepoints(i) - cell_trajectories_lastnonan));

			logp += observed_data.rows() * EvaluateMissingValue(time_offset);
		} else {
			for (int j = 0; j < observed_data.rows(); j++) {
				logp += EvaluateValue(observed_data(j, i), x, 0);
			}
		}
	}

	logp *= weight;

	return true;
}

bool DataLikelihoodTimeCoursePopulationAverage::NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis, bool is_newborn)
{
	if (species_ix == std::numeric_limits<size_t>::max()) {
		return true;
	}

	std::vector<size_t>& our_species_ixs = species_map[species_ix];
	for (size_t i = 0; i < our_species_ixs.size(); i++) {
		size_t our_species_ix = our_species_ixs[i];

		if (entered_mitosis || !include_only_cells_that_went_through_mitosis) {
			Real y = x;
			if (include_only_cells_that_went_through_mitosis) {
				y /= mitotic_population_size;
			} else {
				y /= current_population_size;
			}
			if (!parallel_population_averages.empty()) {
				Real& target = parallel_population_averages[parallel_evaluation_ix](timepoint_ix, our_species_ix);
				if (std::isnan(target)) {
					target = y;
				} else {
					target += y;
				}
			} else {
				Real& target = population_average(timepoint_ix, our_species_ix);
				if (std::isnan(target)) {
					target = y;
				} else {
					target += y;
				}
			}
		}
	}

	return true;
}
