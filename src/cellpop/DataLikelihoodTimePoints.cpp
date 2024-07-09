#include "Utils.h"
#include "DataLikelihoodTimePoints.h"
#include "Experiment.h"
#include "ProbabilityDistributions.h"
#include "../../dependencies/HungarianAlgorithm-master/matching.h"
#include <boost/algorithm/string.hpp>

DataLikelihoodTimePoints::DataLikelihoodTimePoints(size_t parallel_evaluations)
{
}

DataLikelihoodTimePoints::~DataLikelihoodTimePoints()
{
}

bool DataLikelihoodTimePoints::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm)
{
	if (!DataLikelihoodBase::Load(xml_node, experiment, varset, data_file, vm)) {
		return false;
	}

	bool result = true;

	std::string species_name = xml_node.get<std::string>("<xmlattr>.species_name");

	std::string time_dimension_name;
	size_t num_dimensions = 0;
	result &= data_file.GetDimensionName(experiment->GetName(), data_name, 0, time_dimension_name);
	result &= data_file.GetDimensionCount(experiment->GetName(), data_name, &num_dimensions);

	// Check whether there is the right amount of data
	if (num_dimensions != 2) {
		LOGERROR("Need 2 dimensional data");
		return false;
	}

	// Retrieve data from the datafile
	size_t num_timepoints;
	result &= data_file.GetDimensionSize(experiment->GetName(), time_dimension_name, &num_timepoints);
	result &= data_file.GetValues(experiment->GetName(), time_dimension_name, 0, num_timepoints, timepoints);

	size_t num_cells;
	std::string cell_dimension_name;
	result &= data_file.GetDimensionName(experiment->GetName(), data_name, 1, cell_dimension_name);
	result &= data_file.GetDimensionSize(experiment->GetName(), cell_dimension_name, &num_cells);

	observed_data.resize(num_timepoints);
	if (num_dimensions == 2) {
		for (size_t i = 0; i < num_timepoints; i++) {
			VectorReal od;
			result &= data_file.GetValuesDim2(experiment->GetName(), data_name, i, 0, num_cells, od);
			observed_data[i] = od;
		}
	} else {
		LOGERROR("TODO");
		return false;
	}

	if (species_name.find_first_of(';') != std::string::npos) {
		bcm3::tokenize(species_name, species_names, ";");
		for (size_t i = 0; i < species_names.size(); i++) {
			boost::trim(species_names[i]);
		}
	} else {
		species_names.resize(1, species_name);
	}

	if (experiment->GetMaxNumberOfCells() < observed_data.size()) {
		LOGERROR("Maximum number of simulated cells (%zu) in the experiment is not sufficient for the amount of cells in the data (%zu)", experiment->GetMaxNumberOfCells(), observed_data.size());
		return false;
	}
	cell_trajectories.resize(experiment->GetMaxNumberOfCells(), MatrixReal::Constant(num_timepoints, species_names.size(), std::numeric_limits<Real>::quiet_NaN()));

	// Request the experiment to provide the required information during the simulation
	for (size_t j = 0; j < species_names.size(); j++) {
		const std::string& species_name = species_names[j];

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
					for (size_t i = 0; i < num_timepoints; i++) {
						experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, ESynchronizeCellTrajectory::None);
					}
				}
				species_map[species_ix].push_back(j);
			}
		} else if (species_name.find_first_of('/') != std::string::npos) {
			LOGERROR("Division currently not supported for time points data");
			return false;
		} else {
			size_t species_ix = experiment->GetCVodeSpeciesByName(species_names[j]);
			if (species_ix == std::numeric_limits<size_t>::max()) {
				return false;
			}
			if (species_map[species_ix].empty()) {
				for (size_t i = 0; i < num_timepoints; i++) {
					experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, ESynchronizeCellTrajectory::None);
				}
			}
			species_map[species_ix].push_back(j);
		}
	}

	return result;
}

void DataLikelihoodTimePoints::Reset()
{
	for (size_t i = 0; i < cell_trajectories.size(); i++) {
		cell_trajectories[i].setConstant(std::numeric_limits<Real>::quiet_NaN());
	}
}

bool DataLikelihoodTimePoints::Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp)
{
	Real stdev = GetCurrentSTDev(transformed_values, non_sampled_parameters);
	Real data_offset = GetCurrentDataOffset(transformed_values, non_sampled_parameters);
	Real data_scale = GetCurrentDataScale(transformed_values, non_sampled_parameters);

	for (ptrdiff_t ti = 0; ti < timepoints.size(); ti++) {
		int n = std::max((size_t)observed_data[ti].size(), cell_trajectories.size());

		size_t finite_data_count = 0;
		for (size_t i = 0; i < observed_data[ti].size(); i++) {
			if (!std::isnan(observed_data[ti](i))) {
				finite_data_count++;
			}
		}
		if (finite_data_count == 0) {
			continue;
		}

		size_t finite_sim_count = 0;
		for (size_t j = 0; j < cell_trajectories.size(); j++) {
			if (!std::isnan(cell_trajectories[j](ti, 0))) {
				finite_sim_count++;
			}
		}

		if (finite_sim_count < finite_data_count) {
			// Not enough simulated cells at this timepoint; won't be able to match all cells so no need to go further and HA doesn't like it
			logp = -std::numeric_limits<Real>::infinity();
			return true;
		}

		MatrixReal cell_likelihoods(finite_data_count, finite_sim_count);
		cell_likelihoods.setConstant(-std::numeric_limits<Real>::infinity());
		size_t i_ix = 0;
		for (size_t i = 0; i < observed_data[ti].size(); i++) {
			if (std::isnan(observed_data[ti](i))) {
				continue;
			}

			size_t j_ix = 0;
			for (size_t j = 0; j < cell_trajectories.size(); j++) {
				if (std::isnan(cell_trajectories[j](ti, 0))) {
					continue;
				}


				Real cell_logp = 0.0;
				for (int l = 0; l < species_names.size(); l++) {
					Real x = cell_trajectories[j](ti, l);
					Real y = observed_data[ti](i, l);

					if (error_model == ErrorModel::Normal) {
						cell_logp += bcm3::LogPdfNormal(y, x, stdev);
					} else if (error_model == ErrorModel::StudentT4) {
						cell_logp += bcm3::LogPdfTnu4(y, x, stdev);
					} else {
						assert(false);
						cell_logp = std::numeric_limits<Real>::quiet_NaN();
					}
				}
				cell_likelihoods(i_ix, j_ix) = cell_logp;
				j_ix++;
			}
			i_ix++;
		}

		MatrixReal matching;
		if (!HungarianMatching(cell_likelihoods, matching, HUNGARIAN_MATCH_MAX)) {
			LOGERROR("Could not find matching");
			logp = -std::numeric_limits<Real>::infinity();
			return true;
		}
		for (size_t i = 0; i < finite_data_count; i++) {
			for (size_t j = 0; j < finite_sim_count; j++) {
				if (matching(i, j) == 1) {
					logp += cell_likelihoods(i, j);
				}
			}
		}
	}

	logp *= weight;

	return true;
}

bool DataLikelihoodTimePoints::NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis)
{
	if (species_ix == std::numeric_limits<size_t>::max()) {
		return true;
	}

	std::vector<size_t>& our_species_ixs = species_map[species_ix];
	for (size_t i = 0; i < our_species_ixs.size(); i++) {
		size_t our_species_ix = our_species_ixs[i];

		if (cell_ix < cell_trajectories.size()) {
			Real& val = cell_trajectories[cell_ix](timepoint_ix, our_species_ix);
			if (std::isnan(val)) {
				val = x;
			} else {
				val += x;
			}
		}
	}

	return true;
}
