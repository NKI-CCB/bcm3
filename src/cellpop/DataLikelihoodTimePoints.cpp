#include "Utils.h"
#include "DataLikelihoodTimePoints.h"
#include "Experiment.h"
#include "ProbabilityDistributions.h"
#include "VectorUtils.h"
#include "../../dependencies/HungarianAlgorithm-master/matching.h"
#include <boost/algorithm/string.hpp>

DataLikelihoodTimePoints::DataLikelihoodTimePoints(size_t parallel_evaluations)
	: synchronize(ESynchronizeCellTrajectory::None)
	, use_only_nondivided(false)
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

	use_only_nondivided = xml_node.get<bool>("<xmlattr>.use_only_nondivided", false);

	std::string synchronize_str = xml_node.get<std::string>("<xmlattr>.synchronize", "");
	if (synchronize_str == "" || synchronize_str == "none") {
		synchronize = ESynchronizeCellTrajectory::None;
	} else if (synchronize_str == "DNA_replication_start") {
		synchronize = ESynchronizeCellTrajectory::DNAReplicationStart;
	} else if (synchronize_str == "PCNA_gfp_increase") {
		synchronize = ESynchronizeCellTrajectory::PCNA_gfp_increase;
	} else if (synchronize_str == "mitosis" || synchronize_str == "nuclear_envelope_breakdown") {
		synchronize = ESynchronizeCellTrajectory::NuclearEnvelopeBreakdown;
	} else if (synchronize_str == "anaphase" || synchronize_str == "anaphase_onset") {
		synchronize = ESynchronizeCellTrajectory::AnaphaseOnset;
	} else {
		LOGERROR("Synchronization is specified as \"%s\" which is not a recognized synchronization point", synchronize_str.c_str());
		return false;
	}

	std::string use_only_cell_ix_str = vm["cellpop.use_only_cell_ix"].as<std::string>();
	std::vector<std::string> use_only_cell_ix_tokens;
	bcm3::tokenize(use_only_cell_ix_str, use_only_cell_ix_tokens, ",");

	std::string time_dimension_name;
	size_t num_dimensions = 0;
	result &= data_file.GetDimensionName(experiment->GetName(), data_name, 0, time_dimension_name);
	result &= data_file.GetDimensionCount(experiment->GetName(), data_name, &num_dimensions);

	// Check whether there is the right amount of data
	if (num_dimensions != 2 && num_dimensions != 3) {
		LOGERROR("Need 2 or 3 dimensional data");
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
	matched_data.resize(num_timepoints);
	if (num_dimensions == 2) {
		if (use_only_cell_ix_tokens.empty()) {
			for (size_t i = 0; i < num_timepoints; i++) {
				VectorReal od;
				result &= data_file.GetValuesDim2(experiment->GetName(), data_name, i, 0, num_cells, od);
				observed_data[i] = od;
				matched_data[i] = VectorReal::Constant(num_cells, std::numeric_limits<Real>::quiet_NaN());
			}
		} else {
			for (size_t i = 0; i < num_timepoints; i++) {
				observed_data[i].setConstant(use_only_cell_ix_tokens.size(), 1, std::numeric_limits<Real>::quiet_NaN());
				matched_data[i].setConstant(use_only_cell_ix_tokens.size(), 1, std::numeric_limits<Real>::quiet_NaN());
				for (size_t j = 0; j < use_only_cell_ix_tokens.size(); j++) {
					int ix = boost::lexical_cast<int>(use_only_cell_ix_tokens[j]);
					if (ix >= num_cells) {
						LOGERROR("Requested to use cell %d, but data contains only %u cells", ix, num_cells);
						return false;
					} else {
						result &= data_file.GetValue(experiment->GetName(), data_name, i, ix, &observed_data[i](j, 0));
					}
				}
				num_cells = use_only_cell_ix_tokens.size();
			}
		}
	} else if (num_dimensions == 3) {
		size_t num_markers;
		std::string marker_dim_name;
		result &= data_file.GetDimensionName(experiment->GetName(), data_name, 2, marker_dim_name);
		result &= data_file.GetDimensionSize(experiment->GetName(), marker_dim_name, &num_markers);

		if (use_only_cell_ix_tokens.empty()) {
			for (size_t i = 0; i < num_timepoints; i++) {
				observed_data[i] = MatrixReal::Zero(num_cells, num_markers);
				for (size_t j = 0; j < num_markers; j++) {
					VectorReal od;
					result &= data_file.GetValuesDim2(experiment->GetName(), data_name, i, 0, j, num_cells, od);
					observed_data[i].col(j) = od;
				}
				matched_data[i] = MatrixReal::Constant(num_cells, num_markers, std::numeric_limits<Real>::quiet_NaN());
			}
		} else {
			LOGERROR("Not implemented yet");
			return false;
		}
	} else {
		ASSERT(false);
		return false;
	}

	std::string species_name = xml_node.get<std::string>("<xmlattr>.species_name");
	if (species_name.find_first_of(';') != std::string::npos) {
		bcm3::tokenize(species_name, species_names, ";");
		for (size_t i = 0; i < species_names.size(); i++) {
			boost::trim(species_names[i]);
		}
	} else {
		species_names.resize(1, species_name);
	}

	if (experiment->GetMaxNumberOfCells() < num_cells) {
		LOGERROR("Maximum number of simulated cells (%zu) in the experiment is not sufficient for the amount of cells in the data (%zu)", experiment->GetMaxNumberOfCells(), num_cells);
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
						experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, synchronize);
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
					experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, synchronize);
				}
			}
			species_map[species_ix].push_back(j);
		}
	}

	if (synchronize != ESynchronizeCellTrajectory::None) {
		// If we're going to synchronize, make sure we simulate the cells as long as the full duration of the time course
		// (this is needed in case there are negative timepoints)
		Real last_tp = timepoints(timepoints.size() - 1);
		Real full_duration = last_tp - timepoints(0);
		if (full_duration > last_tp) {
			experiment->AddSimulationTimepoints(this, full_duration, num_timepoints + 1, std::numeric_limits<size_t>::max(), synchronize);
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
	std::vector<Real> stdevs(species_names.size());
	std::vector<Real> data_offsets(species_names.size());
	std::vector<Real> data_scales(species_names.size());
	for (size_t i = 0; i < species_names.size(); i++) {
		stdevs[i] = GetCurrentSTDev(transformed_values, non_sampled_parameters, i);
		data_offsets[i] = GetCurrentDataOffset(transformed_values, non_sampled_parameters, i);
		data_scales[i] = GetCurrentDataScale(transformed_values, non_sampled_parameters, i);
	}

	for (ptrdiff_t ti = 0; ti < timepoints.size(); ti++) {
		matched_data[ti].setConstant(std::numeric_limits<Real>::quiet_NaN());

		size_t finite_data_count = 0;
		for (ptrdiff_t i = 0; i < observed_data[ti].rows(); i++) {
			if (observed_data[ti].row(i).array().isFinite().any()) {
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
		Eigen::MatrixXi observed_data_indices(finite_data_count, finite_sim_count);
		Eigen::MatrixXi trajectory_indices(finite_data_count, finite_sim_count);
		cell_likelihoods.setConstant(-std::numeric_limits<Real>::infinity());
		size_t i_ix = 0;
		for (ptrdiff_t i = 0; i < observed_data[ti].rows(); i++) {
			if (!observed_data[ti].row(i).array().isFinite().any()) {
				continue;
			}

			size_t j_ix = 0;
			for (size_t j = 0; j < cell_trajectories.size(); j++) {
				if (std::isnan(cell_trajectories[j](ti, 0))) {
					continue;
				}

				Real cell_logp = 0.0;
				for (int l = 0; l < species_names.size(); l++) {
					Real x = data_offsets[l] + data_scales[l] * cell_trajectories[j](ti, l);
					Real y = observed_data[ti](i, l);
					if (std::isnan(y)) {
						// At least one value of this cell should be non-nan (due to check for any finite in the row; above
						// the other nan's we can ignore as missing values
						continue;
					}

					if (error_model == ErrorModel::Normal) {
						cell_logp += bcm3::LogPdfNormal(y, x, stdevs[l]);
					} else if (error_model == ErrorModel::StudentT4) {
						cell_logp += bcm3::LogPdfTnu4(y, x, stdevs[l]);
					} else {
						assert(false);
						cell_logp = std::numeric_limits<Real>::quiet_NaN();
					}
				}
				cell_likelihoods(i_ix, j_ix) = cell_logp;
				observed_data_indices(i_ix, j_ix) = i;
				trajectory_indices(i_ix, j_ix) = j;
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
#if 0
					int data_ix = observed_data_indices(i, j);
					int traj_ix = trajectory_indices(i, j);
					for (int l = 0; l < species_names.size(); l++) {
						matched_data[ti](data_ix, l) = data_offset + data_scale * cell_trajectories[traj_ix](ti, l);
					}
#endif
				}
			}
		}
#if 1
		for (ptrdiff_t i = 0; i < observed_data[ti].rows(); i++) {
			for (ptrdiff_t l = 0; l < species_names.size(); l++) {
				matched_data[ti](i, l) = data_offsets[l] + data_scales[l] * cell_trajectories[i](ti, l);
			}
		}
#endif
	}

	logp *= weight;

	return true;
}

bool DataLikelihoodTimePoints::NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis, bool is_newborn)
{
	if (species_ix == std::numeric_limits<size_t>::max()) {
		return true;
	}

	if (use_only_nondivided && is_newborn) {
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
