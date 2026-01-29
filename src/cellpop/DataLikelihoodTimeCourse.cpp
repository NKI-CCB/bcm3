#include "Utils.h"
#include "Correlation.h"
#include "DataLikelihoodTimeCourse.h"
#include "Experiment.h"
#include "ProbabilityDistributions.h"
//#include "../../dependencies/HungarianAlgorithm-master/matching.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/algorithm/string.hpp>

DataLikelihoodTimeCourse::DataLikelihoodTimeCourse(size_t parallel_evaluations)
	: synchronize(ESynchronizeCellTrajectory::None)
{
}

DataLikelihoodTimeCourse::~DataLikelihoodTimeCourse()
{
}

bool DataLikelihoodTimeCourse::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm)
{
	if (!DataLikelihoodTimeCourseBase::Load(xml_node, experiment, varset, data_file, vm)) {
		return false;
	}

	bool result = true;

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

	// First check whether we are using only a subset of the cells in the data
	std::string use_only_cell_ix_str = vm["cellpop.use_only_cell_ix"].as<std::string>();
	std::vector<std::string> use_only_cell_ix_tokens;
	bcm3::tokenize(use_only_cell_ix_str, use_only_cell_ix_tokens, ",");

	// Retrieve data from the datafile
	size_t num_dimensions = 0;
	result &= data_file.GetDimensionCount(experiment->GetName(), data_name, &num_dimensions);
	if (num_dimensions == 1) {
		LOG("Time course likelihood for data %s; data has 1 dimension; assuming only 1 cell and no replicates.", data_name.c_str());

		if (use_only_cell_ix_str != "-1") {
			LOGERROR("use_only_cell_ix has been specified, but there is only 1 cell in the data.");
			return false;
		}

		observed_data.resize(1);
		VectorReal od;
		result &= data_file.GetValues(experiment->GetName(), data_name, 0, timepoints.size(), od);
		observed_data[0] = od;
	} else if (num_dimensions == 2) {
		size_t num_cells;
		std::string cell_dimension_name;
		result &= data_file.GetDimensionName(experiment->GetName(), data_name, 1, cell_dimension_name);
		result &= data_file.GetDimensionSize(experiment->GetName(), cell_dimension_name, &num_cells);
		LOG("Time course likelihood for data %s; data has 2 dimension; assuming the second dimension (\"%s\") are the cells (%zu cells).", data_name.c_str(), cell_dimension_name.c_str(), num_cells);

		if (use_only_cell_ix_str == "-1") {
			observed_data.resize(num_cells);
			for (size_t i = 0; i < num_cells; i++) {
				VectorReal od;
				result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, i, timepoints.size(), od);
				observed_data[i] = od;
			}
		} else {
			observed_data.resize(use_only_cell_ix_tokens.size());
			for (size_t i = 0; i < use_only_cell_ix_tokens.size(); i++) {
				int ix = boost::lexical_cast<int>(use_only_cell_ix_tokens[i]);
				if (ix >= num_cells) {
					LOGERROR("Requested to use cell %d, but data contains only %u cells", ix, num_cells);
					return false;
				} else {
					VectorReal od;
					result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, ix, timepoints.size(), od);
					observed_data[i] = od;
				}
			}
		}
	} else if (num_dimensions == 3) {
		size_t num_cells, num_markers;
		std::string cell_dimension_name, marker_dim_name;
		result &= data_file.GetDimensionName(experiment->GetName(), data_name, 1, cell_dimension_name);
		result &= data_file.GetDimensionSize(experiment->GetName(), cell_dimension_name, &num_cells);
		result &= data_file.GetDimensionName(experiment->GetName(), data_name, 2, marker_dim_name);
		result &= data_file.GetDimensionSize(experiment->GetName(), marker_dim_name, &num_markers);
		LOG("Time course likelihood for data %s; data has 3 dimension; assuming the second dimension (\"%s\") are the cells (%zu cells) and the third dimension are different markers.", data_name.c_str(), cell_dimension_name.c_str(), num_cells);

		if (use_only_cell_ix_str == "-1") {
			observed_data.resize(num_cells);
			for (size_t i = 0; i < num_cells; i++) {
				observed_data[i] = MatrixReal::Zero(timepoints.size(), num_markers);
				for (size_t j = 0; j < num_markers; j++) {
					VectorReal od;
					result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, i, j, timepoints.size(), od);
					observed_data[i].col(j) = od;
				}
			}
		} else {
			observed_data.resize(use_only_cell_ix_tokens.size());
			for (size_t i = 0; i < use_only_cell_ix_tokens.size(); i++) {
				int ix = boost::lexical_cast<int>(use_only_cell_ix_tokens[i]);
				if (ix >= num_cells) {
					LOGERROR("Requested to use cell %d, but data contains only %u cells", ix, num_cells);
					return false;
				} else {
					observed_data[i] = MatrixReal::Zero(timepoints.size(), num_markers);
					for (size_t j = 0; j < num_markers; j++) {
						VectorReal od;
						result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, ix, j, timepoints.size(), od);
						observed_data[i].col(j) = od;
					}
				}
			}
		}
	} else {
		LOGERROR("Time course likelihood for data %s; data has %zu dimensions but can only handle 1, 2 or 3 dimensions.", data_name.c_str(), num_dimensions);
		return false;
	}
	
	// Load parent information if available
	if (data_file.VariableExists(experiment->GetName(), "parent")) {
		std::vector<unsigned int> cell_ids(observed_data.size());
		if (use_only_cell_ix_str == "-1") {
			for (unsigned int i = 0; i < observed_data.size(); i++) {
				data_file.GetValue(experiment->GetName(), "cell_id", i, cell_ids.data() + i);
			}
		} else {
			for (unsigned int i = 0; i < observed_data.size(); i++) {
				data_file.GetValue(experiment->GetName(), "cell_id", boost::lexical_cast<int>(use_only_cell_ix_tokens[i]), cell_ids.data() + i);
			}
		}

		observed_children.resize(observed_data.size());
		for (size_t i = 0; i < observed_data.size(); i++) {
			int parent_ix;

			if (use_only_cell_ix_str == "-1") {
				data_file.GetValue(experiment->GetName(), "parent", i, &parent_ix);
			} else {
				data_file.GetValue(experiment->GetName(), "parent", boost::lexical_cast<int>(use_only_cell_ix_tokens[i]), &parent_ix);
			}

			if (parent_ix != std::numeric_limits<int>::min()) {
				auto cii = std::find(cell_ids.begin(), cell_ids.end(), parent_ix);
				if (cii == cell_ids.end()) {
					LOGERROR("Could not find cell %d for parent of cell %d", parent_ix, i);
					return false;
				} else {
					observed_children[cii - cell_ids.begin()].insert(i);
				}
			} else {
				observed_cells_with_no_parents.push_back(i);
			}
		}
	} else {
		for (size_t i = 0; i < observed_data.size(); i++) {
			observed_cells_with_no_parents.push_back(i);
		}
	}

	// Check number of cells and allocate memory for receiving data from the simulator
	if (experiment->GetMaxNumberOfCells() < observed_data.size()) {
		LOGERROR("Maximum number of simulated cells (%zu) in the experiment is not sufficient for the amount of cells in the data (%zu)", experiment->GetMaxNumberOfCells(), observed_data.size());
		return false;
	}
	if (experiment->GetMaxNumberOfCells() > observed_data.size()) {
		// This could fairly easily be supported by fixing the hungarian matching; but estimating variances with more simulated than observed cells leads to 
		// very difficult posteriors in the variance parameters, due to the Sobol variability
		LOGERROR("Simulating more cells (%zu) than there are cells in the data (%zu) - currently not supported.", experiment->GetMaxNumberOfCells(), observed_data.size());
		return false;
	}

	cell_trajectories.resize(experiment->GetMaxNumberOfCells(), MatrixReal::Constant(timepoints.size(), species_names.size(), std::numeric_limits<Real>::quiet_NaN()));
	matched_trajectories.resize(observed_data.size(), MatrixReal::Constant(observed_data[0].rows(), observed_data[0].cols(), std::numeric_limits<Real>::quiet_NaN()));
	trajectory_matching.resize(cell_trajectories.size());
	hungarian_matching_edges.resize(observed_cells_with_no_parents.size() * experiment->GetMaxNumberOfCells());

	result &= RequestSimulationInfo(experiment, ESynchronizeCellTrajectory::None);

	if (synchronize != ESynchronizeCellTrajectory::None) {
		// If we're going to synchronize, make sure we simulate the cells as long as the full duration of the time course
		// (this is needed in case there are negative timepoints)
		Real last_tp = timepoints(timepoints.size() - 1);
		Real full_duration = last_tp - timepoints(0);
		if (full_duration > last_tp) {
			result &= experiment->AddSimulationTimepoints(this, full_duration, timepoints.size() + 1, std::numeric_limits<size_t>::max(), synchronize);
		}
	}

	return result;
}

bool DataLikelihoodTimeCourse::PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
{
	if (!DataLikelihoodTimeCourseBase::PostInitialize(varset, non_sampled_parameter_names)) {
		return false;
	}

	return true;
}

void DataLikelihoodTimeCourse::Reset()
{
	for (size_t i = 0; i < cell_trajectories.size(); i++) {
		cell_trajectories[i].setConstant(std::numeric_limits<Real>::quiet_NaN());
	}
	for (size_t i = 0; i < observed_data.size(); i++) {
		matched_trajectories[i].setConstant(std::numeric_limits<Real>::quiet_NaN());
	}
	simulated_cell_parents.clear();
	simulated_cell_children.clear();
}

bool DataLikelihoodTimeCourse::Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp)
{
	if (!PrepateEvaluation(values, transformed_values, non_sampled_parameters, logp)) {
		return false;
	}

	for (size_t i = 0; i < cell_trajectories.size(); i++) {
		for (int j = 0; j < species_names.size(); j++) {
			cell_trajectories[i].col(j).array() *= data_scales[j];
			cell_trajectories[i].col(j).array() += data_offsets[j];
		}
	}

	if (use_signal_saturation) {
		for (size_t i = 0; i < cell_trajectories.size(); i++) {
			//cell_trajectories[i].array() = saturation_scale / (1.0 + exp(-saturation_hill * cell_trajectories[i].array())) - 0.5 * saturation_scale;
			cell_trajectories[i] *= -1.0; // *= -saturation_hill // the hill coefficient can just be absorbed into the data_scale
			cell_trajectories[i].array() = cell_trajectories[i].array().exp();
			cell_trajectories[i].array() += 1.0;
			cell_trajectories[i].array() = 1.0 / cell_trajectories[i].array();
			cell_trajectories[i] *= saturation_scale;
			cell_trajectories[i].array() -= 0.5 * saturation_scale;
		}
	}

	logp = 0.0;

	// Calculate likelihood of every observed parent cell (recursing into the children) against every simulated cell
	// TODO - relay from the likelihood file how many cells are going to be simulated
	int n = std::max(observed_cells_with_no_parents.size(), simulated_cell_parents.size());
	MatrixReal cell_likelihoods(n, n);
	cell_likelihoods.setZero();

	std::vector< std::vector< std::vector<int> > > matched_hierarchies;
	if (!observed_children.empty()) {
		matched_hierarchies.resize(observed_cells_with_no_parents.size());
		for (size_t i = 0; i < observed_cells_with_no_parents.size(); i++) {
			matched_hierarchies[i].resize(simulated_cell_parents.size());
			for (size_t j = 0; j < simulated_cell_parents.size(); j++) {
				matched_hierarchies[i][j].resize(cell_trajectories.size());
			}
		}
	}

	Real maxdiff = 0.0;
	int edge_count = 0;
	for (size_t i = 0; i < observed_cells_with_no_parents.size(); i++) {
		size_t finite_count = 0;
		size_t observed_cell_ix = observed_cells_with_no_parents[i];
		for (size_t j = 0; j < simulated_cell_parents.size(); j++) {
			if (simulated_cell_parents[j] == std::numeric_limits<size_t>::max()) {
				if (!observed_children.empty()) {
					cell_likelihoods(i, j) = CalculateCellLikelihood(observed_cell_ix, j, stdevs, proportional_stdevs, missing_simulation_time_stdev, &matched_hierarchies[i][j]);
				} else {
					cell_likelihoods(i, j) = CalculateCellLikelihood(observed_cell_ix, j, stdevs, proportional_stdevs, missing_simulation_time_stdev, nullptr);
				}
				if (cell_likelihoods(i, j) != cell_likelihoods(i, j)) {
					logp = -std::numeric_limits<Real>::infinity();
					return false;
				}

				if (cell_likelihoods(i,j) > -std::numeric_limits<Real>::infinity()) {
					finite_count++;
				}
			} else {
				cell_likelihoods(i, j) = -std::numeric_limits<Real>::infinity();
			}

			hungarian_matching_edges[edge_count].left = i;
			hungarian_matching_edges[edge_count].right = j;
			hungarian_matching_edges[edge_count].cost = -cell_likelihoods(i, j);
			edge_count++;
		}
		if (finite_count < observed_cells_with_no_parents.size()) {
			// Too many infinite likelihoods for this observed cell; won't be able to match all cells so no need to go further and HA doesn't like it
			logp = -std::numeric_limits<Real>::infinity();
			return true;
		}
	}

	std::vector<int> matching = hungarianMinimumWeightPerfectMatching(n, hungarian_matching_edges, edge_count);
	if (matching.size() != observed_cells_with_no_parents.size()) {
		logp = -std::numeric_limits<Real>::infinity();
		return true;
	} else {
		for (int i = 0; i < observed_cells_with_no_parents.size(); i++) {
			int j = matching[i];
			if (j == -1) {
				LOGERROR("Could not find matching");
				logp = -std::numeric_limits<Real>::infinity();
				return true;
			} else {
				logp += cell_likelihoods(i, j);

				size_t oi = observed_cells_with_no_parents[i];

				if (optimize_offset_scale) {
					for (int l = 0; l < species_names.size(); l++) {
						Real offset = 0.0, scale = 1.0;
						OptimizeOffsetScale(observed_data[oi], cell_trajectories[j], l, offset, scale);
						matched_trajectories[oi].col(l).array() = offset + scale * cell_trajectories[j].col(l).array();
					}
				} else {
					matched_trajectories[oi] = cell_trajectories[j];
				}

				trajectory_matching[oi] = j;

				if (!observed_children.empty()) {
					for (auto ci = observed_children[oi].begin(); ci != observed_children[oi].end(); ci++) {
						size_t child = matched_hierarchies[i][j][*ci];
						matched_trajectories[*ci] = cell_trajectories[child];
						trajectory_matching[*ci] = child;
					}
				}
			}
		}
	}

	logp *= weight;

	return true;
}

bool DataLikelihoodTimeCourse::NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis, bool is_newborn)
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
				if (use_log_ratio) {
					if (species_ix == log_ratio_numerator_ix[i]) {
						if (val < 1e-16) {
							val = 1e-16;
						}
						val = bcm3::log10(x / val);
					} else {
						if (x < 1e-16) {
							val = bcm3::log10(val / 1e-16);
						} else {
							val = bcm3::log10(val / x);
						}
					}
				} else {
					val += x;
				}
			}
		}
	}

	return true;
}

void DataLikelihoodTimeCourse::NotifyStartingCells(size_t cell_ix)
{
	simulated_cell_parents.resize(cell_ix + 1, std::numeric_limits<size_t>::max());
	simulated_cell_children.resize(cell_ix + 1, std::pair<size_t, size_t>(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
}

void DataLikelihoodTimeCourse::NotifyParents(size_t parent, size_t child)
{
	if (child >= simulated_cell_parents.size()) {
		simulated_cell_parents.resize(child + 1, std::numeric_limits<size_t>::max());
	}
	simulated_cell_parents[child] = parent;

	if (parent >= simulated_cell_children.size()) {
		simulated_cell_children.resize(parent + 1, std::pair<size_t, size_t>(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max()));
	}
	if (simulated_cell_children[parent].first == std::numeric_limits<size_t>::max()) {
		simulated_cell_children[parent].first = child;
	} else if (simulated_cell_children[parent].second == std::numeric_limits<size_t>::max()) {
		simulated_cell_children[parent].second = child;
	} else {
		// Multipolar division?! not implemented
		ASSERT(false);
	}
}

Real DataLikelihoodTimeCourse::CalculateCellLikelihood(size_t observed_cell_ix, size_t simulated_cell_ix, const std::vector<Real>& stdevs, const std::vector<Real>& proportional_stdevs, Real missing_simulation_time_stdev, std::vector<int>* matched_hierarchy)
{
	Real cell_logp = 0.0;

	if (simulated_cell_ix == std::numeric_limits<size_t>::max()) {
		// Apparently the simulated cell did not divide, but we do have observed data, so
		// accumulate missing simulation penalties for every data point, including further recursion
		for (int l = 0; l < species_names.size(); l++) {
			for (int k = 0; k < timepoints.size(); k++) {
				Real y = observed_data[observed_cell_ix](k);
				if (!std::isnan(y)) {
					cell_logp += CalculateMissingValueLikelihood(simulated_cell_ix, k, l, missing_simulation_time_stdev);
				}
			}
		}
		if (!observed_children.empty()) {
			for (std::set<size_t>::const_iterator ci = observed_children[observed_cell_ix].begin(); ci != observed_children[observed_cell_ix].end(); ++ci) {
				cell_logp += CalculateCellLikelihood(*ci, std::numeric_limits<size_t>::max(), stdevs, proportional_stdevs, missing_simulation_time_stdev, matched_hierarchy);
			}
		}
	} else {
		for (int l = 0; l < species_names.size(); l++) {
			Real offset = 0.0;
			Real scale = 1.0;
			if (optimize_offset_scale) {
				OptimizeOffsetScale(observed_data[observed_cell_ix], cell_trajectories[simulated_cell_ix], l, offset, scale);
			}

			Real minus_log_sigma = -log(stdevs[l]);
			Real inv_two_sigma_sq = 1.0 / (2.0 * stdevs[l] * stdevs[l]);

			for (int k = 0; k < timepoints.size(); k++) {
				Real y = observed_data[observed_cell_ix](k, l);
				if (std::isnan(y)) {
					continue;
				}

				Real x = offset + scale * cell_trajectories[simulated_cell_ix](k, l);
				if (std::isnan(x)) {
					cell_logp += CalculateMissingValueLikelihood(simulated_cell_ix, k, l, missing_simulation_time_stdev);
				} else {
					if (error_model == ErrorModel::Normal) {
						Real d = y - x;
						cell_logp += minus_log_sigma - 0.91893853320467274178032973640562 - d * d * inv_two_sigma_sq;
					} else if (error_model == ErrorModel::ProportionalNormal) {
						cell_logp += bcm3::LogPdfNormal(x, y, proportional_stdevs[l] * std::max(x, 0.0));
					} else if (error_model == ErrorModel::AdditiveProportionalNormal) {
						cell_logp += bcm3::LogPdfNormal(x, y, stdevs[l] + proportional_stdevs[l] * std::max(x, 0.0));
					} else if (error_model == ErrorModel::StudentT4) {
						cell_logp += bcm3::LogPdfTnu4(y, x, stdevs[l]);
					} else {
						assert(false);
						cell_logp = std::numeric_limits<Real>::quiet_NaN();
					}
				}
			}

			if (cell_logp == -std::numeric_limits<Real>::infinity()) {
				// Negative infinity - no need to go further
				return cell_logp;
			}
		}

		// Simulation does not have multipolar division, so only need to consider two simulated children, if any
		if (!observed_children.empty() && !observed_children[observed_cell_ix].empty()) {
			if (simulated_cell_ix < simulated_cell_children.size() && simulated_cell_children[simulated_cell_ix].first != std::numeric_limits<size_t>::max()) {
				std::pair<size_t, size_t> simulated_children = simulated_cell_children[simulated_cell_ix];
				MatrixReal children_matching(2, observed_children[observed_cell_ix].size());
				//std::vector< std::vector<double> > children_matching(2, std::vector<double>(observed_children[observed_cell_ix].size()));
				size_t cii = 0;
				size_t finite_count1 = 0, finite_count2 = 0;
				for (std::set<size_t>::const_iterator ci = observed_children[observed_cell_ix].begin(); ci != observed_children[observed_cell_ix].end(); ++ci, ++cii) {
					children_matching(0, cii) = CalculateCellLikelihood(*ci, simulated_children.first, stdevs, proportional_stdevs, missing_simulation_time_stdev, matched_hierarchy);
					children_matching(1, cii) = CalculateCellLikelihood(*ci, simulated_children.second, stdevs, proportional_stdevs, missing_simulation_time_stdev, matched_hierarchy);

					if (children_matching(0, cii) > -std::numeric_limits<Real>::infinity()) {
						finite_count1++;
					}
					if (children_matching(1, cii) > -std::numeric_limits<Real>::infinity()) {
						finite_count2++;
					}
				}

				if (finite_count1 == 0 || finite_count2 == 0) {
					return -std::numeric_limits<Real>::infinity();
				}

				if (observed_children[observed_cell_ix].size() == 1) {
					if (children_matching(0, 0) > children_matching(1, 0)) {
						cell_logp += children_matching(0, 0);
						matched_hierarchy->at(*observed_children[observed_cell_ix].begin()) = simulated_children.first;
					} else {
						cell_logp += children_matching(1, 0);
						matched_hierarchy->at(*observed_children[observed_cell_ix].begin()) = simulated_children.second;
					}
				} else {
#if TODO
					MatrixReal matching;
					if (!HungarianMatching(children_matching, matching, HUNGARIAN_MATCH_MAX)) {
						//LOGERROR("Could not find matching");
						//std::cout << children_matching << std::endl;
						//std::cout << matching << std::endl;
						return -std::numeric_limits<Real>::infinity();
					}
					std::set<size_t>::const_iterator ci = observed_children[observed_cell_ix].begin();
					for (int i = 0; i < matching.rows(); i++, ci++) {
						for (int j = 0; j < matching.cols(); j++) {
							if (matching(i, j) == 1) {
								cell_logp += children_matching(i, j);
								if (j == 0) {
									matched_hierarchy[*ci] = simulated_children.first;
								} else {
									matched_hierarchy[*ci] = simulated_children.second;
								}
							}
						}
					}
#endif
				}
			} else {
				for (std::set<size_t>::const_iterator ci = observed_children[observed_cell_ix].begin(); ci != observed_children[observed_cell_ix].end(); ++ci) {
					cell_logp = CalculateCellLikelihood(*ci, std::numeric_limits<size_t>::max(), stdevs, proportional_stdevs, missing_simulation_time_stdev, matched_hierarchy);
				}
			}
		}
	}

	return cell_logp;
}

Real DataLikelihoodTimeCourse::CalculateMissingValueLikelihood(size_t simulated_cell_ix, int timepoint_ix, int species_ix, Real missing_simulation_time_stdev)
{
	if (simulated_cell_ix == std::numeric_limits<size_t>::max()) {
		return EvaluateMissingValue(timepoints(timepoint_ix));
	} else {
		Real cell_trajectories_firstnonan = timepoints(timepoints.size() - 1);
		for (int m = 0; m < timepoints.size(); m++) {
			if (!std::isnan(cell_trajectories[simulated_cell_ix](m, species_ix))) {
				cell_trajectories_firstnonan = timepoints(m);
				break;
			}
		}
		Real cell_trajectories_lastnonan = timepoints(0);
		for (int m = timepoints.size() - 1; m >= 0; m--) {
			if (!std::isnan(cell_trajectories[simulated_cell_ix](m, species_ix))) {
				cell_trajectories_lastnonan = timepoints(m);
				break;
			}
		}
		Real time_offset = std::min(std::abs(timepoints(timepoint_ix) - cell_trajectories_firstnonan), std::abs(timepoints(timepoint_ix) - cell_trajectories_lastnonan));
		return EvaluateMissingValue(time_offset);
	}
}
