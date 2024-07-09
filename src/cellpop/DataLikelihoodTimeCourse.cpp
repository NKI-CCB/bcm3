#include "Utils.h"
#include "DataLikelihoodTimeCourse.h"
#include "Experiment.h"
#include "ProbabilityDistributions.h"
#include "../../dependencies/HungarianAlgorithm-master/matching.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/algorithm/string.hpp>

DataLikelihoodTimeCourse::DataLikelihoodTimeCourse(size_t parallel_evaluations)
	: use_population_average(true)
	, use_log_ratio(false)
	, include_only_cells_that_went_through_mitosis(false)
	, synchronize(ESynchronizeCellTrajectory::None)
	, fixed_missing_simulation_time_stdev_ix(std::numeric_limits<size_t>::max())
	, fixed_missing_simulation_time_stdev_non_sampled_ix(std::numeric_limits<size_t>::max())
	, fixed_missing_simulation_time_stdev(300.0)
{
	if (parallel_evaluations > 1) {
		parallel_population_averages.resize(parallel_evaluations);
	}
}

DataLikelihoodTimeCourse::~DataLikelihoodTimeCourse()
{
}

bool DataLikelihoodTimeCourse::Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm)
{
	if (!DataLikelihoodBase::Load(xml_node, experiment, varset, data_file, vm)) {
		return false;
	}

	bool result = true;

	std::string species_name = xml_node.get<std::string>("<xmlattr>.species_name");
	use_population_average = xml_node.get<bool>("<xmlattr>.use_population_average", false);
	use_log_ratio = xml_node.get<bool>("<xmlattr>.use_log_ratio", false);
	include_only_cells_that_went_through_mitosis = xml_node.get<bool>("<xmlattr>.include_only_cells_that_went_through_mitosis", false);
	missing_simulation_time_stdev_str = xml_node.get<std::string>("<xmlattr>.missing_simulation_time_stdev", "");

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

	std::string time_dimension_name;
	size_t num_dimensions = 0;
	result &= data_file.GetDimensionName(experiment->GetName(), data_name, 0, time_dimension_name);
	result &= data_file.GetDimensionCount(experiment->GetName(), data_name, &num_dimensions);

	// Retrieve data from the datafile
	size_t num_timepoints;
	result &= data_file.GetDimensionSize(experiment->GetName(), time_dimension_name, &num_timepoints);
	result &= data_file.GetValues(experiment->GetName(), time_dimension_name, 0, num_timepoints, timepoints);

	if (!result) {
		return false;
	}

	std::string use_only_cell_ix_str = vm["cellpop.use_only_cell_ix"].as<std::string>();
	std::vector<std::string> use_only_cell_ix_tokens;
	bcm3::tokenize(use_only_cell_ix_str, use_only_cell_ix_tokens, ",");

	// Check whether there is data for just one cell or an average (one dimension) or a collection of cells (two+ dimensions)
	if (num_dimensions == 1) {
		observed_data.resize(1);
		VectorReal od;
		result &= data_file.GetValues(experiment->GetName(), data_name, 0, num_timepoints, od);
		observed_data[0] = od;
		if (!result) {
			return false;
		}
	} else if (num_dimensions == 2) {
		size_t num_cells;
		std::string cell_dimension_name;
		result &= data_file.GetDimensionName(experiment->GetName(), data_name, 1, cell_dimension_name);
		result &= data_file.GetDimensionSize(experiment->GetName(), cell_dimension_name, &num_cells);

		if (use_only_cell_ix_str == "-1") {
			observed_data.resize(num_cells);
			for (size_t i = 0; i < num_cells; i++) {
				VectorReal od;
				result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, i, num_timepoints, od);
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
					result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, ix, num_timepoints, od);
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

		if (use_only_cell_ix_str == "-1") {
			observed_data.resize(num_cells);
			for (size_t i = 0; i < num_cells; i++) {
				observed_data[i] = MatrixReal::Zero(num_timepoints, num_markers);
				for (size_t j = 0; j < num_markers; j++) {
					VectorReal od;
					result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, i, j, num_timepoints, od);
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
					observed_data[i] = MatrixReal::Zero(num_timepoints, num_markers);
					for (size_t j = 0; j < num_markers; j++) {
						VectorReal od;
						result &= data_file.GetValuesDim1(experiment->GetName(), data_name, 0, ix, j, num_timepoints, od);
						observed_data[i].col(j) = od;
					}
				}
			}
		}
	}

	if (num_dimensions > 1 && use_population_average) {
		LOGERROR("use_population_average=true has been specified, but the data reference includes data for more than 1 cell");
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

	if (use_population_average) {
		population_average.setConstant(num_timepoints, species_names.size(), 0.0);
		for (size_t i = 0; i < parallel_population_averages.size(); i++) {
			parallel_population_averages[i] = population_average;
		}
	} else {
		if (experiment->GetMaxNumberOfCells() < observed_data.size()) {
			LOGERROR("Maximum number of simulated cells (%zu) in the experiment is not sufficient for the amount of cells in the data (%zu)", experiment->GetMaxNumberOfCells(), observed_data.size());
			return false;
		}
		cell_trajectories.resize(experiment->GetMaxNumberOfCells(), MatrixReal::Constant(num_timepoints, species_names.size(), std::numeric_limits<Real>::quiet_NaN()));
		matched_trajectories.resize(observed_data.size(), MatrixReal::Constant(observed_data[0].rows(), observed_data[0].cols(), std::numeric_limits<Real>::quiet_NaN()));
	}

	// Do we have parent information?
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

	// Request the experiment to provide the required information during the simulation
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
					for (size_t i = 0; i < num_timepoints; i++) {
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
					for (size_t i = 0; i < num_timepoints; i++) {
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
				for (size_t i = 0; i < num_timepoints; i++) {
					experiment->AddSimulationTimepoints(this, timepoints(i), i, species_ix, synchronize);
				}
			}
			species_map[species_ix].push_back(j);
		}
	}

	trajectory_matching.resize(cell_trajectories.size());

	if (synchronize != ESynchronizeCellTrajectory::None) {
		// If we're going to synchronize, make sure we simulate the cells as long as the full duration of the time course
		// (this is needed in case there are negative timepoints)
		Real last_tp = timepoints(timepoints.size() - 1);
		Real full_duration = last_tp - timepoints(0);
		if (full_duration > last_tp) {
			experiment->AddSimulationTimepoints(this, full_duration, num_timepoints+1, std::numeric_limits<size_t>::max(), synchronize);
		}
	}

	return result;
}

bool DataLikelihoodTimeCourse::PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names)
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

	return true;
}

void DataLikelihoodTimeCourse::Reset()
{
	for (size_t i = 0; i < cell_trajectories.size(); i++) {
		cell_trajectories[i].setConstant(std::numeric_limits<Real>::quiet_NaN());
	}
	if (!use_population_average) {
		for (size_t i = 0; i < observed_data.size(); i++) {
			matched_trajectories[i].setConstant(std::numeric_limits<Real>::quiet_NaN());
		}
	}
	if (use_population_average) {
		population_average.setConstant(0.0);
		for (size_t i = 0; i < parallel_population_averages.size(); i++) {
			parallel_population_averages[i].setConstant(0.0);
		}
	}
	simulated_cell_parents.clear();
	simulated_cell_children.clear();
}

bool DataLikelihoodTimeCourse::Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp)
{
	Real stdev = GetCurrentSTDev(transformed_values, non_sampled_parameters);
	Real data_offset = GetCurrentDataOffset(transformed_values, non_sampled_parameters);
	Real data_scale = GetCurrentDataScale(transformed_values, non_sampled_parameters);

	Real missing_simulation_time_stdev = 1.0;
	if (fixed_missing_simulation_time_stdev_ix != std::numeric_limits<size_t>::max()) {
		missing_simulation_time_stdev = transformed_values[fixed_missing_simulation_time_stdev_ix];
	} else if (fixed_missing_simulation_time_stdev_non_sampled_ix != std::numeric_limits<size_t>::max()) {
		missing_simulation_time_stdev = non_sampled_parameters[fixed_missing_simulation_time_stdev_non_sampled_ix];
	} else if (!std::isnan(fixed_missing_simulation_time_stdev)) {
		missing_simulation_time_stdev = fixed_missing_simulation_time_stdev;
	}

	if (use_population_average) {
		// If parallel_population_averages.size() >= 1, then add the parallel evaluations
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

		population_average.array() *= data_scale;
		population_average.array() += data_offset;
	} else {
		for (size_t i = 0; i < cell_trajectories.size(); i++) {
			cell_trajectories[i].array() *= data_scale;
			cell_trajectories[i].array() += data_offset;
		}
	}

	logp = 0.0;

	if (use_population_average) {
		ASSERT(observed_data.size() == 1);
		size_t num_missing_points = 0;
		for (int i = 0; i < timepoints.size(); i++) {
			if (!std::isnan(observed_data[0](i))) {
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
					Real time_offset = std::min(std::abs(timepoints(i) - cell_trajectories_firstnonan),
												std::abs(timepoints(i) - cell_trajectories_lastnonan));
					if (error_model == ErrorModel::Normal) {
						logp += bcm3::LogPdfNormal(time_offset, 0, missing_simulation_time_stdev);
					} else if (error_model == ErrorModel::StudentT4) {
						logp += bcm3::LogPdfTnu4(time_offset, 0, missing_simulation_time_stdev);
					} else {
						assert(false);
						logp = std::numeric_limits<Real>::quiet_NaN();
					}
				} else {
					if (error_model == ErrorModel::Normal) {
						logp += bcm3::LogPdfNormal(observed_data[0](i), x, stdev);
					} else if (error_model == ErrorModel::StudentT4) {
						logp += bcm3::LogPdfTnu4(observed_data[0](i), x, stdev);
					} else {
						assert(false);
						logp = std::numeric_limits<Real>::quiet_NaN();
					}
				}
			}
		}
	} else {
		// Calculate likelihood of every observed parent cell (recursing into the children) against every simulated cell
		// TODO - relay from the likelihood file how many cells are going to be simulated
		int n = std::max(observed_cells_with_no_parents.size(), simulated_cell_parents.size());
		MatrixReal cell_likelihoods(n, n);
		cell_likelihoods.setZero();
		//std::vector< std::vector<double> > cell_likelihoods(observed_cells_with_no_parents.size(), std::vector<double>(simulated_cell_parents.size()));
		std::vector< std::vector< std::vector<int> > > matched_hierarchy;
		matched_hierarchy.resize(observed_cells_with_no_parents.size());

		std::string out = "cell_likelihoods:\n";
		Real maxdiff = 0.0;
		for (size_t i = 0; i < observed_cells_with_no_parents.size(); i++) {
			matched_hierarchy[i].resize(simulated_cell_parents.size());

			size_t finite_count = 0;
			size_t observed_cell_ix = observed_cells_with_no_parents[i];
			for (size_t j = 0; j < simulated_cell_parents.size(); j++) {
				matched_hierarchy[i][j].resize(cell_trajectories.size());

				if (simulated_cell_parents[j] == std::numeric_limits<size_t>::max()) {
					cell_likelihoods(i, j) = CalculateCellLikelihood(observed_cell_ix, j, stdev, missing_simulation_time_stdev, matched_hierarchy[i][j]);
					if (cell_likelihoods(i, j) != cell_likelihoods(i, j)) {
						logp = -std::numeric_limits<Real>::infinity();
						return false;
					}
					//cell_likelihoods[i][j] = -CalculateCellLikelihood(observed_cell_ix, j, stdev);
					//if (j > 0) {
					//	maxdiff = std::max(maxdiff, fabs(cell_likelihoods[i][j] - cell_likelihoods[i][0]));
					//}

					out += std::to_string(cell_likelihoods(i,j)) + ",\t";
					if (cell_likelihoods(i,j) > -std::numeric_limits<Real>::infinity()) {
					//if (cell_likelihoods[i][j] > -std::numeric_limits<Real>::infinity()) {
					//if (cell_likelihoods[i][j] < 1e7) {
						finite_count++;
					}
				} else {
					cell_likelihoods(i, j) = -std::numeric_limits<Real>::infinity();
				}
			}
			if (finite_count < observed_cells_with_no_parents.size()) {
				// Too many infinite likelihoods for this observed cell; won't be able to match all cells so no need to go further and HA doesn't like it
				logp = -std::numeric_limits<Real>::infinity();
				return true;
			}
			out += "\n";
		}
		//LOG(out.c_str());

		MatrixReal matching;
		if (!HungarianMatching(cell_likelihoods, matching, HUNGARIAN_MATCH_MAX)) {
			LOGERROR("Could not find matching");
			LOGERROR(out.c_str());
			//std::cout << matching << std::endl;
			logp = -std::numeric_limits<Real>::infinity();
			return true;
		}
		//std::cout << "Likelihoods:" << std::endl;
		//std::cout << cell_likelihoods << std::endl;
		//std::cout << "Matching:" << std::endl;
		//std::cout << matching << std::endl;

		for (int i = 0; i < observed_cells_with_no_parents.size(); i++) {
			for (int j = 0; j < simulated_cell_parents.size(); j++) {
				if (matching(i, j) == 1) {
					logp += cell_likelihoods(i, j);

					size_t oi = observed_cells_with_no_parents[i];
					matched_trajectories[oi] = cell_trajectories[j];
					trajectory_matching[oi] = j;

					if (!observed_children.empty()) {
						for (auto ci = observed_children[oi].begin(); ci != observed_children[oi].end(); ci++) {
							size_t child = matched_hierarchy[i][j][*ci];
							matched_trajectories[*ci] = cell_trajectories[child];
							trajectory_matching[*ci] = child;
						}
					}
				}
			}
		}

		// Fill in the matched data for output storage
	}

	logp *= weight;

	return true;
}

bool DataLikelihoodTimeCourse::NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis)
{
	if (species_ix == std::numeric_limits<size_t>::max()) {
		return true;
	}

	std::vector<size_t>& our_species_ixs = species_map[species_ix];
	for (size_t i = 0; i < our_species_ixs.size(); i++) {
		size_t our_species_ix = our_species_ixs[i];

		if (use_population_average) {
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
		} else {
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

Real DataLikelihoodTimeCourse::CalculateCellLikelihood(size_t observed_cell_ix, size_t simulated_cell_ix, Real stdev, Real missing_simulation_time_stdev, std::vector<int>& matched_hierarchy)
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
				cell_logp += CalculateCellLikelihood(*ci, std::numeric_limits<size_t>::max(), stdev, missing_simulation_time_stdev, matched_hierarchy);
			}
		}
	} else {
		for (int l = 0; l < species_names.size(); l++) {
			for (int k = 0; k < timepoints.size(); k++) {
				Real y = observed_data[observed_cell_ix](k, l);
				if (std::isnan(y)) {
					continue;
				}

				Real x = cell_trajectories[simulated_cell_ix](k, l);
				if (std::isnan(x)) {
					cell_logp += CalculateMissingValueLikelihood(simulated_cell_ix, k, l, missing_simulation_time_stdev);
				} else {
					if (error_model == ErrorModel::Normal) {
						cell_logp += bcm3::LogPdfNormal(y, x, stdev);
					} else if (error_model == ErrorModel::StudentT4) {
						cell_logp += bcm3::LogPdfTnu4(y, x, stdev);
					} else {
						assert(false);
						cell_logp = std::numeric_limits<Real>::quiet_NaN();
					}

					if (cell_logp == -std::numeric_limits<Real>::infinity()) {
						// Infinity - no need to go further, and the HA doesn't like it
						return cell_logp;
					}
				}
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
					children_matching(0, cii) = CalculateCellLikelihood(*ci, simulated_children.first, stdev, missing_simulation_time_stdev, matched_hierarchy);
					children_matching(1, cii) = CalculateCellLikelihood(*ci, simulated_children.second, stdev, missing_simulation_time_stdev, matched_hierarchy);

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
						matched_hierarchy[*observed_children[observed_cell_ix].begin()] = simulated_children.first;
					} else {
						cell_logp += children_matching(1, 0);
						matched_hierarchy[*observed_children[observed_cell_ix].begin()] = simulated_children.second;
					}
				} else {
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
				}
			} else {
				for (std::set<size_t>::const_iterator ci = observed_children[observed_cell_ix].begin(); ci != observed_children[observed_cell_ix].end(); ++ci) {
					cell_logp = CalculateCellLikelihood(*ci, std::numeric_limits<size_t>::max(), stdev, missing_simulation_time_stdev, matched_hierarchy);
				}
			}
		}
	}

	return cell_logp;
}

Real DataLikelihoodTimeCourse::CalculateMissingValueLikelihood(size_t simulated_cell_ix, int timepoint_ix, int species_ix, Real missing_simulation_time_stdev)
{
	if (simulated_cell_ix == std::numeric_limits<size_t>::max()) {
		if (error_model == ErrorModel::Normal) {
			return bcm3::LogPdfNormal(timepoints(timepoint_ix), 0, missing_simulation_time_stdev);
		} else if (error_model == ErrorModel::StudentT4) {
			return bcm3::LogPdfTnu4(timepoints(timepoint_ix), 0, missing_simulation_time_stdev);
		} else {
			assert(false);
			return std::numeric_limits<Real>::quiet_NaN();
		}
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
		Real time_offset = std::min(std::abs(timepoints(timepoint_ix) - cell_trajectories_firstnonan),
			std::abs(timepoints(timepoint_ix) - cell_trajectories_lastnonan));
		if (error_model == ErrorModel::Normal) {
			return bcm3::LogPdfNormal(time_offset, 0, missing_simulation_time_stdev);
		} else if (error_model == ErrorModel::StudentT4) {
			return bcm3::LogPdfTnu4(time_offset, 0, missing_simulation_time_stdev);
		} else {
			assert(false);
			return std::numeric_limits<Real>::quiet_NaN();
		}
	}
}
