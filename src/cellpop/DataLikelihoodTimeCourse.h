#pragma once

#include "DataLikelihoodTimeCourseBase.h"
#include "hungarian.h"

class DataLikelihoodTimeCourse : public DataLikelihoodTimeCourseBase
{
public:
	DataLikelihoodTimeCourse(size_t parallel_evaluations);
	virtual ~DataLikelihoodTimeCourse();

	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names);
	virtual void Reset();
	virtual bool Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp);
	virtual bool NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis, bool is_newborn);
	virtual void NotifyStartingCells(size_t cell_ix);
	virtual void NotifyParents(size_t parent, size_t child);

	inline const std::vector<size_t>& GetTrajectoryMatching() const { return trajectory_matching; }
	inline const VectorReal& GetTimepoints() const { return timepoints; }
	inline const size_t GetNumObservedData() const { return observed_data.size(); }
	inline const MatrixReal& GetObservedData(size_t ix) const { return observed_data[ix]; }
	inline const MatrixReal& GetSimulatedData(size_t ix) const { return matched_trajectories[ix]; }

private:
	Real CalculateCellLikelihood(size_t observed_cell_ix, size_t simulated_cell_ix, const std::vector<Real>& stdevs, const std::vector<Real>& proportional_stdevs, Real missing_simulation_time_stdev, std::vector<int>* matched_hierarchy);
	Real CalculateMissingValueLikelihood(size_t simulated_cell_ix, int timepoint_ix, int species_ix, Real missing_simulation_time_stdev);

	// Data
	std::vector<MatrixReal> observed_data;				// vector of cells, each matrix is timepoint x marker
	std::vector<size_t> observed_cells_with_no_parents;
	std::vector< std::set<size_t> > observed_children;

	// Settings
	ESynchronizeCellTrajectory synchronize;

	// Runtime
	std::vector<MatrixReal> cell_trajectories;
	std::vector<size_t> simulated_cell_parents;
	std::vector< std::pair<size_t, size_t> > simulated_cell_children;
	std::vector<size_t> trajectory_matching;
	std::vector<MatrixReal> matched_trajectories;

	std::vector<WeightedBipartiteEdge> hungarian_matching_edges;
};
