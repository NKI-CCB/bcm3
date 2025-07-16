#pragma once

#include "DataLikelihoodBase.h"

class DataLikelihoodTimeCourse : public DataLikelihoodBase
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
	inline const MatrixReal& GetSimulatedData(size_t ix) const { return use_population_average ? population_average : matched_trajectories[ix]; }

private:
	Real CalculateCellLikelihood(size_t observed_cell_ix, size_t simulated_cell_ix, const std::vector<Real>& stdevs, Real missing_simulation_time_stdev, std::vector<int>& matched_hierarchy);
	Real CalculateMissingValueLikelihood(size_t simulated_cell_ix, int timepoint_ix, int species_ix, Real missing_simulation_time_stdev);
	void OptimizeOffsetScale(MatrixReal& observed, MatrixReal& simulated, int col_ix, Real& offset, Real& scale);

	bool use_population_average;
	bool use_log_ratio;
	bool include_only_cells_that_went_through_mitosis;
	bool population_relative_to_time_average;
	ESynchronizeCellTrajectory synchronize;
	std::string missing_simulation_time_stdev_str;
	size_t fixed_missing_simulation_time_stdev_ix;
	size_t fixed_missing_simulation_time_stdev_non_sampled_ix;
	Real fixed_missing_simulation_time_stdev;

	bool use_signal_saturation;
	bool per_cell_optimize_offset_scale;
	Real saturation_scale;
	std::string saturation_scale_str;
	size_t saturation_scale_ix;

	VectorReal timepoints;
	std::vector<MatrixReal> observed_data;
	std::vector<size_t> observed_cells_with_no_parents;
	std::vector< std::set<size_t> > observed_children;

	MatrixReal population_average;
	std::vector<MatrixReal> parallel_population_averages;

	std::vector<MatrixReal> cell_trajectories;
	std::vector<size_t> simulated_cell_parents;
	std::vector< std::pair<size_t, size_t> > simulated_cell_children;
	std::vector<size_t> trajectory_matching;
	std::vector<MatrixReal> matched_trajectories;

	std::vector<std::string> species_names;
	std::map< size_t, std::vector<size_t> > species_map;
	std::vector<size_t> log_ratio_numerator_ix;
};
