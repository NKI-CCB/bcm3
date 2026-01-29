#pragma once

#include "DataLikelihoodBase.h"
#include "hungarian.h"

class DataLikelihoodTimeCourseBase : public DataLikelihoodBase
{
public:
	DataLikelihoodTimeCourseBase();
	virtual ~DataLikelihoodTimeCourseBase();

	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names);

	inline const VectorReal& GetTimepoints() const { return timepoints; }

protected:
	bool RequestSimulationInfo(Experiment* experiment, ESynchronizeCellTrajectory synchronize);
	bool PrepateEvaluation(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp);
	Real EvaluateValue(Real simulated, Real observed, size_t stdev_ix);
	Real EvaluateMissingValue(Real time_offset);
	void OptimizeOffsetScale(const MatrixReal& observed, const MatrixReal& simulated, int col_ix, Real& offset, Real& scale) const;

	// Data
	VectorReal timepoints;

	// Settings
	bool include_only_cells_that_went_through_mitosis;
	bool use_log_ratio;
	bool optimize_offset_scale;
	Real optimize_offset_min;
	Real optimize_offset_max;
	Real optimize_scale_min;
	Real optimize_scale_max;
	bool use_signal_saturation;
	Real saturation_scale_fixed_value;
	std::string saturation_scale_str;
	size_t saturation_scale_ix;

	std::string missing_simulation_time_stdev_str;
	size_t fixed_missing_simulation_time_stdev_ix;
	size_t fixed_missing_simulation_time_stdev_non_sampled_ix;
	Real fixed_missing_simulation_time_stdev;

	std::vector<std::string> species_names;
	std::map< size_t, std::vector<size_t> > species_map;
	std::vector<size_t> log_ratio_numerator_ix;

	// Runtime
	std::vector<Real> stdevs;
	std::vector<Real> proportional_stdevs;
	std::vector<Real> data_offsets;
	std::vector<Real> data_scales;
	Real saturation_scale;
	Real missing_simulation_time_stdev;
	std::vector<Real> minus_log_sigma;
	std::vector<Real> inv_two_sigma_sq;
};
