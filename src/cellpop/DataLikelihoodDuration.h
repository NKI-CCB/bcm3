#pragma once

#include "DataLikelihoodBase.h"

class DataLikelihoodDuration : public DataLikelihoodBase
{
public:
	DataLikelihoodDuration(size_t parallel_evaluations);
	virtual ~DataLikelihoodDuration();

	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file);
	virtual bool PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names);
	virtual void Reset();
	virtual bool Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp);
	virtual void NotifyDuration(size_t cell_ix, Real duration);

	inline const VectorReal& GetObservedData() const { return observed_durations; }
	inline const VectorReal& GetSimulatedData() const { return matched_durations; }

private:
	VectorReal observed_durations;
	VectorReal simulated_durations;
	VectorReal matched_durations;
	std::vector<int> used_durations;
};
