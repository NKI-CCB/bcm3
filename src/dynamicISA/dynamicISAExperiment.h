#pragma once

#include "VariableSet.h"

class SignalingModel;

class dynamicISAExperiment
{
public:
	dynamicISAExperiment();
	~dynamicISAExperiment();

	bool Load(std::shared_ptr<const bcm3::VariableSet> varset, std::shared_ptr<SignalingModel> model, boost::property_tree::ptree xml_node);
	bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp);

	inline size_t GetNumData() const { return data.size(); }
	inline size_t GetNumConditions() const { return condition_names.size(); }
	inline size_t GetNumReplicates(size_t data_ix) const { return data[data_ix].values.size(); }
	void GetObservedData(size_t data_ix, MatrixReal& out_values) const;
	void GetModeledData(size_t data_ix, VectorReal& out_values) const;
	void GetModeledActivities(MatrixReal& out_values) const;

	inline const std::string& GetName() const { return name; }

private:
	std::string name;
	std::shared_ptr<SignalingModel> model;

	struct ModelConditions
	{
		ModelConditions();

		std::string model_name;
		VectorReal values;
		size_t model_ix;
		size_t parameter_ix;
	};

	std::vector<std::string> condition_names;
	std::vector<ModelConditions> condition_specifications;

	std::vector<size_t> base_param_map;
	std::vector<size_t> decay_param_map;
	std::vector<size_t> strength_param_map;

	VectorReal base_param_values;
	VectorReal decay_param_values;
	VectorReal strength_param_values;
	VectorReal inhib_param_values;
	VectorReal expression_mixing_param_values;
	VectorReal fixed_species_values;
	VectorReal expression_levels;

	struct Data
	{
		std::string model_name;
		std::vector<VectorReal> values;
		size_t model_ix;
		Real timepoint;

		size_t parameter_base_ix;
		size_t parameter_scale_ix;
		size_t parameter_sd_ix;
		size_t parameter_sd_incr_ix;
	};
	std::vector<Data> data;

	VectorReal data_timepoints;
	std::vector< std::tuple<size_t, Real, Real> > treatments;

	MatrixReal modeled_activities;
	std::vector<VectorReal> modeled_data;
};
