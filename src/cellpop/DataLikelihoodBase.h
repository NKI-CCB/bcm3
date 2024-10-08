#pragma once

#include "Experiment.h"
#include "NetCDFDataFile.h"
#include "VariableSet.h"

#include <boost/program_options.hpp>

class DataLikelihoodBase
{
public:
	DataLikelihoodBase();
	virtual ~DataLikelihoodBase() = 0;

	static std::unique_ptr<DataLikelihoodBase> Create(const boost::property_tree::ptree& xml_node, size_t parallel_evaluations);

	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names);
	virtual void Reset() {}
	virtual bool Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp) = 0;

	virtual bool NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis, bool is_newborn) { return true; }
	virtual void NotifyStartingCells(size_t cell_ix) {}
	virtual void NotifyParents(size_t parent, size_t child) {}
	virtual void NotifyDuration(size_t cell_ix, Real duration) {}

protected:
	Real GetCurrentSTDev(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, size_t i = 0);
	Real GetCurrentDataOffset(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, size_t i = 0);
	Real GetCurrentDataScale(const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, size_t i = 0);

	enum class ErrorModel {
		Normal,
		StudentT4,
	};
	std::string data_name;
	Real weight;

	std::string stdev_str;
	std::string offset_str;
	std::string scale_str;

	ErrorModel error_model;
	std::vector<size_t> offset_ix;
	std::vector<size_t> scale_ix;
	std::vector<size_t> stdev_ix;
	std::vector<size_t> non_sampled_offset_ix;
	std::vector<size_t> non_sampled_scale_ix;
	std::vector<size_t> non_sampled_stdev_ix;
	std::vector<Real> fixed_offset_value;
	std::vector<Real> fixed_scale_value;
	std::vector<Real> fixed_stdev_value;

private:
	bool ParseString(std::string str, size_t& var_ix, size_t& non_sampled_var_ix, Real& fixed_value, const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names);
};
