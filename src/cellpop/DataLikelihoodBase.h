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

	virtual bool NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t completed_population_size, size_t parallel_evaluation_ix, bool entered_mitosis) { return true; }
	virtual void NotifyStartingCells(size_t cell_ix) {}
	virtual void NotifyParents(size_t parent, size_t child) {}
	virtual void NotifyDuration(size_t cell_ix, Real duration) {}

protected:
	std::string data_name;
	Real weight;

	std::string stdev_str;
	std::string offset_str;
	std::string scale_str;

	size_t offset_ix;
	size_t scale_ix;
	size_t stdev_ix;
	size_t non_sampled_offset_ix;
	size_t non_sampled_scale_ix;
	size_t non_sampled_stdev_ix;
	Real fixed_offset_value;
	Real fixed_scale_value;
	Real fixed_stdev_value;
	Real fixed_abs_value;
	Real fixed_rel_value;
	Real stdev_multiplication_factor;
};
