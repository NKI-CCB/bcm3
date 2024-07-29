#pragma once

#include "DataLikelihoodBase.h"

class DataLikelihoodTimePoints : public DataLikelihoodBase
{
public:
	DataLikelihoodTimePoints(size_t parallel_evaluations);
	virtual ~DataLikelihoodTimePoints();

	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm);
	virtual void Reset();
	virtual bool Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp);
	virtual bool NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis);

	inline const VectorReal& GetTimepoints() const { return timepoints; }
	inline const MatrixReal& GetObservedData(size_t ix) const { return observed_data[ix]; }
	inline size_t GetNumSimulatedCells() const { return cell_trajectories.size(); }
	inline const MatrixReal& GetSimulatedData(size_t ix) const { return  matched_data[ix]; }

private:
	VectorReal timepoints;
	std::vector<MatrixReal> observed_data;
	std::vector<MatrixReal> cell_trajectories;
	std::vector<MatrixReal> matched_data;
	std::vector<std::string> species_names;
	std::map< size_t, std::vector<size_t> > species_map;
};
