#pragma once

#include "DataLikelihoodTimeCourseBase.h"
#include "hungarian.h"

class DataLikelihoodTimeCoursePopulationAverage : public DataLikelihoodTimeCourseBase
{
public:
	DataLikelihoodTimeCoursePopulationAverage(size_t parallel_evaluations);
	virtual ~DataLikelihoodTimeCoursePopulationAverage();

	virtual bool Load(const boost::property_tree::ptree& xml_node, Experiment* experiment, const bcm3::VariableSet& varset, const bcm3::NetCDFDataFile& data_file, const boost::program_options::variables_map& vm);
	virtual bool PostInitialize(const bcm3::VariableSet& varset, const std::vector<std::string>& non_sampled_parameter_names);
	virtual void Reset();
	virtual bool Evaluate(const VectorReal& values, const VectorReal& transformed_values, const VectorReal& non_sampled_parameters, Real& logp);
	virtual bool NotifySimulatedValue(size_t timepoint_ix, Real x, size_t species_ix, size_t cell_ix, size_t current_population_size, size_t mitotic_population_size, size_t parallel_evaluation_ix, bool entered_mitosis, bool is_newborn);
	virtual void NotifyStartingCells(size_t cell_ix);
	virtual void NotifyParents(size_t parent, size_t child);

	inline const MatrixReal& GetObservedData() const { return observed_data; }
	inline const MatrixReal& GetSimulatedData() const { return population_average; }

private:
	bool relative_to_time_average;

	MatrixReal observed_data;		// replicates x timepoints
	MatrixReal population_average;	// timepoints x species
	std::vector<MatrixReal> parallel_population_averages;
};
