#pragma once

#include "fISAExperiment.h"

#if TODO

class SignalingNetwork;

class fISAExperimentDrugRange : public fISAExperiment
{
public:
	fISAExperimentDrugRange(size_t numthreads);
	virtual ~fISAExperimentDrugRange();
	
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);

private:
	virtual bool LoadTypeSpecificNodes(const boost::property_tree::ptree& xml_node, const bcm3::NetCDFDataFile& datafile);
	virtual bool InitializeParallelData();
	virtual bool ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile);

	struct DataPart : public DataPartBase
	{
		std::vector<MatrixReal> data;
	};
	
	struct ParallelData : public ParallelDataBase {
		std::vector< std::vector<VectorReal> > activities;
		std::vector< std::vector<VectorReal> > modeled_data_values;
	};

	std::vector<DataPart> data;

	std::string drug_species_name;
	size_t drug_model_ix;
	VectorReal drug_concentrations;

	std::vector<ParallelData> parallel_data;
};

#endif
