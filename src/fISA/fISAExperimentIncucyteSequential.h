#pragma once

#include "fISAExperiment.h"

#if TODO

class SignalingNetwork;
class fISAExperimentSingleCondition;

class fISAExperimentIncucyteSequential : public fISAExperiment
{
public:
	fISAExperimentIncucyteSequential(size_t numthreads);
	virtual ~fISAExperimentIncucyteSequential();
	
	virtual bool EvaluateLogProbability(size_t threadix, const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);

private:
	virtual bool LoadTypeSpecificNodes(const boost::property_tree::ptree& xml_node, const bcm3::NetCDFDataFile& datafile);
	virtual bool InitializeParallelData();
	virtual bool ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile);

	struct DataPart : public DataPartBase
	{
		struct tdist {
			Real mup;
			Real mua;
			MatrixReal cov;
			MatrixReal invcov;
		};
		std::vector<tdist> tdists;

		struct gmm {
			Real mup[3];
			Real mua[3];
			MatrixReal cov[3];
			MatrixReal invcov[3];
			Real weights[3];
			Real logncweight[3];
		};
		std::vector<gmm> gmms;
	};
	
	struct ParallelData : public ParallelDataBase {
		std::vector< std::vector<VectorReal> > activities;
		std::vector< std::vector<VectorReal> > modeled_data_values;
	};

	DataPart data;
	MatrixReal response_summary;
	Real data_weight;

	std::string drug_species_name;
	size_t drug_model_ix;
	VectorReal drug_concentrations;

	std::vector<ParallelData> parallel_data;
	fISAExperimentSingleCondition* relative_experiment;
};

#endif
