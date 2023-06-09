#pragma once

#include "fISAExperiment.h"

class SignalingNetwork;
class fISAExperimentSingleCondition;

class fISAExperimentIncucyteSequential : public fISAExperiment
{
public:
	fISAExperimentIncucyteSequential();
	virtual ~fISAExperimentIncucyteSequential();
	
	virtual bool StartEvaluateLogProbability(const VectorReal& values, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);
	virtual bool FinishEvaluateLogProbability(const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);

	virtual size_t GetNumData() const { return drug_concentrations.size() * 2; }
	virtual size_t GetNumReplicates(size_t data_ix) const { return 1; }
	virtual void GetObservedData(size_t data_ix, MatrixReal& out_values) const;
	virtual void GetModeledData(size_t threadix, size_t data_ix, VectorReal& out_values) const;
	virtual void GetModeledActivities(size_t threadix, MatrixReal& out_values) const;

private:
	virtual bool LoadTypeSpecificNodes(const boost::property_tree::ptree& xml_node, const bcm3::NetCDFDataFile& datafile);
	virtual bool InitializeParallelData(size_t evaluation_threads);
	virtual bool ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile);
	bool EvaluateCellLine(void* cell_line_ix, size_t eval_thread_ix);

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

	DataPart data;
	MatrixReal response_summary;

	std::string drug_species_name;
	size_t drug_model_ix;
	VectorReal drug_concentrations;

	std::vector<ParallelDataBase> parallel_data;

	std::vector< std::vector<VectorReal> > activities;
	std::vector< std::vector<VectorReal> > modeled_data_values;

	const fISAExperimentSingleCondition* relative_experiment;
};
