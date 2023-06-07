#pragma once

#include "fISAExperiment.h"

class SignalingNetwork;

class fISAExperimentSingleCondition : public fISAExperiment
{
public:
	fISAExperimentSingleCondition();
	virtual ~fISAExperimentSingleCondition();

	virtual bool StartEvaluateLogProbability(const VectorReal& values, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);
	virtual bool FinishEvaluateLogProbability(const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);

	virtual size_t GetNumData() const { return data.size(); }
	virtual size_t GetNumReplicates(size_t data_ix) const { return data[data_ix].data.size(); }
	virtual void GetObservedData(size_t data_ix, MatrixReal& out_values) const;
	virtual void GetModeledData(size_t threadix, size_t data_ix, VectorReal& out_values) const;
	virtual void GetModeledActivities(size_t threadix, MatrixReal& out_values) const;

protected:
	virtual bool InitializeParallelData(size_t evaluation_threads);
	virtual bool ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile);
	bool EvaluateCellLine(void* cell_line_ix, size_t eval_thread_ix);
	
	struct DataPart : public DataPartBase
	{
		std::vector<VectorReal> data;
	};
	
	struct ParallelData : public ParallelDataBase {
		std::vector<VectorReal> multi_activities;
		std::vector<VectorReal> multi_modeled_data_values;
		VectorReal multisolve_logp;
	};

	std::vector<DataPart> data;

	std::vector<ParallelData> parallel_data;
	std::vector<VectorReal> stored_activities;
	std::vector<VectorReal> stored_modeled_data_values;

	friend class fISAExperimentDrugRange;
	friend class fISAExperimentIncucyteSequential;
};
