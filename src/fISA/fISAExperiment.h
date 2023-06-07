#pragma once

#include <condition_variable>

#include "NetCDFDataFile.h"
#include "VariableSet.h"
#include "Spinlock.h"
#include "TaskManager.h"

class SignalingNetwork;

class fISAExperiment
{
public:
	virtual ~fISAExperiment();
	
	virtual bool StartEvaluateLogProbability(const VectorReal& values, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments) = 0;
	virtual bool FinishEvaluateLogProbability(const VectorReal& values, Real& logp, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments) = 0;

	static std::unique_ptr<fISAExperiment> Create(const boost::property_tree::ptree& xml_node, SignalingNetwork* network, std::shared_ptr<const bcm3::VariableSet> varset, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, size_t numthreads, bcm3::TaskManager* task_manager);
	inline const std::string& GetName() const { return Name; }

	size_t GetNumCellLines() const { return cell_lines.size(); }
	inline const std::string& GetCellLineName(size_t ix) const { return cell_lines[ix]; }
	virtual size_t GetNumData() const { return 0; }
	virtual size_t GetNumReplicates(size_t data_ix) const { return 0; }
	virtual void GetObservedData(size_t data_ix, MatrixReal& out_values) const {}
	virtual void GetModeledData(size_t threadix, size_t data_ix, VectorReal& out_values) const {}
	virtual void GetModeledActivities(size_t threadix, MatrixReal& out_values) const {}

#if 0
	void SetBootstrap(unsigned long rngseed);
	void SetSafeBayesIx(size_t ix);
	void AddSafeBayesCvIx(size_t ix);
#endif

protected:
	enum LikelihoodFunction
	{
		LF_Normal,
		LF_TruncatedNormal,
		LF_StudentT,
		LF_TruncatedStudentT,
		LF_Binomial,
		LF_Beta,
		LF_OrderedProbit,
		LF_OrderedRobit,

		LF_Invalid
	};

	struct ModelConditions
	{
		ModelConditions();

		std::string model_name;
		VectorReal values;
		size_t model_ix;
		size_t parameter_ix;
	};

	struct ExpressionLevel
	{
		ExpressionLevel();

		std::string model_name;
		VectorReal values;
		size_t model_ix;
		size_t base_parameter_ix;
		size_t scale_parameter_ix;
	};

	struct DataPartBase
	{
		DataPartBase();

		std::string model_name;
		std::string data_name;
		std::string base_scale_sd_suffix;

		LikelihoodFunction likelihood_fn;

		bool use_base;
		bool use_scale;
		bool fixed_base;
		bool fixed_sd;
		bool scale_var_with_mean;
		Real weight;
		
		size_t model_ix;
		size_t parameter_base_ix;
		size_t parameter_scale_ix;
		size_t parameter_sd_ix;
		size_t parameter_precision_ix;
		Real fixed_base_value;
		Real fixed_sd_value;
		Real fixed_sd_scale_value;
		Real fixed_precision_value;
		bool data_is_inactive_form;
		bool scale_per_cell_line;

		size_t expression_ix;
		size_t expression_ix2;
		size_t expression_sum_parameter_ix;

		size_t model_ix2;
		size_t data_sum_parameter_ix;
		size_t parameter_sd_base_ix;

		std::vector<size_t> ordered_probit_thresholds;

		std::string relative_to_experiment;
	};
	
	struct ParallelDataBase {
		std::vector<Real> data_logC;
		VectorReal expression;
		std::vector<VectorReal> relative_to_experiment_values;
	};
	
	fISAExperiment();
	bool Load(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);
	bool ParseDataPartBase(const boost::property_tree::ptree& xml_node, DataPartBase& p, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments);
	bool ParseDataFileReference(const std::string& data_name, std::string& data_name_without_index, std::vector<size_t>& data_indices, const bcm3::NetCDFDataFile& datafile);
	void PrepareActivitiesCalculation(VectorReal& activities, VectorReal& expression, const Real* transformed_values, size_t cell_ix) const;

	virtual bool LoadTypeSpecificNodes(const boost::property_tree::ptree& xml_node, const bcm3::NetCDFDataFile& datafile) { return true; }
	virtual bool InitializeParallelData(size_t evaluation_threads) = 0;
	virtual bool ParseDataNode(const boost::property_tree::ptree& xml_node, const std::vector< std::unique_ptr<fISAExperiment> >& other_experiments, const bcm3::NetCDFDataFile& datafile) = 0;

	SignalingNetwork* network;
	std::shared_ptr<const bcm3::VariableSet> varset;
	bcm3::TaskManager* task_manager;

	std::string Name;
	std::string DataFile;

	std::vector<std::string> cell_lines;
	std::vector<ExpressionLevel> expression_levels;
	std::vector<ModelConditions> conditions;
	std::vector<Real> transformed_values;

	VectorReal cell_line_logp;
	std::vector<uint64> evaluation_tasks;

#if 0
	std::vector<size_t> cell_lines_bootstrap_count;
	size_t safe_bayes_ix;
	std::set<size_t> safe_bayes_cv_ix;
#endif
};
