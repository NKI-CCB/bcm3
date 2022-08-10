#pragma once

#include "Likelihood.h"
#include "Prior.h"
#include "ProgressIndicator.h"
#include "RNG.h"
#include "SamplerPTChain.h"
#include "TaskManager.h"

#include <boost/program_options.hpp>

namespace bcm3 {

class SampleHandler
{
public:
	virtual void ReceiveSample(const VectorReal& sample, Real lprior, Real llh, Real temperature) = 0;
};

class SamplerPT
{
public:
	SamplerPT(size_t threads, size_t max_memory_use);
	~SamplerPT();

	void SetVariableSet(std::shared_ptr<VariableSet> varset);
	void SetPrior(std::shared_ptr<Prior> prior);
	void SetLikelihood(std::shared_ptr<Likelihood> likelihood, size_t threadix);
	void SetOutputPath(const std::string& path);
	void SetProgressIndicator(std::shared_ptr<ProgressIndicator> pi);

	bool LoadSettings(const boost::program_options::variables_map& vm);
	bool Initialize();
	void AddSampleHandler(std::shared_ptr<SampleHandler> handler);
	bool Run();

	inline size_t GetNumOutputSamples() const { return num_samples; }
	inline const VectorReal& GetOutputTemperatures() const { return fixed_temperatures; }

	static void AddOptionsDescription(boost::program_options::options_description& pod);

private:
	bool EvaluatePrior(size_t threadix, const VectorReal& values, Real& lprior);
	bool EvaluateLikelihood(size_t threadix, const VectorReal& values, Real& llh);
	bool EvaluatePriorLikelihood(size_t threadix, const VectorReal& values, Real& lprior, Real& llh);

	void DoExchangeMove(size_t sample_ix);
	bool DoMutateMove();
	void EmitSample(size_t sample_ix);
	void LogStatistics();
	void UpdateProgress(Real p, bool force);

	// Target distribution
	std::shared_ptr<VariableSet> varset;
	std::shared_ptr<Prior> prior;
	std::vector< std::shared_ptr<Likelihood> > parallel_likelihoods;

	// Sampling settings
	size_t numthreads;
	VectorReal fixed_temperatures;
	size_t num_samples;
	size_t use_every_nth;
	size_t num_variables;
	std::string proposal_type;
	unsigned long long rngseed;
	std::string output_path;
	bool output_proposal_adaptation;

	Real exchange_probability;
	size_t history_max_samples;
	size_t adapt_proposal_samples;
	size_t adapt_proposal_times;
	size_t adapt_proposal_max_samples;
	Real proposal_scaling_learning_rate;
	size_t proposal_scaling_ema_period;
	size_t stop_proposal_scaling;
	Real target_acceptance_rate;
	Real proposal_t_df;

	int clustered_blocking_nn;
	int clustered_blocking_nn2;
	int clustered_blocking_n_clusters;

	// Run-time variables
	RNG rng;
	VectorReal current_temperatures;
	std::vector< std::unique_ptr<SamplerPTChain> > chains;
	std::vector< std::shared_ptr<SampleHandler> > sample_handlers;
	std::weak_ptr<ProgressIndicator> progress;
	std::atomic<size_t> num_likelihood_evaluations;
	size_t proposal_adaptations_done;
	bool proposal_scaling_adaptations_done;

	std::unique_ptr<TaskManager> task_manager;
	struct AsyncParams {
		RNG rng;
	};
	std::vector<AsyncParams> async;

	friend class SamplerPTChain;
};

}
