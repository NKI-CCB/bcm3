#pragma once

#include "Likelihood.h"
#include "Prior.h"
#include "ProgressIndicator.h"
#include "RNG.h"
#include "SampleHandler.h"
#include "TaskManager.h"

#include <boost/program_options.hpp>

namespace bcm3 {

	class Sampler
	{
	public:
		Sampler(size_t threads, size_t max_memory_use);
		virtual ~Sampler();

		void SetVariableSet(std::shared_ptr<VariableSet> varset);
		void SetPrior(std::shared_ptr<Prior> prior);
		void SetLikelihood(std::shared_ptr<Likelihood> likelihood, size_t threadix);
		void SetOutputPath(const std::string& path);
		void SetProgressIndicator(std::shared_ptr<ProgressIndicator> pi);

		void AddSampleHandler(std::shared_ptr<SampleHandler> handler);

		virtual bool LoadSettings(const boost::program_options::variables_map& vm);
		virtual bool Initialize();

		bool Run();

		inline size_t GetNumOutputSamples() const { return num_samples; }
		inline const VectorReal& GetOutputTemperatures() const { return temperatures; }

		static void AddOptionsDescription(boost::program_options::options_description& pod);

	protected:
		bool EvaluatePrior(size_t threadix, const VectorReal& values, Real& lprior);
		bool EvaluateLikelihood(size_t threadix, const VectorReal& values, Real& llh);
		bool EvaluatePriorLikelihood(size_t threadix, const VectorReal& values, Real& lprior, Real& llh);
		void UpdateProgress(Real p, bool force);

		virtual bool RunImpl() = 0;
		virtual void LogStatistics() {}

		// Target distribution
		std::shared_ptr<VariableSet> varset;
		std::shared_ptr<Prior> prior;
		std::vector< std::shared_ptr<Likelihood> > parallel_likelihoods;

		struct DirichletConstraint
		{
			ptrdiff_t residual_ix;
			std::vector<ptrdiff_t> var_ix;
		};
		std::vector<DirichletConstraint> dirichlet_constraints;

		// Sampling settings
		size_t num_threads;
		size_t num_samples;
		size_t use_every_nth;
		size_t num_variables;
		unsigned long long rngseed;
		std::string output_path;
		VectorReal temperatures;

		// Runtime variables
		std::vector< std::shared_ptr<SampleHandler> > sample_handlers;
		std::weak_ptr<ProgressIndicator> progress;
		int logged_progress;
		RNG rng;
		std::atomic<size_t> num_likelihood_evaluations;

		// Parallelization
		std::unique_ptr<TaskManager> task_manager;
		struct AsyncParams {
			RNG rng;
		};
		std::vector<AsyncParams> async;
	};

}
