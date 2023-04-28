#pragma once

#include "Sampler.h"

namespace bcm3 {

	class SamplerPTChain;

	class SamplerPT : public Sampler
	{
	public:
		SamplerPT(size_t threads, size_t max_memory_use);
		~SamplerPT();

		virtual bool LoadSettings(const boost::program_options::variables_map& vm);
		virtual bool Initialize();

		static void AddOptionsDescription(boost::program_options::options_description& pod);

	protected:
		virtual bool RunImpl();
		virtual void LogStatistics();

		void DoExchangeMove(size_t sample_ix);
		bool DoMutateMove();
		void EmitSample(size_t sample_ix);

		enum class ESwappingScheme {
			StochasticRandom,
			StochasticEvenOdd,
			DeterministicEvenOdd
		};

		// Sampling settings
		std::string blocking_strategy;
		std::string proposal_type;
		bool proposal_transform_to_unbounded;
		bool output_proposal_adaptation;

		ESwappingScheme swapping_scheme;
		Real exchange_probability;
		size_t num_exploration_steps;
		size_t max_history_size;
		size_t adapt_proposal_samples;
		size_t adapt_proposal_times;
		size_t adapt_proposal_max_history_samples;
		Real proposal_scaling_learning_rate;
		size_t proposal_scaling_ema_period;
		size_t stop_proposal_scaling;
		Real target_acceptance_rate;

		// Run-time variables
		std::vector< std::unique_ptr<SamplerPTChain> > chains;
		size_t proposal_adaptations_done;
		bool proposal_scaling_adaptations_done;
		bool previous_swap_even;

		friend class SamplerPTChain;
	};

}
