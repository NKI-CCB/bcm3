#pragma once

#include "BlockingStrategy.h"
#include "Proposal.h"
#include "SampleHistory.h"
#include <boost/dynamic_bitset.hpp>

namespace bcm3 {

	class SamplerPT;

	class SamplerPTChain
	{
	public:
		SamplerPTChain(SamplerPT* sampler);

		bool Initialize(size_t history_size, size_t history_subsampling);

		void SetTemperature(Real temperature, bool is_fixed);
		bool AdaptProposal(size_t thread);
		void AdaptProposalAsync();
		bool AdaptProposalWait();
		bool FindStartingPosition();
		bool MutateMove(size_t thread);
		void MutateMoveAsync();
		bool MutateMoveWait();
		bool ExchangeMove(SamplerPTChain& other);
		void LogStatistics() const;
		void LogProposalInfo() const;

		inline Real GetLogPrior() const { return lprior; }
		inline Real GetLogLikelihood() const { return llh; }
		inline Real GetLogPowerPosterior() const { return lpowerposterior; }
		inline const VectorReal& GetVarValues() const { return current_var_values; }
		inline Real GetTemperature() const { return temperature; }
		inline bool GetIsFixedTemperature() const { return temperature_fixed; }

	private:
		bool AsyncDoMutateMove(void* user, size_t thread);
		bool AsyncDoAdaptProposal(void* user, size_t thread);
		std::unique_ptr<Proposal> CreateProposalInstance(std::vector<ptrdiff_t> variable_indices, RNG& rng);
		bool TestSample(Real new_lpowerposterior, Real log_mh_ratio, RNG& rng);
		void AcceptMutate(bool accept, const VectorReal& nv, Real lprior, Real llh, Real lpowerposterior);
		void AddCurrentSampleToHistory();

		// Settings
		SamplerPT* sampler;
		Real temperature;
		bool temperature_fixed;

		// Run-time variables
		std::unique_ptr<BlockingStrategy> blocking_strategy;
		std::unique_ptr<SampleHistory> sample_history;
		VectorReal current_var_values;
		Real lprior;
		Real llh;
		Real lpowerposterior;

		size_t attempted_mutate;
		size_t attempted_exchange;
		size_t accepted_mutate;
		size_t accepted_exchange;

		struct Block
		{
			std::vector<ptrdiff_t> variable_indices;
			std::unique_ptr<Proposal> proposal;
		};
		std::vector<Block> variable_blocks;

		// For parallelization
		uint64 async_task;
		bool async_burnin;
		std::atomic<int> async_failure;
	};

}
