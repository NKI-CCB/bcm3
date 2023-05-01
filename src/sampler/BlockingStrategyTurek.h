#pragma once

#include "BlockingStrategy.h"

namespace bcm3 {

	// Blocking strategy based on Turek et al, Automated Parameter Blocking for Efficient Markov Chain Monte Carlo Sampling, Bayesian Analysis 2017.
	class BlockingStrategyTurek : public BlockingStrategy
	{
	public:
		virtual std::vector< std::vector<ptrdiff_t> > GetBlocks(const std::unique_ptr<SampleHistory>& sample_history, const std::shared_ptr<SampleHistoryClustering> clustering);
	};

}
