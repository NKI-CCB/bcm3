#pragma once

#include "BlockingStrategy.h"

namespace bcm3 {

	class BlockingStrategyNoBlocking : public BlockingStrategy
	{
	public:
		virtual std::vector< std::vector<ptrdiff_t> > GetBlocks(const std::unique_ptr<SampleHistory>& sample_history);
	};

}
