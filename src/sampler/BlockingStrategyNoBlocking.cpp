#include "Utils.h"
#include "BlockingStrategyNoBlocking.h"

namespace bcm3 {

	std::vector< std::vector<ptrdiff_t> > BlockingStrategyNoBlocking::GetBlocks(const std::unique_ptr<SampleHistory>& sample_history, const std::shared_ptr<SampleHistoryClustering> clustering)
	{
		std::vector< std::vector<ptrdiff_t> > blocks(num_variables);
		for (ptrdiff_t i = 0; i < num_variables; i++) {
			blocks[i].resize(1);
			blocks[i][0] = i;
		}
		return blocks;
	}

}
