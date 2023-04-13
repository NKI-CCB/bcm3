#include "Utils.h"
#include "BlockingStrategyOneBlock.h"

namespace bcm3 {

	std::vector< std::vector<ptrdiff_t> > BlockingStrategyOneBlock::GetBlocks(const std::unique_ptr<SampleHistory>& sample_history)
	{
		std::vector< std::vector<ptrdiff_t> > blocks(1);
		blocks[0].resize(num_variables);
		for (ptrdiff_t i = 0; i < num_variables; i++) {
			blocks[0][i] = i;
		}
		return blocks;
	}

}
