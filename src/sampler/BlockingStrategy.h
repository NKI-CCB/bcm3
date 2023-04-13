#pragma once

namespace bcm3 {

	class SampleHistory;

	class BlockingStrategy
	{
	public:
		virtual void Initialize(size_t num_variables);
		virtual std::vector< std::vector<ptrdiff_t> > GetBlocks(const std::unique_ptr<SampleHistory>& sample_history) = 0;

	protected:
		size_t num_variables;
	};
}
