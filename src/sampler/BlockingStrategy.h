#pragma once

namespace bcm3 {

	class SampleHistory;
	class SampleHistoryClustering;

	class BlockingStrategy
	{
	public:
		virtual void Initialize(size_t num_variables);
		virtual bool UsesClustering();

		virtual std::vector< std::vector<ptrdiff_t> > GetBlocks(const std::unique_ptr<SampleHistory>& sample_history, const std::shared_ptr<SampleHistoryClustering> clustering) = 0;

	protected:
		size_t num_variables;
	};
}
