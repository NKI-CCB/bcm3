#pragma once

#include "Proposal.h"

namespace bcm3 {

	class ProposalParametricMixture : public Proposal
	{
	public:
		ProposalParametricMixture();
		virtual ~ProposalParametricMixture();

		virtual VectorReal GetNewSample(const VectorReal& current_position) = 0;
		virtual Real GetProposalProbability(const VectorReal& sample) = 0;

	protected:
		virtual void InitializeImpl(const std::unique_ptr<SampleHistory>& sample_history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices);
	};

}
