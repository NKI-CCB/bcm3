#include "Utils.h"
#include "ProposalParametricMixture.h"

namespace bcm3 {

	ProposalParametricMixture::ProposalParametricMixture()
	{
	}

	ProposalParametricMixture::~ProposalParametricMixture()
	{
	}

	VectorReal ProposalParametricMixture::GetNewSample(const VectorReal& current_position)
	{
		VectorReal transformed = Transform(current_position);
	}

	Real ProposalParametricMixture::GetProposalProbability(const VectorReal& sample)
	{
		VectorReal transformed = Transform(sample);
	}

	void ProposalParametricMixture::InitializeImpl(const std::unique_ptr<SampleHistory>& sample_history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices)
	{

	}
}
