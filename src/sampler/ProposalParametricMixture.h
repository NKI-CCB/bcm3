#pragma once

#include "Proposal.h"

namespace bcm3 {

	class GMM;

	class ProposalParametricMixture : public Proposal
	{
	public:
		ProposalParametricMixture();
		virtual ~ProposalParametricMixture();

		virtual void Update(RNG& rng);
		virtual void NotifyAccepted(bool accepted);
		virtual void LogInfo() const;

	protected:
		virtual bool InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info);
		virtual void GetNewSampleImpl(const VectorReal& current_position, VectorReal& new_position, Real& log_mh_ratio, RNG& rng);

		// Settings
		Real t_dof;

		// Runtime variables
		std::shared_ptr<GMM> gmm;
		VectorReal scales;
		VectorReal acceptance_rate_emas;
		ptrdiff_t selected_component;
	};

}
