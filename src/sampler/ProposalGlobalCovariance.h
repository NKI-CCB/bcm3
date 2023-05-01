#pragma once

#include "Proposal.h"

namespace bcm3 {

	class ProposalGlobalCovariance : public Proposal
	{
	public:
		ProposalGlobalCovariance();
		virtual ~ProposalGlobalCovariance();

		virtual void LogInfo() const;

	protected:
		virtual bool InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info);
		virtual void GetNewSampleImpl(const VectorReal& current_position, VectorReal& new_position, Real& log_mh_ratio, RNG& rng);

		// Settings
		Real t_dof;

		// Runtime variables
		MatrixReal covariance;
		Eigen::LLT<MatrixReal> covariance_llt;
		MatrixReal covariance_decomp;
		Real logC;
	};

}
