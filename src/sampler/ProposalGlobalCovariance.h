#pragma once

#include "Proposal.h"

namespace bcm3 {

	class ProposalGlobalCovariance : public Proposal
	{
	public:
		ProposalGlobalCovariance();
		virtual ~ProposalGlobalCovariance();

		virtual void GetNewSample(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, RNG& rng);
		virtual Real CalculateMHRatio(const VectorReal& current_position, ptrdiff_t curpos_cluster_assignment, const VectorReal& new_position, ptrdiff_t newpos_cluster_assignment);

		virtual void LogInfo() const;

	protected:
		virtual bool InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info);

		// Settings
		Real t_dof;

		// Runtime variables
		MatrixReal covariance;
		Eigen::LLT<MatrixReal> covariance_llt;
		MatrixReal covariance_decomp;
		Real logC;
	};

}
