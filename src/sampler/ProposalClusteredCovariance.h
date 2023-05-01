#pragma once

#include "Proposal.h"

namespace bcm3 {

	class GMM;

	class ProposalClusteredCovariance : public Proposal
	{
	public:
		ProposalClusteredCovariance();
		virtual ~ProposalClusteredCovariance();

		virtual bool UsesClustering();
		virtual void Update(RNG& rng);
		virtual void NotifyAccepted(bool accepted);
		virtual void LogInfo() const;
		virtual void WriteToFile(const std::string& fn, std::string adaptation_group);

	protected:
		virtual bool InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info);
		virtual void GetNewSampleImpl(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, Real& log_mh_ratio, RNG& rng);

		// Settings
		Real t_dof;

		// Runtime variables
		std::shared_ptr<GMM> gmm; // We're not actually using a full Gaussian mixture; we're just using the class as a convenient way to store multiple Gaussians
		VectorReal scales;
		VectorReal acceptance_rate_emas;
		ptrdiff_t selected_component;
	};

}