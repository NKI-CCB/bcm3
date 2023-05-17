#pragma once

#include "Proposal.h"

namespace bcm3 {

	class GMM;

	class ProposalClusteredCovariance : public Proposal
	{
	public:
		ProposalClusteredCovariance();
		virtual ~ProposalClusteredCovariance();

		virtual void GetNewSample(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, RNG& rng);
		virtual Real CalculateMHRatio(const VectorReal& current_position, ptrdiff_t curpos_cluster_assignment, const VectorReal& new_position, ptrdiff_t newpos_cluster_assignment);

		virtual bool UsesClustering();
		virtual void Update(RNG& rng);
		virtual void NotifyAccepted(bool accepted);
		virtual void LogInfo() const;
		virtual void WriteToFile(const std::string& fn, std::string adaptation_group, std::vector<ptrdiff_t>& variable_indices);

	protected:
		virtual bool InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info);

		// Runtime variables
		std::shared_ptr<GMM> gmm; // We're not actually using a full Gaussian mixture; we're just using the class as a convenient way to store multiple Gaussians
		VectorReal scales;
		VectorReal acceptance_rate_emas;
		ptrdiff_t selected_component;
	};

}
