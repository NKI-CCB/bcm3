#pragma once

#include "ProposalGaussianMixture.h"

namespace bcm3 {

	class GMM;

	class ProposalGaussianMixtureFitInR : public ProposalGaussianMixture
	{
	public:
		ProposalGaussianMixtureFitInR();
		virtual ~ProposalGaussianMixtureFitInR();

	protected:
		virtual bool InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info);

		static std::mutex netcdf_mutex;
	};

}
