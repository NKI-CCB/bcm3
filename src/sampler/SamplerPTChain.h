#pragma once

#include "GMM.h"
#include <boost/dynamic_bitset.hpp>

namespace bcm3 {

class SamplerPT;

class SamplerPTChain
{
public:
	SamplerPTChain(SamplerPT* sampler);

	bool Initialize(size_t history_size, size_t history_subsampling);

	void SetTemperature(Real temperature, bool is_fixed);
	bool AdaptProposal(size_t thread);
	void AdaptProposalAsync();
	bool AdaptProposalWait();
	bool FindStartingPosition();
	bool MutateMove(size_t thread);
	void MutateMoveAsync();
	bool MutateMoveWait();
	bool ExchangeMove(SamplerPTChain& other);
	void LogStatistics() const;
	void LogProposalInfo() const;

	inline Real GetLogPrior() const { return lprior; }
	inline Real GetLogLikelihood() const { return llh; }
	inline Real GetLogPowerPosterior() const { return lpowerposterior; }
	inline const VectorReal& GetVarValues() const { return current_var_values; }
	inline Real GetTemperature() const { return temperature; }
	inline bool GetIsFixedTemperature() const { return temperature_fixed; }

private:
	bool AsyncDoMutateMove(void* user, size_t thread);
	bool AsyncDoAdaptProposal(void* user, size_t thread);
	bool AdaptProposalGMM(size_t thread);
	bool AdaptProposalBlocked(size_t thread);
	bool AdaptProposalClusteredBlocked(size_t thread);
	ptrdiff_t GetSampleCluster(const VectorReal& x);
	void UpdateScale(Real& scale, Real current_acceptance_rate_ema, RNG& rng) const;
	Real GetTDistScaleFactor(RNG& rng) const;
	bool TestSample(Real new_lpowerposterior, Real log_mh_ratio, RNG& rng);
	void AcceptMutate(bool accept, const VectorReal& nv, Real lprior, Real llh, Real lpowerposterior, Real& current_acceptance_rate_ema);
	void AddCurrentSampleToHistory();

	enum struct ProposalType {
		Blocked,
		GMM,
		ClusteredBlocked
	};

	// Settings
	SamplerPT* sampler;
	Real temperature;
	bool temperature_fixed;
	ProposalType proposal_type;

	// Run-time variables
	VectorReal current_var_values;
	Real lprior;
	Real llh;
	Real lpowerposterior;

	size_t attempted_mutate;
	size_t attempted_exchange;
	size_t accepted_mutate;
	size_t accepted_exchange;

	Eigen::MatrixXf sample_history;
	ptrdiff_t sample_history_n;
	ptrdiff_t sample_history_n_s;
	ptrdiff_t sample_history_subsampling;

	// For blocked sampler
	struct Block
	{
		std::vector<int> variable_indices;
		MatrixReal cov;
		Eigen::LLT<MatrixReal> covariance_llt;
		MatrixReal covariance_decomp;
		Real logC;
		Real sigma;
		Real scale;
		Real current_acceptance_rate_ema;
	};
	std::vector<Block> sampler_blocks;

	// For GMM sampler
	std::shared_ptr<GMM> gmm_proposal;
	VectorReal scales;
	VectorReal acceptance_rate_ema;

	// For clustered blocked sampler
	struct Cluster
	{
		std::vector<Block> blocks;
	};
	MatrixReal clustered_blocking_scaled_samples;
	VectorReal clustered_blocking_variable_scaling;
	VectorReal clustered_blocking_sample_scale;
	MatrixReal clustered_blocking_spectral_decomposition;
	MatrixReal clustered_blocking_kmeans_centroids;
	std::vector<Cluster> clustered_blocking_blocks;
	std::vector< std::vector<ptrdiff_t> > clustered_blocking_nearest_neighbors;
	std::vector< boost::dynamic_bitset<> > clustered_blocking_nearest_neighbors_bitset;
	ptrdiff_t current_cluster_assignment;
	size_t adaptation_iter;
	VectorReal get_cluster_buffer;


	// For parallelization
	uint64 async_task;
	bool async_burnin;
	std::atomic<int> async_failure;
};

}
