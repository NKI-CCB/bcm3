#include "Utils.h"
#include "BlockingStrategy.h"
#include "BlockingStrategyClusteredTurek.h"
#include "BlockingStrategyNoBlocking.h"
#include "BlockingStrategyOneBlock.h"
#include "BlockingStrategyTurek.h"
#include "GMM.h"
#include "NetCDFBundler.h"
#include "Prior.h"
#include "Proposal.h"
#include "ProposalClusteredCovariance.h"
#include "ProposalGlobalCovariance.h"
#include "ProposalParametricMixture.h"
#include "SampleHistory.h"
#include "SampleHistoryClustering.h"
#include "SamplerPT.h"
#include "SamplerPTChain.h"
#include "SummaryStats.h"
#include "VectorUtils.h"

#include <fstream>
#include <boost/filesystem.hpp>

namespace bcm3 {

	SamplerPTChain::SamplerPTChain(SamplerPT* sampler)
		: sampler(sampler)
		, temperature(std::numeric_limits<Real>::quiet_NaN())
		, temperature_fixed(true)
		, lprior(-std::numeric_limits<Real>::infinity())
		, llh(-std::numeric_limits<Real>::infinity())
		, lpowerposterior(-std::numeric_limits<Real>::infinity())
		, adaptation_iteration(0)
		, current_cluster_assignment(-1)
		, attempted_mutate(0)
		, attempted_exchange(0)
		, accepted_mutate(0)
		, accepted_exchange(0)
		, async_task(0)
		, async_burnin(false)
	{
	}

	SamplerPTChain::~SamplerPTChain()
	{
	}

	bool SamplerPTChain::Initialize(size_t history_size, size_t history_subsampling)
	{
		if (sampler->num_variables > std::numeric_limits<int>::max()) {
			LOGERROR("Too many variables");
			return false;
		}

		lprior = -std::numeric_limits<Real>::infinity();
		llh = -std::numeric_limits<Real>::infinity();
		lpowerposterior = -std::numeric_limits<Real>::infinity();
		current_var_values = VectorReal::Zero(sampler->num_variables);

		attempted_mutate = 0;
		attempted_exchange = 0;
		accepted_mutate = 0;
		accepted_exchange = 0;

		if (sampler->blocking_strategy == "one_block") {
			blocking_strategy = std::make_unique<BlockingStrategyOneBlock>();
		} else if (sampler->blocking_strategy == "no_blocking") {
			blocking_strategy = std::make_unique<BlockingStrategyNoBlocking>();
		} else if (sampler->blocking_strategy == "Turek") {
			blocking_strategy = std::make_unique<BlockingStrategyTurek>();
		} else if (sampler->blocking_strategy == "clustered_autoblock") {
			blocking_strategy = std::make_unique<BlockingStrategyClusteredTurek>();
		} else {
			LOGERROR("Unknown blocking strategy \"%s\"", sampler->blocking_strategy.c_str());
			return false;
		}
		blocking_strategy->Initialize(sampler->num_variables);
		
		sample_history = std::make_unique<SampleHistory>();
		sample_history->Initialize(sampler->num_variables, history_size, history_subsampling);

		if (blocking_strategy->UsesClustering() && temperature != 0.0) {
			sample_history_clustering = std::make_shared<SampleHistoryClustering>(sampler->adapt_proposal_max_history_samples);
			current_cluster_assignment = -1;
		}

		if (!AdaptProposal(0)) {
			return false;
		}

		if (sampler->output_proposal_adaptation && temperature == sampler->temperatures.tail(1)(0)) {
			std::string proposal_output_fn = sampler->output_path + "sampler_adaptation.nc";
			if (boost::filesystem::exists(proposal_output_fn)) {
				boost::filesystem::remove(proposal_output_fn);
			}
		}

		return true;
	}

	void SamplerPTChain::SetTemperature(Real temperature, bool is_fixed)
	{
		this->temperature = temperature;
		temperature_fixed = is_fixed;
		lpowerposterior = lprior + temperature * llh;
	}

	bool SamplerPTChain::AdaptProposal(size_t thread)
	{
		bool log_info = (temperature == sampler->temperatures.tail(1)(0)) ? true : false;

		// No adaptation necessary for a chain at temperature 0 as samples are just drawn from the prior
		if (temperature == 0.0) {
			return true;
		}

		// First cluster the history if the blocking strategy or proposal wants it
		if (sample_history_clustering && adaptation_iteration > 0) {
			size_t discard_first_samples = 0;
			if (adaptation_iteration == 1) {
				discard_first_samples = 100;
			}
			if (!sample_history_clustering->Cluster(sample_history, discard_first_samples, sampler->async[thread].rng, log_info)) {
				return false;
			}
		}

		// Update variable blocks according to the blocking strategy
		std::vector< std::vector<ptrdiff_t> > blocks = blocking_strategy->GetBlocks(sample_history, sample_history_clustering);
		if (log_info) {
			LOG("Blocking strategy returned %zu block(s)", blocks.size());
		}

		// Update the proposal for every variable block
		variable_blocks.clear();
		variable_blocks.resize(blocks.size());
		for (ptrdiff_t i = 0; i < variable_blocks.size(); i++) {
			if (log_info) {
				LOG(" Block %zd", i+1);
			}
			variable_blocks[i].variable_indices = blocks[i];
			variable_blocks[i].proposal = CreateProposalInstance(blocks[i], sampler->async[thread].rng);
		}

		if (sampler->output_proposal_adaptation && temperature == sampler->temperatures.tail(1)(0)) {
			std::string proposal_output_fn = sampler->output_path + "sampler_adaptation.nc";
			for (ptrdiff_t i = 0; i < variable_blocks.size(); i++) {
				variable_blocks[i].proposal->WriteToFile(proposal_output_fn, std::string("adapt") + std::to_string(adaptation_iteration) + std::string("_block") + std::to_string(i + 1));
			}
		}
		adaptation_iteration++;

		// Discard the sample history up to this point
		sample_history->Reset();

		return true;
	}

	void SamplerPTChain::AdaptProposalAsync()
	{
		async_failure = false;
		TTask task = boost::bind(&SamplerPTChain::AsyncDoAdaptProposal, this, boost::placeholders::_1, boost::placeholders::_2);
		async_task = sampler->task_manager->AddTask(task, NULL);
	}

	bool SamplerPTChain::AdaptProposalWait()
	{
		sampler->task_manager->WaitTask(async_task);
		return async_failure == 0;
	}

	bool SamplerPTChain::FindStartingPosition()
	{
		llh = -std::numeric_limits<Real>::infinity();
		lprior = -std::numeric_limits<Real>::infinity();
		lpowerposterior = -std::numeric_limits<Real>::infinity();

		static const size_t MAX_TRIES = 2000;
		for (size_t i = 0; i < MAX_TRIES; i++) {
			if (!sampler->prior->Sample(current_var_values, &sampler->rng)) {
				LOGERROR("Sampling from prior failed");
				return false;
			}

			if (!sampler->EvaluatePriorLikelihood(0, current_var_values, lprior, llh)) {
				return false;
			}
			lpowerposterior = lprior + temperature * llh;
			if (lpowerposterior > -std::numeric_limits<Real>::infinity()) {
				break;
			}
		}

		if (std::isnan(lpowerposterior) || lpowerposterior == -std::numeric_limits<Real>::infinity()) {
			LOGERROR("Could not find starting position with power posterior != inf after %u tries", MAX_TRIES);
			return false;
		} else {
			return true;
		}
	}

	bool SamplerPTChain::MutateMove(size_t thread)
	{
		RNG& rng = sampler->async[thread].rng;

		if (temperature == 0.0) {
			// Just sample from the prior
			if (!sampler->prior->Sample(current_var_values, &sampler->rng)) {
				LOGERROR("Sampling from prior failed");
				return false;
			}

			if (!sampler->EvaluatePriorLikelihood(thread, current_var_values, lprior, llh)) {
				return false;
			}
			if (llh == -std::numeric_limits<Real>::infinity()) {
				// 0 * -inf is undefined; but a negative-infinity log likelihood is usually a computational problem, so assume it's
				// an extremely large negative value in which case the multiplication with a 0 temperature should still give 0
				lpowerposterior = lprior;
			} else {
				lpowerposterior = lprior + temperature * llh;
			}

			attempted_mutate++;
			accepted_mutate++;
		} else {
			if (sample_history_clustering) {
				// Check that cluster assignment is up to date
				if (current_cluster_assignment == -1) {
					current_cluster_assignment = sample_history_clustering->GetSampleCluster(current_var_values);
				}
			}

			// Sample each block separately
			for (ptrdiff_t i = 0; i < variable_blocks.size(); i++) {
				Block& b = variable_blocks[i];

				// Give the proposal a chance to update based on previous acceptance rates
				b.proposal->Update(rng);

				// Retrieve a new sample from the proposal distribution for the variables in this block
				VectorReal block_current_position(b.variable_indices.size());
				for (ptrdiff_t i = 0; i < b.variable_indices.size(); i++) {
					block_current_position(i) = current_var_values(b.variable_indices[i]);
				}
				VectorReal block_new_position;
				b.proposal->GetNewSample(block_current_position, current_cluster_assignment, block_new_position, rng);

				// Construct the full new vector of variable values
				VectorReal new_var_values = current_var_values;
				for (ptrdiff_t i = 0; i < b.variable_indices.size(); i++) {
					new_var_values(b.variable_indices[i]) = block_new_position(i);
				}

				// Evaluate prior, likelihood & posterior
				Real new_lprior = -std::numeric_limits<Real>::infinity(), new_llh = -std::numeric_limits<Real>::infinity();
				if (!sampler->EvaluatePriorLikelihood(thread, new_var_values, new_lprior, new_llh)) {
					return false;
				}
				Real new_lpowerposterior = new_lprior + temperature * new_llh;

				// Calculate Metropolis-Hastings ratio
				ptrdiff_t new_cluster_assignment = -1;
				if (sample_history_clustering) {
					new_cluster_assignment = sample_history_clustering->GetSampleCluster(new_var_values);
				}
				Real log_mh_ratio = b.proposal->CalculateMHRatio(block_current_position, current_cluster_assignment, block_new_position, new_cluster_assignment);

				bool accept = TestSample(new_lpowerposterior, log_mh_ratio, rng);

				if (accept) {
					accepted_mutate++;

					current_var_values = new_var_values;
					current_cluster_assignment = new_cluster_assignment;
					lprior = new_lprior;
					llh = new_llh;
					lpowerposterior = new_lpowerposterior;
				}

				b.proposal->NotifyAccepted(accept);
			}

			sample_history->AddSample(current_var_values);
		}

		return true;
	}

	void SamplerPTChain::MutateMoveAsync()
	{
		async_failure = false;
		TTask task = boost::bind(&SamplerPTChain::AsyncDoMutateMove, this, boost::placeholders::_1, boost::placeholders::_2);
		async_task = sampler->task_manager->AddTask(task, NULL);
	}

	bool SamplerPTChain::MutateMoveWait()
	{
		sampler->task_manager->WaitTask(async_task);
		return async_failure == 0;
	}

	bool SamplerPTChain::ExchangeMove(SamplerPTChain& other)
	{
		attempted_exchange++;

		SamplerPTChain& chain1 = *this;
		SamplerPTChain& chain2 = other;
	
		Real proposed_lpowerposterior1, proposed_lpowerposterior2;
		if (chain1.temperature == 0.0) {
			// This is to handle the case where chain2.llh == -std::numeric_limits<Real>::infinity();
			proposed_lpowerposterior1 = chain2.lprior;
		} else {
			proposed_lpowerposterior1 = chain1.temperature * chain2.llh + chain2.lprior;
		}
		if (chain2.temperature == 0.0) {
			// This is to handle the case where chain1.llh == -std::numeric_limits<Real>::infinity();
			proposed_lpowerposterior2 = chain1.lprior;
		} else {
			proposed_lpowerposterior2 = chain2.temperature * chain1.llh + chain1.lprior;
		}

		bool swap = false;
		Real transition_probability = (proposed_lpowerposterior1 + proposed_lpowerposterior2) - (chain1.lpowerposterior + chain2.lpowerposterior);
		transition_probability = exp(transition_probability);
		transition_probability = std::min((Real)1.0, transition_probability);

		Real alpha = sampler->rng.GetReal();
		if (alpha < transition_probability) {
			swap = true;
		}

		if (swap) {
			if (chain1.llh == -std::numeric_limits<Real>::infinity() || chain2.llh == -std::numeric_limits<Real>::infinity()) {
				printf("Swapped %g/%g between %g and %g!\n", chain1.llh, chain2.llh, chain1.temperature, chain2.temperature);
			}

			accepted_exchange++;
			std::swap(chain1.current_var_values, chain2.current_var_values);
			std::swap(chain1.llh, chain2.llh);
			std::swap(chain1.lprior, chain2.lprior);
			chain1.current_cluster_assignment = -1;
			chain2.current_cluster_assignment = -1;
			chain1.lpowerposterior = proposed_lpowerposterior1;
			chain2.lpowerposterior = proposed_lpowerposterior2;
		}

		sample_history->AddSample(current_var_values);
		other.sample_history->AddSample(current_var_values);
		return swap;
	}

	void SamplerPTChain::LogStatistics() const
	{
		LOG("%11.7f | %12.5f | %14.5f",
			temperature,
			accepted_mutate / (Real)attempted_mutate,
			accepted_exchange / (Real)attempted_exchange);
	}

	void SamplerPTChain::LogProposalInfo() const
	{
		LOG("Proposal information:");
		for (size_t i = 0; i < variable_blocks.size(); i++) {
			const Block& b = variable_blocks[i];

			std::string ixs = std::to_string(b.variable_indices[0]);
			for (size_t vi = 1; vi < b.variable_indices.size(); vi++) {
				ixs += std::string(",") + std::to_string(b.variable_indices[vi]);
			}

			LOG(" Block %u: %s", i + 1, ixs.c_str());

			b.proposal->LogInfo();
		}
	}

	bool SamplerPTChain::AsyncDoMutateMove(void* user, size_t thread)
	{
		bool result = MutateMove(thread);
		if (!result) {
			async_failure = 1;
		}
		return result;
	}

	bool SamplerPTChain::AsyncDoAdaptProposal(void* user, size_t thread)
	{
		bool result = AdaptProposal(thread);
		if (!result) {
			async_failure = 1;
		}
		return result;
	}

	std::unique_ptr<Proposal> SamplerPTChain::CreateProposalInstance(std::vector<ptrdiff_t> variable_indices, RNG& rng)
	{
		std::unique_ptr<Proposal> proposal;
		if (sampler->proposal_type == "global_covariance") {
			proposal = std::make_unique<ProposalGlobalCovariance>();
		} else if (sampler->proposal_type == "parametric_mixture") {
			proposal = std::make_unique<ProposalParametricMixture>();
		} else if (sampler->proposal_type == "clustered_covariance") {
			proposal = std::make_unique<ProposalClusteredCovariance>();
		} else {
			LOGERROR("Unknown proposal type \"%s\"", sampler->proposal_type.c_str());
			return proposal;
		}

		bool log_info = (temperature == sampler->temperatures.tail(1)(0)) ? true : false;
		if (!proposal->Initialize(*sample_history, sample_history_clustering, sampler->adapt_proposal_max_history_samples, sampler->proposal_transform_to_unbounded, sampler->prior, variable_indices, rng, log_info)) {
			LOGERROR("  Proposal initialization failed.");
			proposal.reset();
		}

		return proposal;
	}

	bool SamplerPTChain::TestSample(Real new_lpowerposterior, Real log_mh_ratio, RNG& rng)
	{
		attempted_mutate++;
		if (new_lpowerposterior > -std::numeric_limits<Real>::infinity()) {
			Real transition_probability = new_lpowerposterior - lpowerposterior;
			transition_probability = exp(transition_probability + log_mh_ratio);
			transition_probability = std::min((Real)1.0, transition_probability);

			Real alpha = rng.GetReal();
			if (alpha < transition_probability) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

}
