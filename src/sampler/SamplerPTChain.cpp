#include "Utils.h"
#include "BlockingStrategy.h"
#include "BlockingStrategyNoBlocking.h"
#include "BlockingStrategyOneBlock.h"
#include "BlockingStrategyTurek.h"
#include "GMM.h"
#include "NetCDFBundler.h"
#include "Prior.h"
#include "Proposal.h"
#include "ProposalGlobalCovariance.h"
#include "SampleHistory.h"
#include "SamplerPT.h"
#include "SamplerPTChain.h"
#include "SummaryStats.h"
#include "VectorUtils.h"

#include <fstream>

namespace bcm3 {

SamplerPTChain::SamplerPTChain(SamplerPT* sampler)
	: sampler(sampler)
	, temperature(std::numeric_limits<Real>::quiet_NaN())
	, lprior(-std::numeric_limits<Real>::infinity())
	, llh(-std::numeric_limits<Real>::infinity())
	, lpowerposterior(-std::numeric_limits<Real>::infinity())
	, attempted_mutate(0)
	, attempted_exchange(0)
	, accepted_mutate(0)
	, accepted_exchange(0)
	, async_task(0)
	, async_burnin(false)
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
	} else {
		LOGERROR("Unknown blocking strategy \"%s\"", sampler->blocking_strategy.c_str());
		return false;
	}
	blocking_strategy->Initialize(sampler->num_variables);

#if 0
	{
		proposal_type = ProposalType::Blocked;
	} else if (sampler->proposal_type == "GMM") {
		proposal_type = ProposalType::GMM;
	} else if (sampler->proposal_type == "clustered") {
		proposal_type = ProposalType::ClusteredBlocked;
	} else {
		LOGERROR("Unknown proposal type \"%s\"", sampler->proposal_type.c_str());
		return false;
	}
#endif

	sample_history = std::make_unique<SampleHistory>();
	sample_history->Initialize(sampler->num_variables, history_size, history_subsampling);

	if (!AdaptProposal(0)) {
		return false;
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
	// No adaptation necessary for a chain at temperature 0 as samples are just drawn from the prior
	if (temperature == 0.0) {
		return true;
	}

	// First update blocks according to the blocking strategy
	std::vector< std::vector<ptrdiff_t> > blocks = blocking_strategy->GetBlocks(sample_history);

	// Update the proposal for every variable block
	variable_blocks.clear();
	variable_blocks.resize(blocks.size());
	for (ptrdiff_t i = 0; i < variable_blocks.size(); i++) {
		variable_blocks[i].variable_indices = blocks[i];
		variable_blocks[i].proposal = CreateProposalInstance(blocks[i]);
	}

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
		// Sample each block separately
		for (ptrdiff_t i = 0; i < variable_blocks.size(); i++) {
			Block& b = variable_blocks[i];

			// Give the proposal a chance to update based on previous acceptance rates
			b.proposal->Update(rng);

			// Retrieve a new sample from the proposal distribution
			VectorReal current_position(b.variable_indices.size());
			for (ptrdiff_t i = 0; i < b.variable_indices.size(); i++) {
				current_position(i) = current_var_values(b.variable_indices[i]);
			}
			VectorReal new_position;
			Real log_mh_ratio;
			b.proposal->GetNewSample(current_position, new_position, log_mh_ratio, rng);

			// Evaluate prior, likelihood & posterior
			Real new_lprior = -std::numeric_limits<Real>::infinity(), new_llh = -std::numeric_limits<Real>::infinity();
			if (!sampler->EvaluatePriorLikelihood(thread, new_position, new_lprior, new_llh)) {
				return false;
			}
			Real new_lpowerposterior = new_lprior + temperature * new_llh;

			bool accept = TestSample(new_lpowerposterior, log_mh_ratio, rng);
			AcceptMutate(accept, new_position, new_lprior, new_llh, new_lpowerposterior);
			b.proposal->NotifyAccepted(accept);
		}
		AddCurrentSampleToHistory();
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
		chain1.lpowerposterior = proposed_lpowerposterior1;
		chain2.lpowerposterior = proposed_lpowerposterior2;
	}

	AddCurrentSampleToHistory();
	other.AddCurrentSampleToHistory();
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
#if 0
	switch (proposal_type) {
	case ProposalType::Blocked:
		LOG(" Blocked proposal:");
		for (size_t i = 0; i < sampler_blocks.size(); i++) {
			const Block& b = sampler_blocks[i];

			std::string ixs = std::to_string(b.variable_indices[0]);
			for (size_t vi = 1; vi < b.variable_indices.size(); vi++) {
				ixs += std::string(",") + std::to_string(b.variable_indices[vi]);
			}

			LOG("  block %u: %8.5f - %s", i + 1, b.scale, ixs.c_str());
		}
		break;

	case ProposalType::GMM:
		LOG(" GMM proposal with %u components", gmm_proposal->GetNumComponents());
		for (size_t i = 0; i < gmm_proposal->GetNumComponents(); i++) {
			LOG("  Component %u; scale=%8.5f, weight=%8.5f, condition number=%6g", i + 1, scales(i), gmm_proposal->GetWeights()(i), 1.0 / gmm_proposal->GetCovarianceDecomp(i).rcond());
		}
		break;

	case ProposalType::ClusteredBlocked:
		LOG(" Clustered blocked proposal:");
		LOG("  Variable blocking:");
		for (size_t i = 0; i < clustered_blocking_blocks[0].blocks.size(); i++) {
			const Block& b = clustered_blocking_blocks[0].blocks[i];

			std::string ixs = std::to_string(b.variable_indices[0]);
			for (size_t vi = 1; vi < b.variable_indices.size(); vi++) {
				ixs += std::string(",") + std::to_string(b.variable_indices[vi]);
			}

			LOG("   block %u: %s", i + 1, ixs.c_str());
		}
		LOG("  Clusters:");
		for (size_t j = 0; j < clustered_blocking_blocks.size(); j++) {
			LOG("   cluster %u", j + 1);
			for (size_t i = 0; i < clustered_blocking_blocks[j].blocks.size(); i++) {
				const Block& b = clustered_blocking_blocks[j].blocks[i];
				LOG("    block %u: scale: %8.5f \tacceptance: %8.5f", i + 1, b.scale, b.current_acceptance_rate_ema);
			}
		}
		break;

	default:
		LOGERROR("Invalid proposal type");
		break;
	}
#endif
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

std::unique_ptr<Proposal> SamplerPTChain::CreateProposalInstance(std::vector<ptrdiff_t> variable_indices)
{
	std::unique_ptr<Proposal> proposal;
	if (sampler->proposal_type == "global_covariance") {
		proposal = std::make_unique<ProposalGlobalCovariance>();
	} else {
		LOGERROR("Unknown proposal type \"%s\"", sampler->proposal_type.c_str());
		return proposal;
	}

	proposal->Initialize(sample_history, sampler->prior, variable_indices);
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

void SamplerPTChain::AcceptMutate(bool accept, const VectorReal& nv, Real lprior, Real llh, Real lpowerposterior)
{
	if (accept) {
		accepted_mutate++;

		current_var_values = nv;
		this->lprior = lprior;
		this->llh = llh;
		this->lpowerposterior = lpowerposterior;
	}
}

void SamplerPTChain::AddCurrentSampleToHistory()
{
	sample_history->AddSample(current_var_values);
}

}
