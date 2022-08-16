#include "Utils.h"
#include "Clustering.h"
#include "GMM.h"
#include "NetCDFBundler.h"
#include "Prior.h"
#include "SamplerPT.h"
#include "SamplerPTChain.h"
#include "SummaryStats.h"
#include "VectorUtils.h"

#include <fstream>

namespace bcm3 {

inline void reflect_bounds(Real& x, Real lb, Real ub)
{
	while (1) {
		if (x < lb) {
			x = lb + (lb - x);
		} else if (x > ub) {
			x = ub - (x - ub);
		} else {
			break;
		}
	}
}

SamplerPTChain::SamplerPTChain(SamplerPT* sampler)
	: sampler(sampler)
	, temperature(std::numeric_limits<Real>::quiet_NaN())
	, proposal_type(ProposalType::ClusteredBlocked)
	, lprior(-std::numeric_limits<Real>::infinity())
	, llh(-std::numeric_limits<Real>::infinity())
	, lpowerposterior(-std::numeric_limits<Real>::infinity())
	, attempted_mutate(0)
	, attempted_exchange(0)
	, accepted_mutate(0)
	, accepted_exchange(0)
	, sample_history_n(0)
	, sample_history_n_s(0)
	, sample_history_subsampling(1)
	, current_cluster_assignment(0)
	, adaptation_iter(0)
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

	sample_history.resize(sampler->num_variables, history_size);
	sample_history_n = 0;
	sample_history_n_s = 0;
	sample_history_subsampling = history_subsampling;

	if (sampler->proposal_type == "autoblock") {
		proposal_type = ProposalType::Blocked;
	} else if (sampler->proposal_type == "GMM") {
		proposal_type = ProposalType::GMM;
	} else if (sampler->proposal_type == "clustered") {
		proposal_type = ProposalType::ClusteredBlocked;
	} else {
		LOGERROR("Unknown proposal type \"%s\"", sampler->proposal_type.c_str());
		return false;
	}

	switch (proposal_type) {
		case ProposalType::Blocked:
		{
			// Start with scalar sampler
			sampler_blocks.resize(sampler->num_variables);
			for (size_t i = 0; i < sampler->num_variables; i++) {
				Real prior_var;
				if (!sampler->prior->EvaluateMarginalVariance(i, prior_var)) {
					// TODO - somehow come up with something reasonable?
					prior_var = 1.0;
				}
				sampler_blocks[i].variable_indices.resize(1, (int)i);
				sampler_blocks[i].sigma = sqrt(prior_var);
				sampler_blocks[i].scale = 1.0;
				sampler_blocks[i].current_acceptance_rate_ema = sampler->target_acceptance_rate;
			}
			break;
		}

		case ProposalType::ClusteredBlocked:
		{
			// Start with scalar sampler
			clustered_blocking_blocks.resize(1);
			clustered_blocking_blocks[0].blocks.resize(sampler->num_variables);
			for (size_t i = 0; i < sampler->num_variables; i++) {
				Real prior_var;
				if (!sampler->prior->EvaluateMarginalVariance(i, prior_var)) {
					// TODO - somehow come up with something reasonable?
					prior_var = 1.0;
				}
				clustered_blocking_blocks[0].blocks[i].variable_indices.resize(1, (int)i);
				clustered_blocking_blocks[0].blocks[i].sigma = sqrt(prior_var);
				clustered_blocking_blocks[0].blocks[i].scale = 1.0;
				clustered_blocking_blocks[0].blocks[i].current_acceptance_rate_ema = sampler->target_acceptance_rate;
			}
			current_cluster_assignment = 0;
			adaptation_iter = 0;

			clustered_blocking_variable_scaling.setOnes(sampler->num_variables);
			break;
		}

		case ProposalType::GMM:
		{
			// Start with vector sampler
			sampler_blocks.resize(1);
			sampler_blocks[0].variable_indices.resize(sampler->num_variables);
			sampler_blocks[0].covariance_decomp.setConstant(sampler->num_variables, sampler->num_variables, 0);
			MatrixReal cov(sampler->num_variables, sampler->num_variables);
			for (size_t i = 0; i < sampler->num_variables; i++) {
				Real prior_var;
				if (!sampler->prior->EvaluateMarginalVariance(i, prior_var)) {
					// TODO - somehow come up with something reasonable?
					prior_var = 1.0;
				}
				sampler_blocks[0].variable_indices[i] = (int)i;
				sampler_blocks[0].covariance_decomp(i, i) = sqrt(prior_var);
			}
			sampler_blocks[0].scale = 1.0 / (Real)sampler->num_variables;
			sampler_blocks[0].current_acceptance_rate_ema = sampler->target_acceptance_rate;
			break;
		}

		default:
			LOGERROR("Invalid proposal type");
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
	// No adaptation necessary for a chain at 0 temperature as samples are just drawn from the prior
	if (temperature == 0.0) {
		return true;
	}

	switch (proposal_type) {
		case ProposalType::Blocked:
			if (!AdaptProposalBlocked(thread)) {
				return false;
			}
			break;

		case ProposalType::GMM:
			if (!AdaptProposalGMM(thread)) {
				return false;
			}
			break;

		case ProposalType::ClusteredBlocked:
			if (!AdaptProposalClusteredBlocked(thread)) {
				return false;
			}
			break;

		default:
			LOGERROR("Invalid proposal type");
			return false;
	}

	// Reset sample history
	sample_history_n = 0;
	sample_history_n_s = 0;

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
		if (proposal_type == ProposalType::ClusteredBlocked) {
			current_cluster_assignment = 0;
		}

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
			// 0 * -inf is undefined; but the negative-infinity is really a computational problem and should rather be some extremely large negative value
			// so the multiplication with a 0 temperature should still give 0
			lpowerposterior = lprior;
		} else {
			lpowerposterior = lprior + temperature * llh;
		}

		attempted_mutate++;
		accepted_mutate++;
	} else {
		switch (proposal_type) {
		case ProposalType::Blocked:
		{	
			// Sample each block separately
			for (size_t i = 0; i < sampler_blocks.size(); i++) {
				Block& b = sampler_blocks[i];

				UpdateScale(b.scale, b.current_acceptance_rate_ema, rng);
				Real use_scale = b.scale * GetTDistScaleFactor(rng);

				// Propose new parameters
				VectorReal new_values = current_var_values;
				if (b.variable_indices.size() == 1) {
					// Scalar
					int ix = b.variable_indices[0];
					Real x = rng.GetNormal(current_var_values(ix), use_scale * b.sigma);
					reflect_bounds(x, sampler->prior->GetLowerBound(ix), sampler->prior->GetUpperBound(ix));
					new_values(ix) = x;
				} else {
					// Multivariate normal random value
					VectorReal x = rng.GetMultivariateUnitNormal(b.variable_indices.size());
					x = b.covariance_decomp * x;
					for (size_t vi = 0; vi < b.variable_indices.size(); vi++) {
						int ix = b.variable_indices[vi];
						x(vi) *= use_scale;
						x(vi) += current_var_values(ix);
						reflect_bounds(x(vi), sampler->prior->GetLowerBound(ix), sampler->prior->GetUpperBound(ix));
						new_values(ix) = x(vi);
					}
				}

				// Evaluate prior, likelihood & posterior
				Real new_lprior = -std::numeric_limits<Real>::infinity(), new_llh = -std::numeric_limits<Real>::infinity();
				if (!sampler->EvaluatePriorLikelihood(thread, new_values, new_lprior, new_llh)) {
					return false;
				}
				Real new_lpowerposterior = new_lprior + temperature * new_llh;

				bool accept = TestSample(new_lpowerposterior, 0.0, rng);
				AcceptMutate(accept, new_values, new_lprior, new_llh, new_lpowerposterior, b.current_acceptance_rate_ema);
			}
			AddCurrentSampleToHistory();
			break;
		}

		case ProposalType::GMM:
		{
			VectorReal fwd_resp = gmm_proposal->CalculateResponsibilities(current_var_values);
			unsigned int component = rng.Sample(fwd_resp);

			UpdateScale(scales(component), acceptance_rate_ema(component), rng);
			Real use_scale = scales(component) * GetTDistScaleFactor(rng);

			// Propose new parameters
			VectorReal x = rng.GetMultivariateUnitNormal(current_var_values.size());
			x *= use_scale;
			x = gmm_proposal->GetCovarianceDecomp(component).matrixL() * x;
			x += current_var_values;
			for (size_t i = 0; i < current_var_values.size(); i++) {
				reflect_bounds(x(i), sampler->prior->GetLowerBound(i), sampler->prior->GetUpperBound(i));
			}

			// Evaluate prior, likelihood & posterior
			Real new_lprior = -std::numeric_limits<Real>::infinity(), new_llh = -std::numeric_limits<Real>::infinity();
			if (!sampler->EvaluatePriorLikelihood(thread, x, new_lprior, new_llh)) {
				return false;
			}
			Real new_lpowerposterior = new_lprior + temperature * new_llh;

			// Calculate Metropolis-Hastings ratio
			VectorReal rev_resp = gmm_proposal->CalculateResponsibilities(x);
			Real fwd_logp = 0;
			Real rev_logp = 0;
			VectorReal v = (x - current_var_values) / use_scale;
			for (size_t i = 0; i < gmm_proposal->GetNumComponents(); i++) {
				VectorReal s = gmm_proposal->GetCovarianceDecomp(i).matrixL().solve(v);
				Real p = gmm_proposal->GetLogC(i) - 0.5 * s.dot(s) + log(fwd_resp(i));
				fwd_logp = bcm3::logsum(fwd_logp, p);

				s = gmm_proposal->GetCovarianceDecomp(i).matrixL().solve(-v);
				p = gmm_proposal->GetLogC(i) - 0.5 * s.dot(s) + log(rev_resp(i));
				rev_logp = bcm3::logsum(rev_logp, p);
			}

			bool accept = TestSample(new_lpowerposterior, rev_logp - fwd_logp, rng);
			AcceptMutate(accept, x, new_lprior, new_llh, new_lpowerposterior, acceptance_rate_ema(component));
			AddCurrentSampleToHistory();
			break;
		}

		case ProposalType::ClusteredBlocked:
		{
			// Sample each block separately
			Cluster& c = clustered_blocking_blocks[current_cluster_assignment];
			for (size_t i = 0; i < c.blocks.size(); i++) {
				Block& b = c.blocks[i];

				UpdateScale(b.scale, b.current_acceptance_rate_ema, rng);
				Real use_scale = b.scale;// *GetTDistScaleFactor(rng);

				// Propose new parameters
				VectorReal new_values = current_var_values;
				if (b.variable_indices.size() == 1) {
					// Scalar
					int ix = b.variable_indices[0];
					Real x = rng.GetNormal(current_var_values(ix), use_scale * b.sigma);
					reflect_bounds(x, sampler->prior->GetLowerBound(ix), sampler->prior->GetUpperBound(ix));
					new_values(ix) = x;
				} else {
					// Multivariate normal random value
					VectorReal x = rng.GetMultivariateUnitNormal(b.variable_indices.size());
					x = b.covariance_decomp * x;
					for (size_t vi = 0; vi < b.variable_indices.size(); vi++) {
						int ix = b.variable_indices[vi];
						x(vi) *= use_scale;
						x(vi) += current_var_values(ix);
						reflect_bounds(x(vi), sampler->prior->GetLowerBound(ix), sampler->prior->GetUpperBound(ix));
						new_values(ix) = x(vi);
					}
				}

				// Evaluate prior, likelihood & posterior
				Real new_lprior = -std::numeric_limits<Real>::infinity(), new_llh = -std::numeric_limits<Real>::infinity();
				if (!sampler->EvaluatePriorLikelihood(thread, new_values, new_lprior, new_llh)) {
					return false;
				}
				Real new_lpowerposterior = new_lprior + temperature * new_llh;

				// Find cluster of proposed point
				ptrdiff_t new_cluster_assignment = GetSampleCluster(new_values);
				Real log_mh_ratio;
				if (new_cluster_assignment != current_cluster_assignment) {
					Block& newb = clustered_blocking_blocks[new_cluster_assignment].blocks[i];

					Real log_fwd_p, log_bwd_p;
					if (b.variable_indices.size() == 1) {
						int ix = b.variable_indices[0];
						Real v = new_values(ix) - current_var_values(ix);

						log_fwd_p = b.logC - 0.5 * v * v / square(b.scale * b.sigma);
						log_bwd_p = newb.logC - 0.5 * v * v / square(newb.scale * newb.sigma);
					} else {
						VectorReal v(b.variable_indices.size());
						for (size_t vi = 0; vi < b.variable_indices.size(); vi++) {
							v(vi) = new_values(b.variable_indices[vi]) - current_var_values(b.variable_indices[vi]);
						}
						v /= b.scale;
						b.covariance_llt.matrixL().solveInPlace(v);
						log_fwd_p = b.logC - 0.5 * v.dot(v);

						for (size_t vi = 0; vi < b.variable_indices.size(); vi++) {
							v(vi) = current_var_values(b.variable_indices[vi]) - new_values(b.variable_indices[vi]);
						}
						v /= newb.scale;
						newb.covariance_llt.matrixL().solveInPlace(v);
						log_bwd_p = newb.logC - 0.5 * v.dot(v);
					}

					log_mh_ratio = log_bwd_p - log_fwd_p;
				} else {
					// Proposal is symmetrical within cluster
					log_mh_ratio = 0.0;
				}

				bool accept = TestSample(new_lpowerposterior, log_mh_ratio, rng);
				AcceptMutate(accept, new_values, new_lprior, new_llh, new_lpowerposterior, b.current_acceptance_rate_ema);

				if (accept) {
					current_cluster_assignment = new_cluster_assignment;
				}
			}
			AddCurrentSampleToHistory();
			break;
		}
		break;

		default:
			LOGERROR("Invalid proposal type");
			return false;
		}
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

bool SamplerPTChain::AdaptProposalGMM(size_t thread)
{
	RNG& rng = sampler->async[thread].rng;

	ptrdiff_t n_history_samples = std::min((ptrdiff_t)sample_history.cols(), sample_history_n);
	MatrixReal sample_history_d = sample_history.block(0, 0, sampler->num_variables, n_history_samples).cast<Real>();

	// Calculate effective sample size
	VectorReal ess(sampler->num_variables);
	for (size_t i = 0; i < sampler->num_variables; i++) {
		Real rho_t = 0;

		Real mu = mean((VectorReal)sample_history_d.row(i));
		Real sigmaSq = var((VectorReal)sample_history_d.row(i));

		int lag_max = std::max(5, (int)(10 * log10((Real)n_history_samples)));
		for (int lag = 1; lag < lag_max; lag++) {
			Real c = acf((VectorReal)sample_history_d.row(i), (int)lag, mu, sigmaSq);
			rho_t += c;
		}

		ess(i) = (Real)n_history_samples / (1.0 + 2.0 * rho_t);
	}
	Real min_ess = ess.minCoeff();
	if (temperature == 1.0) {
		LOG("Fitting GMMs...");
		LOG("Number of samples in history: %u", n_history_samples);
		LOG("Minimum effective sample size: %g", min_ess);
	}

	// Fit GMMs for increasing number of components and select the one with the lowest AIC
	Real best_aic = std::numeric_limits<Real>::infinity();
	gmm_proposal.reset();
	static const size_t num_components[7] = { 1, 2, 3, 5, 8, 13, 21 };
	for (size_t i = 0; i < 7; i++) {
		std::shared_ptr<GMM> gmm = std::make_shared<GMM>();
		if (min_ess < num_components[i] * sampler->num_variables * 2) {
			if (temperature == 1.0) {
				LOG("GMM num_components=%2d - not enough effective samples", num_components[i], gmm->GetAIC());
			}
		} else if (gmm->Fit(sample_history, n_history_samples, num_components[i], rng)) {
			if (temperature == 1.0) {
				LOG("GMM num_components=%2d - AIC=%.6g", num_components[i], gmm->GetAIC());

#if 0
				for (size_t j = 0; j < gmm->GetNumComponents(); j++) {
					Eigen::SelfAdjointEigenSolver<MatrixReal> eig;
					eig.compute(gmm->GetCovariance(j));
					std::stringstream sstream;
					sstream << eig.eigenvalues().transpose();
					LOG(" Component %u; covariance eigenvalues: %s", j, sstream.str().c_str());
				}
#endif
			}

			if (gmm->GetAIC() < best_aic) {
				gmm_proposal = gmm;
				best_aic = gmm->GetAIC();
			}
		} else {
			if (temperature == 1.0) {
				LOG("GMM num_components=%2d failed", i + 1);
			}
		}
	}

	if (gmm_proposal) {
		// Start sampling at a scale of 1
		scales.setOnes(gmm_proposal->GetNumComponents());
		acceptance_rate_ema.setConstant(gmm_proposal->GetNumComponents(), sampler->target_acceptance_rate);
	}

	return true;
}

bool SamplerPTChain::AdaptProposalBlocked(size_t thread)
{
	ptrdiff_t n_history_samples = std::min((ptrdiff_t)sample_history.cols(), sample_history_n);
	MatrixReal sample_history_d = sample_history.block(0, 0, sampler->num_variables, n_history_samples).cast<Real>();

	MatrixReal empirical_covariance = cov(sample_history_d.transpose());
	MatrixReal empirical_correlation = cor(sample_history_d.transpose());

	// Calculate distance for clustering
	MatrixReal distance = MatrixReal::Zero(empirical_correlation.cols(), empirical_correlation.cols());
	for (ptrdiff_t i = 0; i < empirical_correlation.rows(); i++) {
		for (ptrdiff_t j = 0; j < empirical_correlation.cols(); j++) {
			distance(i, j) = 1.0 - fabs(empirical_correlation(i, j));
		}
	}

	std::vector< std::set<size_t> > clusters;
	bcm3::TreeCluster(distance, clusters, 0.5);
	sampler_blocks.clear();
	sampler_blocks.resize(clusters.size());
	for (size_t i = 0; i < clusters.size(); i++) {
		for (std::set<size_t>::iterator vi = clusters[i].begin(); vi != clusters[i].end(); ++vi) {
			sampler_blocks[i].variable_indices.push_back((int)*vi);
		}

		if (clusters[i].size() > 1) {
			MatrixReal cov = MatrixReal::Zero(clusters[i].size(), clusters[i].size());
			for (size_t x = 0; x < clusters[i].size(); x++) {
				for (size_t y = 0; y < clusters[i].size(); y++) {
					cov(x, y) = empirical_covariance(sampler_blocks[i].variable_indices[x], sampler_blocks[i].variable_indices[y]);
				}
			}
			sampler_blocks[i].covariance_llt = cov.llt();
			sampler_blocks[i].covariance_decomp = sampler_blocks[i].covariance_llt.matrixL();
			Real det = 0.0;
			for (size_t j = 0; j < clusters[i].size(); j++) {
				det += log(sampler_blocks[i].covariance_decomp(j, j));
			}
			sampler_blocks[i].logC = -det - 0.5 * clusters[i].size() * log(2.0 * M_PI);
		} else {
			int ix = sampler_blocks[i].variable_indices[0];

			// Make sure the proposal variance is not too small
			Real var = empirical_covariance(ix,ix);
			Real prior_var;
			if (!sampler->prior->EvaluateMarginalVariance(ix, prior_var)) {
				// TODO - somehow come up with something reasonable?
				prior_var = 1.0;
			}
			var = std::max(var, (Real)1e-6 * prior_var);
			sampler_blocks[i].sigma = sqrt(var);
		}

		sampler_blocks[i].scale = 1.0;
		sampler_blocks[i].current_acceptance_rate_ema = sampler->target_acceptance_rate;
	}

	return true;
}

bool SamplerPTChain::AdaptProposalClusteredBlocked(size_t thread)
{
	ptrdiff_t n = std::min((ptrdiff_t)sample_history.cols(), sample_history_n);
	if (n < 1) {
		LOGERROR("Insufficient history samples");
		return false;
	}
	LOG("Clustered blocking adaptation - samples in history: %d", n);

	NetCDFBundler update_info_output;
	bool output_update_info = false;
	std::string output_update_info_group = std::string("info") + std::to_string(adaptation_iter);
	if (sampler->output_proposal_adaptation && temperature == 1.0) {
		if (update_info_output.Open(sampler->output_path + "sampler_adaptation.nc")) {
			update_info_output.AddGroup(output_update_info_group);
			output_update_info = true;
		}
	}

	// Calculate standard deviation over entire history
	MatrixReal sample_history_d = sample_history.block(0, 0, sampler->num_variables, n).cast<Real>();
	clustered_blocking_variable_scaling = rowSd(sample_history_d);

	for (ptrdiff_t i = 0; i < clustered_blocking_variable_scaling.size(); i++) {
		if (clustered_blocking_variable_scaling(i) <= 0.0 || std::isinf(clustered_blocking_variable_scaling(i)) || std::isnan(clustered_blocking_variable_scaling(i))) {
			LOGERROR("0, infinite or nan sd in sample history for variable %d", i);
			clustered_blocking_variable_scaling.setOnes(sampler->num_variables);
			return false;
		}
	}

	// Retrieve all unique samples
	std::vector<ptrdiff_t> use_sample_ix;
	for (ptrdiff_t i = 1; i < n; i++) {
		bool duplicate = false;
		for (ptrdiff_t j = 0; j < i; j++) {
			Eigen::VectorXf v = sample_history.col(i) - sample_history.col(j);
			float d = v.dot(v);
			if (d < std::numeric_limits<float>::epsilon()) {
				duplicate = true;
			}
		}
		if (!duplicate) {
			if (adaptation_iter == 0 && i < 100) {
				// Throw away very first samples as burnin, even for the adaptation
			} else {
				use_sample_ix.push_back(i);
			}
		}
	}
	if (use_sample_ix.size() < (ptrdiff_t)sampler->clustered_blocking_nn2 + 1) {
		LOGERROR("Insufficient history samples to do spectral clustering");
		clustered_blocking_variable_scaling.setOnes(sampler->num_variables);
		return false;
	} else {
		if (temperature == 1.0) {
			LOG("Clustered blocking adaptation - unique samples: %u", use_sample_ix.size());
		}
		size_t max_samples_for_clustering = sampler->adapt_proposal_max_samples * sampler->num_variables;
		if (use_sample_ix.size() > max_samples_for_clustering) {
			while (use_sample_ix.size() > max_samples_for_clustering) {
				unsigned int drop_sample = sampler->async[thread].rng.GetUnsignedInt(use_sample_ix.size()-1);
				use_sample_ix.erase(use_sample_ix.begin() + drop_sample);
			}
			if (temperature == 1.0) {
				LOG("Clustered blocking adaptation - downsampling to %u samples for spectral clustering", max_samples_for_clustering);
			}
		}
	}
	clustered_blocking_scaled_samples = MatrixReal::Zero(use_sample_ix.size(), sampler->num_variables);
	for (size_t i = 0; i < use_sample_ix.size(); i++) {
		for (size_t j = 0; j < sampler->num_variables; j++) {
			clustered_blocking_scaled_samples(i, j) = sample_history_d(j, use_sample_ix[i]) / clustered_blocking_variable_scaling(j);
		}
	}
	n = clustered_blocking_scaled_samples.rows();

	if (output_update_info) {
		update_info_output.AddMatrix(output_update_info_group, "clustering_input_samples", clustered_blocking_scaled_samples);
		update_info_output.AddVector(output_update_info_group, "clustering_input_sample_scaling", clustered_blocking_variable_scaling);
	}

	// Calculate kernel matrix
	MatrixReal K(n, n);
	VectorReal dists(n);
	clustered_blocking_sample_scale = VectorReal::Ones(n);
	clustered_blocking_nearest_neighbors.resize(n);
	for (ptrdiff_t si = 0; si < n; si++) {
		for (ptrdiff_t sj = 0; sj < n; sj++) {
			VectorReal v = clustered_blocking_scaled_samples.row(si) - clustered_blocking_scaled_samples.row(sj);
			dists(sj) = v.dot(v);
			K(si, sj) = dists(sj);
		}

		std::vector<ptrdiff_t> ranks = rank(dists);
		auto it = std::find(ranks.begin(), ranks.end(), sampler->clustered_blocking_nn + 1);
		ASSERT(it != ranks.end());

		clustered_blocking_sample_scale(si) = sqrt(dists(it - ranks.begin()));
		clustered_blocking_nearest_neighbors[si].resize(sampler->clustered_blocking_nn2);
		for (int i = 1; i < sampler->clustered_blocking_nn2 + 1; i++) {
			auto it = std::find(ranks.begin(), ranks.end(), i + 1);
			ASSERT(it != ranks.end());
			clustered_blocking_nearest_neighbors[si][i - 1] = it - ranks.begin();
		}
	}

	for (ptrdiff_t si = 1; si < n; si++) {
		for (ptrdiff_t sj = 0; sj < si; sj++) {
			Real cnns = 0.0;
			for (int i = 0; i < sampler->clustered_blocking_nn2; i++) {
				for (int j = 0; j < sampler->clustered_blocking_nn2; j++) {
					if (clustered_blocking_nearest_neighbors[si][i] == clustered_blocking_nearest_neighbors[sj][j]) {
						cnns += 1.0;
						break;
					}
				}
			}

			Real x = K(si, sj) / (clustered_blocking_sample_scale(si) * clustered_blocking_sample_scale(sj) * (cnns + 1.0));
			K(si, sj) = exp(-x);
			K(sj, si) = K(si, sj);
		}
	}
	K.diagonal().array() = 0.0;

	if (output_update_info) {
		update_info_output.AddMatrix(output_update_info_group, "K", K);
	}

	Eigen::DiagonalMatrix<Real, Eigen::Dynamic> D;
	D.diagonal() = rowSum(K).cwiseSqrt().cwiseInverse();
	MatrixReal L = D * K * D;

	Eigen::SelfAdjointEigenSolver<MatrixReal> eig;
	eig.compute(L);
	//if (temperature == 1.0) {
	//	for (int i = 0; i < D.cols(); i++) {
	//		LOG("Clustered blocking adaptation - eigenvalue %d: %g", i, eig.eigenvalues()(i));
	//	}
	//}

	MatrixReal Y(n, sampler->clustered_blocking_n_clusters);
	for (int i = 0; i < sampler->clustered_blocking_n_clusters; i++) {
		Y.col(i) = eig.eigenvectors().col(n - i - 1);
	}
	for (ptrdiff_t i = 0; i < n; i++) {
		Real d = Y.row(i).dot(Y.row(i));
		Y.row(i) *= rsqrt(d);
	}
	clustered_blocking_spectral_decomposition = Y;

	if (output_update_info) {
		update_info_output.AddMatrix(output_update_info_group, "Y", Y);
	}

	std::vector<ptrdiff_t> assignment;
	NaiveKMeans(Y, sampler->clustered_blocking_n_clusters, 10, 100, clustered_blocking_kmeans_centroids, assignment, sampler->async[thread].rng);

	if (output_update_info) {
		std::vector<int> assignment_int(assignment.size());
		for (ptrdiff_t i = 0; i < assignment.size(); i++) {
			ASSERT(assignment[i] <= std::numeric_limits<int>::max());
			assignment_int[i] = (int)assignment[i];
		}
		update_info_output.AddVector(output_update_info_group, "assignment", assignment_int);
	}

	// Could perhaps do slightly better here and assign the whole history rather than the subsample;
	// then calculate the covariances over the whole history

	// Calculate covariance in every cluster
	MatrixReal max_abs_corr = MatrixReal::Constant(sampler->num_variables, sampler->num_variables, 0.0);
	std::vector< MatrixReal > covariances(sampler->clustered_blocking_n_clusters);
	for (ptrdiff_t i = 0; i < sampler->clustered_blocking_n_clusters; i++) {
		size_t count = 0;
		for (ptrdiff_t j = 0; j < n; j++) {
			if (assignment[j] == i) {
				count++;
			}
		}

		MatrixReal selected(count, sampler->num_variables);
		ptrdiff_t k = 0;
		for (ptrdiff_t j = 0; j < n; j++) {
			if (assignment[j] == i) {
				selected.row(k++) = sample_history_d.col(use_sample_ix[j]);
			}
		}

		covariances[i] = cov(selected);
		MatrixReal corr = cor(selected);
		max_abs_corr.array() = max_abs_corr.array().max(corr.array().abs());
	}

	// Calculate distance for clustering based on the maximum correlation found in any of the clusters
	MatrixReal distance = MatrixReal::Zero(max_abs_corr.cols(), max_abs_corr.cols());
	for (ptrdiff_t i = 0; i < max_abs_corr.rows(); i++) {
		for (ptrdiff_t j = 0; j < max_abs_corr.cols(); j++) {
			distance(i, j) = 1.0 - max_abs_corr(i, j);
		}
	}

	// Cluster the variables
	std::vector< std::set<size_t> > variable_clusters;
	bcm3::TreeCluster(distance, variable_clusters, 0.5);

	// Use the same variable clustering for all the sample clusters, but use only the samples belonging
	// to a cluster to calculate the covariance for this block of variables in this cluster
	clustered_blocking_blocks.resize(sampler->clustered_blocking_n_clusters);
	for (ptrdiff_t ci = 0; ci < sampler->clustered_blocking_n_clusters; ci++) {
		Cluster& c = clustered_blocking_blocks[ci];
		c.blocks.clear();
		c.blocks.resize(variable_clusters.size());

		std::vector<int> block_assignment(sampler->num_variables);
		MatrixReal joint_covariance_for_output = MatrixReal::Zero(sampler->num_variables, sampler->num_variables);

		for (size_t i = 0; i < variable_clusters.size(); i++) {
			for (std::set<size_t>::iterator vi = variable_clusters[i].begin(); vi != variable_clusters[i].end(); ++vi) {
				c.blocks[i].variable_indices.push_back((int)*vi);
			}

			if (output_update_info && ci == 0) {
				update_info_output.AddVector(output_update_info_group, "block" + std::to_string(i) + std::string("_variables"), c.blocks[i].variable_indices);
			}

			if (variable_clusters[i].size() > 1) {
				for (std::vector<int>::iterator bi = c.blocks[i].variable_indices.begin(); bi != c.blocks[i].variable_indices.end(); bi++) {
					block_assignment[*bi] = i;
				}

				MatrixReal cov = MatrixReal::Zero(variable_clusters[i].size(), variable_clusters[i].size());
				for (size_t x = 0; x < variable_clusters[i].size(); x++) {
					for (size_t y = 0; y < variable_clusters[i].size(); y++) {
						cov(x, y) = covariances[ci](c.blocks[i].variable_indices[x], c.blocks[i].variable_indices[y]);
						joint_covariance_for_output(c.blocks[i].variable_indices[x], c.blocks[i].variable_indices[y]) = cov(x, y);
					}
				}
				c.blocks[i].cov = cov;
				c.blocks[i].covariance_llt = cov.llt();
				c.blocks[i].covariance_decomp = c.blocks[i].covariance_llt.matrixL();

				Real det = 0.0;
				for (size_t j = 0; j < variable_clusters[i].size(); j++) {
					det += log(c.blocks[i].covariance_decomp(j, j));
				}
				c.blocks[i].logC = -det - 0.5 * variable_clusters[i].size() * log(2.0 * M_PI);
			} else {
				int ix = c.blocks[i].variable_indices[0];

				// Make sure the proposal variance is not too small
				Real var = covariances[ci](ix, ix);
				Real prior_var;
				if (!sampler->prior->EvaluateMarginalVariance(ix, prior_var)) {
					// TODO - somehow come up with something reasonable?
					prior_var = 1.0;
				}
				var = std::max(var, (Real)1e-6 * prior_var);
				c.blocks[i].sigma = sqrt(var);

				block_assignment[ix] = i;
				joint_covariance_for_output(ix, ix) = var;
			}

			c.blocks[i].scale = 1.0;
			c.blocks[i].current_acceptance_rate_ema = sampler->target_acceptance_rate;
		}

		if (output_update_info) {
			update_info_output.AddVector(output_update_info_group, "cluster" + std::to_string(ci) + std::string("_block_assignment"), block_assignment);
			update_info_output.AddMatrix(output_update_info_group, "cluster" + std::to_string(ci) + std::string("_covariance"), joint_covariance_for_output);
		}
	}

	current_cluster_assignment = GetSampleCluster(current_var_values);

	if (output_update_info) {
		update_info_output.Close();
	}
	adaptation_iter++;
	return true;
}

ptrdiff_t SamplerPTChain::GetSampleCluster(const VectorReal& x)
{
	ptrdiff_t n = clustered_blocking_scaled_samples.rows();
	if (n <= 0) {
		ASSERT(clustered_blocking_blocks.size() == 1);
		return 0;
	}
	
	VectorReal y = x.cwiseQuotient(clustered_blocking_variable_scaling);
	VectorReal dists(n);
	for (ptrdiff_t i = 0; i < n; i++) {
		VectorReal v = clustered_blocking_scaled_samples.row(i).transpose() - y;
		dists(i) = v.dot(v);
	}

	std::vector<ptrdiff_t> ranks = rank(dists);

	auto it = std::find(ranks.begin(), ranks.end(), sampler->clustered_blocking_nn + 1);
	ASSERT(it != ranks.end());
	Real scale = sqrt(dists(it - ranks.begin()));

	std::vector<ptrdiff_t> nearest_neighbors(sampler->clustered_blocking_nn2);
	for (ptrdiff_t i = 1; i < sampler->clustered_blocking_nn2 + 1; i++) {
		auto it = std::find(ranks.begin(), ranks.end(), i + 1);
		ASSERT(it != ranks.end());
		nearest_neighbors[i - 1] = it - ranks.begin();
	}

	VectorReal B = VectorReal::Zero(n);
	for (ptrdiff_t si = 0; si < n; si++) {
		Real cnns = 0.0;
		for (ptrdiff_t i = 0; i < sampler->clustered_blocking_nn2; i++) {
			for (ptrdiff_t j = 0; j < sampler->clustered_blocking_nn2; j++) {
				if (clustered_blocking_nearest_neighbors[si][i] == nearest_neighbors[j]) {
					cnns += 1.0;
					break;
				}
			}
		}

		B(si) = exp(-dists(si) / (scale * clustered_blocking_sample_scale(si) * (cnns + 1.0)));
	}

	VectorReal f = B.transpose() * clustered_blocking_spectral_decomposition;

	Real max_dist = -std::numeric_limits<Real>::infinity();
	ptrdiff_t max_ix = 0;
	for (ptrdiff_t i = 0; i < sampler->clustered_blocking_n_clusters; i++) {
		//VectorReal v = f.transpose() - clustered_blocking_kmeans_centroids.row(i);
		Real d = f.dot(clustered_blocking_kmeans_centroids.row(i));

		if (d > max_dist) {
			max_dist = d;
			max_ix = i;
		}
	}

	return max_ix;
}

void SamplerPTChain::UpdateScale(Real& scale, Real current_acceptance_rate_ema, RNG& rng) const
{
	if (sampler->adapt_proposal_samples > 0 && !sampler->proposal_scaling_adaptations_done) {
		Real learn_rate = 1.0 + rng.GetReal() * sampler->proposal_scaling_learning_rate;
		if (current_acceptance_rate_ema < 0.952381 * sampler->target_acceptance_rate) {
			scale /= learn_rate;
			scale = std::max(scale, (Real)1e-4);
		} else if (current_acceptance_rate_ema > 1.05 * sampler->target_acceptance_rate) {
			scale *= learn_rate;
			scale = std::min(scale, (Real)10.0);
		}
	}
}

Real SamplerPTChain::GetTDistScaleFactor(RNG& rng) const
{
	if (sampler->proposal_t_df > 0.0) {
		Real w = rng.GetGamma(0.5 * sampler->proposal_t_df, 0.5 * sampler->proposal_t_df);
		return bcm3::rsqrt(w);
	} else {
		return 1.0;
	}
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

void SamplerPTChain::AcceptMutate(bool accept, const VectorReal& nv, Real lprior, Real llh, Real lpowerposterior, Real& current_acceptance_rate_ema)
{
	const Real ema_alpha = 2.0 / (sampler->proposal_scaling_ema_period + 1);
	if (accept) {
		accepted_mutate++;

		current_var_values = nv;
		this->lprior = lprior;
		this->llh = llh;
		this->lpowerposterior = lpowerposterior;

		current_acceptance_rate_ema += (1.0 - current_acceptance_rate_ema) * ema_alpha;
	} else {
		current_acceptance_rate_ema += (0.0 - current_acceptance_rate_ema) * ema_alpha;
	}
}

void SamplerPTChain::AddCurrentSampleToHistory()
{
	sample_history_n_s++;
	if (sample_history_n_s == sample_history_subsampling) {
		size_t sample_history_ix = sample_history_n;
		if (sample_history_ix >= sample_history.cols()) {
			// Wrap around and start filling from the start again
			sample_history_ix = sample_history_ix % sample_history.cols();
		}
		sample_history.col(sample_history_ix) = current_var_values.cast<float>();
		/*
		if (temperature == 1.0) {
			LOG("Adding sample at history ix %zd: %g,%g,%g,%g,%g,%g",
				sample_history_ix,
				current_var_values(0),
				current_var_values(1),
				current_var_values(2),
				current_var_values(3),
				current_var_values(4),
				current_var_values(5));
		}
		*/
		sample_history_n++;
		sample_history_n_s = 0;
	}
}

}
