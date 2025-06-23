#include "Utils.h"
#include "GMM.h"
#include "NetCDFBundler.h"
#include "Prior.h"
#include "ProposalGaussianMixture.h"
#include "SummaryStats.h"

namespace bcm3 {

	ProposalGaussianMixture::ProposalGaussianMixture(bool select_with_adjusted_AIC)
		: selected_component(-1)
		, select_with_adjusted_AIC(select_with_adjusted_AIC)
	{
	}

	ProposalGaussianMixture::~ProposalGaussianMixture()
	{
	}

	void ProposalGaussianMixture::GetNewSample(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, RNG& rng)
	{
		VectorReal fwd_resp = gmm->CalculateResponsibilities(current_position);
		selected_component = rng.Sample(fwd_resp);

		Real t_scale;
		if (t_dof > 0.0) {
			Real w = rng.GetGamma(0.5 * t_dof, 0.5 * t_dof);
			t_scale = bcm3::rsqrt(w);
		} else {
			t_scale = 1.0;
		}

		// Propose new parameters
		VectorReal x = rng.GetMultivariateUnitNormal(num_variables);
		x = gmm->GetCovarianceDecomp(selected_component).matrixL() * x;
		x *= t_scale * scales(selected_component);
		new_position = x + current_position;

		for (ptrdiff_t i = 0; i < num_variables; i++) {
			new_position(i) = ReflectOnBounds(new_position(i), variable_bounds[i].lower, variable_bounds[i].upper);
		}
	}

	Real ProposalGaussianMixture::CalculateMHRatio(const VectorReal& current_position, ptrdiff_t curpos_cluster_assignment, const VectorReal& new_position, ptrdiff_t newpos_cluster_assignment)
	{
		VectorReal fwd_resp = gmm->CalculateResponsibilities(current_position);
		VectorReal rev_resp = gmm->CalculateResponsibilities(new_position);
		Real fwd_logp = -std::numeric_limits<Real>::infinity();
		Real rev_logp = -std::numeric_limits<Real>::infinity();
		for (ptrdiff_t i = 0; i < gmm->GetNumComponents(); i++) {
			VectorReal v = (new_position - current_position) / (scales(i));

			VectorReal s = gmm->GetCovarianceDecomp(i).matrixL().solve(v);
			Real p = -log(square(scales(i))) + gmm->GetLogC(i) - 0.5 * s.dot(s) + log(fwd_resp(i));
			fwd_logp = bcm3::logsum(fwd_logp, p);

			s = gmm->GetCovarianceDecomp(i).matrixL().solve(-v);
			p = -log(square(scales(i))) + gmm->GetLogC(i) - 0.5 * s.dot(s) + log(rev_resp(i));
			rev_logp = bcm3::logsum(rev_logp, p);
		}
		Real log_mh_ratio = rev_logp - fwd_logp;
		return log_mh_ratio;
	}

	void ProposalGaussianMixture::Update(RNG& rng)
	{
		if (selected_component == -1) {
			// Very first update call; can't update yet
		} else {
			// Only update the scale of the component that was previously selected for a sample.
			Real& scale = scales(selected_component);
			const Real& acceptance_rate_ema = acceptance_rate_emas(selected_component);

			Real learn_rate = 1.0 + rng.GetReal() * scaling_learning_rate * gmm->GetNumComponents();
			if (acceptance_rate_ema < target_acceptance_rate / (1.0 - scaling_learning_rate)) {
				scale /= learn_rate;
				scale = std::max(scale, (Real)1e-4);
			} else if (acceptance_rate_ema > (1 + scaling_learning_rate) * target_acceptance_rate) {
				scale *= learn_rate;
				scale = std::min(scale, (Real)10.0);
			}
		}

		// We fully override the base Proposal scaling adaptation; so no need to call the base Proposal::Update
	}

	void ProposalGaussianMixture::NotifyAccepted(bool accepted)
	{
		const Real ema_alpha = 2.0 / (scaling_ema_period + 1);

		Real& acceptance_rate_ema = acceptance_rate_emas(selected_component);
		if (accepted) {
			acceptance_rate_ema += (1.0 - acceptance_rate_ema) * ema_alpha;
		} else {
			acceptance_rate_ema += (0.0 - acceptance_rate_ema) * ema_alpha;
		}

		// We fully override the base Proposal scaling adaptation; so no need to call the base Proposal::NotifyAccepted
	}

	void ProposalGaussianMixture::LogInfo() const
	{
		LOG("  GMM proposal with %u components", gmm->GetNumComponents());
		for (ptrdiff_t i = 0; i < gmm->GetNumComponents(); i++) {
			LOG("    Component %u; scale=%8.5f, weight=%8.5f, condition number=%6g", i + 1, scales(i), gmm->GetWeights()(i), 1.0 / gmm->GetCovarianceDecomp(i).rcond());
		}
	}

	void ProposalGaussianMixture::WriteToFile(const std::string& fn, std::string adaptation_group, std::vector<ptrdiff_t>& variable_indices)
	{
		NetCDFBundler update_info_output;
		if (update_info_output.Open(fn)) {
			update_info_output.AddGroup(adaptation_group);
			update_info_output.AddVector(adaptation_group, "gmm_weights", gmm->GetWeights());
			for (ptrdiff_t i = 0; i < gmm->GetNumComponents(); i++) {
				update_info_output.AddVector(adaptation_group, "cluster" + std::to_string(i) + std::string("_mean"), gmm->GetMean(i));
				update_info_output.AddMatrix(adaptation_group, "cluster" + std::to_string(i) + std::string("_covariance"), gmm->GetCovariance(i));
			}

			if (transform_to_unbounded) {
				MatrixReal bounds(variable_bounds.size(), 2);
				for (ptrdiff_t i = 0; i < variable_bounds.size(); i++) {
					bounds(i, 0) = variable_bounds[i].lower;
					bounds(i, 1) = variable_bounds[i].upper;
				}
				update_info_output.AddMatrix(adaptation_group, "transform_bounds", bounds);
			}

			update_info_output.Close();
		}
	}

	bool ProposalGaussianMixture::InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info)
	{
		gmm.reset();

		if (history.rows() >= 2) {
			// Calculate effective sample size
			size_t num_history_samples = history.rows();
			VectorReal ess(num_variables);
			for (size_t i = 0; i < num_variables; i++) {
				Real rho_t = 0;

				Real mu = mean((VectorReal)history.col(i));
				Real sigmaSq = var((VectorReal)history.col(i));

				int lag_max = std::max(5, (int)(10 * log10((Real)num_history_samples)));
				for (int lag = 1; lag < lag_max; lag++) {
					Real c = acf((VectorReal)history.col(i), (int)lag, mu, sigmaSq);
					rho_t += c;
				}

				ess(i) = (Real)num_history_samples / (1.0 + 2.0 * rho_t);
			}
			Real min_ess = ess.minCoeff();
			Real AIC_adjust_factor = min_ess / num_history_samples;

			if (log_info) {
				LOG("Fitting GMMs...");
				LOG("Minimum effective sample size: %g out of %zu samples; AIC_adjust_factor = %g", min_ess, num_history_samples, AIC_adjust_factor);
			}

			// Fit GMMs for increasing number of components and select the one with the lowest AIC
			Real best_aic = std::numeric_limits<Real>::infinity();
			static const size_t num_components[7] = { 1, 2, 3, 4, 5, 8, 13 };
			for (size_t i = 0; i < 7; i++) {
				std::shared_ptr<GMM> test_gmm_k = std::make_shared<GMM>();
				if (min_ess < num_components[i] * (1 + std::min(num_variables/2, (size_t)10))) {
					if (log_info) {
						LOG("GMM num_components=%2zu - not enough effective samples", num_components[i], test_gmm_k->GetAIC());
					}
				} else if (test_gmm_k->Fit(history, num_history_samples, num_components[i], rng, num_history_samples / min_ess, log_info)) {
					Real nparam = 0.5 * test_gmm_k->GetAIC() + test_gmm_k->GetLogLikelihood();
					Real adjusted_AIC = 2.0 * nparam - 2.0 * AIC_adjust_factor * test_gmm_k->GetLogLikelihood();
					if (log_info) {
						LOG("GMM num_components=%2zu - AIC=%.6g, adjusted AIC=%.6g", num_components[i], test_gmm_k->GetAIC(), adjusted_AIC);
					}

					if (select_with_adjusted_AIC) {
						if (adjusted_AIC < best_aic) {
							gmm = test_gmm_k;
							best_aic = test_gmm_k->GetAIC();
						}
					} else {
						if (test_gmm_k->GetAIC() < best_aic) {
							gmm = test_gmm_k;
							best_aic = test_gmm_k->GetAIC();
						}
					}
				} else {
					if (log_info) {
						LOG("GMM num_components=%2zu failed", num_components[i]);
					}
				}
			}

#if 1
			if (log_info) {
				if (gmm) {
					LOG("Selected GMM with %zu components (selection through full AIC)", gmm->GetNumComponents());
					for (ptrdiff_t i = 0; i < gmm->GetNumComponents(); i++) {
						std::stringstream str;
						str << gmm->GetMean(i).transpose();
						LOG("Mean component %zd: %s", i, str.str().c_str());

						LOG("Covariance component %zd:", i);
						for (ptrdiff_t j = 0; j < gmm->GetCovariance(i).rows(); j++) {
							std::stringstream().swap(str);
							str << gmm->GetCovariance(i).row(j);
							LOG("%s", str.str().c_str());
						}
					}
				} else {
					LOG("Unable to fit even a single Gaussian, resorting to a single Gaussian with covariance based on the prior");
				}
			}
#endif
		}

		if (!gmm) {
			// If the GMM fitting failed, or if we didn't have any samples to work with because it's the first initialization,
			// we start with vector sampler with a single component with diagonal covariance equal to prior variance
			VectorReal mean(num_variables);
			MatrixReal covariance(num_variables, num_variables);

			covariance.setZero();
			for (size_t i = 0; i < variable_indices.size(); i++) {
				Real prior_mean;
				if (!prior->EvaluateMarginalMean(variable_indices[i], prior_mean)) {
					// TODO - somehow come up with something reasonable?
					prior_mean = 0.0;
				}

				Real prior_var;
				if (!prior->EvaluateMarginalVariance(variable_indices[i], prior_var)) {
					// TODO - somehow come up with something reasonable?
					prior_var = 1.0;
				}

				mean(i) = prior_mean;
				covariance(i, i) = prior_var;
			}

			std::vector<VectorReal> means(1, mean);
			std::vector<MatrixReal> covs(1, covariance);
			VectorReal weights;
			weights.setConstant(1, 1);

			gmm = std::make_shared<GMM>();
			if (!gmm->Set(means, covs, weights)) {
				LOGERROR("Creation of GMM failed");
				return false;
			}
		}

		scales.setConstant(gmm->GetNumComponents(), 2.38 / sqrt(num_variables));
		acceptance_rate_emas.setConstant(gmm->GetNumComponents(), target_acceptance_rate);

		return true;
	}

}
