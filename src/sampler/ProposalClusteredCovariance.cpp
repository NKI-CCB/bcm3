#include "Utils.h"
#include "GMM.h"
#include "NetCDFBundler.h"
#include "Prior.h"
#include "ProposalClusteredCovariance.h"
#include "SampleHistoryClustering.h"
#include "SummaryStats.h"
#include "checks.h"

namespace bcm3 {

	ProposalClusteredCovariance::ProposalClusteredCovariance()
		: t_dof(0.0)
		, selected_component(-1)
	{
	}

	ProposalClusteredCovariance::~ProposalClusteredCovariance()
	{
	}

	Real Proposal::CalculateMHRatio(const VectorReal& current_position, ptrdiff_t curpos_cluster_assignment, const VectorReal& new_position, ptrdiff_t newpos_cluster_assignment)
	{
		return 0.0;
	}

	void ProposalClusteredCovariance::GetNewSample(const VectorReal& current_position, ptrdiff_t history_cluster_assignment, VectorReal& new_position, RNG& rng)
	{
		if (history_cluster_assignment == -1) {
			// There wasn't any history to assign to, in which case we should still have a GMM with 1 component
			ASSERT(gmm->GetNumComponents() == 1);
			selected_component = 0;
		} else {
			ASSERT(history_cluster_assignment < gmm->GetNumComponents());
			selected_component = history_cluster_assignment;
		}

		Real t_scale;
		if (t_dof > 0.0) {
			Real w = rng.GetGamma(0.5 * t_dof, 0.5 * t_dof);
			t_scale = bcm3::rsqrt(w);
		} else {
			t_scale = 1.0;
		}

		// Multivariate normal random value
		VectorReal x = rng.GetMultivariateUnitNormal(num_variables);
		x = gmm->GetCovarianceDecomp(selected_component).matrixL() * x;
		x *= t_scale * scales(selected_component);
		new_position = x + current_position;

		if (!transform_to_unbounded) {
			for (ptrdiff_t i = 0; i < num_variables; i++) {
				new_position(i) = ReflectOnBounds(new_position(i), variable_bounds[i].lower, variable_bounds[i].upper);
			}
		}
	}

	Real ProposalClusteredCovariance::CalculateMHRatio(const VectorReal& current_position, ptrdiff_t curpos_cluster_assignment, const VectorReal& new_position, ptrdiff_t newpos_cluster_assignment)
	{
		if (curpos_cluster_assignment == -1) {
			return 0.0;
		}
		ASSERT(curpos_cluster_assignment < gmm->GetNumComponents());
		ASSERT(newpos_cluster_assignment < gmm->GetNumComponents());

		Real log_mh_ratio;
		if (newpos_cluster_assignment != curpos_cluster_assignment) {
			VectorReal v = new_position - current_position;
			v /= scales(curpos_cluster_assignment);
			gmm->GetCovarianceDecomp(curpos_cluster_assignment).matrixL().solveInPlace(v);
			Real log_fwd_p = gmm->GetLogC(curpos_cluster_assignment) - 0.5 * v.dot(v);

			v = current_position - new_position;
			v /= scales(newpos_cluster_assignment);
			gmm->GetCovarianceDecomp(newpos_cluster_assignment).matrixL().solveInPlace(v);
			Real log_bwd_p = gmm->GetLogC(newpos_cluster_assignment) - 0.5 * v.dot(v);

			log_mh_ratio = log_bwd_p - log_fwd_p;
		} else {
			// Proposal is symmetrical within cluster
			log_mh_ratio = 0.0;
		}
		return log_mh_ratio;
	}

	bool ProposalClusteredCovariance::UsesClustering()
	{
		return true;
	}

	void ProposalClusteredCovariance::Update(RNG& rng)
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

	void ProposalClusteredCovariance::NotifyAccepted(bool accepted)
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

	void ProposalClusteredCovariance::LogInfo() const
	{
		LOG("  Clustered covariance proposal with %u components", gmm->GetNumComponents());
		for (ptrdiff_t i = 0; i < gmm->GetNumComponents(); i++) {
			LOG("   Component %u; scale=%8.5f, weight=%8.5f, condition number=%6g", i + 1, scales(i), gmm->GetWeights()(i), 1.0 / gmm->GetCovarianceDecomp(i).rcond());
		}
	}

	void ProposalClusteredCovariance::WriteToFile(const std::string& fn, std::string adaptation_group, std::vector<ptrdiff_t>& variable_indices)
	{
		NetCDFBundler update_info_output;
		if (update_info_output.Open(fn)) {
			update_info_output.AddGroup(adaptation_group);
			
			std::vector<int> indices(variable_indices.size());
			for (int i = 0; i < variable_indices.size(); i++) {
				ASSERT(variable_indices[i] < std::numeric_limits<int>::max());
				indices[i] = (int)variable_indices[i];
			}
			update_info_output.AddVector(adaptation_group, "variable_indices", indices);

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

	bool ProposalClusteredCovariance::InitializeImpl(const MatrixReal& history, std::shared_ptr<Prior> prior, std::vector<ptrdiff_t>& variable_indices, RNG& rng, bool log_info)
	{
		if (history.rows() < 2) {
			// Start with vector sampler with a single component with diagonal covariance equal to prior variance
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
				return false;
			}
		} else {
			std::vector<VectorReal> means(sample_history_clustering->GetNumClusters(), VectorReal::Zero(num_variables));
			std::vector<MatrixReal> covs(sample_history_clustering->GetNumClusters(), MatrixReal::Identity(num_variables, num_variables));
			for (ptrdiff_t ci = 0; ci < sample_history_clustering->GetNumClusters(); ci++) {
				std::vector<ptrdiff_t> sample_ix = sample_history_clustering->GetSamplesFromCluster(ci);
				MatrixReal cluster_samples = history(sample_ix, Eigen::all);

				means[ci] = colMean(cluster_samples);
				covs[ci] = cov(cluster_samples);
				covs[ci].diagonal().array() += 1e-8;

				if (log_info) {
					LOG("  Cluster %zd based on %zu samples", ci, sample_ix.size());
				}

				if (!is_positive_semi_definite(covs[ci])) {
					LOGWARNING("  Cluster %zd - covariance matrix is not semi positive definite!", ci);
				}
			}

			VectorReal weights;
			weights.setConstant(sample_history_clustering->GetNumClusters(), 1.0 / sample_history_clustering->GetNumClusters());
			gmm = std::make_shared<GMM>();
			if (!gmm->Set(means, covs, weights)) {
				return false;
			}
		}

		scales.setConstant(gmm->GetNumComponents(), 2.38 / sqrt(num_variables));
		acceptance_rate_emas.setConstant(gmm->GetNumComponents(), target_acceptance_rate);

		return true;
	}

}
