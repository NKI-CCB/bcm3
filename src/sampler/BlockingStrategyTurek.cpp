#include "Utils.h"
#include "BlockingStrategyTurek.h"
#include "Clustering.h"
#include "SampleHistory.h"

namespace bcm3 {

	std::vector< std::vector<ptrdiff_t> > BlockingStrategyTurek::GetBlocks(const std::unique_ptr<SampleHistory>& sample_history)
	{
		MatrixReal empirical_correlation = sample_history->GetEmpiricalCorrelation();

		// Calculate distance for clustering
		MatrixReal distance = MatrixReal::Zero(empirical_correlation.cols(), empirical_correlation.cols());
		for (ptrdiff_t i = 0; i < empirical_correlation.rows(); i++) {
			for (ptrdiff_t j = 0; j < empirical_correlation.cols(); j++) {
				distance(i, j) = 1.0 - fabs(empirical_correlation(i, j));
			}
		}

		std::vector< std::set<size_t> > clusters;
		bcm3::TreeCluster(distance, clusters, 0.5);
		std::vector< std::vector<ptrdiff_t> > blocks(clusters.size());
		for (size_t i = 0; i < clusters.size(); i++) {
			for (std::set<size_t>::iterator vi = clusters[i].begin(); vi != clusters[i].end(); ++vi) {
				blocks[i].push_back((int)*vi);
			}
#if 0
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
				Real var = empirical_covariance(ix, ix);
				Real prior_var;
				if (!sampler->prior->EvaluateMarginalVariance(ix, prior_var)) {
					// TODO - somehow come up with something reasonable?
					prior_var = 1.0;
				}
				var = std::max(var, (Real)1e-6 * prior_var);
				sampler_blocks[i].sigma = sqrt(var);
			}

			blocks[i].scale = 1.0;
			blocks[i].current_acceptance_rate_ema = sampler->target_acceptance_rate;
#endif
		}

		return blocks;
	}

}
