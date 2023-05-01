#include "Utils.h"
#include "BlockingStrategyClusteredTurek.h"
#include "Clustering.h"
#include "SampleHistory.h"
#include "SampleHistoryClustering.h"
#include "SummaryStats.h"

namespace bcm3 {

	bool BlockingStrategyClusteredTurek::UsesClustering()
	{
		return true;
	}

	std::vector< std::vector<ptrdiff_t> > BlockingStrategyClusteredTurek::GetBlocks(const std::unique_ptr<SampleHistory>& sample_history, const std::shared_ptr<SampleHistoryClustering> clustering)
	{
		ASSERT(clustering);

		std::vector< std::vector<ptrdiff_t> > blocks;

		if (sample_history->GetSampleCount() > 2) {
			size_t n = clustering->GetNumClusteredSamples();

			// Calculate covariance in every sample cluster
			MatrixReal max_abs_corr = MatrixReal::Constant(num_variables, num_variables, 0.0);
			for (ptrdiff_t i = 0; i < clustering->GetNumClusters(); i++) {
				size_t count = 0;
				for (ptrdiff_t j = 0; j < n; j++) {
					if (clustering->GetClusterAssignment(j) == i) {
						count++;
					}
				}

				MatrixReal selected(count, num_variables);
				ptrdiff_t k = 0;
				for (ptrdiff_t j = 0; j < n; j++) {
					if (clustering->GetClusterAssignment(j) == i) {
						selected.row(k++) = sample_history->GetHistorySample(clustering->GetHistorySampleIx(j)).cast<Real>();
					}
				}

				MatrixReal corr = cor(selected);
				max_abs_corr.array() = max_abs_corr.array().max(corr.array().abs());
			}

			// Calculate distance for clustering based on the maximum correlation found in any of the sample clusters
			MatrixReal distance = MatrixReal::Zero(max_abs_corr.cols(), max_abs_corr.cols());
			for (ptrdiff_t i = 0; i < max_abs_corr.rows(); i++) {
				for (ptrdiff_t j = 0; j < max_abs_corr.cols(); j++) {
					distance(i, j) = 1.0 - max_abs_corr(i, j);
				}
			}

			// Cluster the variables
			std::vector< std::set<size_t> > variable_clusters;
			bcm3::TreeCluster(distance, variable_clusters, 0.5);

			// Create blocks based on the clustered variables
			blocks.resize(variable_clusters.size());
			for (size_t i = 0; i < variable_clusters.size(); i++) {
				for (std::set<size_t>::iterator vi = variable_clusters[i].begin(); vi != variable_clusters[i].end(); ++vi) {
					blocks[i].push_back((int)*vi);
				}
			}
		} else {
			blocks.resize(num_variables);
			for (ptrdiff_t i = 0; i < num_variables; i++) {
				blocks[i].resize(1);
				blocks[i][0] = i;
			}
		}

		return blocks;
	}

}
