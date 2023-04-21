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
		}

		return blocks;
	}

}
