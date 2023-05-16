#pragma once

#include <boost/dynamic_bitset.hpp>

namespace bcm3 {

	class SampleHistory;

	class SampleHistoryClustering
	{
	public:
		SampleHistoryClustering(size_t max_samples, size_t nn, size_t nn2, size_t num_clusters);
		~SampleHistoryClustering();

		bool Cluster(const std::unique_ptr<SampleHistory>& sample_history, size_t discard_first_samples, RNG& rng, bool log_info, const std::string& output_path);
		bool AssignAllHistorySamples(std::vector<ptrdiff_t> history_sample_ix, MatrixReal& samples);

		ptrdiff_t GetSampleCluster(const VectorReal& sample) const;
		std::vector<ptrdiff_t> GetSamplesFromCluster(ptrdiff_t cluster_ix) const;
		std::vector<ptrdiff_t> GetAllSamplesAssignedToCluster(ptrdiff_t cluster_ix) const;
		inline size_t GetNumClusters() const { return num_clusters; }
		inline size_t GetNumClusteredSamples() const { return used_sample_ix.size(); }
		inline ptrdiff_t GetHistorySampleIx(ptrdiff_t i) const { return used_sample_ix[i]; }
		inline ptrdiff_t GetClusterAssignment(ptrdiff_t i) const { return cluster_assignment[i]; }

	private:
		// Settings
		size_t max_samples;
		size_t density_aware_kernel_nn;
		size_t density_aware_kernel_nn2;
		size_t num_clusters;

		// Runtime variables
		std::vector<ptrdiff_t> used_sample_ix;
		MatrixReal scaled_samples;
		VectorReal variable_scaling;
		VectorReal density_aware_kernel_sample_scale;
		std::vector< std::vector<ptrdiff_t> > density_aware_kernel_nearest_neighbors;
		std::vector< boost::dynamic_bitset<> > density_aware_kernel_nearest_neighbors_bitset;
		MatrixReal spectral_decomposition;
		MatrixReal spectral_kmeans_centroids;
		std::vector<ptrdiff_t> cluster_assignment;
		int clustering_iter;

		std::vector<ptrdiff_t> all_sample_cluster_assignment;
	};

}
