#include "Utils.h"
#include "Clustering.h"
#include "RNG.h"
#include "SampleHistory.h"
#include "SampleHistoryClustering.h"
#include "SummaryStats.h"
#include "VectorUtils.h"

namespace bcm3 {

	SampleHistoryClustering::SampleHistoryClustering(size_t max_samples)
		: max_samples(max_samples)
		, density_aware_kernel_nn(3)
		, density_aware_kernel_nn2(7)
		, num_clusters(4)
	{
	}

	SampleHistoryClustering::~SampleHistoryClustering()
	{
	}

	bool SampleHistoryClustering::Cluster(const std::unique_ptr<SampleHistory>& sample_history, size_t discard_first_samples, RNG& rng, bool log_info)
	{
		size_t num_variables = sample_history->samples.rows();

		ptrdiff_t n = sample_history->GetSampleCount();
		if (n < 1) {
			return false;
		}
		if (log_info) {
			LOG("Sample history clustering - samples in history: %d", n);
		}

#if 0
		NetCDFBundler update_info_output;
		bool output_update_info = false;
		std::string output_update_info_group = std::string("info") + std::to_string(adaptation_iter);
		if (sampler->output_proposal_adaptation && temperature == sampler->current_temperatures.tail(1)(0)) {
			if (update_info_output.Open(sampler->output_path + "sampler_adaptation.nc")) {
				update_info_output.AddGroup(output_update_info_group);
				output_update_info = true;
			}
		}
#endif

		// Calculate standard deviation over entire history
		// TODO - can do this with the float matrix as input
		MatrixReal sample_history_d = sample_history->samples.block(0, 0, num_variables, n).cast<Real>();
		variable_scaling = rowSd(sample_history_d);

		for (ptrdiff_t i = 0; i < variable_scaling.size(); i++) {
			if (variable_scaling(i) <= 0.0 || std::isinf(variable_scaling(i)) || std::isnan(variable_scaling(i))) {
				LOGERROR("0, infinite or nan sd in sample history for variable %d", i);
				variable_scaling.setOnes(num_variables);
				return false;
			}
		}

		// Retrieve all unique samples
		used_sample_ix.clear();
		for (ptrdiff_t i = 0; i < n; i++) {
			bool duplicate = false;
			for (ptrdiff_t j = 0; j < i; j++) {
				Eigen::VectorXf v = sample_history->samples.col(i) - sample_history->samples.col(j);
				float d = v.dot(v);
				if (d < std::numeric_limits<float>::epsilon()) {
					duplicate = true;
				}
			}
			if (!duplicate) {
				if (i < discard_first_samples) {
					// Throw away very first samples as burnin, even for the adaptation
				} else {
					used_sample_ix.push_back(i);
				}
			}
		}

		if (used_sample_ix.size() < (ptrdiff_t)density_aware_kernel_nn2 + 1) {
			LOGERROR("Insufficient history samples to do spectral clustering");
			variable_scaling.setOnes(num_variables);
			return false;
		} else {
			if (log_info) {
				LOG("Sample history clustering - unique samples: %u", used_sample_ix.size());
			}
			if (used_sample_ix.size() > max_samples) {
				while (used_sample_ix.size() > max_samples) {
					unsigned int drop_sample = rng.GetUnsignedInt(used_sample_ix.size() - 1);
					used_sample_ix.erase(used_sample_ix.begin() + drop_sample);
				}
				if (log_info) {
					LOG("Sample history clustering - downsampled to %u samples for spectral clustering", max_samples);
				}
			}
		}


		scaled_samples = MatrixReal::Zero(used_sample_ix.size(), num_variables);
		for (ptrdiff_t i = 0; i < used_sample_ix.size(); i++) {
			for (ptrdiff_t j = 0; j < num_variables; j++) {
				scaled_samples(i, j) = sample_history_d(j, used_sample_ix[i]) / variable_scaling(j);
			}
		}
		n = scaled_samples.rows();

#if 0
		if (output_update_info) {
			update_info_output.AddMatrix(output_update_info_group, "clustering_input_samples", clustered_blocking_scaled_samples);
			update_info_output.AddVector(output_update_info_group, "clustering_input_sample_scaling", clustered_blocking_variable_scaling);
		}
#endif

		// Calculate kernel matrix
		MatrixReal K(n, n);
		VectorReal dists(n);
		density_aware_kernel_sample_scale = VectorReal::Ones(n);
		density_aware_kernel_nearest_neighbors.resize(n);
		density_aware_kernel_nearest_neighbors_bitset.resize(n);
		for (ptrdiff_t si = 0; si < n; si++) {
			for (ptrdiff_t sj = 0; sj < n; sj++) {
				VectorReal v = scaled_samples.row(si) - scaled_samples.row(sj);
				dists(sj) = v.dot(v);
				K(si, sj) = dists(sj);
			}

			std::vector<ptrdiff_t> ranks = rank(dists);
			auto it = std::find(ranks.begin(), ranks.end(), density_aware_kernel_nn + 1);
			ASSERT(it != ranks.end());

			density_aware_kernel_sample_scale(si) = sqrt(dists(it - ranks.begin()));
			density_aware_kernel_nearest_neighbors[si].resize(density_aware_kernel_nn2);
			density_aware_kernel_nearest_neighbors_bitset[si].resize(n);
			for (int i = 1; i < density_aware_kernel_nn2 + 1; i++) {
				auto it = std::find(ranks.begin(), ranks.end(), i + 1);
				ASSERT(it != ranks.end());
				density_aware_kernel_nearest_neighbors[si][i - 1] = it - ranks.begin();
				density_aware_kernel_nearest_neighbors_bitset[si][it - ranks.begin()] = true;
			}
		}

		for (ptrdiff_t si = 1; si < n; si++) {
			for (ptrdiff_t sj = 0; sj < si; sj++) {
				Real cnns = 0.0;
				for (ptrdiff_t i = 0; i < density_aware_kernel_nn2; i++) {
					cnns += (Real)density_aware_kernel_nearest_neighbors_bitset[si][density_aware_kernel_nearest_neighbors_bitset[sj][i]];
				}

				Real val = -K(si, sj) / (density_aware_kernel_sample_scale(si) * density_aware_kernel_sample_scale(sj) * (cnns + 1.0));
				K(si, sj) = approx_exp(val);
				K(sj, si) = K(si, sj);
			}
		}
		K.diagonal().array() = 0.0;

#if 0
		if (output_update_info) {
			update_info_output.AddMatrix(output_update_info_group, "K", K);
		}
#endif

		// Calculate graph Laplacian
		Eigen::DiagonalMatrix<Real, Eigen::Dynamic> D;
		D.diagonal() = rowSum(K).cwiseSqrt().cwiseInverse();
		MatrixReal L = D * K * D;

		Eigen::SelfAdjointEigenSolver<MatrixReal> eig;
		eig.compute(L);

		MatrixReal Y(n, num_clusters);
		for (int i = 0; i < num_clusters; i++) {
			Y.col(i) = eig.eigenvectors().col(n - i - 1);
		}
		for (ptrdiff_t i = 0; i < n; i++) {
			Real d = Y.row(i).dot(Y.row(i));
			Y.row(i) *= rsqrt(d);
		}
		spectral_decomposition = Y;

#if 0
		if (output_update_info) {
			update_info_output.AddMatrix(output_update_info_group, "Y", Y);
		}
#endif

		NaiveKMeans(Y, num_clusters, 10, 100, spectral_kmeans_centroids, cluster_assignment, rng);

#if 0
		if (output_update_info) {
			std::vector<int> assignment_int(assignment.size());
			for (ptrdiff_t i = 0; i < assignment.size(); i++) {
				ASSERT(assignment[i] <= std::numeric_limits<int>::max());
				assignment_int[i] = (int)assignment[i];
			}
			update_info_output.AddVector(output_update_info_group, "assignment", assignment_int);
		}
#endif

		return true;
	}
	
	ptrdiff_t SampleHistoryClustering::GetSampleCluster(const VectorReal& sample) const
	{
		ptrdiff_t n = scaled_samples.rows();
		if (n == 0) {
			return -1;
		}

		int nearest_neighbours_needed = std::max(density_aware_kernel_nn, density_aware_kernel_nn2);
		std::vector<ptrdiff_t> nearest_neighbors(nearest_neighbours_needed);
		VectorReal nearest_neighbors_dists = VectorReal::Constant(nearest_neighbours_needed, std::numeric_limits<Real>::infinity());

		VectorReal y = sample.cwiseQuotient(variable_scaling);
		VectorReal dists(n);
		VectorReal buffer;
		for (ptrdiff_t i = 0; i < n; i++) {
			VectorReal x = scaled_samples.row(i);
			x -= y; // Very weird, but if you do the - y directly from the scaled_samples.row() you get errors.
			Real dist = x.dot(x);
			dists(i) = dist;
			for (int j = 0; j < nearest_neighbours_needed; j++) {
				if (dist < nearest_neighbors_dists[j]) {
					// Move the rest backward
					for (int k = nearest_neighbours_needed - 1; k > j; k--) {
						nearest_neighbors[k] = nearest_neighbors[k - 1];
						nearest_neighbors_dists(k) = nearest_neighbors_dists(k - 1);
					}

					// Insert this point
					nearest_neighbors[j] = i;
					nearest_neighbors_dists(j) = dist;
					break;
				}
			}
		}

		// The nearest neighbors list now does not include self, so do not add 1 like was done during adaptation
		Real scale = sqrt(nearest_neighbors_dists[density_aware_kernel_nn]);

		VectorReal B = VectorReal::Zero(n);
		for (ptrdiff_t si = 0; si < n; si++) {
			Real cnns = 0.0;
			for (ptrdiff_t i = 0; i < density_aware_kernel_nn2; i++) {
				cnns += (Real)density_aware_kernel_nearest_neighbors_bitset[si][nearest_neighbors[i]];
			}

			B(si) = exp(-dists(si) / (scale * density_aware_kernel_sample_scale(si) * (cnns + 1.0)));
		}

		VectorReal f = B.transpose() * spectral_decomposition;

		Real max_dist = -std::numeric_limits<Real>::infinity();
		ptrdiff_t max_ix = 0;
		for (ptrdiff_t i = 0; i < num_clusters; i++) {
			Real d = f.dot(spectral_kmeans_centroids.row(i));

			if (d > max_dist) {
				max_dist = d;
				max_ix = i;
			}
		}

		return max_ix;
	}

	std::vector<ptrdiff_t> SampleHistoryClustering::GetSamplesFromCluster(ptrdiff_t cluster_ix) const
	{
		std::vector<ptrdiff_t> samples;
		for (ptrdiff_t i = 0; i < cluster_assignment.size(); i++) {
			if (cluster_assignment[i] == cluster_ix) {
				samples.push_back(used_sample_ix[i]);
			}
		}
		return samples;
	}

}
