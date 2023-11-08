#include "Utils.h"
#include "Clustering.h"
#include "RNG.h"
#include "VectorUtils.h"

extern "C" {
#include "../../dependencies/cluster-1.52a/src/cluster.h"
}

namespace bcm3 {

size_t find_cluster(Node* nodelist, std::vector< std::set<size_t> >& clusters, int i)
{
	int el = i;
	while (el < 0) {
		size_t ix = (-el) - 1;
		if (nodelist[ix].left >= 0) {
			el = nodelist[ix].left;
			break;
		}
		if (nodelist[ix].right >= 0) {
			el = nodelist[ix].right;
			break;
		}
		el = nodelist[ix].right;
	}
	for (size_t ci = 0; ci < clusters.size(); ci++) {
		if (clusters[ci].find(el) != clusters[ci].end()) {
			return ci;
		}
	}
	return std::numeric_limits<size_t>::max();
}

bool TreeCluster(const MatrixReal& distance, std::vector< std::set<size_t> >& clusters, Real cut_height)
{
	if (distance.cols() < 1) {
		return false;
	}

#if 0
	MatrixReal dcopy = distance;
	std::vector<Real*> pointers;
	pointers.resize(dcopy.cols(), NULL);
	for (size_t i = 0; i < (size_t)dcopy.cols(); i++) {
		pointers[i] = dcopy.data() + i * dcopy.rows();
	}
#else
	Eigen::MatrixXd dcopy = distance.cast<double>();
	std::vector<double*> pointers;
	pointers.resize(dcopy.cols(), NULL);
	for (size_t i = 0; i < (size_t)dcopy.cols(); i++) {
		pointers[i] = dcopy.data() + i * dcopy.rows();
	}
#endif

	Node* nodelist = treecluster((int)dcopy.cols(), (int)dcopy.cols(), NULL, NULL, NULL, 0, 0, 'm', &pointers[0]);
	
	clusters.clear();

	for (int i = 0; i < dcopy.cols() - 1; i++) {
		if (nodelist[i].distance < cut_height) {
			if (nodelist[i].left >= 0) {
				if (nodelist[i].right >= 0) {
					// New cluster
					clusters.push_back(std::set<size_t>());
					clusters.back().insert(nodelist[i].left);
					clusters.back().insert(nodelist[i].right);
				} else {
					// Add left to cluster containing right
					size_t ci = find_cluster(nodelist, clusters, nodelist[i].right);
					if (ci == std::numeric_limits<size_t>::max()) {
						LOGERROR("Error in clustering");
						return false;
					}
					clusters[ci].insert(nodelist[i].left);
				}
			} else {
				size_t ci = find_cluster(nodelist, clusters, nodelist[i].left);
				if (nodelist[i].right >= 0) {
					// Add right to cluster containing left
					clusters[ci].insert(nodelist[i].right);
					if (ci == std::numeric_limits<size_t>::max()) {
						LOGERROR("Error in clustering");
						return false;
					}
				} else {
					// Merge right into left
					size_t ci2 = find_cluster(nodelist, clusters, nodelist[i].right);
					clusters[ci].insert(clusters[ci2].begin(), clusters[ci2].end());
					clusters.erase(clusters.begin() + ci2);
				}
			}
		} else {
			if (nodelist[i].left >= 0) {
				clusters.push_back(std::set<size_t>());
				clusters.back().insert(nodelist[i].left);
			}
			if (nodelist[i].right >= 0) {
				clusters.push_back(std::set<size_t>());
				clusters.back().insert(nodelist[i].right);
			}
		}
	}

	return true;
}

void KMeans(int k, size_t nsamples, size_t dim, const Real** data, const Real* weights, int npass, int* assignment, Real* error)
{
	int ifound;
	std::vector<int*> mask(nsamples);
	for (size_t i = 0; i < nsamples; i++) {
		mask[i] = new int[dim];
		for (size_t j = 0; j < dim; j++) {
			mask[i][j] = 1;
		}
	}
	//TODO - single precision support
	kcluster(k, (int)nsamples, (int)dim, const_cast<Real**>(data), mask.data(), const_cast<Real*>(weights), 0, npass, 'a', 'e', assignment, error, &ifound);
	for (size_t i = 0; i < nsamples; i++) {
		delete[] mask[i];
	}
}

bool NaiveKMeans(const MatrixReal& samples, int k, int npass, int niter, MatrixReal& centroids, std::vector<ptrdiff_t>& assignment, RNG& rng)
{
	ptrdiff_t n = samples.rows();
	ptrdiff_t D = samples.cols();

	if (k <= 0) {
		LOGERROR("K-means error: Invalid k: %d", k);
		return false;
	}
	if (k >= n) {
		LOGERROR("K-means error: k larger or equal to n: k=%d n=%d", k, n);
		return false;
	}

	Real previous_distance = std::numeric_limits<Real>::infinity();
	MatrixReal stored_centroids;
	std::vector<ptrdiff_t> stored_assignment;

	for (int passi = 0; passi < npass; passi++) {
		// Random assignment
		assignment.resize(n);
		std::vector<ptrdiff_t> counts(k);
		for (ptrdiff_t i = 0; i < n; i++) {
			assignment[i] = rng.GetUnsignedInt((unsigned int)(k-1));
			counts[assignment[i]]++;
		}

		// Ensure we have one sample in every cluster
		for (int i = 0; i < k; i++) {
			// Since we have k < n, this should always terminate at some point
			while (counts[i] == 0) {
				// Select a random sample to re-assign to this cluster
				ptrdiff_t new_assignment = rng.GetUnsignedInt((unsigned int)(n-1));
				if (counts[new_assignment] > 1) {
					counts[new_assignment]--;
					assignment[new_assignment] = i;
					counts[i]++;
					break;
				}
			}
		}

		Real total_distance = 0.0;
		for (int iteri = 0; iteri < niter; iteri++) {
			// Calculate centroids
			centroids.setZero(k, D);
			VectorReal inc = VectorReal::Zero(k);
			for (ptrdiff_t i = 0; i < n; i++) {
				int cluster = assignment[i];
				inc(cluster) += 1.0;
				centroids.row(cluster) += (samples.row(i) - centroids.row(cluster)) / inc(cluster);
			}

			// Re-assign to closest centroid
			total_distance = 0.0;
			for (ptrdiff_t i = 0; i < n; i++) {
				if (counts[assignment[i]] > 1) {
					VectorReal dists(k);
					for (ptrdiff_t j = 0; j < k; j++) {
						VectorReal v = samples.row(i) - centroids.row(j);
						dists(j) = v.dot(v);
					}
					ptrdiff_t minc = minindex(dists);
					if (minc < 0) {
						LOGERROR("Error in K-means clustering; probably NaN samples?");
						return false;
					} else {
						if (minc != assignment[i]) {
							counts[assignment[i]]--;
							assignment[i] = minc;
							counts[minc]++;
						}
						total_distance += dists(minc);
					}
				} else {
					// Never make a cluster empty
					VectorReal v = samples.row(i) - centroids.row(assignment[i]);
					total_distance += v.dot(v);
				}
			}
		}

		if (total_distance < previous_distance) {
			previous_distance = total_distance;
			stored_centroids = centroids;
			stored_assignment = assignment;
		}
	}

	centroids = stored_centroids;
	assignment = stored_assignment;
	return true;
}


void KMeans(const MatrixReal& samples, int k, int npass, std::vector<int>& assignment, RNG* rng)
{
	ASSERT(samples.IsRowMajor == 0);

	int ifound;
	int nsamples = samples.cols();
	int dim = samples.rows();

	std::vector<int*> mask(nsamples);
	std::vector<double*> sample_p(nsamples);
	double* weights = new double[nsamples];
	for (size_t i = 0; i < nsamples; i++) {
		mask[i] = new int[dim];
		sample_p[i] = const_cast<double*>(samples.data()) + dim;
		for (size_t j = 0; j < dim; j++) {
			mask[i][j] = 1;
		}
		weights[i] = 1.0;
	}

	//TODO - single precision support
	Real error = std::numeric_limits<Real>::quiet_NaN();
	kcluster(k, nsamples, dim, sample_p.data(), mask.data(), weights, 0, npass, 'a', 'e', assignment.data(), &error, &ifound);
	for (size_t i = 0; i < nsamples; i++) {
		delete[] mask[i];
	}
	delete[] weights;
}

}
