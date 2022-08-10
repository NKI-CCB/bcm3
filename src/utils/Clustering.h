#pragma once

namespace bcm3 {

class RNG;

bool TreeCluster(const MatrixReal& distance, std::vector< std::set<size_t> >& clusters, Real cut_height);
void KMeans(int k, size_t nsamples, size_t dim, const Real** data, const Real* weights, int npass, int* assignment, Real* error);
bool NaiveKMeans(const MatrixReal& samples, int k, int npass, int niter, MatrixReal& centroids, std::vector<ptrdiff_t>& assignment, RNG& rng);

}
