#pragma once

#include "RNG.h"

namespace bcm3 {

class GMM
{
public:
	bool Set(const std::vector<VectorReal>& means, const std::vector<MatrixReal>& covariances, VectorReal weights);
	bool Fit(const Eigen::MatrixXf& samples, size_t num_samples, size_t num_components, RNG& rng, Real ess_factor = 1.0);

	Real LogPdf(const VectorReal& x);
	VectorReal CalculateResponsibilities(const VectorReal& x);

	size_t GetNumComponents() const { return components.size(); }
	const VectorReal& GetWeights() const { return weights; }
	const VectorReal& GetMean(size_t component_ix) const { return components[component_ix].mean; }
	const MatrixReal& GetCovariance(size_t component_ix) const { return components[component_ix].covariance; }
	const Eigen::LLT<MatrixReal>& GetCovarianceDecomp(size_t component_ix) const { return components[component_ix].covariance_llt; }
	Real GetLogC(size_t component_ix) const { return components[component_ix].logC; }
	Real GetAIC() const { return aic; }

private:
	bool KMeanspp(const Eigen::MatrixXf& samples, size_t num_samples, size_t num_components, RNG& rng, MatrixReal& sample_weights);
	void CalculateMeanCovariance(const Eigen::MatrixXf& samples, size_t num_samples, const VectorReal& weights, VectorReal& mean, MatrixReal& covariance, Real ess_factor);
	void CalculateMeanCovarianceNERCOME(const Eigen::MatrixXf& samples, size_t num_samples, const VectorReal& weights, VectorReal& mean, MatrixReal& covariance, RNG& rng);
	void EM_maximization(const Eigen::MatrixXf& samples, size_t num_samples, const MatrixReal& sample_weights, Real ess_factor);
	bool EM_expectation(const Eigen::MatrixXf& samples, size_t num_samples, MatrixReal& sample_weights, Real& logl);
	Real LogPdfMVN(const VectorReal& x, const VectorReal& mu, const Eigen::LLT<MatrixReal>& cov_decomposition, Real logC);

	struct Component {
		VectorReal mean;
		MatrixReal covariance;
		Eigen::LLT<MatrixReal> covariance_llt;
		Real logC;
	};

	VectorReal weights;
	std::vector<Component> components;
	Real aic;
};

}
