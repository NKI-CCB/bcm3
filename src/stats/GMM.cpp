#include "Utils.h"
#include "checks.h"
#include "GMM.h"
#include "mvn.h"

namespace bcm3 {

	bool GMM::Set(const std::vector<VectorReal>& means, const std::vector<MatrixReal>& covariances, VectorReal weights)
	{
		ASSERT(means.size() == covariances.size());
		ASSERT(means.size() == weights.size());

		if (means.empty()) {
			return false;
		}

		ASSERT(means[0].size() == covariances[0].rows());
		ASSERT(means[0].size() == covariances[0].cols());

		components.resize(means.size());
		for (ptrdiff_t i = 0; i < means.size(); i++) {
			ASSERT(is_positive_semi_definite(covariances[i]));

			components[i].mean = means[i];
			components[i].covariance = covariances[i];
			components[i].covariance_llt.compute(components[i].covariance);
			if (components[i].covariance_llt.info() != Eigen::Success) {
				return false;
			}

			Real det = 0.0;
			for (size_t j = 0; j < covariances[i].rows(); j++) {
				det += log(components[i].covariance_llt.matrixL()(j, j));
			}
			components[i].logC = -det - 0.5 * covariances[i].rows() * log(2.0 * M_PI);
		}

		this->weights = weights;
		return true;
	}

	bool GMM::Fit(const MatrixReal& samples, size_t num_samples, size_t num_components, RNG& rng, Real ess_factor)
	{
		const size_t maxsteps = 100;
		const Real logl_epsilon = 1e-5;
		size_t D = samples.cols();

		Real logl;
		bool singular = false;
		bool converged = false;
		if (num_components == 1) {
			VectorReal responsibilities = VectorReal::Ones(num_samples);
			components.resize(1);
			CalculateMeanCovariance(samples, num_samples, responsibilities, components[0].mean, components[0].covariance, ess_factor);

			components[0].covariance_llt.compute(components[0].covariance);
			if (components[0].covariance_llt.info() != Eigen::Success) {
				return false;
			}

			Real det = 0.0;
			for (size_t i = 0; i < D; i++) {
				det += log(components[0].covariance_llt.matrixL()(i, i));
			}
			components[0].logC = -det - 0.5 * D * log(2.0 * M_PI);

			logl = 0.0;
			for (size_t j = 0; j < num_samples; j++) {
				Real p = LogPdfMVN(samples.row(j), components[0].mean, components[0].covariance_llt, components[0].logC);
				logl += p;
			}
			weights = VectorReal::Ones(1);

			converged = true;
		} else {
			if (num_samples < 2.0 * D * num_components) {
				// Need at least a bit more than p * K samples
				// Each component needs at least p samples for the regularization of the covariance estimation to work
				return false;
			}

			// Initial allocation of samples to components by k-means++
			MatrixReal responsibilities;
			if (!KMeanspp(samples, num_samples, num_components, rng, responsibilities)) {
				return false;
			}
			for (size_t i = 0; i < components.size(); i++) {
				CalculateMeanCovariance(samples, num_samples, responsibilities.col(i), components[i].mean, components[i].covariance, ess_factor);
			}
			weights.setConstant(num_components, 1.0 / num_components);

			// EM
			//EM_maximization(samples, num_samples, responsibilities);
			Real prev_logl = -std::numeric_limits<Real>::infinity();
			for (size_t i = 0; i < maxsteps; i++) {
				if (!EM_expectation(samples, num_samples, responsibilities, logl)) {
					singular = true;
					break;
				}

				if (logl < prev_logl) {
					// Some numerical issues or singularity?
					converged = true;
					break;
				} else if (logl - prev_logl < logl * logl_epsilon) {
					converged = true;
					break;
				} else {
					prev_logl = logl;
				}

				EM_maximization(samples, num_samples, responsibilities, ess_factor);
			}
		}

		size_t nparam = num_components * (D + D * (D + 1) / 2) + num_components - 1;
		aic = 2 * nparam - 2 * logl;
		return !singular;
	}

	Real GMM::LogPdf(const VectorReal& x)
	{
		Real logp = -std::numeric_limits<Real>::infinity();
		for (size_t i = 0; i < components.size(); i++) {
			VectorReal v = x - components[i].mean;
			components[i].covariance_llt.matrixL().solveInPlace(v);
			Real component_logp = components[i].logC - 0.5 * v.dot(v) + log(weights(i));
			logp = bcm3::logsum(logp, component_logp);
		}
		return logp;
	}

	VectorReal GMM::CalculateResponsibilities(const VectorReal& x)
	{
		VectorReal probs = VectorReal(components.size());
		for (size_t i = 0; i < components.size(); i++) {
			probs(i) = LogPdfMVN(x, components[i].mean, components[i].covariance_llt, components[i].logC) + log(weights(i));
		}

		// Make sure the components with the highest responsibility has a value that doesn't under/overflow after exponentation
		Real lsum = bcm3::logsum(probs);
		probs.array() -= lsum;

		// Transform and normalize
		probs = probs.array().exp();
		return probs / probs.sum();
	}

	bool GMM::KMeanspp(const MatrixReal& samples, size_t num_samples, size_t num_components, RNG& rng, MatrixReal& responsibilities)
	{
		if (num_components < 2) {
			return false;
		}

		components.resize(num_components);
		unsigned int ix = rng.GetUnsignedInt((unsigned int)num_samples-1);
		components[0].mean = samples.row(ix);
		std::set<unsigned int> used_samples;
		used_samples.insert(ix);
	
		for (size_t i = 1; i < num_components; i++) {
			VectorReal mindistsq = VectorReal::Constant(num_samples, 0);
			Real total = 0.0;
			for (size_t j = 0; j < num_samples; j++) {
				if (used_samples.find(j) != used_samples.end()) {
					continue;
				}

				VectorReal v = samples.row(j);
				Real mindistsqsearch = std::numeric_limits<Real>::max();
				for (size_t l = 0; l < i; l++) {
					VectorReal d = v - components[l].mean;
					Real distsq = d.dot(d);
					mindistsqsearch = std::min(mindistsqsearch, distsq);
				}
				mindistsq(j) = mindistsqsearch;
				total += mindistsq(j);
			}

			mindistsq /= total;
			unsigned int newix = rng.Sample(mindistsq);
			components[i].mean = samples.row(newix);
			used_samples.insert(newix);
		}

		responsibilities = MatrixReal::Zero(num_samples, num_components);
		for (size_t i = 0; i < num_samples; i++) {
			Real mindist = std::numeric_limits<Real>::max();
			size_t whichmin = std::numeric_limits<size_t>::max();
			const VectorReal& v = samples.row(i);
			for (size_t j = 0; j < num_components; j++) {
				VectorReal d = v - components[j].mean;
				Real distsq = d.dot(d);
				if (distsq < mindist) {
					mindist = distsq;
					whichmin = j;
				}
			}

			responsibilities(i, whichmin) = 1;
		}

		return true;
	}

	void GMM::CalculateMeanCovariance(const MatrixReal& samples, size_t num_samples, const VectorReal& responsibilities, VectorReal& mean, MatrixReal& covariance, Real ess_factor)
	{
		mean.setZero(samples.cols());
		covariance.setZero(samples.cols(), samples.cols());

		VectorReal x = VectorReal(samples.cols());
		VectorReal d = VectorReal(samples.cols());
		VectorReal d2 = VectorReal(samples.cols());

		Real wsum = 0;
		for (size_t i = 0; i < num_samples; i++) {
			x = samples.row(i);
			Real w = responsibilities(i);
			if (w >= std::numeric_limits<Real>::epsilon()) {
				wsum += w;
				d = x - mean;
				mean += (w / wsum) * d;
				d2 = x - mean;
	#if 0
				for (size_t j = 0; j < samples.rows(); j++) {
					for (size_t k = j; k < samples.rows(); k++) {
						covariance(j, k) += w * d(j) * d2(k);
						covariance(k, j) = covariance(j, k);
					}
				}
	#else
				covariance += w * d * d2.transpose();
	#endif
			}
		}

		if (wsum < 2.0) {
			// Too few samples; can't do anything
			covariance.setIdentity();
			return;
		}

		// Get the final covariance estimate
		covariance /= wsum;

		// Regularization
		Real n_eff = wsum / ess_factor;
		if (n_eff < samples.cols()) {
			// Less samples than variables; have to settle for a diagonal matrix
			VectorReal tmp = covariance.diagonal();
			covariance = tmp.asDiagonal();
		} else {
			// Stein's minimax shrinkage
			Eigen::SelfAdjointEigenSolver<MatrixReal> eig;
			eig.compute(covariance);
			VectorReal shrunk_eigval = eig.eigenvalues();
			for (size_t i = 0; i < shrunk_eigval.size(); i++) {
				shrunk_eigval((shrunk_eigval.size() - 1) - i) *= n_eff / (n_eff + samples.cols() + 1 - 2 * i);
			}
			covariance = eig.eigenvectors() * shrunk_eigval.asDiagonal() * eig.eigenvectors().transpose();
			covariance.diagonal().array() += 1e-8;
		}
	}

	void GMM::EM_maximization(const MatrixReal& samples, size_t num_samples, const MatrixReal& responsibilities, Real ess_factor)
	{
		for (size_t i = 0; i < components.size(); i++) {
			weights(i) = responsibilities.col(i).sum() / (Real)num_samples;
			CalculateMeanCovariance(samples, num_samples, responsibilities.col(i), components[i].mean, components[i].covariance, ess_factor);
		}
	}

	bool GMM::EM_expectation(const MatrixReal& samples, size_t num_samples, MatrixReal& responsibilities, Real& logl)
	{
		VectorReal sample_logl = VectorReal::Constant(num_samples, -std::numeric_limits<Real>::infinity());
		VectorReal v = VectorReal(samples.cols());

		for (size_t i = 0; i < components.size(); i++) {
			components[i].covariance_llt.compute(components[i].covariance);
			if (components[i].covariance_llt.info() != Eigen::Success) {
				return false;
			}

			Real det = 0.0;
			for (size_t j = 0; j < samples.cols(); j++) {
				det += log(components[i].covariance_llt.matrixL()(j, j));
			}
			Real logC = -det - 0.5 * samples.cols() * log(2.0 * M_PI);
			components[i].logC = logC;

			Real log_weight = log(weights(i));
			for (size_t j = 0; j < num_samples; j++) {
				//Real p = LogPdfMVN(samples.col(j), components[i].mean, components[i].covariance_llt, logC) + log(weights(i));
				v = samples.row(j);
				v -= components[i].mean;
				components[i].covariance_llt.matrixL().solveInPlace(v);
				Real p = logC - 0.5 * v.dot(v) + log_weight;

				responsibilities(j, i) = exp(p);
				sample_logl(j) = bcm3::logsum(sample_logl(j), p);
			}
		}
		logl = sample_logl.sum();

		for (size_t i = 0; i < num_samples; i++) {
			Real total_weight = responsibilities.row(i).sum();
			if (total_weight == 0) {
				// Responsibilities for all components are below double precision limits,
				// impossible to compare so just give equal responsibility to all components
				responsibilities.row(i).setConstant(1.0 / components.size());
			} else {
				responsibilities.row(i) /= total_weight;
			}
		}

		return true;
	}

	Real GMM::LogPdfMVN(const VectorReal& x, const VectorReal& mu, const Eigen::LLT<MatrixReal>& cov_decomposition, Real logC)
	{
		VectorReal v = x - mu;
		cov_decomposition.matrixL().solveInPlace(v);
		Real logp = logC - 0.5 * v.dot(v);
		return logp;
	}

}
