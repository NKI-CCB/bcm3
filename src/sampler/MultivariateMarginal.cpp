#include "Utils.h"
#include "MultivariateMarginal.h"
#include "ProbabilityDistributions.h"
#include "RNG.h"
#include "VectorUtils.h"

#include <boost/math/special_functions/gamma.hpp>

namespace bcm3 {

	MultivariateMarginal::MultivariateMarginal()
		: Type(D_Invalid)
		, N(0)
		, eta(1.0)
		, D(0)
		, log_normalization_constant(0.0)
	{
	}

	MultivariateMarginal::~MultivariateMarginal()
	{

	}

	bool MultivariateMarginal::CreateDirichlet()
	{
		Type = D_Dirichlet;
		alpha.resize(0);
		N = 0;
		return true;
	}

	bool MultivariateMarginal::CreateLKJCholesky(Real eta, size_t D)
	{
		if (eta <= 0.0) {
			LOGERROR("eta should be greater than 0");
			return false;
		}

		Type = D_LKJCholesky;
		this->eta = eta;
		this->D = D;
		N = 0;
		return true;
	}

	bool MultivariateMarginal::AddVariable(Real param, size_t& ix)
	{
		if (Type == D_Dirichlet) {
			ix = alpha.size();
			alpha.conservativeResize(alpha.size() + 1);
			alpha[ix] = param;
			N++;
		} else if (Type == D_LKJCholesky) {
			ix = N;
			N++;
		} else {
			LOGERROR("Invalid distribution type");
			return false;
		}

		return true;
	}

	bool MultivariateMarginal::Initialize()
	{
		if (N == 0) {
			LOGERROR("No variables have been added to the multivariate marginal");
			return false;
		}

		if (Type == D_Dirichlet) {
			ASSERT(N == alpha.size());

			Real sum = 0.0;
			Real lprod = 0.0;
			for (int i = 0; i < alpha.size(); i++) {
				sum += alpha[i];
				lprod += boost::math::lgamma(alpha[i]);
			}
			log_normalization_constant = boost::math::lgamma(sum) - lprod;
		} else if (Type == D_LKJCholesky) {
			if (N != D * (D - 1) / 2) {
				LOGERROR("Incorrect number of variables (%zu) for an LKJ distribution of dimension %zu", N, D);
				return false;
			}
		} else {
			LOGERROR("Invalid distribution type");
			return false;
		}

		return true;
	}

	bool MultivariateMarginal::Sample(Real* values, RNG* rng) const
	{
		if (Type == D_Dirichlet) {
			Real sum = 0.0;
			for (int i = 0; i < N; i++) {
				values[i] = rng->GetGamma(alpha[i], 1.0);
				sum += values[i];
			}

			Real inv = 1.0 / sum;
			for (int i = 0; i < N; i++) {
				values[i] *= inv;
			}
		} else if (Type == D_LKJCholesky) {
			Real beta = eta + 0.5 * (D - 1);
			int ix = 0;
			for (int i = 0; i < D - 1; i++) {
				beta -= 0.5;
				for (int j = i + 1; j < D; j++) {
					ASSERT(ix < N);
					values[ix++] = 2.0 * rng->GetBeta(beta, beta) - 1;
				}
			}
		} else {
			LOGERROR("Invalid distribution type");
			return false;
		}

		return true;
	}

	bool MultivariateMarginal::EvaluateLogPDF(Real& logp, const Real* values) const
	{
		if (Type == D_Dirichlet) {
			Real sum = 0.0;
			for (int i = 0; i < N; i++) {
				if (values[i] < 0.0) {
					logp = -std::numeric_limits<Real>::infinity();
					return true;
				} else if (values[i] > 1.0) {
					logp = -std::numeric_limits<Real>::infinity();
					return true;
				} else {
					sum += values[i];
				}
			}
			if (std::abs(sum - 1.0) > 1e-15) {
				logp = -std::numeric_limits<Real>::infinity();
				return true;
			} else {
				logp = 0.0;
				for (int i = 0; i < N; i++) {
					logp += (alpha[i] - 1) * log(values[i]);
				}
				logp += log_normalization_constant;
				return true;
			}
		} else if (Type == D_LKJCholesky) {
			Real beta = eta + 0.5 * (D - 1);
			int ix = 0;
			logp = 0.0;
			for (int i = 0; i < D - 1; i++) {
				beta -= 0.5;
				for (int j = i + 1; j < D; j++) {
					ASSERT(ix < N);
					logp += bcm3::LogPdfBeta(0.5 * (values[ix++] + 1), beta, beta) - 0.5;
				}
			}
		} else {
			LOGERROR("Invalid distribution type");
			return false;
		}

		return true;
	}
	
	bool MultivariateMarginal::EvaluateMarginalMean(size_t i, Real& mean) const
	{
		if (i >= N) {
			LOGERROR("Out of bounds");
			return false;
		}

		if (Type == D_Dirichlet) {
			Real sum_alpha = alpha.sum();
			mean = alpha(i) / sum_alpha;
		} else if (Type == D_LKJCholesky) {
			// Marginal mean is Beta rescaled to [-1,1]
			// Since the two beta parameters are the same, the distribution is symmetric and the mean is 0
			mean = 0.0;
		} else {
			LOGERROR("Invalid distribution type");
			return false;
		}

		return true;
	}

	bool MultivariateMarginal::EvaluateMarginalVariance(size_t i, Real& var) const
	{
		if (i >= N) {
			LOGERROR("Out of bounds");
			return false;
		}

		if (Type == D_Dirichlet) {
			Real sum_alpha = alpha.sum();
			var = alpha(i) * (sum_alpha - alpha(i)) / (sum_alpha * sum_alpha * (sum_alpha + 1.0));
		} else if (Type == D_LKJCholesky) {
			// Marginal mean is Beta rescaled to [-1,1]
			Real alpha = eta - 0.5 * (D - 2);
			var = bcm3::square(alpha) / (bcm3::square(2 * alpha) * (2 * alpha + 1));
		} else {
			LOGERROR("Invalid distribution type");
			return false;
		}

		return true;
	}

	Real MultivariateMarginal::GetLowerBound(size_t i) const
	{
		if (i >= N) {
			LOGERROR("Out of bounds");
			return std::numeric_limits<Real>::quiet_NaN();
		}

		if (Type == D_Dirichlet) {
			return 0.0;
		} else if (Type == D_LKJCholesky) {
			return -1.0;
		} else {
			LOGERROR("Invalid distribution type");
			return std::numeric_limits<Real>::quiet_NaN();
		}
	}

	Real MultivariateMarginal::GetUpperBound(size_t i) const
	{
		if (i >= N) {
			LOGERROR("Out of bounds");
			return std::numeric_limits<Real>::quiet_NaN();
		}

		if (Type == D_Dirichlet) {
			return 1.0;
		} else if (Type == D_LKJCholesky) {
			return 1.0;
		} else {
			LOGERROR("Invalid distribution type");
			return std::numeric_limits<Real>::quiet_NaN();
		}
	}

}
