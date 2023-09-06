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

	bool MultivariateMarginal::AddVariable(Real param, size_t& ix)
	{
		if (Type == D_Dirichlet) {
			ix = alpha.size();
			alpha.conservativeResize(alpha.size() + 1);
			alpha[ix] = param;
			N++;
		} else {
			LOGERROR("Invalid distribution type");
			return false;
		}

		return true;
	}

	bool MultivariateMarginal::Initialize()
	{
		if (Type == D_Dirichlet) {
			ASSERT(N == alpha.size());

			Real sum = 0.0;
			Real lprod = 0.0;
			for (int i = 0; i < alpha.size(); i++) {
				sum += alpha[i];
				lprod += boost::math::lgamma(alpha[i]);
			}
			log_normalization_constant = boost::math::lgamma(sum) - lprod;

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
		} else {
			LOGERROR("Invalid distribution type");
			return std::numeric_limits<Real>::quiet_NaN();
		}
	}

}
