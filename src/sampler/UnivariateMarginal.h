#pragma once

namespace bcm3 {

class RNG;

class UnivariateMarginal
{
public:
	UnivariateMarginal();

	bool Initialize(boost::property_tree::ptree pt);

	static std::unique_ptr<UnivariateMarginal> CreateUniform(Real lower, Real upper);
	static std::unique_ptr<UnivariateMarginal> CreateNormal(Real mu, Real sigma);
	static std::unique_ptr<UnivariateMarginal> CreateExponential(Real lambda);
	static std::unique_ptr<UnivariateMarginal> CreateGamma(Real k, Real theta);
	static std::unique_ptr<UnivariateMarginal> CreateBeta(Real a, Real b);
	static std::unique_ptr<UnivariateMarginal> CreateHalfCauchy(Real scale);
	static std::unique_ptr<UnivariateMarginal> CreateBetaPrime(Real a, Real b, Real scale);

	bool Sample(Real& value, RNG* rng) const;
	bool EvaluatePDF(Real& p, const Real value) const;
	bool EvaluateLogPDF(Real& logp, const Real value) const;
	bool EvaluateLogDerivative(Real& d, const Real value) const;
	bool EvaluateLog2ndDerivative(Real& d, const Real value) const;
	bool EvaluateMean(Real& mean) const;
	bool EvaluateVariance(Real& var) const;
	bool EvaluateCDF(Real& p, const Real value) const;
	bool EvaluateQuantile(Real& value, const Real p) const;
	Real GetLowerBound() const;
	Real GetUpperBound() const;

private:
	enum EDistribution
	{
		D_Invalid,
		D_Uniform,
		D_Normal,
		D_Exponential,
		D_Gamma,
		D_Beta,
		D_HalfCauchy,
		D_BetaPrime,
		D_ExponentialMix,
	};
	
	EDistribution PriorType;

	Real mu;			// Normal
	Real sigma;			// Normal
	Real a;				// Uniform, Beta, BetaPrime
	Real b;				// Uniform, Beta, BetaPrime
	Real lambda;		// Exponential
	Real k;				// Gamma
	Real theta;			// Gamma
	Real scale;			// HalfCauchy, BetaPrime
	Real mixing_coeff;
	Real lambda2;

	// Optimization
	Real OneOverTwoSigmaSq;
	Real LogOneOverSqrtTwoPiSigmaSq;
};

}
