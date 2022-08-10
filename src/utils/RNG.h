#pragma once

#include <boost/math/special_functions/log1p.hpp>
#include <random>

namespace bcm3 {

class RNG
{
public:
	RNG(unsigned long long seed = RNG::GetTimeBasedSeed());
	~RNG();

	void Seed(unsigned long long seed);

			unsigned int	GetUnsignedInt();
			unsigned int	GetUnsignedInt(unsigned int max);
			unsigned int    Sample(const VectorReal& probabilities);
			Real			GetReal();											//< Generates uniform random number in interval [0,1)
	inline	Real			GetUniform(const Real lower, const Real upper);		//< Generates uniform random number in interval [lower,upper)
	inline	Real			GetExponential(const Real mu);
			Real			GetGamma(const Real k, const Real theta);
			Real			GetNormal(const Real mu, const Real sigma);
			VectorReal		GetMultivariateUnitNormal(size_t D);
	inline  Real			GetBernoulli(const Real p);
	inline	Real			GetBeta(const Real a, const Real b);
	inline	Real			GetCauchy(const Real scale);
	inline	Real			GetHalfCauchy(const Real scale);

	static unsigned long long GetTimeBasedSeed();

private:
#if 1
	std::ranlux48 rl;
	double buffer[11];
	void fillbuffer();
#else
	ranluxpp* rl;
#endif
};

inline Real RNG::GetUniform(const Real lower, const Real upper)
{
	ASSERT(lower <= upper);
	return lower + GetReal() * (upper - lower);
}

inline Real RNG::GetExponential(const Real mu)
{
	const Real u = GetReal();
	return -mu * boost::math::log1p(-u);
}

inline Real RNG::GetBernoulli(const Real p)
{
	if (GetReal() < p) {
		return 1.0;
	} else {
		return 0.0;
	}
}

inline Real RNG::GetBeta(const Real a, const Real b)
{
	Real x1 = GetGamma(a, 1.0);
	Real x2 = GetGamma(b, 1.0);
	return x1 / (x1 + x2);
}

inline Real RNG::GetCauchy(const Real scale)
{
	Real x1 = GetReal();
	return scale * tan((Real)M_PI * (x1 - (Real)0.5));
}

inline Real RNG::GetHalfCauchy(const Real scale)
{
	return fabs(GetCauchy(scale));
}

}
